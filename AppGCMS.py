import streamlit as st
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import io

# --- PAGE CONFIG ---
st.set_page_config(page_title="Forensic Lab Suite PRO", layout="wide")

# --- SHARED SESSION STATE ---
if 'unknowns_results' not in st.session_state:
    st.session_state['unknowns_results'] = []
if 'curve_calculated' not in st.session_state:
    st.session_state['curve_calculated'] = False

# --- TABS SETUP ---
tab1, tab2 = st.tabs(["📉 Quantification (GC-FID/MS)", "🔬 Identification (MS Isotopes)"])

# ==========================================
# TAB 1: CHROMATOGRAPHY (QUANTIFICATION)
# ==========================================
with tab1:
    st.header("1. Method Parameters")
    
    col1, col2 = st.columns(2)
    with col1:
        stock_unit = st.text_input("Stock concentration unit (e.g., mg/L)", value="mg/L")
        target_unit = st.text_input("Standard concentration unit (e.g., ug/mL)", value="ug/mL")
        vol_unit = st.text_input("Volume unit (e.g., mL)", value="mL")
    with col2:
        C1_raw = st.number_input(f"Stock Concentration ({stock_unit})", value=100.0)
        V2 = st.number_input(f"Final Standard Volume ({vol_unit})", value=10.0)
        c2_input = st.text_input("Standard Concentrations (comma separated)", value="0.1, 0.5, 1.0, 2.0, 5.0")

    st.markdown("### Internal Standard (IS)")
    use_is = st.checkbox("Use Internal Standard (IS)", value=True)
    if use_is:
        col_is1, col_is2 = st.columns(2)
        with col_is1:
            is_vol = st.number_input("IS Volume added to vial", value=10.0)
        with col_is2:
            is_unit = st.text_input("IS volume unit (e.g., uL)", value="uL")
    else:
        is_vol = 0.0
        is_unit = ""

    # --- Unit Conversion Logic ---
    unit_factors = {'ug/ml': 1.0, 'mg/l': 1.0, 'ppm': 1.0, 'mg/ml': 1000.0, 'g/l': 1000.0, 'ug/l': 0.001, 'ng/ml': 0.001, 'ppb': 0.001}
    s_unit_clean = stock_unit.lower().replace(" ", "")
    t_unit_clean = target_unit.lower().replace(" ", "")
    multiplier = unit_factors.get(s_unit_clean, 1.0) / unit_factors.get(t_unit_clean, 1.0)
    C1_converted = C1_raw * multiplier

    c2_list = [float(x.strip()) for x in c2_input.split(",") if x.strip()]
    v1_list = [(c * V2) / C1_converted for c in c2_list]

    # --- Pipetting Instructions ---
    st.header("2. Pipetting Instructions")
    pipette_dict = {
        f"Conc ({target_unit})": c2_list,
        f"Stock Volume ({vol_unit})": [round(v, 4) for v in v1_list]
    }
    if use_is:
        pipette_dict[f"Add IS ({is_unit})"] = [is_vol] * len(c2_list)
    pipette_dict[f"Fill up to ({vol_unit})"] = [V2] * len(c2_list)
    st.table(pd.DataFrame(pipette_dict))

    # --- Data Entry ---
    st.header("3. Chromatographic Data")
    if use_is:
        entry_df = pd.DataFrame({
            f"Standard ({target_unit})": c2_list, 
            "Peak Area (Analyte)": [0.0] * len(c2_list),
            "Peak Area (IS)": [0.0] * len(c2_list)
        })
    else:
        entry_df = pd.DataFrame({
            f"Standard ({target_unit})": c2_list, 
            "Peak Area": [0.0] * len(c2_list)
        })

    edited_df = st.data_editor(entry_df, use_container_width=True)

    if st.button("Calculate Calibration Curve", type="primary"):
        if use_is:
            std_areas = edited_df["Peak Area (Analyte)"].tolist()
            is_areas = edited_df["Peak Area (IS)"].tolist()
            y_vals = [s / i if i > 0 else 0 for s, i in zip(std_areas, is_areas)]
            st.session_state['std_areas_data'] = std_areas
            st.session_state['is_areas_data'] = is_areas
            st.session_state['ratios_data'] = y_vals
        else:
            y_vals = edited_df["Peak Area"].tolist()
            st.session_state['y_vals_data'] = y_vals

        slope, intercept, r_value, p_value, std_err = stats.linregress(c2_list, y_vals)
        st.session_state['slope'] = slope
        st.session_state['intercept'] = intercept
        st.session_state['r2'] = r_value**2
        st.session_state['curve_has_is'] = use_is 
        st.session_state['c2_list_data'] = c2_list
        st.session_state['y_plot_data'] = y_vals
        st.session_state['curve_calculated'] = True

    # --- Curve Visualization ---
    if st.session_state.get('curve_calculated', False):
        st.success(f"**Equation:** y = {st.session_state['slope']:.4f}x + {st.session_state['intercept']:.4f} | **R²** = {st.session_state['r2']:.4f}")
        
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.scatter(st.session_state['c2_list_data'], st.session_state['y_plot_data'], color='red', label='Calibration Points')
        ax.plot(st.session_state['c2_list_data'], [st.session_state['slope']*x + st.session_state['intercept'] for x in st.session_state['c2_list_data']], color='blue', label='Trendline')
        ax.set_xlabel(f"Concentration ({target_unit})")
        ax.set_ylabel("Ratio (Analyte/IS)" if st.session_state['curve_has_is'] else "Peak Area")
        ax.legend()
        st.pyplot(fig)

    # --- Sample Analysis ---
    st.header("4. Sample Analysis")
    if st.session_state.get('curve_calculated', False):
        cu1, cu2, cu3 = st.columns(3)
        with cu1: unk_name = st.text_input("Sample Name:", value="Sample 1")
        with cu2: unk_area = st.number_input("Analyte Peak Area:", value=0.0)
        with cu3: unk_is_area = st.number_input("IS Peak Area (in Sample):", value=1.0) if st.session_state['curve_has_is'] else 1.0

        if st.button("Add to Report"):
            y_for_calc = unk_area / unk_is_area if st.session_state['curve_has_is'] else unk_area
            res = (y_for_calc - st.session_state['intercept']) / st.session_state['slope']
            st.session_state['unknowns_results'].append({
                "Name": unk_name,
                "Area": unk_area, 
                "Ratio": round(y_for_calc, 4),
                f"Result ({target_unit})": round(res, 4),
                "Time": datetime.now().strftime('%H:%M:%S')
            })
        
        if st.session_state['unknowns_results']:
            st.table(pd.DataFrame(st.session_state['unknowns_results']))
            
            # --- REPORT EXPORT ---
            st.divider()
            if st.button("Clear Data"):
                st.session_state['unknowns_results'] = []
                st.rerun()

# ==========================================
# TAB 2: MS IDENTIFICATION (ISOTOPES)
# ==========================================
with tab2:
    st.header("MS Isotope Analysis & Forensic ID")
    
    # NEW: Retention Time for Forensic Correlation
    col_id1, col_id2 = st.columns(2)
    with col_id1:
        ms_rt = st.number_input("Observed Retention Time (RT) [min]:", value=0.0, step=0.01)
    with col_id2:
        mz_m = st.number_input("m/z of peak M (Molecular Ion):", min_value=1.0, value=146.0, step=1.0)

    st.divider()
    st.subheader("Peak Heights (Raw Detector Data)")
    h1, h2, h3, h4, h5 = st.columns(5)
    with h1: peak_m = st.number_input("M", value=8200.0, format="%.1f")
    with h2: peak_m1 = st.number_input("M+1 (C)", value=550.0, format="%.1f")
    with h3: peak_m2 = st.number_input("M+2", value=5400.0, format="%.1f") 
    with h4: peak_m4 = st.number_input("M+4", value=920.0, format="%.1f") 
    with h5: peak_m6 = st.number_input("M+6", value=0.0, format="%.1f")

    if st.button("Analyze Spectrum 🚀"):
        m_norm = 100.0
        m1_norm = (peak_m1 / peak_m) * 100.0 if peak_m > 0 else 0.0
        m2_norm = (peak_m2 / peak_m) * 100.0 if peak_m > 0 else 0.0
        m4_norm = (peak_m4 / peak_m) * 100.0 if peak_m > 0 else 0.0
        m6_norm = (peak_m6 / peak_m) * 100.0 if peak_m > 0 else 0.0
        
        st.divider()
        st.subheader("🔍 Identification Results")
        st.info(f"**RT Match:** Data correlated to peak at **{ms_rt} min**.")

        # 1. Carbon Estimation
        carbon_count = round(m1_norm / 1.1)
        st.write(f"**Estimated Carbon atoms:** {carbon_count}")
        
        # 2. Halogen Logic
        cl_count, br_count, s_count = 0, 0, 0
        tol = 10.0 
        
        if abs(m2_norm - 100.0) < tol and abs(m4_norm - 33.0) < tol:
            st.warning("🚨 **DETECTED: 3 CHLORINE ATOMS (Cl3)**")
            cl_count = 3
        elif abs(m2_norm - 200.0) < 15.0 and abs(m4_norm - 100.0) < tol:
            st.error("🚨 **DETECTED: 2 BROMINE ATOMS (Br2)**")
            br_count = 2
        elif abs(m2_norm - 133.0) < 15.0 and abs(m4_norm - 33.0) < tol:
            st.error("🚨 **DETECTED: 1 CHLORINE + 1 BROMINE**")
            cl_count, br_count = 1, 1
        elif abs(m2_norm - 66.0) < tol and abs(m4_norm - 11.0) < 5.0:
            st.warning("🚨 **DETECTED: 2 CHLORINE ATOMS (Cl2)**")
            cl_count = 2
        elif abs(m2_norm - 100.0) < tol:
            st.error("🚨 **DETECTED: 1 BROMINE ATOM (Br)**")
            br_count = 1
        elif abs(m2_norm - 33.0) < tol:
            st.warning("🚨 **DETECTED: 1 CHLORINE ATOM (Cl)**")
            cl_count = 1
        elif abs(m2_norm - 4.4) < 1.0:
            st.success("⚠️ **DETECTED: SULFUR (S)**")
            s_count = 1

        # 3. Formula Assembly & Nitrogen Rule
        st.divider()
        n_rule_hint = "Odd mass: Likely odd number of Nitrogens." if int(mz_m) % 2 != 0 else "Even mass."
        st.markdown(f"**Nitrogen Rule:** {n_rule_hint}")
        
        identified_mass = (carbon_count*12) + (cl_count*35) + (br_count*79) + (s_count*32)
        missing = int(mz_m - identified_mass)
        
        if missing >= 0:
            formula = f"C{carbon_count}"
            if cl_count > 0: formula += f" Cl{cl_count}"
            if br_count > 0: formula += f" Br{br_count}"
            if s_count > 0: formula += f" S{s_count}"
            st.success(f"### Proposed Formula Skeleton: {formula} (Residue: {missing} Da)")
        
        # 4. Visualization
        chart_data = pd.DataFrame({
            "m/z Offset": ["M", "M+1", "M+2", "M+4", "M+6"],
            "Intensity [%]": [m_norm, m1_norm, m2_norm, m4_norm, m6_norm]
        }).set_index("m/z Offset")
        st.bar_chart(chart_data)

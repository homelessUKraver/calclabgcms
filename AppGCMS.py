Python
import streamlit as st
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import io

# --- PAGE CONFIG ---
st.set_page_config(page_title="Forensic Lab Suite PRO", layout="wide")

# --- INITIALIZATION ---
if 'unknowns_results' not in st.session_state:
    st.session_state['unknowns_results'] = []
if 'ms_id_results' not in st.session_state:
    st.session_state['ms_id_results'] = []
if 'curve_calculated' not in st.session_state:
    st.session_state['curve_calculated'] = False

# --- TABS SETUP ---
tab1, tab2 = st.tabs(["Quantification (GC)", "Identification (MS Isotopes)"])

# ==========================================
# TAB 1: CHROMATOGRAPHY (QUANTIFICATION)
# ==========================================
with tab1:
    st.header("1. Method Parameters")
    col1, col2 = st.columns(2)
    with col1:
        stock_unit = st.text_input("Stock concentration unit", value="mg/L")
        target_unit = st.text_input("Standard concentration unit", value="ug/mL")
        vol_unit = st.text_input("Volume unit", value="mL")
    with col2:
        C1_raw = st.number_input(f"Stock Concentration ({stock_unit})", value=100.0)
        V2 = st.number_input(f"Final Standard Volume ({vol_unit})", value=10.0)
        c2_input = st.text_input("Standard Concentrations (comma separated)", value="0.1, 0.5, 1.0, 2.0, 5.0")

    use_is = st.checkbox("Use Internal Standard (IS)", value=True, key="is_check")
    if use_is:
        is_vol = st.number_input("IS Volume added", value=10.0)
        is_unit = st.text_input("IS unit", value="uL")
    else:
        is_vol, is_unit = 0.0, ""

    # Unit Conversion Logic
    unit_factors = {'ug/ml': 1.0, 'mg/l': 1.0, 'ppm': 1.0, 'mg/ml': 1000.0, 'g/l': 1000.0, 'ug/l': 0.001}
    multiplier = unit_factors.get(stock_unit.lower(), 1.0) / unit_factors.get(target_unit.lower(), 1.0)
    C1_converted = C1_raw * multiplier
    c2_list = [float(x.strip()) for x in c2_input.split(",") if x.strip()]
    v1_list = [(c * V2) / C1_converted for c in c2_list]

    st.header("2. Pipetting & Data")
    st.table(pd.DataFrame({f"Conc ({target_unit})": c2_list, f"Stock Vol ({vol_unit})": [round(v, 4) for v in v1_list], "Fill up to": V2}))

    entry_df = pd.DataFrame({f"Standard ({target_unit})": c2_list, "Analyte Area": 0.0, "IS Area": 1.0}) if use_is else pd.DataFrame({f"Standard ({target_unit})": c2_list, "Area": 0.0})
    edited_df = st.data_editor(entry_df, use_container_width=True)

    if st.button("Calculate Calibration Curve"):
        y_vals = (edited_df["Analyte Area"] / edited_df["IS Area"]).tolist() if use_is else edited_df["Area"].tolist()
        slope, intercept, r_value, _, _ = stats.linregress(c2_list, y_vals)
        st.session_state.update({'slope': slope, 'intercept': intercept, 'r2': r_value**2, 'curve_calculated': True, 'c2_list': c2_list, 'y_vals': y_vals})
        if use_is: 
            st.session_state['std_areas'] = edited_df["Analyte Area"].tolist()
            st.session_state['is_areas'] = edited_df["IS Area"].tolist()

    if st.session_state.get('curve_calculated'):
        st.success(f"y = {st.session_state['slope']:.4f}x + {st.session_state['intercept']:.4f} (R² = {st.session_state['r2']:.4f})")
        
        st.header("3. Sample Analysis")
        col_u1, col_u2, col_u3 = st.columns(3)
        with col_u1: u_name = st.text_input("Sample Name", value="Sample 1")
        with col_u2: u_area = st.number_input("Analyte Area", value=0.0, key="u_a")
        with col_u3: u_is = st.number_input("IS Area", value=1.0) if use_is else 1.0

        if st.button("Add to Report"):
            ratio = u_area / u_is
            res = (ratio - st.session_state['intercept']) / st.session_state['slope']
            st.session_state['unknowns_results'].append({"Name": u_name, "Area": u_area, "IS_Area": u_is, "Ratio": round(ratio, 4), f"Result ({target_unit})": round(res, 4)})

        if st.session_state['unknowns_results']:
            st.table(pd.DataFrame(st.session_state['unknowns_results']))

            # --- REPORT GENERATOR GC ---
            def gen_gc_report():
                output = io.StringIO()
                output.write(f"--- GC QUANT REPORT ---\nDate: {datetime.now()}\n")
                output.write(f"Equation: y={st.session_state['slope']}x+{st.session_state['intercept']}; R2={st.session_state['r2']}\n\n")
                pd.DataFrame(st.session_state['unknowns_results']).to_csv(output, index=False, sep=';')
                return output.getvalue()

            st.download_button("Download GC Report", gen_gc_report(), "GC_Report.csv", "text/csv")

# ==========================================
# TAB 2: MS IDENTIFICATION (ISOTOPES)
# ==========================================
with tab2:
    st.header("MS Isotope Analysis")
    col_m1, col_m2, col_m3 = st.columns(3)
    with col_m1: ms_name = st.text_input("Substance Name", value="Unknown")
    with col_m2: ms_rt = st.number_input("RT [min]", value=0.0)
    with col_m3: mz_m = st.number_input("m/z M", value=146.0)

    st.subheader("Peak Heights")
    h_cols = st.columns(5)
    p_m = h_cols[0].number_input("M", value=1000.0)
    p_m1 = h_cols[1].number_input("M+1", value=110.0)
    p_m2 = h_cols[2].number_input("M+2", value=0.0)
    p_m4 = h_cols[3].number_input("M+4", value=0.0)
    p_m6 = h_cols[4].number_input("M+6", value=0.0)

    if st.button("Analyze & Save MS Result"):
        m1_n = (p_m1 / p_m) * 100
        m2_n = (p_m2 / p_m) * 100
        c_count = round(m1_n / 1.1)
        
        # Halogen Logic Simplified
        hal = "None"
        if 30 < m2_n < 35: hal = "1 Chlorine (Cl)"
        elif 90 < m2_n < 110: hal = "1 Bromine (Br)"
        elif 60 < m2_n < 70: hal = "2 Chlorine (Cl2)"
        
        skeleton = f"C{c_count} " + (hal if hal != "None" else "")
        residue = mz_m - (c_count * 12) - (35 if "Cl" in hal else 0) - (79 if "Br" in hal else 0)

        st.session_state['ms_id_results'].append({
            "Name": ms_name, "RT": ms_rt, "m/z M": mz_m, 
            "C_atoms": c_count, "Halogens": hal, 
            "Skeleton": skeleton, "Residue": residue
        })

    if st.session_state['ms_id_results']:
        st.subheader("MS Results Table")
        st.table(pd.DataFrame(st.session_state['ms_id_results']))

        # --- REPORT GENERATOR MS ---
        def gen_ms_report():
            output = io.StringIO()
            output.write(f"--- MS IDENTIFICATION REPORT ---\nDate: {datetime.now()}\n\n")
            pd.DataFrame(st.session_state['ms_id_results']).to_csv(output, index=False, sep=';')
            return output.getvalue()

        st.download_button("Download MS Report", gen_ms_report(), "MS_Report.csv", "text/csv")

# --- CLEAR ALL ---
if st.sidebar.button("🗑️ Reset All Data"):
    for key in st.session_state.keys(): del st.session_state[key]
    st.rerun()

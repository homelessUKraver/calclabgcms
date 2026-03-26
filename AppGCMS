import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(page_title="Forensic Lab Suite", page_icon="⚖️", layout="wide")

# --- STYLE ---
st.markdown("""
    <style>
    .main { background-color: #f5f7f9; }
    .stMetric { background-color: #ffffff; padding: 15px; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
    </style>
    """, unsafe_view_safe_base64=True)

tab1, tab2 = st.tabs(["📉 QUANTIFICATION (GC-FID/MS)", "🔬 IDENTIFICATION (MS Isotope Tracker)"])

# --- TAB 1: FULL GC QUANTIFICATION ENGINE ---
with tab1:
    st.header("Quantitative Analysis & Calibration")
    
    with st.expander("🛠️ Calibration Curve Setup (Internal Standard Mode)", expanded=True):
        col_input1, col_input2 = st.columns(2)
        with col_input1:
            concs_raw = st.text_input("Standard Concentrations [mg/mL]:", "5, 10, 20, 40, 60")
        with col_input2:
            ratios_raw = st.text_input("Area Ratios (Area Target / Area IS):", "0.45, 0.92, 1.85, 3.65, 5.50")
        
        try:
            x_cal = np.fromstring(concs_raw, sep=',')
            y_cal = np.fromstring(ratios_raw, sep=',')
            if len(x_cal) == len(y_cal) and len(x_cal) > 1:
                slope, intercept = np.polyfit(x_cal, y_cal, 1)
                r_squared = np.corrcoef(x_cal, y_cal)[0,1]**2
                st.success(f"Curve Ready: **y = {slope:.4f}x + {intercept:.4f}** (R² = {r_squared:.4f})")
            else:
                st.warning("Awaiting valid data points...")
        except:
            st.error("Check data format (use commas to separate numbers)")

    st.subheader("🧪 Sample Analysis")
    c1, c2, c3 = st.columns(3)
    with c1:
        sample_area = st.number_input("Analyte Peak Area (Sample):", value=5200.0)
    with c2:
        is_area = st.number_input("Internal Standard Area (Sample):", value=1500.0)
    with c3:
        dilution = st.number_input("Dilution Factor:", value=1.0)

    if st.button("Calculate Concentration", key="calc_btn"):
        sample_ratio = sample_area / is_area
        # x = (y - b) / m
        conc_result = ((sample_ratio - intercept) / slope) * dilution
        
        st.divider()
        res_col1, res_col2 = st.columns(2)
        res_col1.metric("Calculated Concentration", f"{conc_result:.3f} mg/mL")
        res_col2.metric("Area Ratio", f"{sample_ratio:.3f}")

# --- TAB 2: MS IDENTIFICATION ENGINE ---
with tab2:
    st.header("Substance Identification & Isotope Pattern")
    
    # Forensic context: RT + Mass
    col_info1, col_info2 = st.columns(2)
    with col_info1:
        obs_rt = st.number_input("Observed Retention Time (RT) [min]:", value=8.25, step=0.01)
    with col_info2:
        target_mz = st.number_input("Molecular Ion (m/z) of peak M:", value=197.0)

    st.divider()
    st.subheader("Isotope Cluster Data")
    
    # Peak Heights for Isotopes
    h1, h2, h3, h4, h5 = st.columns(5)
    with h1: m_h = st.number_input("Height M", value=5200.0)
    with h2: m1_h = st.number_input("Height M+1", value=620.0)
    with h3: m2_h = st.number_input("Height M+2", value=1750.0)
    with h4: m4_h = st.number_input("Height M+4", value=0.0)
    with h5: m6_h = st.number_input("Height M+6", value=0.0)

    if st.button("Run Forensic Analysis 🚀"):
        # 1. Normalization
        m1_pct = (m1_h / m_h) * 100
        m2_pct = (m2_h / m_h) * 100
        m4_pct = (m4_h / m_h) * 100
        
        # 2. Logic Engines
        c_atoms = round(m1_pct / 1.1)
        
        st.info(f"**RT Analysis:** Peak detected at **{obs_rt} min**.")
        st.write(f"**Carbon Estimate:** Found approx. **{c_atoms}** Carbon atoms.")

        # Halogen Detection Logic
        cl_found, br_found = 0, 0
        if 25 < m2_pct < 45 and m4_pct < 5:
            st.warning("🚨 **IDENTIFIED: 1 Chlorine atom (Cl)**")
            cl_found = 1
        elif 90 < m2_pct < 110:
            st.error("🚨 **IDENTIFIED: 1 Bromine atom (Br)**")
            br_found = 1
        elif 60 < m2_pct < 80:
            st.warning("🚨 **IDENTIFIED: 2 Chlorine atoms (Cl2)**")
            cl_found = 2
        # ... (more logic can be added here)

        # 3. Nitrogen Rule & Formula Assembly
        st.divider()
        st.subheader("🧩 Molecular Formula Assembly")
        
        # Nitrogen Rule
        if int(target_mz) % 2 != 0:
            st.info("💡 **Nitrogen Rule:** Odd mass detected. Likely contains an odd number of Nitrogen atoms.")
            n_count = 1
        else:
            n_count = 0

        # Assembly (simplified for forensic hint)
        known_mass = (c_atoms * 12) + (cl_found * 35) + (br_found * 79) + (n_count * 14)
        missing = int(target_mz - known_mass)
        
        # Display Formula Hint
        formula = f"C{c_atoms}"
        if n_count > 0: formula += f" N{n_count}"
        if cl_found > 0: formula += f" Cl{cl_found}"
        if br_found > 0: formula += f" Br{br_found}"
        if missing > 0: formula += f" (H/O residue: {missing} Da)"
        
        st.success(f"**Proposed Skeleton:** {formula}")

        # 4. Viz
        st.subheader("Isotope Visualization")
        mz_axis = [target_mz, target_mz+1, target_mz+2, target_mz+3, target_mz+4]
        height_axis = [100.0, m1_pct, m2_pct, 0.0, m4_pct]
        chart_df = pd.DataFrame({"m/z": mz_axis, "Intensity [%]": height_axis}).set_index("m/z")
        st.bar_chart(chart_df)

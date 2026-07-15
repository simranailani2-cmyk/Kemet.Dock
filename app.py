import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

st.set_page_config(
    page_title="Kemet Dock",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for a minimalist, modern layout
st.markdown(
    """
    <style>
    .reportview-container {
        background: #f8f9fa;
    }
    .sidebar .sidebar-content {
        background: #ffffff;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("Kemet Dock")
st.markdown("### Molecular Docking Portal")

@st.cache_data
def load_data():
    return pd.read_csv("kemet_data.csv")

df = load_data()

st.header("Phase 1: Search Database")

col1, col2 = st.columns(2)

with col1:
    search_by = st.selectbox("Search by", ["Common Name", "Protein Target"])

with col2:
    if search_by == "Common Name":
        options = df["Common Name"].unique().tolist()
    else:
        options = df["Protein Target"].unique().tolist()

    selected_option = st.selectbox(f"Select {search_by}", options)

# Filter dataframe based on selection
if search_by == "Common Name":
    selected_data = df[df["Common Name"] == selected_option]
else:
    selected_data = df[df["Protein Target"] == selected_option]

def get_image_base64(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            return base64.b64encode(buffered.getvalue()).decode()
    except:
        pass
    return ""

st.markdown("---")

if not selected_data.empty:
    for idx, row in selected_data.iterrows():
        # Get molecule image
        img_b64 = get_image_base64(row['SMILES'])
        img_tag = f'<img src="data:image/png;base64,{img_b64}" style="width: 100%; border-radius: 8px;"/>' if img_b64 else ''

        card_html = f"""
        <div style="
            background: white;
            padding: 24px;
            border-radius: 12px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            margin-bottom: 24px;
            border: 1px solid #e1e4e8;">

            <h2 style="color: #2c3e50; margin-top: 0;">{row['Common Name']}
                <span style="font-size: 0.6em; color: #7f8c8d; font-style: italic;">({row['Botanical Name']})</span>
            </h2>

            <div style="display: flex; flex-wrap: wrap; gap: 24px;">
                <!-- Left Column: Molecule Image -->
                <div style="flex: 1; min-width: 250px; text-align: center;">
                    <div style="background: #f8f9fa; padding: 16px; border-radius: 8px; border: 1px solid #e1e4e8;">
                        {img_tag}
                        <p style="margin-top: 12px; font-weight: bold; color: #34495e;">{row['Active Phytochemical']}</p>
                    </div>
                </div>

                <!-- Right Column: Details -->
                <div style="flex: 2; min-width: 300px;">
                    <div style="margin-bottom: 16px;">
                        <h4 style="color: #e67e22; border-bottom: 2px solid #f39c12; padding-bottom: 4px;">Historical Context</h4>
                        <p><strong>Ancient Symbols:</strong> <span style="font-size: 1.2em;">{row['Ancient Symbols']}</span></p>
                        <p><strong>Claim Transliteration:</strong> <em>{row['Claim Transliteration']}</em></p>
                        <p><strong>Papyrus Source:</strong> {row['Papyrus Source']}</p>
                        <p><strong>Ancient Claim:</strong> {row['Ancient Claim']}</p>
                    </div>

                    <div>
                        <h4 style="color: #2980b9; border-bottom: 2px solid #3498db; padding-bottom: 4px;">Modern Data</h4>
                        <p><strong>Family:</strong> {row['Family']}</p>
                        <p><strong>Modern Use:</strong> {row['Modern Use']}</p>
                        <p><strong>Protein Target:</strong> {row['Protein Target']}</p>
                        <p><strong>PDB ID:</strong> {row['PDB ID']}</p>
                        <p style="word-break: break-all;"><strong>SMILES:</strong> {row['SMILES']}</p>
                    </div>
                </div>
            </div>
        </div>
        """

        st.markdown(card_html, unsafe_allow_html=True)

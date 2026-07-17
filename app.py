import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO
from docking_engine import fetch_receptor, prepare_ligand, smart_cavity_finder, run_vina_docking

import parse_pdbqt
import adme_profiler
import logging


def clear_interaction_state():
    for key in list(st.session_state.keys()):
        if key.startswith('interactions_df_') or key.startswith('highlight_atoms_'):
            st.session_state[key] = None

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
    search_by = st.selectbox("Search by", ["Common Name", "Protein Target"], on_change=clear_interaction_state)

with col2:
    if search_by == "Common Name":
        options = df["Common Name"].unique().tolist()
    else:
        options = df["Protein Target"].unique().tolist()

    selected_option = st.selectbox(f"Select {search_by}", options, on_change=clear_interaction_state)

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
                        <p><strong>Active Phytochemical:</strong> {row['Active Phytochemical']}</p>
                        <p><strong>Protein Target:</strong> {row['Protein Target']}</p>
                        <p><strong>PDB ID:</strong> {row['PDB ID']}</p>
                        <p style="word-break: break-all;"><strong>SMILES:</strong> {row['SMILES']}</p>
                    </div>
                </div>
            </div>
        </div>
        """


        components.html(card_html, height=600, scrolling=True)

        # Add a setup expander for Grid Box controls
        with st.expander(f"Docking Setup for {row['Common Name']}", expanded=True):
            pdb_id = row['PDB ID']
            smiles = row['SMILES']

            if st.button(f"Initialize Receptor {pdb_id}", key=f"init_{idx}"):
                with st.spinner(f"Fetching receptor {pdb_id} and calculating cavity..."):
                    output_pdb = f"{pdb_id}.pdb"
                    receptor_pdbqt = fetch_receptor(pdb_id, output_pdb)
                    if receptor_pdbqt:
                        center, dims = smart_cavity_finder(output_pdb)
                        st.session_state[f'center_{idx}'] = center.tolist() if hasattr(center, 'tolist') else list(center)
                        st.session_state[f'dims_{idx}'] = dims.tolist() if hasattr(dims, 'tolist') else list(dims)
                        st.session_state[f'rec_pdbqt_{idx}'] = receptor_pdbqt
                        st.session_state[f'setup_done_{idx}'] = True
                    else:
                        st.error("Failed to fetch receptor.")

            if st.session_state.get(f'setup_done_{idx}', False):
                st.markdown("**Grid Box Parameters**")

                center = st.session_state[f'center_{idx}']
                dims = st.session_state[f'dims_{idx}']

                col_cx, col_cy, col_cz = st.columns(3)
                cx = col_cx.number_input("Center X", value=float(center[0]), format="%.3f", key=f"cx_{idx}")
                cy = col_cy.number_input("Center Y", value=float(center[1]), format="%.3f", key=f"cy_{idx}")
                cz = col_cz.number_input("Center Z", value=float(center[2]), format="%.3f", key=f"cz_{idx}")

                col_sx, col_sy, col_sz = st.columns(3)
                sx = col_sx.number_input("Size X", value=float(dims[0]), format="%.3f", key=f"sx_{idx}")
                sy = col_sy.number_input("Size Y", value=float(dims[1]), format="%.3f", key=f"sy_{idx}")
                sz = col_sz.number_input("Size Z", value=float(dims[2]), format="%.3f", key=f"sz_{idx}")

                if st.button(f"Run Vina Docking", key=f"dock_{idx}"):
                    with st.spinner("Preparing docking pipeline..."):
                        receptor_pdbqt = st.session_state[f'rec_pdbqt_{idx}']

                        st.write(f"Preparing ligand...")
                        ligand_pdbqt, uff_delta = prepare_ligand(smiles, "ligand.pdbqt")

                        if receptor_pdbqt and ligand_pdbqt:
                            st.write(f"Running AutoDock Vina...")
                            vina_output = run_vina_docking(receptor_pdbqt, ligand_pdbqt, [cx, cy, cz], [sx, sy, sz])

                            st.success("Docking complete!")

                            # Parse vina_output
                            lines = vina_output.split('\n')
                            data = []
                            parsing = False
                            for line in lines:
                                if line.startswith('-----+------------+----------+----------'):
                                    parsing = True
                                    continue
                                if parsing:
                                    parts = line.split()
                                    if len(parts) == 4 and parts[0].isdigit():
                                        try:
                                            mode = int(parts[0])
                                            affinity = float(parts[1])
                                            rmsd_ub = float(parts[3])
                                            data.append({
                                                'Binding Mode': mode,
                                                'Affinity (kcal/mol)': affinity,
                                                'RMSD': rmsd_ub
                                            })
                                        except ValueError:
                                            pass
                                    elif len(parts) == 0 or 'Writing' in line:
                                        break

                            if data:
                                st.dataframe(pd.DataFrame(data), hide_index=True, use_container_width=False)
                                st.session_state[f'docking_data_{idx}'] = data
                                st.session_state[f'docking_done_{idx}'] = True
                                st.session_state[f'pdb_id_{idx}'] = pdb_id
                                st.session_state[f'smiles_{idx}'] = smiles
                                st.session_state[f'uff_delta_{idx}'] = uff_delta
                            else:
                                st.error("Could not parse Vina output.")
                        else:
                            st.error("Failed to prepare receptor or ligand.")


        if st.session_state.get(f'docking_done_{idx}', False):
            st.markdown("### Select Docking Pose")
            data = st.session_state[f'docking_data_{idx}']
            options = [f"Mode {d['Binding Mode']} (Affinity: {d['Affinity (kcal/mol)']} kcal/mol)" for d in data]
            selected_mode_str = st.selectbox("Generated Poses", options, key=f"pose_select_{idx}")
            selected_idx = options.index(selected_mode_str)
            selected_mode_data = data[selected_idx]

            poses = parse_pdbqt.extract_poses("ligand_out.pdbqt")

            if poses and selected_idx < len(poses):
                selected_pose_str = poses[selected_idx]

                # Metric Card
                uff_delta = st.session_state.get(f'uff_delta_{idx}', 0.0)
                receptor_pdbqt = f"{st.session_state[f'pdb_id_{idx}']}.pdbqt"
                interactions_data = parse_pdbqt.calc_interactions(selected_pose_str.split('\n'), receptor_pdbqt)

                interactions_df = pd.DataFrame(interactions_data) if interactions_data is not None else pd.DataFrame()

                if not interactions_df.empty and "Receptor Residue" in interactions_df.columns:
                    interacting_res = list(interactions_df["Receptor Residue"].unique())
                else:
                    interacting_res = []

                col1, col2, col3 = st.columns(3)
                col1.metric("Pose Affinity", f"{selected_mode_data['Affinity (kcal/mol)']} kcal/mol")
                col2.metric("UFF Minimization Delta", f"{uff_delta:.2f} kcal/mol")
                col3.metric("Interacting Residues", str(len(interacting_res)))

                st.write(f"Interacting receptor residues: {', '.join(interacting_res) if interacting_res else 'None'}")

                st.markdown("### Interaction Analysis")
                if not interactions_df.empty:
                    int_col1, int_col2 = st.columns([1, 1])
                    with int_col1:
                        st.dataframe(interactions_df, hide_index=True)
                    with int_col2:
                        try:
                            mol = Chem.MolFromSmiles(st.session_state[f'smiles_{idx}'])
                            if mol:
                                highlight_atoms = []
                                for atom_str in interactions_df["Ligand Atom"].unique():
                                    try:
                                        atom_idx = int(atom_str.split(" ")[1]) - 1
                                        if atom_idx < mol.GetNumAtoms():
                                            highlight_atoms.append(atom_idx)
                                        else:
                                            logging.warning(f"Atom index {atom_idx} out of bounds for molecule.")
                                    except ValueError:
                                        pass
                                st.session_state[f'interactions_df_{idx}'] = interactions_df
                                st.session_state[f'highlight_atoms_{idx}'] = highlight_atoms
                                img = Draw.MolToImage(mol, highlightAtoms=highlight_atoms, size=(400, 400))
                                st.image(img)
                        except Exception as e:
                            st.write("2D interaction map unavailable.")
                else:
                    st.info('No significant interactions found for this phytochemical.')
                    try:
                        mol = Chem.MolFromSmiles(st.session_state[f'smiles_{idx}'])
                        if mol:
                            img = Draw.MolToImage(mol, size=(400, 400))
                            st.image(img)
                    except Exception as e:
                        pass

                # 3Dmol.js rendering
                st.markdown("### 3D Interaction Viewer")

                style_options = ['cartoon', 'sphere', 'stick', 'line', 'cross']
                style_col1, style_col2 = st.columns(2)
                with style_col1:
                    receptor_style = st.selectbox("Receptor Style", style_options, index=0, key=f"rec_style_{idx}")
                with style_col2:
                    ligand_style = st.selectbox("Ligand Style", style_options, index=2, key=f"lig_style_{idx}")

                show_surface = st.checkbox("Show Pocket Cavity Mesh", value=True, key=f"surf_{idx}")

                with open(receptor_pdbqt, 'r') as f:
                    receptor_data = f.read()

                viewer_html = f'''
                <div id="container-{idx}" style="height: 500px !important; width: 100% !important; display: block; position: relative;" class="viewer_3Dmoljs"
                     data-backgroundcolor="0xffffff" data-style="stick"></div>
                <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                <script>
                    var initViewer_{idx} = setInterval(function() {{
                        if (typeof $3Dmol !== 'undefined') {{
                            clearInterval(initViewer_{idx});
                            var viewer = $3Dmol.createViewer("container-{idx}", {{defaultcolors: $3Dmol.rasmolElementColors}});
                            var receptor_data = `{receptor_data}`;
                            var ligand_data = `{selected_pose_str}`;

                            viewer.addModel(receptor_data, "pdb");
                            viewer.setStyle({{model: 0}}, {{{receptor_style}: {{color: 'spectrum'}} }});

                            if ({'true' if show_surface else 'false'}) {{
                                viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.8, color: 'white'}}, {{model: 0}});
                            }}

                            viewer.addModel(ligand_data, "pdb");
                            viewer.setStyle({{model: 1}}, {{{ligand_style}: {{colorscheme: 'greenCarbon'}} }});

                            viewer.zoomTo();
                            viewer.render();
                        }}
                    }}, 100);
                </script>
                '''
                components.html(viewer_html, height=550)

                # Phase 2 and 3: ADMET & Design
                st.markdown("---")
                st.header("ADMET & Design")
                smiles = st.session_state[f'smiles_{idx}']
                variants = adme_profiler.generate_variants(smiles)

                orig_adme = adme_profiler.get_admet(smiles)
                adme_data = []
                if orig_adme:
                    orig_adme["Molecule"] = "Original Phytochemical"
                    adme_data.append(orig_adme)

                if variants:
                    var_adme = adme_profiler.get_admet(variants[0])
                    if var_adme:
                        var_adme["Molecule"] = "Redesign Variant 1"
                        adme_data.append(var_adme)

                if adme_data:
                    df_adme = pd.DataFrame(adme_data)
                    cols = ['Molecule'] + [c for c in df_adme.columns if c != 'Molecule']
                    df_adme = df_adme[cols]
                    st.dataframe(df_adme, hide_index=True)

                    # Side-by-side visual comparison
                    comp_col1, comp_col2 = st.columns(2)
                    with comp_col1:
                        mol1 = Chem.MolFromSmiles(smiles)
                        if mol1:
                            st.image(Draw.MolToImage(mol1, size=(400, 400)), caption='Original Phytochemical')
                    with comp_col2:
                        if variants:
                            mol2 = Chem.MolFromSmiles(variants[0])
                            if mol2:
                                st.image(Draw.MolToImage(mol2, size=(400, 400)), caption='Redesign Variant 1')

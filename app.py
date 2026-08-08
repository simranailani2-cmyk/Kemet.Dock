import os
import urllib.parse
import json
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


def generate_html_report(plant_name, smiles, receptor_name, pdb_id, df_adme_html, interactions_df_html, interactions_var_df_html=None):
    variant_section = ""
    if interactions_var_df_html:
        variant_section = f'''
        <div class="data-section">
            <h2>Redesign Variant Bond Information</h2>
            {interactions_var_df_html}
        </div>
        '''

    html = f'''
    <html>
    <head>
    <style>
        body {{
            font-family: sans-serif;
            background-color: #f5f5dc; /* Egyptian sand/brown color palette */
            color: #5c4033; /* Dark brown text */
            padding: 20px;
        }}
        h1 {{
            color: #8b4513; /* SaddleBrown */
            text-align: center;
        }}
        .data-section {{
            background-color: #ffffff;
            border: 1px solid #d2b48c; /* Tan border */
            padding: 15px;
            margin-bottom: 20px;
            border-radius: 5px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }}
        th, td {{
            border: 1px solid #d2b48c;
            padding: 8px;
            text-align: left;
        }}
        th {{
            background-color: #f0e68c; /* Khaki */
        }}
        .footer {{
            text-align: right;
            font-size: small;
            color: gray;
            margin-top: 30px;
        }}
    </style>
    </head>
    <body>
        <h1>Kemet Dock: Egyptian Phytochemical Analysis</h1>

        <div class="data-section">
            <h2>Session Information</h2>
            <p><strong>Active Phytochemical Name:</strong> {plant_name}</p>
            <p><strong>SMILES String:</strong> {smiles}</p>
            <p><strong>Receptor Name:</strong> {receptor_name}</p>
            <p><strong>PDB ID:</strong> {pdb_id}</p>
        </div>

        <div class="data-section">
            <h2>Bond Information</h2>
            {interactions_df_html}
        </div>
        {variant_section}

        <div class="data-section">
            <h2>ADMET Properties</h2>
            {df_adme_html}
        </div>

        <div class="footer">
            by Simran Ailani
        </div>
    </body>
    </html>
    '''
    return html

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


# Unified Header: Full Title block on left, Cat and Pharaoh on the right, bottom-aligned
col_title, col_spacer, col_cat, col_phar = st.columns(
    [6, 1, 1, 1],
    vertical_alignment="bottom"
)

with col_title:
    # Main Header Image and Title
    st.title("𓆎𓏏𓈇 Kemet Dock 𓆎𓏏𓈇")
    st.markdown("### Molecular Docking Portal")

# col_spacer remains empty to push the images to the right side

with col_cat:
    # Cat on the left side of the Pharaoh
    st.image("https://upload.wikimedia.org/wikipedia/commons/a/ae/Egyptian_Cat.svg", width=80)

with col_phar:
    # Pharaoh on the far right
    st.image("https://upload.wikimedia.org/wikipedia/commons/a/ab/Pharaoh_in_war.svg", width=80)
# Sidebar Images - Deities
st.sidebar.image("https://upload.wikimedia.org/wikipedia/commons/7/71/Seshat.svg", width=150)
st.sidebar.markdown("<br>", unsafe_allow_html=True)
st.sidebar.image("https://upload.wikimedia.org/wikipedia/commons/c/c3/Thoth.svg", width=150)

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


        components.iframe(f'data:text/html;charset=utf-8,{urllib.parse.quote(card_html)}', height=600, scrolling=True)


        if st.button("Reset Environment", key=f"reset_{idx}"):
            keys_to_clear = [
                f'setup_done_{idx}', f'center_{idx}', f'dims_{idx}', f'rec_pdbqt_{idx}',
                f'docking_done_{idx}', f'docking_data_{idx}', f'pdb_id_{idx}',
                f'smiles_{idx}', f'uff_delta_{idx}', f'interactions_df_{idx}',
                f'highlight_atoms_{idx}'
            ]
            for k in keys_to_clear:
                if k in st.session_state:
                    del st.session_state[k]

            temp_files = [
                'ligand.pdbqt',
                'ligand_out.pdbqt',
                f"{row['PDB ID']}.pdb",
                f"{row['PDB ID']}.pdbqt"
            ]
            for tf in temp_files:
                if os.path.exists(tf):
                    try:
                        os.remove(tf)
                    except OSError:
                        pass

            st.rerun()

        st.header("Phase 2: Receptor & Grid Setup")

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

                st.header("Phase 3: Docking Simulation & Results")

                if st.button(f"Run Vina Docking", key=f"dock_{idx}"):
                    with st.spinner("Preparing docking pipeline..."):
                        receptor_pdbqt = st.session_state[f'rec_pdbqt_{idx}']

                        st.write(f"Preparing ligand...")
                        ligand_pdbqt, uff_delta = prepare_ligand(smiles, "ligand.pdbqt")

                        if receptor_pdbqt and ligand_pdbqt:
                            st.write(f"Running AutoDock Vina...")
                            progress_bar = st.progress(0, text="Starting docking...")
                            import subprocess
                            import os

                            vina_path = os.path.abspath("./vina") if os.path.exists("./vina") else "vina"
                            vina_command = [
                                vina_path, "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt),
                                "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                                "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
                                "--exhaustiveness", "16", "--out", "docking_poses.pdbqt"
                            ]

                            process = subprocess.Popen(vina_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                            output_log, progress_count = [], 0

                            while True:
                                char = process.stdout.read(1).decode("utf-8", errors="ignore")
                                if not char: break
                                output_log.append(char)
                                if char == '*':
                                    progress_count += 1
                                    if progress_bar:
                                        pct = min(100, int((progress_count / 50) * 100))
                                        progress_bar.progress(pct / 100.0, text=f"Exploring binding modes... {pct}%")

                            process.wait()
                            if process.returncode != 0:
                                st.error("Engine failed!")

                            vina_output = "".join(output_log)
                            progress_bar.empty()

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
                    interacting_res_data = interactions_df.to_dict("records")
                    interacting_res = list(interactions_df["Receptor Residue"].unique())
                else:
                    interacting_res_data = []
                    interacting_res = []

                col1, col2, col3 = st.columns(3)
                col1.metric("Pose Affinity", f"{selected_mode_data['Affinity (kcal/mol)']} kcal/mol")
                col2.metric("UFF Minimization Delta", f"{uff_delta:.2f} kcal/mol")
                col3.metric("Interacting Residues", str(len(interacting_res)))

                st.write(f"Interacting receptor residues: {', '.join(interacting_res) if interacting_res else 'None'}")

                st.markdown("### Interaction Analysis")
                if not interactions_df.empty:
                    st.dataframe(interactions_df, hide_index=True)
                else:
                    st.info('No significant interactions found for this phytochemical.')

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


                                var interactionsData = {json.dumps(interacting_res_data)};
                                var interactingRes = {json.dumps(interacting_res)};
                                interactionsData.forEach(function(interaction) {{
                                    var res = interaction["Receptor Residue"];
                                    var rx = interaction["Receptor XYZ"][0];
                                    var ry = interaction["Receptor XYZ"][1];
                                    var rz = interaction["Receptor XYZ"][2];
                                    var lx = interaction["Ligand XYZ"][0];
                                    var ly = interaction["Ligand XYZ"][1];
                                    var lz = interaction["Ligand XYZ"][2];

                                    // Add dashed line
                                    viewer.addCylinder({{
                                        start: {{x: lx, y: ly, z: lz}},
                                        end: {{x: rx, y: ry, z: rz}},
                                        radius: 0.1,
                                        color: 'yellow',
                                        dashed: true,
                                        fromCap: 1,
                                        toCap: 1
                                    }});

                                    // Anchor label to receptor atom
                                    viewer.addLabel(res, {{
                                        position: {{x: rx, y: ry, z: rz}},
                                        fontColor: 'white',
                                        backgroundColor: 'black',
                                        backgroundOpacity: 0.5,
                                        fontsize: 12
                                    }});
                                }});


                            viewer.zoomTo();
                            viewer.render();
                        }}
                    }}, 100);
                </script>
                '''
                components.iframe(f'data:text/html;charset=utf-8,{urllib.parse.quote(viewer_html)}', height=550)

                smiles = st.session_state[f'smiles_{idx}']
                variants = adme_profiler.generate_variants(smiles)

                if variants:
                    st.markdown("---")
                    st.header("Phase 4: Redesign Variant Docking")

                    selected_variant = st.selectbox("Select Redesign Variant", variants, key=f"variant_select_{idx}")

                    if st.button(f"Initiate Docking for Redesigned Variant", key=f"dock_var_{idx}"):
                        try:
                            center = st.session_state[f'center_{idx}']
                            dims = st.session_state[f'dims_{idx}']
                            cx = st.session_state.get(f"cx_{idx}", center[0])
                            cy = st.session_state.get(f"cy_{idx}", center[1])
                            cz = st.session_state.get(f"cz_{idx}", center[2])
                            sx = st.session_state.get(f"sx_{idx}", dims[0])
                            sy = st.session_state.get(f"sy_{idx}", dims[1])
                            sz = st.session_state.get(f"sz_{idx}", dims[2])

                            ligand_var_pdbqt, uff_delta_var = prepare_ligand(selected_variant, "ligand_var.pdbqt")
                            progress_bar_var = st.progress(0, text="Starting docking for variant...")
                            vina_output_var = run_vina_docking(receptor_pdbqt, ligand_var_pdbqt, [cx, cy, cz], [sx, sy, sz], progress_bar=progress_bar_var)
                            progress_bar_var.empty()

                            lines_var = vina_output_var.split('\n')
                            data_var = []
                            parsing_var = False
                            for line in lines_var:
                                if line.startswith('-----+------------+----------+----------'):
                                    parsing_var = True
                                    continue
                                if parsing_var:
                                    parts_var = line.split()
                                    if len(parts_var) == 4 and parts_var[0].isdigit():
                                        try:
                                            mode_var = int(parts_var[0])
                                            affinity_var = float(parts_var[1])
                                            rmsd_ub_var = float(parts_var[3])
                                            data_var.append({
                                                'Binding Mode': mode_var,
                                                'Affinity (kcal/mol)': affinity_var,
                                                'RMSD': rmsd_ub_var
                                            })
                                        except ValueError:
                                            pass
                                    elif len(parts_var) == 0 or 'Writing' in line:
                                        break

                            if data_var:
                                st.session_state[f'docking_var_data_{idx}'] = data_var
                                st.session_state[f'docking_var_done_{idx}'] = True
                                st.session_state[f'uff_delta_var_{idx}'] = uff_delta_var
                                st.success("Variant docking complete!")
                            else:
                                st.error("Could not parse Vina output for variant.")
                        except Exception as e:
                            st.error(f"Error during variant docking: {e}")

                    if st.session_state.get(f'docking_var_done_{idx}', False):
                        data_var = st.session_state[f'docking_var_data_{idx}']
                        options_var = [f"Mode {d['Binding Mode']} (Affinity: {d['Affinity (kcal/mol)']} kcal/mol)" for d in data_var]
                        selected_mode_str_var = st.selectbox("Variant Generated Poses", options_var, key=f"pose_select_var_{idx}")
                        selected_idx_var = options_var.index(selected_mode_str_var)
                        selected_mode_data_var = data_var[selected_idx_var]

                        poses_var = parse_pdbqt.extract_poses("ligand_var_out.pdbqt")
                        if poses_var and selected_idx_var < len(poses_var):
                            selected_pose_str_var = poses_var[selected_idx_var]
                            uff_delta_var = st.session_state.get(f'uff_delta_var_{idx}', 0.0)
                            interactions_data_var = parse_pdbqt.calc_interactions(selected_pose_str_var.split('\n'), receptor_pdbqt)
                            interactions_df_var = pd.DataFrame(interactions_data_var) if interactions_data_var is not None else pd.DataFrame()
                            st.session_state[f'interactions_var_df_{idx}'] = interactions_df_var

                            if not interactions_df_var.empty and "Receptor Residue" in interactions_df_var.columns:
                                interacting_res_var_data = interactions_df_var.to_dict("records")
                                interacting_res_var = list(interactions_df_var["Receptor Residue"].unique())
                            else:
                                interacting_res_var_data = []
                                interacting_res_var = []

                            col1_var, col2_var, col3_var = st.columns(3)
                            col1_var.metric("Pose Affinity", f"{selected_mode_data_var['Affinity (kcal/mol)']} kcal/mol")
                            col2_var.metric("UFF Minimization Delta", f"{uff_delta_var:.2f} kcal/mol")
                            col3_var.metric("Interacting Residues", str(len(interacting_res_var)))

                            st.write(f"Interacting receptor residues: {', '.join(interacting_res_var) if interacting_res_var else 'None'}")

                            if not interactions_df_var.empty:
                                int_col1_var, int_col2_var = st.columns([1, 1])
                                with int_col1_var:
                                    st.dataframe(interactions_df_var, hide_index=True)
                            else:
                                st.info("No significant interactions found for the redesigned variant.")

                            with open(receptor_pdbqt, 'r') as f:
                                receptor_data_var = f.read()

                            viewer_html_var = f'''
                            <div id="container-var-{idx}" style="height: 500px !important; width: 100% !important; display: block; position: relative;" class="viewer_3Dmoljs"
                                 data-backgroundcolor="0xffffff" data-style="stick"></div>
                            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                            <script>
                                var initViewer_var_{idx} = setInterval(function() {{
                                    if (typeof $3Dmol !== 'undefined') {{
                                        clearInterval(initViewer_var_{idx});
                                        var viewer = $3Dmol.createViewer("container-var-{idx}", {{defaultcolors: $3Dmol.rasmolElementColors}});
                                        var receptor_data = `{receptor_data_var}`;
                                        var ligand_data = `{selected_pose_str_var}`;

                                        viewer.addModel(receptor_data, "pdb");
                                        viewer.setStyle({{model: 0}}, {{{receptor_style}: {{color: 'spectrum'}} }});

                                        if ({'true' if show_surface else 'false'}) {{
                                            viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.8, color: 'white'}}, {{model: 0}});
                                        }}

                                        viewer.addModel(ligand_data, "pdb");
                                        viewer.setStyle({{model: 1}}, {{{ligand_style}: {{colorscheme: 'greenCarbon'}} }});


                                        var interactionsData = {json.dumps(interacting_res_var_data)};
                                        var interactingRes = {json.dumps(interacting_res_var)};
                                        interactionsData.forEach(function(interaction) {{
                                            var res = interaction["Receptor Residue"];
                                            var rx = interaction["Receptor XYZ"][0];
                                            var ry = interaction["Receptor XYZ"][1];
                                            var rz = interaction["Receptor XYZ"][2];
                                            var lx = interaction["Ligand XYZ"][0];
                                            var ly = interaction["Ligand XYZ"][1];
                                            var lz = interaction["Ligand XYZ"][2];

                                            // Add dashed line
                                            viewer.addCylinder({{
                                                start: {{x: lx, y: ly, z: lz}},
                                                end: {{x: rx, y: ry, z: rz}},
                                                radius: 0.1,
                                                color: 'yellow',
                                                dashed: true,
                                                fromCap: 1,
                                                toCap: 1
                                            }});

                                            // Anchor label to receptor atom
                                            viewer.addLabel(res, {{
                                                position: {{x: rx, y: ry, z: rz}},
                                                fontColor: 'white',
                                                backgroundColor: 'black',
                                                backgroundOpacity: 0.5,
                                                fontsize: 12
                                            }});
                                        }});


                                        viewer.zoomTo();
                                        viewer.render();
                                    }}
                                }}, 100);
                            </script>
                            '''
                            components.iframe(f'data:text/html;charset=utf-8,{urllib.parse.quote(viewer_html_var)}', height=550)

                # Phase 3 relocated: ADMET & Design
                if st.session_state.get(f'docking_var_done_{idx}', False):
                    st.markdown("---")
                    st.header("Phase 5: ADMET & Design")

                    orig_adme = adme_profiler.get_admet(smiles)
                    adme_data = []
                    if orig_adme:
                        orig_adme["Molecule"] = "Original Phytochemical"
                        adme_data.append(orig_adme)

                    if variants:
                        for i, var_smiles in enumerate(variants):
                            var_adme = adme_profiler.get_admet(var_smiles)
                            if var_adme:
                                var_adme["Molecule"] = f"Redesign Variant {i+1}"
                                adme_data.append(var_adme)

                    if adme_data:
                        df_adme = pd.DataFrame(adme_data)
                        cols = ['Molecule'] + [c for c in df_adme.columns if c != 'Molecule' and c != 'Violation Details']
                        display_df = df_adme[cols]
                        st.dataframe(display_df, hide_index=True)

                        for data in adme_data:
                            violation_details = data.get("Violation Details", [])
                            if violation_details:
                                details_str = ", ".join(violation_details)
                                st.warning(f"⚠️ {data['Molecule']} Violations: {details_str}")

                # --- REPORT DOWNLOAD AND DEVELOPER SIGNATURE ---
                st.markdown("---")

                # Ensure ADMET DataFrame exists before converting
                df_adme_html = df_adme.to_html(index=False) if 'df_adme' in locals() and not df_adme.empty else "<p>No ADMET properties available.</p>"

                # Ensure interactions DataFrame exists before converting
                interactions_df = st.session_state.get(f'interactions_df_{idx}', None)
                interactions_df_html = interactions_df.to_html(index=False) if interactions_df is not None and not interactions_df.empty else "<p>No significant interactions found.</p>"

                interactions_var_df = st.session_state.get(f'interactions_var_df_{idx}', None)
                interactions_var_df_html = interactions_var_df.to_html(index=False) if interactions_var_df is not None and not interactions_var_df.empty else None

                report_html = generate_html_report(
                    plant_name=row['Active Phytochemical'],
                    smiles=st.session_state.get(f'smiles_{idx}', row['SMILES']),
                    receptor_name=row['Protein Target'],
                    pdb_id=row['PDB ID'],
                    df_adme_html=df_adme_html,
                    interactions_df_html=interactions_df_html,
                    interactions_var_df_html=interactions_var_df_html
                )

                col1, col2 = st.columns([1, 1])
                with col1:
                    filename = f"kemet_dock_{row['Common Name'].replace(' ', '_')}_report.html"
                    st.download_button(
                        label='Download Report',
                        data=report_html,
                        file_name=filename,
                        mime='text/html',
                        key=f"download_btn_{idx}"
                    )

                with col2:
                    st.markdown("<div style='text-align: right; color: gray; font-size: small;'>by Simran Ailani</div>", unsafe_allow_html=True)

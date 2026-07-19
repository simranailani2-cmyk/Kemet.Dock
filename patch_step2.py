import re

with open("app.py", "r") as f:
    content = f.read()

# We need to insert the layout at the end of the Phase 3 block, within the loop over selected_data.iterrows().
# The block currently ends with:
#                             if mol2:
#                                 st.image(Draw.MolToImage(mol2, size=(400, 400)), caption='Redesign Variant 1')
# Let's insert our logic right after that.

insertion_code = """
                    # --- REPORT DOWNLOAD AND DEVELOPER SIGNATURE ---
                    st.markdown("---")

                    # Ensure ADMET DataFrame exists before converting
                    df_adme_html = df_adme.to_html(index=False) if 'df_adme' in locals() and not df_adme.empty else "<p>No ADMET properties available.</p>"

                    # Ensure interactions DataFrame exists before converting
                    interactions_df = st.session_state.get(f'interactions_df_{idx}', None)
                    interactions_df_html = interactions_df.to_html(index=False) if interactions_df is not None and not interactions_df.empty else "<p>No significant interactions found.</p>"

                    report_html = generate_html_report(
                        plant_name=row['Active Phytochemical'],
                        smiles=st.session_state.get(f'smiles_{idx}', row['SMILES']),
                        receptor_name=row['Protein Target'],
                        pdb_id=row['PDB ID'],
                        df_adme_html=df_adme_html,
                        interactions_df_html=interactions_df_html
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
"""

# Be very careful matching the exact indentation.
target_lines = """
                            if mol2:
                                st.image(Draw.MolToImage(mol2, size=(400, 400)), caption='Redesign Variant 1')"""

if "REPORT DOWNLOAD AND DEVELOPER SIGNATURE" not in content:
    content = content.replace(target_lines, target_lines + insertion_code)
    with open("app.py", "w") as f:
        f.write(content)
    print("Code inserted successfully.")
else:
    print("Code already inserted.")

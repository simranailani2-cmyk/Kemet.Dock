import re

with open("app.py", "r") as f:
    content = f.read()

func_code = """
def generate_html_report(plant_name, smiles, receptor_name, pdb_id, df_adme_html, interactions_df_html):
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

"""

if "def generate_html_report" not in content:
    # Insert it right before the UI layout begins
    target_str = "st.set_page_config"
    content = content.replace(target_str, func_code + target_str)

    with open("app.py", "w") as f:
        f.write(content)
    print("Function added successfully.")
else:
    print("Function already exists.")

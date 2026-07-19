import re

with open("app.py", "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "st.header(\"ADMET & Design\")" in line:
        print(f"ADMET Header found at line {i+1}")
    if "selected_data.iterrows()" in line:
        print(f"Iterrows found at line {i+1}")

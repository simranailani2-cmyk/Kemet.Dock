import sys

def patch_file():
    with open('app.py', 'r') as f:
        lines = f.readlines()

    start_idx = 126 # Index for line 127
    end_idx = 139 # Index for line 139

    # We replace from start_idx to end_idx (inclusive, which means we slice up to end_idx)
    # Actually, python slice [start_idx:end_idx] means lines[126:139] which is lines 127 to 139

    replacement_code = """
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
"""

    new_lines = lines[:start_idx] + [replacement_code] + lines[end_idx:]

    with open('app.py', 'w') as f:
        f.writelines(new_lines)

if __name__ == '__main__':
    patch_file()

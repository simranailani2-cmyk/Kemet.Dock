import math

def extract_poses(pdbqt_file):
    poses = []
    current_pose = []
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    current_pose = []
                elif line.startswith('ENDMDL'):
                    poses.append("".join(current_pose))
                else:
                    if not line.startswith('ENDROOT') and not line.startswith('TORSDOF') and not line.startswith('ROOT'):
                        current_pose.append(line)
        return poses
    except FileNotFoundError:
        return []

def calc_interactions(ligand_lines, receptor_pdbqt):
    # Load receptor coordinates
    receptor_atoms = []
    with open(receptor_pdbqt, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_name = line[17:20].strip()
                    res_num = line[22:26].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    receptor_atoms.append((f"{res_name}{res_num}", x, y, z))
                except ValueError:
                    pass

    # Extract ligand coordinates
    ligand_atoms = []
    for line in ligand_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ligand_atoms.append((x, y, z))
            except ValueError:
                pass

    # Find interacting residues (< 4 Angstroms) and construct dicts
    interactions = []
    seen = set()
    for i, (lx, ly, lz) in enumerate(ligand_atoms):
        for res, rx, ry, rz in receptor_atoms:
            dist = math.sqrt((lx-rx)**2 + (ly-ry)**2 + (lz-rz)**2)
            if dist <= 4.0:
                key = (res, i)
                if key not in seen:
                    if dist < 3.5:
                        bond_type = "Polar/H-Bond"
                    else:
                        bond_type = "Hydrophobic"

                    interactions.append({
                        "Receptor Residue": res,
                        "Ligand Atom": f"Atom {i+1}",
                        "Distance": round(dist, 2),
                        "Bond Type": bond_type
                    })
                    seen.add(key)

    # Sort interactions by distance
    interactions.sort(key=lambda x: x["Distance"])
    return interactions

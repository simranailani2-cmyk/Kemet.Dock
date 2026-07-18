import os
import urllib.request
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def download_vina_executable():
    url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
    vina_path = "./vina"
    if not os.path.exists(vina_path):
        urllib.request.urlretrieve(url, vina_path)
        os.chmod(vina_path, 0o755)
    return vina_path

def fetch_receptor(pdb_id, output_pdb):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, output_pdb)
    except Exception as e:
        print(f"Error downloading PDB {pdb_id}: {e}")
        return None

    with open(output_pdb, 'r') as f:
        pdb_str = f.read()

    lines = []
    for line in pdb_str.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            element = line[76:78].strip()
            if not element:
                element = line[12:16].strip()[0]
            ad_type = element
            if element == 'H': ad_type = 'HD'
            elif element == 'C': ad_type = 'C'
            elif element == 'N': ad_type = 'N'
            elif element == 'O': ad_type = 'OA'
            elif element == 'S': ad_type = 'SA'
            elif element == 'P': ad_type = 'P'
            elif element == 'CL': ad_type = 'Cl'
            elif element == 'BR': ad_type = 'Br'
            elif element == 'F': ad_type = 'F'
            elif element == 'I': ad_type = 'I'
            elif element == 'FE': ad_type = 'Fe'
            elif element == 'ZN': ad_type = 'Zn'
            elif element == 'CA': ad_type = 'Ca'
            elif element == 'MG': ad_type = 'Mg'
            elif element == 'MN': ad_type = 'Mn'

            new_line = line[:66].ljust(66) + "    +0.000 " + ad_type.ljust(2)
            lines.append(new_line)

    pdbqt_output = output_pdb.replace('.pdb', '.pdbqt')
    with open(pdbqt_output, 'w') as f:
        f.write("\n".join(lines))

    return pdbqt_output

def prepare_ligand(smiles, output_pdbqt):
    uff_delta = 0.0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, 0.0

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    # Memory-safe UFF minimization
    if mol.GetNumAtoms() <= 4000:
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            if ff:
                e1 = ff.CalcEnergy()
                ff.Minimize()
                e2 = ff.CalcEnergy()
                uff_delta = e1 - e2
        except Exception as e:
            print(f"UFF optimization failed: {e}")
            pass

    # MMFF optimization
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        print(f"MMFF optimization failed: {e}")
        pass

    pdb_block = Chem.MolToPDBBlock(mol)

    lines = ["ROOT"]
    for line in pdb_block.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            element = line[76:78].strip()
            if not element:
                element = line[12:16].strip()[0]
            ad_type = element
            if element == 'H': ad_type = 'HD'
            elif element == 'C': ad_type = 'C'
            elif element == 'N': ad_type = 'N'
            elif element == 'O': ad_type = 'OA'
            elif element == 'S': ad_type = 'SA'
            elif element == 'P': ad_type = 'P'
            elif element == 'CL': ad_type = 'Cl'
            elif element == 'BR': ad_type = 'Br'
            elif element == 'F': ad_type = 'F'
            elif element == 'I': ad_type = 'I'
            elif element == 'FE': ad_type = 'Fe'
            elif element == 'ZN': ad_type = 'Zn'
            elif element == 'CA': ad_type = 'Ca'
            elif element == 'MG': ad_type = 'Mg'
            elif element == 'MN': ad_type = 'Mn'

            new_line = line[:66].ljust(66) + "    +0.000 " + ad_type.ljust(2)
            lines.append(new_line)
    lines.append("ENDROOT")
    lines.append("TORSDOF 0")

    with open(output_pdbqt, 'w') as f:
        f.write("\n".join(lines))
    return output_pdbqt, uff_delta

def smart_cavity_finder(pdb_file):
    try:
        mol = Chem.MolFromPDBFile(pdb_file, sanitize=False)
        if mol:
            conf = mol.GetConformer()
            coords = conf.GetPositions()
            if len(coords) > 0:
                center = np.mean(coords, axis=0)
                mins = np.min(coords, axis=0)
                maxs = np.max(coords, axis=0)
                dims = maxs - mins
                dims = np.clip(dims, a_min=10.0, a_max=60.0)
                return [float(c) for c in center], [float(d) for d in dims]
    except Exception as e:
        print(f"Error parsing PDB for cavity: {e}")

    return [0.0, 0.0, 0.0], [20.0, 20.0, 20.0]

def run_vina_docking(receptor_pdbqt, ligand_pdbqt, center, dims):
    vina_path = download_vina_executable()

    cmd = [
        vina_path,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(dims[0]),
        "--size_y", str(dims[1]),
        "--size_z", str(dims[2]),
        "--exhaustiveness", "1"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res.stdout

import re
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def generate_variants(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return []

    variants = set()

    # Cleavage/redesign rules (SMARTS -> replacement SMILES)
    rules = [
        ('[OX2H]', 'F'),
        ('[OX2H]', 'Cl'),
        ('[OX2H]', 'N'),
        ('[NX3H2]', 'O'),
        ('[NX3H2]', 'F'),
        ('C(=O)[OX2H1]', 'C(=O)N'),
        ('C(=O)[OX2H1]', 'C(=O)OC'),

    ]

    for smarts, rep_smiles in rules:
        pattern = Chem.MolFromSmarts(smarts)
        replacement = Chem.MolFromSmiles(rep_smiles)
        if pattern and replacement and mol.HasSubstructMatch(pattern):
            try:
                res = Chem.ReplaceSubstructs(mol, pattern, replacement)
                for r in res:
                    try:
                        Chem.SanitizeMol(r)
                        variants.add(Chem.MolToSmiles(r))
                    except:
                        pass
            except:
                pass

    return list(variants)

def get_admet(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    ri = mol.GetRingInfo()
    max_ring_size = max([len(r) for r in ri.AtomRings()]) if ri.AtomRings() else 0

    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1

    return {
        "Molecular Weight": round(mw, 2),
        "LogP": round(logp, 2),
        "TPSA": round(tpsa, 2),
        "H-Bond Donors": hbd,
        "H-Bond Acceptors": hba,
        "Max Ring Size": max_ring_size,
        "Lipinski Violations": violations
    }

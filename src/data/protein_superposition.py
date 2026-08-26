import pymol
from pymol import cmd

# pymol.finish_launching(['pymol', '-cq'])

# cmd.fetch('6wdz', 'target_protein')
# cmd.fetch('6we0', 'mobile_protein')

# rms_result = cmd.super('mobile_protein', 'target_protein', object="local_align_obj")

# print(f"Alignment RMSD: {rms_result[0]:.3f} Å")
# print(f"Aligned Atoms: {rms_result[1]}")

# cmd.save('mobile_aligned.pdb', 'mobile_protein')


cmd.fetch('6wdz', 'mobile_protein')
cmd.fetch('6we0', 'target_protein')

rms_result = cmd.super('mobile_protein', 'target_protein', object="local_align_obj")

# sDBM for WDV occurs from residue 90-99
query = "target_protein and chain A and resi 90-99 and name CA"

# Build index -> residue metadata lookup
atom_info = {}

for obj in ("target_protein", "mobile_protein"):
    model = cmd.get_model(f"{obj} and name CA")
    for atom in model.atom:
        atom_info[(obj, atom.index)] = {
            "chain": atom.chain,
            "resi": atom.resi,
            "resn": atom.resn,
        }

# Find structural-alignment partners of the query residues
partners = []

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
    "VAL": "V",
}

def one_letter(resn):
    return AA3_TO_1.get(resn.upper(), "X")

for atom_a, atom_b in cmd.get_raw_alignment("local_align_obj"):
    if atom_a[0] == "target_protein":
        target_atom, mobile_atom = atom_a, atom_b
    else:
        target_atom, mobile_atom = atom_b, atom_a

    # Check whether this target atom is in the requested region
    target_index = target_atom[1]
    if cmd.count_atoms(f"({query}) and index {target_index}"):
        t = atom_info[target_atom]
        m = atom_info[mobile_atom]
        partners.append((t, m))

for target, mobile in partners:
    print(
        f"Target: {one_letter(target['resn'])}"
        f"-> Mobile: {one_letter(mobile['resn'])}"
    )
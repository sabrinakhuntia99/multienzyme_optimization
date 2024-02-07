# Extract amino acid sequence
def extract_aa_sequence(structure):
    amino_acid_sequence = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is an amino acid (exclude water molecules, ligands, etc.)
                if residue.get_id()[0] == " " and residue.get_resname() not in ["HOH", "H2O", "WAT"]:
                    amino_acid_sequence += residue.get_resname()

    return amino_acid_sequence

def three_to_one(aa_sequence):
    three_to_one_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    one_letter_sequence = ""
    for i in range(0, len(aa_sequence), 3):
        three_aa = aa_sequence[i:i + 3]
        one_aa = three_to_one_dict.get(three_aa)
        if one_aa:
            one_letter_sequence += one_aa
        else:
            one_letter_sequence += 'X'  # Placeholder for unknown amino acids
    return one_letter_sequence

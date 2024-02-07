from __future__ import division
from Bio.PDB import *
from data_extraction import extract_aa_sequence, three_to_one
from proteo_sim import cleave_with_enzyme

p = PDBParser()

#  AlphaFold structure for test
structure = p.get_structure('O60361',
                            r"C:\Users\Sabrina\PycharmProjects\structural_proteomics\venv\AF-O60361-F1-model_v4.pdb")

# Extract one-letter AA sequence
aa_sequence =(extract_aa_sequence(structure))
sequence = three_to_one(aa_sequence)
print(sequence)

# Define protease cleavage sites
enzymes = {
    "Trypsin": "KR",
    "Lys-C": "K"
}

# Choose the enzyme
enzyme_name = "Trypsin"
enzyme_cleavage_site = enzymes[enzyme_name]

# Perform cleavage
cleaved_peptides = cleave_with_enzyme(sequence, enzyme_cleavage_site)

# Print the cleaved peptides
print("Cleaved peptides using", enzyme_name + ":")
for peptide in cleaved_peptides:
    print(peptide)

from __future__ import division
from Bio.PDB import *
from data_extraction import extract_aa_sequence, three_to_one
from proteo_sim import cleave_with_enzyme, count_detectable_peptides
import matplotlib.pyplot as plt

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
    "trypsin": "KR",
    "lys-c": "K",
    "rasp-n": "D",
    "ide-z": "KR",
    "arg-c": "R",
    "asp-n": "D",
    "chymotrypsin": "FYWL",
    "elastase": "GASV",
    "glu-c": "E",
    "pepsin": "FWY",
    "protev-plus": "ENLYFQG"
}

# Choose the enzyme
selected_enzymes = ["chymotrypsin"]

# Initialize list to store cleaved peptides
cleaved_peptides = []

# Perform cleavage with selected enzymes
for enzyme_name in selected_enzymes:
    enzyme_cleavage_site = enzymes[enzyme_name]
    cleaved_peptides += cleave_with_enzyme(sequence, enzyme_cleavage_site)

# Print cleaved peptides
print("Cleaved peptides using", selected_enzymes, ":")
for peptide in cleaved_peptides:
    print(peptide)

# Number of peptides that satisfy the mass spectrometer’s detectability conditions
# (peptide length of 7-30 amino acids, ideally 15)


# Print the count of detectable peptides
detectable_peptides = count_detectable_peptides(cleaved_peptides)
print("Number of peptides that satisfy the mass spectrometer's detectability conditions:", detectable_peptides)

# Number of proteins that can be detected

# Distribution of peptide lengths

# Collect peptide lengths
peptide_lengths = [len(peptide) for peptide in cleaved_peptides]

# Plot histogram
plt.hist(peptide_lengths, bins=range(min(peptide_lengths), max(peptide_lengths) + 1), edgecolor='black')

# Highlight bar for peptides with 15 amino acids length
plt.axvline(x=15, color='red', linestyle='--', linewidth=2)

plt.xlabel('Peptide Length')
plt.ylabel('Frequency')
plt.title('Distribution of Digested Peptide Lengths')
plt.grid(True)
plt.show()

# Average size of digested peptides
average_size = sum(len(peptide) for peptide in cleaved_peptides) / len(cleaved_peptides)
print("Average size of digested peptides:", average_size)

# Maximum size of undigested peptides
max_size = max(len(peptide) for peptide in cleaved_peptides)
print("Maximum size of undigested peptides:", max_size)

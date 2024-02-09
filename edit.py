from __future__ import division
from Bio.PDB import *
from data_extraction import extract_aa_sequence, three_to_one
from proteo_sim import cleave_with_enzyme, count_detectable_peptides
import matplotlib.pyplot as plt
from itertools import combinations

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

# Initialize variables to track the best combination
best_detectable_peptides = 0
best_average_size = 0
best_max_size = 0
best_enzymes_used = []

# Iterate over all enzymes and combinations of enzymes
for enzyme_count in range(1, len(enzymes) + 1):
    for selected_enzymes in combinations(enzymes.keys(), enzyme_count):
        # Initialize list to store cleaved peptides for this combination of enzymes
        cleaved_peptides = []
        # Perform cleavage with selected enzymes
        for enzyme_name in selected_enzymes:
            enzyme_cleavage_site = enzymes[enzyme_name]
            cleaved_peptides += cleave_with_enzyme(sequence, enzyme_cleavage_site)
        # Count the number of detectable peptides
        detectable_peptides = count_detectable_peptides(cleaved_peptides)
        # Calculate the average size of digested peptides
        average_size = sum(len(peptide) for peptide in cleaved_peptides) / len(cleaved_peptides)
        # Calculate the maximum size of undigested peptides
        max_size = max(len(peptide) for peptide in cleaved_peptides)
        # Update the best combination if this one is better
        if detectable_peptides > best_detectable_peptides \
                or (detectable_peptides == best_detectable_peptides and average_size > best_average_size) \
                or (detectable_peptides == best_detectable_peptides and average_size == best_average_size and max_size > best_max_size):
            best_detectable_peptides = detectable_peptides
            best_average_size = average_size
            best_max_size = max_size
            best_enzymes_used = selected_enzymes

# Print results for the best combination
print("Best combination of enzymes:", best_enzymes_used)
print("Number of peptides that satisfy the mass spectrometer's detectability conditions:", best_detectable_peptides)
print("Average size of digested peptides:", best_average_size)
print("Maximum size of undigested peptides:", best_max_size)

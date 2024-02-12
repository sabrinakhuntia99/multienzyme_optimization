from __future__ import division
from Bio.PDB import *
from data_extraction import extract_aa_sequence, three_to_one
from proteo_sim import cleave_with_enzyme, count_detectable_peptides
import matplotlib.pyplot as plt

# Initialize PDB parser
p = PDBParser()

# List of UniProt IDs
uni_prot_ids = [
    'Q03692', 'Q5QPC7', 'Q5QPC8', 'C9JMN2', 'H7C381', 'P12107', 'A2AAS7', 'P13942', 'E7ES50',
    'E7EX21', 'Q5TAT6', 'J3QT75', 'J3QT83', 'Q05707', 'Q4G0W3', 'P39059', 'A6NCT7', 'H7BY97', 'H7BZL8',
    'H7C3F0', 'Q07092', 'A2A2Y8', 'Q9UMD9', 'H7BXV5', 'P39060', 'Q14993', 'Q5JVU1', 'I3L3H7', 'P02452',
    'P08123', 'B7ZBI4', 'B7ZBI5', 'Q9P218', 'A6PVD9', 'F5GZK2', 'H0YDH6', 'Q96P44', 'H0YAX7', 'Q8NFW1',
    'Q86Y22', 'E9PNK8', 'F8WDM8', 'Q17RW2', 'A8MWQ5', 'D6R8Y2', 'E9PNV9', 'Q9BXS0', 'Q96A83', 'Q5T1U7',
    'Q8IZC6', 'H7BZU0', 'H7C3P2', 'Q2UY09', 'P02458', 'P02461', 'P02462', 'A2A352', 'P08572', 'H7BXM4',
    'Q01955', 'P53420', 'H0Y9R8', 'P29400', 'B4DZ39', 'Q14031', 'P20908', 'P05997', 'P25940', 'P12109',
    'C9JH44', 'H7C0M5', 'P12110', 'C9JNG9', 'A8TX70', 'H0Y393', 'H0Y9T2', 'A6NMZ7', 'C9JBL3', 'C9JTN9',
    'P27658', 'E9PP49', 'P25067', 'P20849', 'B1AKJ3', 'H0Y409', 'Q14055', 'Q14050', 'Q4VXW1'
]

# Get user input for enzyme selection
enzyme_names = input("Enter the enzyme names separated by commas (e.g., trypsin,chymotrypsin): ").lower().split(',')

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

# Initialize list to store cleaved peptides
cleaved_peptides_trypsin = []

# Perform cleavage with trypsin for all UniProt IDs
for uni_prot_id in uni_prot_ids:
    # Construct the full PDB file path
    pdb_file_path = f"C:\\Users\\Sabrina\\PycharmProjects\\multienzyme_optimization\\venv\\AF-{uni_prot_id}-F1-model_v4.pdb"

    # Parse the PDB file
    structure = p.get_structure(uni_prot_id, pdb_file_path)

    # Extract one-letter AA sequence
    aa_sequence = extract_aa_sequence(structure)

    # Perform cleavage with trypsin
    cleaved_peptides_trypsin += cleave_with_enzyme(aa_sequence, enzymes["trypsin"])

# Initialize list to store cleaved peptides
cleaved_peptides_chymotrypsin = []

# Perform cleavage with chymotrypsin for all peptides cleaved with trypsin
for peptide in cleaved_peptides_trypsin:
    cleaved_peptides_chymotrypsin += cleave_with_enzyme(peptide, enzymes["chymotrypsin"])

# Combine all peptides
all_cleaved_peptides = cleaved_peptides_trypsin + cleaved_peptides_chymotrypsin

# Print the count of detectable peptides
detectable_peptides = count_detectable_peptides(all_cleaved_peptides)
print("Number of peptides that satisfy the mass spectrometer's detectability conditions:", detectable_peptides)

# Calculate theoretical average detected sequence coverage
total_collagen_length = sum(len(extract_aa_sequence(p.get_structure(uni_prot_id, f"C:\\Users\\Sabrina\\PycharmProjects\\multienzyme_optimization\\venv\\AF-{uni_prot_id}-F1-model_v4.pdb"))) for uni_prot_id in uni_prot_ids)

# Calculate total length of detectable peptides
total_detectable_length = sum(len(peptide) for peptide in all_cleaved_peptides if 7 <= len(peptide) <= 17)

# Calculate coverage
coverage = total_detectable_length / total_collagen_length * 100  # Convert to percentage

# Print coverage
print("Theoretical average detected sequence coverage by trypsin and chymotrypsin:", coverage)

# Distribution of peptide lengths
peptide_lengths = [len(peptide) for peptide in all_cleaved_peptides]

# Plot histogram
plt.hist(peptide_lengths, bins=range(min(peptide_lengths), max(peptide_lengths) + 1), edgecolor='black')
plt.xlabel('Peptide Length')
plt.ylabel('Frequency')
plt.title('Distribution of Digested Peptide Lengths')
plt.grid(True)
plt.show()

# Average size of digested peptides
average_size = sum(len(peptide) for peptide in all_cleaved_peptides) / len(all_cleaved_peptides)
print("Average size of digested peptides:", average_size)

# Maximum size of undigested peptides
max_size = max(len(peptide) for peptide in all_cleaved_peptides)
print("Maximum size of undigested peptides:", max_size)

# Count of proteins with less than two peptides
proteins_less_than_two_peptides = sum(1 for protein_id in uni_prot_ids if sum(1 for peptide in all_cleaved_peptides if peptide.startswith(protein_id)) < 2)
print("Number of proteins not detected with at least two peptides:", proteins_less_than_two_peptides)

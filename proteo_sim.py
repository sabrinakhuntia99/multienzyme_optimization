def cleave_with_enzyme(sequence, enzymes):
    cleavage_sites = [0]  # Start with the beginning of the sequence

    # Find cleavage sites for each enzyme
    for enzyme in enzymes:
        for i in range(len(sequence)):
            if sequence[i] in enzyme:
                cleavage_sites.append(i + 1)  # Add 1 because we want the position after the cleavage site

    # Sort and remove duplicates from cleavage sites
    cleavage_sites = sorted(set(cleavage_sites))

    # Generate cleaved peptides
    cleaved_peptides = [sequence[cleavage_sites[i]:cleavage_sites[i + 1]] for i in range(len(cleavage_sites) - 1)]

    # Add the last peptide
    cleaved_peptides.append(sequence[cleavage_sites[-1]:])

    return cleaved_peptides

def count_detectable_peptides(cleaved_peptides):
    # Initialize counter for peptides that satisfy the detectability conditions
    detectable_peptides = 0

    # Iterate through cleaved peptides and count those that satisfy the conditions
    for peptide in cleaved_peptides:
        peptide_length = len(peptide)
        if 7 <= peptide_length <= 30:
            detectable_peptides += 1

    return detectable_peptides


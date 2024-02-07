def cleave_with_enzyme(sequence, enzyme):
    cleaved_peptides = []
    cleavage_sites = [0]  # Start with the beginning of the sequence

    # Find cleavage sites
    for i in range(len(sequence)):
        if sequence[i] in enzyme:
            cleavage_sites.append(i + 1)  # Add 1 because we want the position after the cleavage site

    # Generate cleaved peptides
    for i in range(len(cleavage_sites) - 1):
        start = cleavage_sites[i]
        end = cleavage_sites[i + 1]
        cleaved_peptides.append(sequence[start:end])

    # Add the last peptide
    cleaved_peptides.append(sequence[cleavage_sites[-1]:])

    return cleaved_peptides
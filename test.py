import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read the TSV file into a pandas DataFrame
file_path = r'C:\Users\Sabrina\PycharmProjects\multienzyme_optimization\promega_enzymes.tsv'
data = pd.read_csv(file_path, sep='\t')

# Filter the DataFrame for enzymes with cleavage site type 'T'
terminal_enzymes = data[data['CLEAVE_SITE_TYPE'] == 'T']

# Initialize an empty correlation matrix for enzymes
enzymes = terminal_enzymes['ENZYME'].unique()
correlation_matrix = pd.DataFrame(index=enzymes, columns=enzymes)

# Iterate through enzymes and calculate correlations based on amino acids
for enzyme1, amino_acid_list1 in terminal_enzymes[['ENZYME', 'AMINO_ACIDS']].values:
    if pd.notnull(amino_acid_list1):
        acids1 = [acid.strip() for acid in amino_acid_list1.split('-')]
        for enzyme2, amino_acid_list2 in terminal_enzymes[['ENZYME', 'AMINO_ACIDS']].values:
            if pd.notnull(amino_acid_list2):
                acids2 = [acid.strip() for acid in amino_acid_list2.split('-')]
                correlation = len(set(acids1).intersection(acids2)) / len(set(acids1).union(acids2))
                correlation_matrix.loc[enzyme1, enzyme2] = correlation
                correlation_matrix.loc[enzyme2, enzyme1] = correlation

# Plotting the heatmap with a green gradient
plt.figure(figsize=(12, 10))
sns.heatmap(correlation_matrix.astype(float), annot=True, cmap='Greens', fmt='.2f')
plt.title('Correlation between Enzymes Based on Amino Acids for Terminal Cleavage')
plt.xlabel('Enzymes')
plt.ylabel('Enzymes')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()


# Initialize an empty correlation matrix for enzymes based on pH
enzymes = terminal_enzymes['ENZYME'].unique()
ph_matrix = pd.DataFrame(index=enzymes, columns=enzymes)

# Iterate through enzymes and calculate correlations based on pH
for enzyme1, ph1 in terminal_enzymes[['ENZYME', 'PH']].values:
    if pd.notnull(ph1):
        for enzyme2, ph2 in terminal_enzymes[['ENZYME', 'PH']].values:
            if pd.notnull(ph2):
                # Handling pH ranges by taking the midpoint if present
                ph1_values = [float(val) for val in ph1.split('-')]
                ph2_values = [float(val) for val in ph2.split('-')]
                ph1_avg = sum(ph1_values) / len(ph1_values)
                ph2_avg = sum(ph2_values) / len(ph2_values)
                correlation = abs(ph1_avg - ph2_avg)  # Using absolute difference for correlation
                ph_matrix.loc[enzyme1, enzyme2] = correlation
                ph_matrix.loc[enzyme2, enzyme1] = correlation

# Plotting the heatmap with a blue gradient (inverted) based on pH difference
plt.figure(figsize=(12, 10))
sns.heatmap(ph_matrix.astype(float), annot=True, cmap='Blues_r', fmt='.2f')  # '_r' suffix for reversed colormap
plt.title('Difference in pH Values between Enzymes for Terminal Cleavage')
plt.xlabel('Enzymes')
plt.ylabel('Enzymes')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()
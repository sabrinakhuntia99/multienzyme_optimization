import pandas as pd
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv(r"C:\Users\Sabrina\PycharmProjects\multienzyme_optimization\enzymes_protease.csv")

# Extract relevant columns
enzyme_data = data[['ENZYME', 'CLEAVAGE_TERMINAL', 'CLEAVAGE_RESIDUES']]

# Pivot the table to get a matrix of enzymes and their cleavage sites
pivot_table = pd.pivot_table(enzyme_data, index='ENZYME', columns=['CLEAVAGE_TERMINAL', 'CLEAVAGE_RESIDUES'],
                             aggfunc=lambda x: 1 if len(x) > 0 else 0, fill_value=0)

# Calculate cosine similarity between enzymes
similarity_matrix = cosine_similarity(pivot_table)

# Convert similarity matrix to DataFrame
similarity_df = pd.DataFrame(similarity_matrix, columns=pivot_table.index, index=pivot_table.index)

# Plot the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(similarity_df, cmap='viridis', annot=True, fmt=".2f")
plt.title('Enzyme Compatibility Based on Cleavage Site')
plt.xlabel('Enzymes')
plt.ylabel('Enzymes')
plt.tight_layout()
plt.show()

# Extract enzyme families
enzyme_families = data['CLASS']

# Count the frequency of each unique enzyme family
enzyme_family_freq = enzyme_families.value_counts()

# Plotting the frequency bar graph
plt.figure(figsize=(10, 6))
enzyme_family_freq.plot(kind='bar', color='skyblue')
plt.xlabel('Enzyme Family')
plt.ylabel('Number of Proteases')
plt.title('Enzyme Family Frequency')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

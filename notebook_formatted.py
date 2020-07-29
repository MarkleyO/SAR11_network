#%%

# Importing libraries, and reading CSV files into the dataFrame format.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import re
import random
import math

SAR11_network_edges_raw_df = pd.read_csv("fullPubique.7.mcl.cyto.txt default edge.csv")
SAR11_network_nodes_raw_df = pd.read_csv("fullPubique.7.mcl.cyto.txt default node.csv")
SAR11_PSM_raw_df = pd.read_excel("HTCC1062.PSM.counts.xlsx", header=5) # contains coordinates, strand, description
SAR11_modules_raw_df = pd.read_excel("Modules.xlsx")

ecoli_network_nodes_raw_df = pd.read_csv("ecoli.7.mcl.cyto.txt default node.csv")
ecoli_network_edges_raw_df = pd.read_csv("ecoli.7.mcl.cyto.txt default edge.csv")

# Once files are imported, Edges are examined.
# In order to continue, we want to combine the Connected Nodes (in locus_tag format); Start Coordinates;
# Status; Coefficient; Distance into a single DataFrame

#%%
# SAR11 edges in the 'SAR11_network_edges_raw_df' file are listed in 'C134_XXXX_* (interacts with) C134_XXXX_*' format
# under the 'name' column. This is converted into the 'SAR11_XXXX' format and separated into a 'Source' and
# 'Target' column. These two columns, along with the correlation coefficient ('Column 3') and Up/Downregulation status
# ('status') are saved into the DataFrame 'SAR11_edge_list'.

# Pattern with which to strip data from is identified
pattern = '(.*)\s+\(interacts with\)\s+(.*)'
# Gathering sources
SAR11_edge_sources = [re.search(pattern, edge).group(1) for edge in SAR11_network_edges_raw_df['name']]
# Gathering targets
SAR11_edge_targets = [re.search(pattern, edge).group(2) for edge in SAR11_network_edges_raw_df['name']]
# Dict to convert into locus format. Generated from columns in 'SAR11_modules_raw_df' file.
graph_name_to_locus_tag_dict = dict(zip(SAR11_modules_raw_df['name (in graph)'],
                                        SAR11_modules_raw_df['locus_tag']))

# Storing data into the DataFrame
SAR11_edge_list = pd.DataFrame({'Sources': [graph_name_to_locus_tag_dict[source] for source in SAR11_edge_sources],
                               'Targets': [graph_name_to_locus_tag_dict[target] for target in SAR11_edge_targets]})
SAR11_edge_list['Correlation'] = SAR11_network_edges_raw_df['Column 3']
SAR11_edge_list['Status'] = SAR11_network_edges_raw_df['status']

# %%
# Initial ordering of Genes in genome is as observed and provided by data. This is placed into dictionary format.
SAR11_transcription_start_sites_dict = dict(zip(SAR11_PSM_raw_df['locus_tag'], SAR11_PSM_raw_df['start_coord']))
tss = list(SAR11_transcription_start_sites_dict.values())

collected_distances = []

# Edge transcription locations are then retrieved from dictionary and appended to DataFrame
SAR11_edge_list['Source Loc'] = [SAR11_transcription_start_sites_dict.get(source) for source in SAR11_edge_list['Sources']]
SAR11_edge_list['Target Loc'] = [SAR11_transcription_start_sites_dict.get(target) for target in SAR11_edge_list['Targets']]

# Base pair distance between two locations is calculated and appended to the DataFrame. Distances will be the lesser
# of a) the absolute value of the difference between the two, or b) the circumference minus three distances. Total
# genome distance is considered to be 1,308,584. Distances of 0 are a result of edges between SAR11_0932, and
# SAR11_0932-1. Distances are then appended to the DataFrame.

distances = []
for i in range(len(SAR11_edge_list)):
    difference = abs(SAR11_edge_list['Source Loc'][i] - SAR11_edge_list['Target Loc'][i])
    alternate_distance = 1308584 - difference
    distances.append(min(difference, alternate_distance))

SAR11_edge_list['Distance'] = distances
collected_distances.append(list(SAR11_edge_list['Distance']))

# A Plot is made to illustrate this distribution.
plt.figure()
a = sns.jointplot(SAR11_edge_list['Distance'], SAR11_edge_list['Correlation'], kind='scatter', s=1)
a.set_axis_labels("Base Pair Distance on Genome", "Edge Correlation Coefficient")
plt.title("SAR11 Observed Gene Locations")

# This process is then repeated below, only with randomised transcription start site positions, and recalculated
# distances.

# %%
random.shuffle(tss)
randomised_dict_1 = dict(zip(SAR11_transcription_start_sites_dict.keys(), tss))

SAR11_edge_list['Source Loc'] = [randomised_dict_1.get(source) for source in SAR11_edge_list['Sources']]
SAR11_edge_list['Target Loc'] = [randomised_dict_1.get(target) for target in SAR11_edge_list['Targets']]

distances = []
for i in range(len(SAR11_edge_list)):
    difference = abs(SAR11_edge_list['Source Loc'][i] - SAR11_edge_list['Target Loc'][i])
    alternate_distance = 1308584 - difference
    distances.append(min(difference, alternate_distance))

SAR11_edge_list['Distance'] = distances
collected_distances.append(list(SAR11_edge_list['Distance']))

plt.figure()
a = sns.jointplot(SAR11_edge_list['Distance'], SAR11_edge_list['Correlation'], kind='scatter', s=1)
a.set_axis_labels("Base Pair Distance on Genome", "Edge Correlation Coefficient")
plt.title("SAR11 Random Dist 1")

# %%
random.shuffle(tss)
randomised_dict_2 = dict(zip(SAR11_transcription_start_sites_dict.keys(), tss))

SAR11_edge_list['Source Loc'] = [randomised_dict_2.get(source) for source in SAR11_edge_list['Sources']]
SAR11_edge_list['Target Loc'] = [randomised_dict_2.get(target) for target in SAR11_edge_list['Targets']]

distances = []
for i in range(len(SAR11_edge_list)):
    difference = abs(SAR11_edge_list['Source Loc'][i] - SAR11_edge_list['Target Loc'][i])
    alternate_distance = 1308584 - difference
    distances.append(min(difference, alternate_distance))

SAR11_edge_list['Distance'] = distances
collected_distances.append(list(SAR11_edge_list['Distance']))

plt.figure()
a = sns.jointplot(SAR11_edge_list['Distance'], SAR11_edge_list['Correlation'], kind='scatter', s=1)
a.set_axis_labels("Base Pair Distance on Genome", "Edge Correlation Coefficient")
plt.title("SAR11 Random Dist 2")

# %%
random.shuffle(tss)
randomised_dict_3 = dict(zip(SAR11_transcription_start_sites_dict.keys(), tss))

SAR11_edge_list['Source Loc'] = [randomised_dict_3.get(source) for source in SAR11_edge_list['Sources']]
SAR11_edge_list['Target Loc'] = [randomised_dict_3.get(target) for target in SAR11_edge_list['Targets']]

distances = []
for i in range(len(SAR11_edge_list)):
    difference = abs(SAR11_edge_list['Source Loc'][i] - SAR11_edge_list['Target Loc'][i])
    alternate_distance = 1308584 - difference
    distances.append(min(difference, alternate_distance))

SAR11_edge_list['Distance'] = distances
collected_distances.append(list(SAR11_edge_list['Distance']))

plt.figure()
a = sns.jointplot(SAR11_edge_list['Distance'], SAR11_edge_list['Correlation'], kind='scatter', s=1)
a.set_axis_labels("Base Pair Distance on Genome", "Edge Correlation Coefficient")
plt.title("SAR11 Random Dist 3")

# %%
plt.show()

#%%
print("Comparing SAR11 Distance Distribution to Random Distributions:")
print(stats.ks_2samp(collected_distances[0], collected_distances[1]))
print(stats.ks_2samp(collected_distances[0], collected_distances[2]))
print(stats.ks_2samp(collected_distances[0], collected_distances[3]))
print("\nComparing Random Distributions Amongst Each Other:")
print(stats.ks_2samp(collected_distances[1], collected_distances[2]))
print(stats.ks_2samp(collected_distances[1], collected_distances[3]))
print(stats.ks_2samp(collected_distances[2], collected_distances[3]))

for i in range(len(collected_distances)):
    collected_distances[i] = [dist for dist in collected_distances[i] if dist > 0]

print("\nComparing SAR11 Distance Distribution to Random Distributions:")
print(stats.ks_2samp(collected_distances[0], collected_distances[1]))
print(stats.ks_2samp(collected_distances[0], collected_distances[2]))
print(stats.ks_2samp(collected_distances[0], collected_distances[3]))
print("\nComparing Random Distributions Amongst Each Other:")
print(stats.ks_2samp(collected_distances[1], collected_distances[2]))
print(stats.ks_2samp(collected_distances[1], collected_distances[3]))
print(stats.ks_2samp(collected_distances[2], collected_distances[3]))


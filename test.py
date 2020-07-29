import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

w=10
h=10
fig=plt.figure(figsize=(8, 8))
columns = 4
rows = 5

for i in range(1, columns*rows +1):
    img = np.random.randint(10, size=(h,w))
    fig.add_subplot(rows, columns, i)
    plt.imshow(img)
plt.show()


a = np.random.uniform(size=4999)
b = np.random.uniform(size=4999)
c = np.random.normal(size=4999)
d = np.random.uniform(size=100)
print(stats.kstest(a, stats.uniform.cdf))
print(stats.kstest(b, stats.uniform.cdf))
print(stats.kstest(c, stats.uniform.cdf))
print(stats.ks_2samp(a, b))
print(stats.ks_2samp(a, c))
print(stats.ks_2samp(b, d))

# raw_df_SAR11_nodes = pd.read_csv("fullPubique.7.mcl.cyto.txt default node.csv")
#
# SAR11_categories = pd.read_excel("HTCC1062.PSM.counts.xlsx", header=5, usecols=['locus_tag', 'category'])
#
# cluster_dict = SAR11_categories.set_index('locus_tag').to_dict()['category']
# SAR11_formatted_data = np.empty((984, 5))
#
# for row in raw_df_SAR11_nodes.iterrows():
#     # selecting corresponding row from formatted data
#     formatted_data_row = SAR11_formatted_data[row[0]]
#     formatted_data_row[0] = int(row[1][12][6:])
#     formatted_data_row[1] = int(row[1][5])
#     formatted_data_row[2] = row[1][1]
#     formatted_data_row[3] = row[1][2]
#     formatted_data_row[4] = cluster_dict.get(row[1][12])
#
# SAR11_df = pd.DataFrame(SAR11_formatted_data, columns=['SAR11 Gene', 'Degree', 'Betweenness Centrality',
#                                                        'Closeness Centrality', 'Cluster'])
# print(SAR11_df)
#
# grid = sns.JointGrid(x='SAR11 Gene', y='Degree', data=SAR11_df)
# g = grid.plot_joint(sns.scatterplot, hue='Cluster', data=SAR11_df)
# sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] == 1, 'SAR11 Gene'], ax=g.ax_marg_x, legend=False, shade=True)
# sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] > 0, 'SAR11 Gene'], ax=g.ax_marg_x, legend=False, shade=True)
# # sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] == 3, 'SAR11 Gene'], ax=g.ax_marg_x, legend=False, shade=True)
# # sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] == 4, 'SAR11 Gene'], ax=g.ax_marg_x, legend=False, shade=True)
# sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] == 1, 'Degree'], ax=g.ax_marg_y, vertical=True, legend=True, shade=True)
# sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] > 0, 'Degree'], ax=g.ax_marg_y, vertical=True, legend=True, shade=True)
# # sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] == 3, 'Degree'], ax=g.ax_marg_y, vertical=True, legend=True)
# # sns.kdeplot(SAR11_df.loc[SAR11_df['Cluster'] == 4, 'Degree'], ax=g.ax_marg_y, vertical=True, legend=True)
#
# g2 = sns.jointplot("SAR11 Gene", "Degree", data=SAR11_df, kind="reg")

# tips = sns.load_dataset('tips')
# print(tips)
# grid = sns.JointGrid(x='total_bill', y='tip', data=tips)
#
# g = grid.plot_joint(sns.scatterplot, hue='smoker', data=tips)
# sns.kdeplot(tips.loc[tips['smoker']=='Yes', 'total_bill'], ax=g.ax_marg_x, legend=False)
# sns.kdeplot(tips.loc[tips['smoker']=='No', 'total_bill'], ax=g.ax_marg_x, legend=False)
# sns.kdeplot(tips.loc[tips['smoker']=='Yes', 'tip'], ax=g.ax_marg_y, vertical=True, legend=False)
# sns.kdeplot(tips.loc[tips['smoker']=='No', 'tip'], ax=g.ax_marg_y, vertical=True, legend=False)

plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def get_circular_distance(source_index, target_index, circumference):
    gene_locations = sorted([source_index, target_index])  # high 1375, low 8
    return min((gene_locations[1] - gene_locations[0]), ((circumference - gene_locations[1]) + gene_locations[0]))


def read_SAR11_edges(raw_df_SAR11):
    counts7 = 0
    counts2 = 0
    SAR11_formatted_data = np.empty((8662, 6))
    for row in raw_df_SAR11.iterrows():
        # selecting corresponding row from formatted data
        formatted_data_row = SAR11_formatted_data[row[0]]
        # Setting location in genome
        spliced_name_str = row[1][2].split(" ")
        formatted_data_row[0] = int(spliced_name_str[0][5:9])
        formatted_data_row[1] = int(spliced_name_str[3][5:9])
        # Setting the distance between genes
        formatted_data_row[2] = get_circular_distance(formatted_data_row[0], formatted_data_row[1], 1375)
        if get_circular_distance(formatted_data_row[0], formatted_data_row[1], 1375) <= 7:
            counts7 += 1
            if get_circular_distance(formatted_data_row[0], formatted_data_row[1], 1375) <= 2:
                counts2 += 1
        # Setting weight of edge
        formatted_data_row[3] = row[1][0]
        # Setting Edge Betweenness
        formatted_data_row[4] = row[1][1]
        # Setting status of edge (1 is up, 0 is down)
        if row[1][3] == "up":
            formatted_data_row[5] = 1
        else:
            formatted_data_row[5] = 0

        SAR11_df = pd.DataFrame(SAR11_formatted_data, columns=['SAR11 Source Gene', 'SAR11 Target Gene', 'Distance',
                                                               'Edge Weight', 'Edge Betweenness', 'Status'])
    print("Sar counts = ", counts7, counts2)
    return SAR11_df


def read_ecoli_edges(raw_df_ecoli):
    counts7 = 0
    counts2 = 0
    ecoli_formatted_data = np.empty((50846, 6))
    for row in raw_df_ecoli.iterrows():
        # selecting corresponding row from formatted data
        formatted_data_row = ecoli_formatted_data[row[0]]
        # Setting location in genome
        spliced_name_str = row[1][2].split(" ")
        spliced_source_str = spliced_name_str[0].split("_")
        spliced_target_str = spliced_name_str[3].split("_")
        formatted_data_row[0] = int(spliced_source_str[1][1:])
        formatted_data_row[1] = int(spliced_target_str[1][1:])

        # Setting the distance between genes
        formatted_data_row[2] = get_circular_distance(formatted_data_row[0], formatted_data_row[1], 4624)
        if get_circular_distance(formatted_data_row[0], formatted_data_row[1], 4624) <= 23:
            counts7 += 1
            if get_circular_distance(formatted_data_row[0], formatted_data_row[1], 4624) <= 2:
                counts2 += 1
        # Setting weight of edge
        formatted_data_row[3] = np.float64(row[1][0])
        # Setting Edge Betweenness
        formatted_data_row[4] = np.float64(row[1][1])
        # Setting status of edge (1 is up, 0 is down)
        if row[1][3] == "up":
            formatted_data_row[5] = 1
        else:
            formatted_data_row[5] = 0

    ecoli_df = pd.DataFrame(ecoli_formatted_data, columns=['Ecoli Source Gene', 'Ecoli Target Gene', 'Distance',
                                                           'Edge Weight', 'Edge Betweenness', 'Status'])
    print("ecoli counts = ", counts7, counts2)
    return ecoli_df


def read_SAR11_nodes(raw_df_SAR11_nodes, col1):
    cluster_dict = col1.set_index('locus_tag').to_dict()['category']
    SAR11_formatted_data = np.empty((984, 5))

    for row in raw_df_SAR11_nodes.iterrows():
        # selecting corresponding row from formatted data
        formatted_data_row = SAR11_formatted_data[row[0]]
        formatted_data_row[0] = int(row[1][12][6:])
        formatted_data_row[1] = int(row[1][5])
        formatted_data_row[2] = row[1][1]
        formatted_data_row[3] = row[1][2]
        formatted_data_row[4] = cluster_dict.get(row[1][12])

    SAR11_df = pd.DataFrame(SAR11_formatted_data, columns=['SAR11 Gene', 'Degree', 'Betweenness Centrality',
                                                           'Closeness Centrality', 'Cluster'])

    return SAR11_df


def read_ecoli_nodes(raw_df_ecoli_nodes):
    ecoli_formatted_data = np.empty((3401, 4))
    for row in raw_df_ecoli_nodes.iterrows():
        # selecting corresponding row from formatted data
        formatted_data_row = ecoli_formatted_data[row[0]]
        spliced_str = row[1][8].split("_")
        formatted_data_row[0] = int(spliced_str[1][1:])
        formatted_data_row[1] = int(row[1][4])
        formatted_data_row[2] = row[1][1]
        formatted_data_row[3] = row[1][2]

    ecoli_df = pd.DataFrame(ecoli_formatted_data, columns=['E.Coli Gene', 'Degree', 'Betweenness Centrality',
                                                           'Closeness Centrality'])

    return ecoli_df


def main():
    if(1):
        raw_df_SAR11 = pd.read_csv("fullPubique.7.mcl.cyto.txt default edge.csv", header=0,
                                   usecols=["Column 3", "name", "status", "EdgeBetweenness"])

        raw_df_ecoli = pd.read_csv("ecoli.7.mcl.cyto.txt default edge.csv", header=0,
                                   usecols=["Column 3", "name", "status", "EdgeBetweenness"])

        # raw_df_SAR11_nodes = pd.read_csv("fullPubique.7.mcl.cyto.txt default node.csv")
        #
        # raw_df_ecoli_nodes = pd.read_csv("ecoli.7.mcl.cyto.txt default node.csv")
        #
        # SAR11_categories = pd.read_excel("HTCC1062.PSM.counts.xlsx", header=5, usecols=['locus_tag', 'category'])
        # SAR11_clusters = pd.read_excel("Modules.xlsx")

        # SAR11_nodes = read_SAR11_nodes(raw_df_SAR11_nodes, SAR11_categories)
        # print("\n\n---------SAR11 Nodes---------")
        # print(SAR11_nodes)
        # print("MINIMUMS:\n", SAR11_nodes.min())
        # print("MAXIMUMS:\n", SAR11_nodes.max())

        # ecoli_nodes = read_ecoli_nodes(raw_df_ecoli_nodes)
        # print("\n\n---------E.Coli Nodes---------")
        # print(ecoli_nodes)
        # print("MINIMUMS:\n", ecoli_nodes.min())
        # print("MAXIMUMS:\n", ecoli_nodes.max())

        SAR11_edges = read_SAR11_edges(raw_df_SAR11)
        # print("\n\n---------SAR11 Edges---------")
        # print(SAR11_edges)
        # print("MINIMUMS:\n", SAR11_edges.min())
        # print("MAXIMUMS:\n", SAR11_edges.max())

        ecoli_edges = read_ecoli_edges(raw_df_ecoli)
        # print("\n\n---------E.Coli Edges---------")
        # print(ecoli_edges)
        # print("MINIMUMS:\n", ecoli_edges.min())
        # print("MAXIMUMS:\n", ecoli_edges.max())

    # plt.figure()
    # s_dist_edge_weight = SAR11_edges.plot.scatter(x='Distance',
    #                                               y='Edge Weight',
    #                                               c='Status',
    #                                               s=1,
    #                                               colormap='cool')
    #
    # plt.figure()
    # e_dist_edge_weight = ecoli_edges.plot.scatter(x='Distance',
    #                                               y='Edge Weight',
    #                                               c='Status',
    #                                               s=1,
    #                                               colormap='cool')

    # KDE Plots
    plt.figure()
    a = sns.jointplot(SAR11_edges['Distance'], SAR11_edges['Edge Weight'], shade=True, kind='kde')
    a.set_axis_labels("Gene ID # Distance on Genome", "Edge Correlation Coefficient")
    plt.title("SAR11")
    # plt.ylabel("Edge Correlation Coefficient")
    # plt.xlabel("Gene ID # Distance on Genome")
    plt.figure()
    upregulated = SAR11_edges[SAR11_edges['Status'] == 1]
    downregulated = SAR11_edges[SAR11_edges['Status'] == 0]
    sar11_kdeplot = sns.kdeplot(data=upregulated['Distance'], data2=upregulated['Edge Weight'], cmap='Reds',
                                shade=True)
    plt.figure()
    sar11_kdeplot = sns.kdeplot(data=downregulated['Distance'], data2=downregulated['Edge Weight'], cmap='Blues',
                                shade=True)

    plt.figure()
    b = sns.jointplot(ecoli_edges['Distance'], ecoli_edges['Edge Weight'], shade=True, kind='kde')
    b.set_axis_labels("Gene ID # Distance on Genome", "Edge Correlation Coefficient")
    plt.title("E.Coli")
    # plt.ylabel("Edge Correlation Coefficient")
    # plt.xlabel("Gene ID # Distance on Genome")
    plt.figure()
    upregulated = ecoli_edges[ecoli_edges['Status'] == 1]
    downregulated = ecoli_edges[ecoli_edges['Status'] == 0]
    ecoli_kdeplot = sns.kdeplot(data=upregulated['Distance'], data2=upregulated['Edge Weight'], cmap='Reds',
                                shade=True)
    plt.figure()
    ecoli_kdeplot = sns.kdeplot(data=downregulated['Distance'], data2=downregulated['Edge Weight'], cmap='Blues',
                                shade=True)

    # APPLYING Freedman-Diaconis rule, Bin Width = 2(IQR(x) / Cube Root of n)
    sar11_q75, sar11_q25 = np.percentile(SAR11_edges['Distance'], [25, 75])
    sar11_iqr = sar11_q75 - sar11_q25
    ecoli_q75, ecoli_q25 = np.percentile(ecoli_edges['Distance'], [25, 75])
    ecoli_iqr = ecoli_q75 - ecoli_q25
    sar11_bin_width = int(abs(np.round(2 * (sar11_iqr / (np.cbrt(len(SAR11_edges['Distance'])))))))  # 36
    ecoli_bin_width = int(abs(np.round(2 * (ecoli_iqr / (np.cbrt(len(ecoli_edges['Distance'])))))))  # 42

    # Creating Histograms
    plt.figure()
    s_edge_dist_hist = SAR11_edges.hist(column='Distance',
                                        bins=int(round(SAR11_edges['Distance'].max() / sar11_bin_width)))
    concatenated_data_set = [x for x in SAR11_edges['Distance'] if x > sar11_bin_width]
    sar11_edge_subset = concatenated_data_set
    # plt.figure()
    # plt.hist(concatenated_data_set, bins=18)
    plt.title("SAR11 Distances Histogram")
    plt.ylabel("Edge Counts")
    plt.xlabel("Distance in Bins of Width 36")
    print("Kolmogorov - Smirnov Test SAR11 Subset vs. Uniform:")
    print(stats.kstest(sar11_edge_subset, stats.uniform(loc=sar11_bin_width, scale=SAR11_edges['Distance'].max()).cdf))
    print()

    plt.figure()
    e_edge_dist_hist = ecoli_edges.hist(column='Distance',
                                        bins=int(round(ecoli_edges['Distance'].max() / ecoli_bin_width)))
    concatenated_data_set = [x for x in ecoli_edges['Distance'] if (SAR11_edges['Distance'].max() > x > sar11_bin_width)]
    ecoli_edge_subset = concatenated_data_set
    concatenated_data_set = [x for x in ecoli_edges['Distance'] if (ecoli_edges['Distance'].max() > x > sar11_bin_width)]
    ecoli_edge_subset_2 = concatenated_data_set
    # plt.figure()
    # plt.hist(ecoli_edge_subset, bins=18)
    plt.title("E.Coli Distances Histogram")
    plt.ylabel("Edge Counts")
    plt.xlabel("Distance in Bins of Width 42")

    rand_a = np.random.uniform(size=4999)
    rand_b = np.random.uniform(size=4999)

    print("Kolmogorov - Smirnov Test E.Coli Subset 1 vs. Uniform:")
    print(stats.kstest(ecoli_edge_subset, stats.uniform(loc=sar11_bin_width, scale=SAR11_edges['Distance'].max()).cdf))
    print()
    print("Kolmogorov - Smirnov Test E.Coli Subset 2 vs. Uniform:")
    print(stats.kstest(ecoli_edge_subset_2, stats.uniform(loc=sar11_bin_width, scale=ecoli_edges['Distance'].max()).cdf))
    print()
    print("Kolmogorov - Smirnov Test SAR11 vs E.Coli: \n", stats.ks_2samp(sar11_edge_subset, ecoli_edge_subset))
    print()
    print("Kolmogorov - Smirnov Test Random Uniform 1 vs Random Uniform 2:\n", stats.ks_2samp(rand_a, rand_b))
    print()
    print("Kolmogorov - Smirnov Test SAR11 vs E.Coli: \n", stats.ks_2samp(SAR11_edges['Edge Weight'],
                                                                          ecoli_edges['Edge Weight']))

    if(1):
        print()
        # s_dist_edge_weight = SAR11_edges.plot.scatter(x='Distance',
        #                                               y='Edge Weight',
        #                                               c='Status',
        #                                               s=1,
        #                                               colormap='cool')
        #
        # plt.figure()
        # e_dist_edge_weight = ecoli_edges.plot.scatter(x='Distance',
        #                                               y='Edge Weight',
        #                                               c='Status',
        #                                               s=1,
        #                                               colormap='cool')
        #
        # e_dist_edge_weight = ecoli_edges.plot.scatter(x='Distance',
        #                                               y='Edge Weight',
        #                                               c='Edge Betweenness',
        #                                               s=1,
        #                                               colormap='plasma')
        #
        #
        # for row in SAR11_edges.iterrows():
        #     if row[1][5] == 0:
        #         row[1][3] *= -1
        #
        # figure_6 = SAR11_edges.plot.scatter(x='Distance',
        #                                     y='Edge Weight',
        #                                     s=1)
        #
        # for row in ecoli_edges.iterrows():
        #     if row[1][5] == 0:
        #         row[1][3] *= -1
        #
        # figure_7 = ecoli_edges.plot.scatter(x='Distance',
        #                                     y='Edge Weight',
        #                                     s=1)
        #
        #
        # SAR11_nodes = SAR11_nodes.sort_values(by='SAR11 Gene')
        # sar_grid = sns.JointGrid(x='SAR11 Gene', y='Degree', data=SAR11_nodes)
        # sar_gene_deg = sar_grid.plot_joint(sns.lineplot, color='darkslateblue')
        # sns.regplot(SAR11_nodes.loc[SAR11_nodes['Cluster'] == 1, 'SAR11 Gene'],
        #             SAR11_nodes.loc[SAR11_nodes['Cluster'] == 1, 'Degree'], data=SAR11_nodes, scatter_kws={'s': 1},
        #             ax=sar_gene_deg.ax_marg_x, color='mediumorchid') # the core
        # sns.regplot(SAR11_nodes.loc[SAR11_nodes['Cluster'] > 1, 'SAR11 Gene'],
        #             SAR11_nodes.loc[SAR11_nodes['Cluster'] > 1, 'Degree'], data=SAR11_nodes, scatter_kws={'s': 1},
        #             ax=sar_gene_deg.ax_marg_x, color='moccasin')  # the core
        # sns.kdeplot(SAR11_nodes.loc[SAR11_nodes['Cluster'] == 1, 'Degree'], ax=sar_gene_deg.ax_marg_y,
        #             vertical=True, legend=False, shade=True, color='mediumorchid') # the core
        # sns.kdeplot(SAR11_nodes.loc[SAR11_nodes['Cluster'] > 0, 'Degree'], ax=sar_gene_deg.ax_marg_y,
        #             vertical=True, legend=False, shade=True, color='moccasin')
        # plt.fill_between(SAR11_nodes['SAR11 Gene'].values, SAR11_nodes['Degree'].values, facecolor='darkslateblue')
        #
        # ecoli_nodes = ecoli_nodes.sort_values(by='E.Coli Gene')
        # ecoli_grid = sns.JointGrid(x='E.Coli Gene', y='Degree', data=ecoli_nodes)
        # ecoli_gene_deg = ecoli_grid.plot_joint(sns.lineplot, color='darkslateblue')
        # sns.regplot(ecoli_nodes.loc[ecoli_nodes['E.Coli Gene'] > -1, 'E.Coli Gene'],
        #             ecoli_nodes.loc[ecoli_nodes['E.Coli Gene'] > -1, 'Degree'],
        #             ax=ecoli_gene_deg.ax_marg_x, scatter_kws={'s': 1}, color='mediumorchid')
        # sns.kdeplot(ecoli_nodes.loc[ecoli_nodes['Degree'] > 0, 'Degree'],
        #             ax=ecoli_gene_deg.ax_marg_y, legend=False, shade=True, color='mediumorchid', vertical=True)
        # sns.distplot(ecoli_nodes.loc[ecoli_nodes['Degree'] > -1, 'Degree'],
        #              ax=ecoli_gene_deg.ax_marg_y, vertical=True, color='mediumorchid', bins=100)
        # plt.fill_between(ecoli_nodes['E.Coli Gene'].values, ecoli_nodes['Degree'].values, facecolor='darkslateblue')
        #
        # sns.kdeplot(data=SAR11_edges['Edge Weight'], data2=SAR11_edges['Distance'])

    plt.show()


if __name__ == "__main__":
    main()

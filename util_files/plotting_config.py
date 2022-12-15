

# these are the features we will use for plotting and ML
features = ['count', 'variant_of_concern', 'variant', 'count_proportion_to_max',
            'count_proportion_to_sum', 'degree_centrality', 'closeness_centrality',
            'betweenness_centrality', 'current_flow_betweenness_centrality', 'nfd',
            'fdc', 'formanCurvature', 'cluster_coeff']

save_dir = "plots/"
boxplot_save_dir = save_dir + "boxplots/"
violinplot_save_dir = save_dir + "violinplots/"
pairplot_save_dir = save_dir + "pairplots/"
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from util_files import plotting_config


def graph_nodes_to_df(graph):
    # create a dataframe of the nodes in the graph with their attributes
    df = pd.DataFrame.from_dict(dict(graph.nodes(data=True)), orient='index')
    return df



def boxplot(df, x_feature=None, y_feature=None, stratifier=None):
    # create a boxplot of the y_feature stratified by the x_feature
    sns.boxplot(x=x_feature, y=y_feature, data=df, hue=stratifier)
    plt.savefig(plotting_config.boxplot_save_dir + x_feature + "_" + y_feature + ".png")
    plt.clf()

def violinplot(df, x_feature=None, y_feature=None, stratifier=None):
    # create a violinplot of the y_feature stratified by the x_feature
    sns.violinplot(x=x_feature, y=y_feature, data=df, hue=stratifier)
    plt.savefig(plotting_config.violinplot_save_dir + x_feature + "_" + y_feature + ".png")
    plt.clf()
def main():
    G = nx.read_gml("gml_files/graph_with_voc_normalized_attrs.gml")
    # print a list of feature names for a node
    print(list(G.nodes(data=True)["B.1.1.130"].keys()))
    df = graph_nodes_to_df(G)
    for feature in plotting_config.features:
        boxplot(df, x_feature="variant_of_concern", y_feature=feature)
        violinplot(df, x_feature="variant_of_concern", y_feature=feature)
    # create a pairplot of the features
    sns.pairplot(df, hue="variant_of_concern")
    plt.savefig(plotting_config.pairplot_save_dir + "pairplot.png")
    print(df)




if __name__ == "__main__":
    main()
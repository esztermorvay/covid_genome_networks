import os

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import network_stats
from util import get_all_files_in_dir_as_list, gml_dir
from util_files import plotting_config


def graph_nodes_to_df(graph):
    # create a dataframe of the nodes in the graph with their attributes
    df = pd.DataFrame.from_dict(dict(graph.nodes(data=True)), orient='index')
    return df



def boxplot(df,x_feature=None, y_feature=None,  save_name="", stratifier=None):
    # create a boxplot of the y_feature stratified by the x_feature
    sns.boxplot(x=x_feature, y=y_feature, data=df, hue=stratifier)
    plt.tight_layout()
    plt.savefig(plotting_config.boxplot_save_dir + save_name + "/"  +  x_feature + "_" + y_feature + ".png")
    plt.clf()

def violinplot(df, x_feature=None, y_feature=None,  save_name="", stratifier=None):
    # create a violinplot of the y_feature stratified by the x_feature
    sns.violinplot(x=x_feature, y=y_feature, data=df, hue=stratifier)
    plt.tight_layout()
    plt.savefig(plotting_config.violinplot_save_dir + save_name + "/" + x_feature + "_" + y_feature + ".png")
    plt.clf()


def main():
    gml_files = get_all_files_in_dir_as_list(gml_dir)
    for gml_file in gml_files:
        ext = gml_file.split(".")[-1]
        if gml_file == "community" or ext != "gml" or "graph_with_voc_normalized" not in gml_file:
            continue

        # make directories within the subdirectories if they dont exist yet
        # check if directory exists
        print(gml_file)
        # if not os.path.exists(plotting_config.boxplot_save_dir + gml_file):
        #     os.mkdir(plotting_config.boxplot_save_dir+gml_file)
        #     os.mkdir(plotting_config.violinplot_save_dir+gml_file)
        #     os.mkdir(plotting_config.pairplot_save_dir+gml_file)
        # load the graph
        G = nx.read_gml(gml_dir + "/" + gml_file)
        df = graph_nodes_to_df(G)
        print(network_stats.get_t_scores(df))


        # # print a list of feature names for a node
        # print(list(G.nodes(data=True)["B.1.1.130"].keys()))
        # df = graph_nodes_to_df(G)
        # for feature in plotting_config.features:
        #     try:
        #         boxplot(df, x_feature="variant_of_concern", y_feature=feature, save_name=gml_file)
        #         violinplot(df, x_feature="variant_of_concern", y_feature=feature, save_name=gml_file)
        #     except:
        #         print("Could not plot feature: " + feature)
        # # create a pairplot of the features
        # sns.pairplot(df, hue="variant_of_concern")
        # plt.savefig(plotting_config.pairplot_save_dir + gml_file + "/pairplot.png")
        # plt.clf()
        # print(df)




if __name__ == "__main__":
    main()
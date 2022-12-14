import sys

import networkx as nx
import util
import pandas as pd

from graph_visualizer import draw_graph
from set_attr_util import *
# TODO normalize
# TODO invert
def get_vars_of_concern(csv_file):
    """
    return a dictionary of the variants of concern by pango lineage to name from the csv file
    :param csv_file:
    :return:
    """
   # read in the csv
    df = pd.read_csv(csv_file)
    # drop all columns from the df except Pango Lineage and WHO Label
    print(df.columns)

    # get the pango lineage and the name of the variant of concern
    # remove the Nextstrain Clade and Other column from the df
    df = df.drop(columns=['Nextstrain Clade', 'Other'])
    df = df.dropna()
    print(df)

    #remove the \xa0 character from the values in the Pango Lineage column
    df["Pango Lineage"] = df["Pango Lineage"].str.replace("\xa0", "")
   # create a dictionary of the "Pango Lineage" to the "WHO Label" column
    pango_to_name = dict(zip(df["Pango Lineage"], df["WHO Label"]))
    return pango_to_name

def normalize_edges(graph):
    """
    normalize the edge weights in the graph using linear normalization
    :param graph:
    :return:
    """
    # get the max weight
    max_weight = 0
    # set min weigh tto the max value of an int
    min_weight = sys.maxsize
    for edge in graph.edges:
        if graph.edges[edge]['weight'] > max_weight:
            max_weight = graph.edges[edge]['weight']
        if graph.edges[edge]['weight'] < min_weight:
            min_weight = graph.edges[edge]['weight']
    # normalize the weights
    for edge in graph.edges:
        graph.edges[edge]['weight'] = (graph.edges[edge]['weight'] - min_weight  +1 ) / (max_weight - min_weight + 1)



def normalize_counts(graph):
    """
    normalize the counts of the nodes in the graph
    :param graph:
    :return:
    """
    sum = 0
    # get the max count
    max_count = 0
    for node in graph.nodes:
        if graph.nodes[node]['count'] > max_count:
            max_count = graph.nodes[node]['count']
        sum += graph.nodes[node]['count']
    # normalize the counts
    for node in graph.nodes:
        graph.nodes[node]['count_proportion_to_max'] = graph.nodes[node]['count'] / max_count
        graph.nodes[node]['count_proportion_to_sum'] = graph.nodes[node]['count'] / sum

def add_variant_of_concern(graph, pango_to_name):
    """
    add the variant of concern attribute to the graph
    :param graph:
    :param pango_to_name:
    :return:
    """
    # iterate through the nodes in the graph
    for node in graph.nodes:
        graph.nodes[node]['variant_of_concern'] = False
        graph.nodes[node]['variant'] = "None"
        for lineage in pango_to_name:
            if lineage == node or lineage+"." in node:
                graph.nodes[node]['variant_of_concern'] = True
                graph.nodes[node]['variant'] = pango_to_name[lineage]
                break


def main():
    graph_name = "graph_raw.gml"
    graph = nx.read_gml(graph_name)
    pango_to_name = get_vars_of_concern(util.util_dir + "/variants_of_concern.csv")
    add_variant_of_concern(graph, pango_to_name)
    nx.write_gml(graph, "gml_files/graph_with_voc.gml")
    normalize_edges(graph)
    normalize_counts(graph)
    set_distances(graph)
    set_degree_centrality(graph)
    set_closeness_centrality(graph)
    nx.write_gml(graph, "gml_files/graph_with_voc_normalized.gml")
    # thresholding(graph, 0.9)
    set_betweenness_centralities(graph)
    set_current_flow_betweenness_centrality(graph)
    node_dimension(graph, weight="weight")
    FDC(graph, weight="weight")
    # count the number of edges with weight 0
    count = 0
    for edge in graph.edges:
        if graph.edges[edge]['weight'] == 0:
            print(edge, graph.edges[edge])
            count += 1
    print(count)
    graph = set_forman(graph)
    set_cluster_coeffs(graph)
    #TODO community detection
    graph = set_orc(graph)
    # FIXME ORC must be done on mac bc fork() method doesn't work on windows
    nx.write_gml(graph, "gml_files/graph_with_voc_normalized_attrs.gml")
    print("edges:", len(graph.edges))
    print('done')
    # thresholding(graph, 0.9)
    graph2 = set_inverse_weights(graph)
    threshold_values = [.01, .03, .05, .1, 1]
    for threshold in threshold_values:
        thresholding_low(graph2, threshold)
        # do community detection here
        graph2 = set_orc(graph2, "orc_inversed")

        nx.write_gml(graph2, "gml_files/graph_inversed_thresholded" + str(threshold) + ".gml")

        draw_graph(graph2, "graph_inversed_thresholded" + str(threshold), "variant_of_concern")

if __name__ == "__main__":
    main()
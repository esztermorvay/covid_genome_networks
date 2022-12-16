import community as community_louvain
from GraphRicciCurvature.OllivierRicci import OllivierRicci

import networkx as nx

from graph_visualizer import draw_graph
from util import get_all_files_in_dir_as_list, gml_dir

import pandas as pd


def set_louvain_communities(G):
    partition = community_louvain.best_partition(G)
    nx.set_node_attributes(G, partition, 'louvain_community')
    return G

def set_ricci_flow_community(G):
    orc = OllivierRicci(G, alpha=0.5, verbose="TRACE")
    orc.compute_ricci_flow(iterations=10)
    cc = orc.ricci_community()
    print(cc)

def main():
    # iterate through all the gml files in the gml_files directory
    gml_files = get_all_files_in_dir_as_list(gml_dir)
    for gml_file in gml_files:
        # load the graph
        G = nx.read_gml(gml_dir + "/" + gml_file)
        # set the louvain communities
        G = set_louvain_communities(G)
        # set the ricci flow communities
        # set_ricci_flow_community(G)
        draw_graph(G, "communities" + gml_file, "louvain_community")
        nx.write_gml(G, gml_dir + "/" + gml_file)
        variants = []
        communities = []
        # iterate through nodes in graph
        for node in G.nodes:
            variants.append(G.nodes[node]["variant"])
            communities.append(G.nodes[node]["louvain_community"])
        # create a dataframe with the variant and community
        df = pd.DataFrame({"variant": variants, "community": communities})
        # make sure all columns are printed
        pd.set_option('display.max_columns', None)
        # print the crosstab
        print(pd.crosstab(df["variant"], df["community"]))


if __name__ == "__main__":
    main()
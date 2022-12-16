import community as community_louvain
from GraphRicciCurvature.OllivierRicci import OllivierRicci

import networkx as nx

from util import get_all_files_in_dir_as_list, gml_dir


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
        set_ricci_flow_community(G)
        nx.write_gml(G, gml_dir + "/" + gml_file)

if __name__ == "__main__":
    main()
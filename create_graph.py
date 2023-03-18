import json
import os
import networkx as nx

scores_dirs= ["scores", "garbage_collection/scores", "garbage_collection2/scores", "garbage_collection3/scores"]
counts_dirs = ["counts", "garbage_collection/counts", "garbage_collection2/counts", "garbage_collection3/counts"]

def create_network_counts():
    # create a dict of all the counts
    counts = {}
    for counts_dir in counts_dirs:
        for counts_file in os.listdir(counts_dir):
            if counts_file.endswith(".json"):
                with open(counts_dir + "/" + counts_file, "r") as f:
                    try:
                        data = json.load(f)
                        for key in data:
                            counts[key] = data[key]
                    except:
                        continue
    # save the counts
    with open("counts.json", "w") as f:
        json.dump(counts, f, indent=4)

    # create a networkx network with the counts as node attributes
    G = nx.Graph()
    for key in counts:
        G.add_node(key, count=counts[key])
    return G

def add_edges(G):
    # add edges to the network
    for scores_dir in scores_dirs:
        for scores_file in os.listdir(scores_dir):
            if scores_file.endswith(".json"):
                with open(scores_dir + "/" + scores_file, "r") as f:
                    try:
                        data = json.load(f)
                        for key in data:
                            node1 = key.split("_")[0]
                            node2 = key.split("_")[1]
                            G.add_edge(node1, node2, weight=data[key])
                    except:
                        continue

    return G
def main():
    G = create_network_counts()
    G = add_edges(G)
    print(G.nodes.data())
    # print edges
    print(G.edges.data())
    nx.write_gml(G, "full_graph/gml_files/raw_graph.gml", stringizer=str)





# if name == main
if __name__ == "__main__":
    main()
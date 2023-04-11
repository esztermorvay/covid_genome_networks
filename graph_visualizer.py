import networkx as nx
from sklearn import preprocessing
import matplotlib.pyplot as plt
def draw_graph(G, save_name, clustering_label="variant"):
    """
    A helper function to draw a nx graph with community.
    """
    plt.clf()
    complex_list = nx.get_node_attributes(G, clustering_label)

    le = preprocessing.LabelEncoder()
    node_color = le.fit_transform(list(complex_list.values()))
    # draw with labels
    # label_dict
    nx.draw_spring(G,nodelist=G.nodes(),
                   node_color=node_color,
                   cmap=plt.cm.rainbow,
                   alpha=0.8)
    plt.legend()
    plt.savefig("full_graph/graphs/" + save_name + ".png")
    # plt.show()


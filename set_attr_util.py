import pandas as pd
import networkx as nx
import numpy as np
# import seaborn as sns
import math
import copy
# import community as community_louvain
from GraphRicciCurvature.OllivierRicci import OllivierRicci
from GraphRicciCurvature.FormanRicci import FormanRicci
from collections import Counter
import scipy.stats as stats
import statistics
import os
#TODO orc
#TODO forman
#TODO nfd
#TODO fdc
#TODO community detection





def set_orc(graph, name="orc"):
    """
    set the orc of each node in the graph
    :param graph:
    :return:

    """
    # get the orc of each node
    orc = OllivierRicci(graph, alpha=0.5, verbose="TRACE", cache_maxsize=2000000)
    orc.compute_ricci_curvature()
    # orc_dict = orc.get_ricci_curvature()
    # nx.set_node_attributes(graph, orc_dict, name)
    return orc.G.copy()

def set_forman(graph):
    """
    set the forman of each node in the graph
    :param graph:
    :return:
    """
    # get the forman of each node
    forman = FormanRicci(graph)
    forman.compute_ricci_curvature()
    # forman_dict = forman.get_ricci_curvature()
    G = forman.G
    return G

def set_cluster_coeffs(graph):
    cluster_coeffs = nx.clustering(graph, weight="weight")
    nx.set_node_attributes(graph, cluster_coeffs, 'cluster_coeff')
    return graph

def set_betweenness_centralities(graph):
    """
    set the betweenness centrality of each node in the graph
    :param graph:
    :return:
    """
    # get the betweenness centrality of each node
    betweenness_centrality = nx.betweenness_centrality(graph, weight='weight')
    nx.set_node_attributes(graph, betweenness_centrality, 'betweenness_centrality')
    return graph

def set_current_flow_betweenness_centrality(graph):
    """
    set the current flow betweenness centrality of each node in the graph
    :param graph:
    :return:
    """
    # get the current flow betweenness centrality of each node
    current_flow_betweenness_centrality = nx.current_flow_betweenness_centrality(graph, weight='weight')
    nx.set_node_attributes(graph, current_flow_betweenness_centrality, 'current_flow_betweenness_centrality')
    return graph

def set_distances(g):
    g_distance_dict = {(e1, e2): 1 - weight for e1, e2, weight in g.edges(data='weight')}
    nx.set_edge_attributes(g, g_distance_dict, 'distance')
    return g

def set_inverse_weights(g):
    g2 = g.copy()
    for edge in g.edges:
        g2.edges[edge]['weight'] = 1 - g2.edges[edge]['weight']
    return g2

def set_degree_centrality(g, weighted=False):
    if weighted:
        degree_sums = dict(g.degree(weight="weight"))
        # get the max degree
        max_degree = max(degree_sums.values())
        # get the degree centrality of each node
        degree_centrality = {node: degree / max_degree for node, degree in degree_sums.items()}
        # set the attribute
        nx.set_node_attributes(g, degree_centrality, 'degree_centrality')
    else:
        degree_centrality = nx.degree_centrality(g)
        nx.set_node_attributes(g, degree_centrality, 'degree_centrality')
    return g

def set_closeness_centrality(g):
    closeness_centrality = nx.closeness_centrality(g, distance='distance')
    nx.set_node_attributes(g, closeness_centrality, 'closeness_centrality')
    return g

def node_dimension(G, weight=None):
    node_dimension = {}
    for node in G.nodes():
        grow = []
        r_g = []
        num_g = []
        num_nodes = 0
        if weight == None:
            spl = nx.single_source_shortest_path_length(G, node)
        else:
            spl = nx.single_source_dijkstra_path_length(G, node)
        for s in spl.values():
            if s > 0:
                grow.append(s)
        grow.sort()
        num = Counter(grow)
        for i, j in num.items():
            num_nodes += j
            if i > 0:
                # if np.log(num_nodes) < 0.95*np.log(G.number_of_nodes()):
                r_g.append(i)
                num_g.append(num_nodes)
                if np.log(num_nodes) > 0.9 * np.log(G.number_of_nodes()):
                    break
        x = np.log(r_g)
        y = np.log(num_g)
        #         if len(r_g) < 3:
        #             print("local",node)

        if len(r_g) > 1:
            try:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            except:
                print(node)
            node_dimension[node] = slope
        else:
            node_dimension[node] = 0
    nx.set_node_attributes(G, node_dimension, 'nfd')
    return G

def FDC(G, weight=None):
    node_dimension = {}
    Gn = G.number_of_nodes()
    for node in G.nodes():
        grow = []
        r_g = []
        num_g = []
        num_nodes = 0
        deg = 0
        if weight == None:
            spl = nx.single_source_shortest_path_length(G, node)
        else:
            spl = nx.single_source_dijkstra_path_length(G, node)
        for s in spl.values():
            if s > 0:
                grow.append(s)
        grow.sort()
        num = Counter(grow)
        deg = list(num.values())[0]
        for i, j in num.items():
            num_nodes += j
            if i > 0:
                # if np.log(num_nodes) < 0.95*np.log(G.number_of_nodes()):
                r_g.append(i)
                num_g.append(num_nodes)
                if np.log(num_nodes) > 0.9 * np.log(Gn):
                    break
        # if len(r_g) < 3:
        # print("local",node,len(r_g))
        if len(r_g) > 1:
            try:
                slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(r_g), np.log(num_g))
            except:
                print(node)
            node_dimension[node] = slope / np.log(Gn / deg)
        else:
            node_dimension[node] = 0
    nx.set_node_attributes(G, node_dimension, 'fdc')
    return G

def thresholding(g, threshold):
    for edge in g.edges:
        if g.edges[edge]['weight'] < threshold:
            g.remove_edge(edge[0], edge[1])
    return g

def thresholding_low(g, upper_bound):
    """only keep edges lower than a certain threshold"""
    for edge in g.edges:
        if g.edges[edge]['weight'] > upper_bound:
            g.remove_edge(edge[0], edge[1])
    return g
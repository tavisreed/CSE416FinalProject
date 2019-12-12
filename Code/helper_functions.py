import pickle
import pypdb as pyp
from bioservices import *
import networkx as nx
import collections
import numpy as np
import scipy
import operator
import random
import community
import matplotlib.pyplot as plt
import os


# This file contains functions used in the main file or in the mining file. Some functions are no longer used,
# but have been left here just in case they are needed once more

# Here we create a class for Protien Structures
class Structure:
    def __init__(self, name):
        self.name = name
        self.ligands = set()
        self.subunits = set()
        self.functions = set()
        self.processes = set()

    def Add_Ligand(self, ligand_list):
        self.ligands = set(ligand_list)

    def Add_Subunits(self, subunits):
        self.subunits = set(subunits)

    def Add_Functions(self, function):
        self.functions = set(function)

    def Add_Processes(self, process):
        self.processes = set(process)


# Here we create a function for retrieving a subunit's functions or proccesses
def get_FP(string, type):
    list = []
    # Go through all the lines if the string
    for line in string.splitlines():

        if line.find(type) != -1:
            for item in line.split(";"):
                if item.find(type) != -1:
                    list.append(item.split(":")[1])
    return list


# get list of ligands for a structure
def get_ligands_and_uniprot(the_id):
    lig_list = []
    if pyp.get_ligands(the_id)['ligandInfo'] is not None:  # check if this is a dictionary
        # print("here")
        list_ligands = pyp.get_ligands(the_id)['ligandInfo']['ligand']

        if type(list_ligands) == list:  # multiple ligands -- list
            for index, lig in enumerate(list_ligands):
                # print(list_ligands[index]['@chemicalID'])
                lig_list.append(list_ligands[index]['@chemicalID'])

        else:  # one ligand -- dictionary
            lig_list.append(list_ligands['@chemicalID'])

        # get uniportkb ids of subunits
        u = UniProt()
        u.TIMEOUT = 1000
        res = u.search(the_id, columns='id', frmt='tab')
        if type(res) == int or res == None:
            print("Int or None")
            return None, None, [], []

        res = res.split()

        try:
            res.pop(0)
        except:
            return None, None, [], []

        # for residue in res:
        #
        functions = []
        processes = []
        # For each subunit find functions and processes
        if len(res) > 10:
            return None, None, [], []
        for subunit in res:
            info = u.retrieve(subunit, frmt='txt')
            temp_functions = get_FP(info, "F:")
            functions = functions + temp_functions
            temp_processes = get_FP(info, "P:")
            processes = processes + temp_processes
        functions = set(functions)
        processes = set(processes)
        # print(functions)
        # print(processes)
        # print(lig_list)
        return res, lig_list, functions, processes
    else:
        # print("no ligands")
        return None, None, [], []


# A function to save a dictionary
def saveDict(dictonary, filename):
    with open(filename + '.pkl', 'wb') as pkl_file:
        pickle.dump(dictonary, pkl_file, protocol=2)


# A function to read a dictionary
def readDict(filename, var_name):
    with open(filename + '.pkl', 'rb') as pkl_file:
        var_name = pickle.load(pkl_file)
    return var_name


# Create an edge in the graph between a property node and a protein node
def create_Edge(struct, G, property):
    prop_list = get_property(struct, property)
    for prop in prop_list:
        G.add_node(prop, type=property)
        G.add_node(struct.name, type='protein')
        G.add_edge(struct.name, prop)


# Get a given property list
def get_property(struct, property):
    if property == 'ligands':
        return struct.ligands
    elif property == 'subunits':
        return struct.subunits
    elif property == 'functions':
        return struct.functions
    elif property == 'processes':
        return struct.processes


# Get the average of a property in the dictonary
def get_mean_property(dictonary, property):
    mean = 0
    for (id, struct) in dictonary.items():
        mean += len(get_property(struct, property))
    return mean / len(dictonary)


# Get all unique types of a property
def get_all_property(dictonary, property):
    temp_list = []
    for (id, struct) in dictonary.items():
        temp_list = temp_list + list(get_property(struct, property))

    return set(temp_list)


# Get the degree distribution and the expected degree
def degree_dist(G):
    degree_vals_raw = list(nx.degree(G))
    degree_vals_raw[0][1]
    degree_vals = []
    for num in degree_vals_raw:
        degree_vals.append(num[1])

    degree_vals = sorted(degree_vals)
    iters = collections.Counter(degree_vals)
    x_axis = list(iters.keys())
    y_axis = list(iters.values())
    y_axis_norm = np.array(y_axis) / len(G.nodes())

    expected_degree_val = np.dot(list(x_axis), y_axis_norm)
    return x_axis, y_axis, expected_degree_val


# A simple sorter
def btw_sorter(elem):
    # Return the second element in the element tuple
    return elem[1]


# Calculate the modularity score of a graph or subgraph
def ModularityScore(G):
    # If the graph is empty, just return an error string
    try:
        # Get the adjaceny matrix
        A = nx.adjacency_matrix(G).todense()
    except:
        return "Empty_Graph"
    A = np.array(A)
    # Get the number of edges
    m = len(list(G.edges()))
    # Avoid divding by zero for a graph with no edegs
    if m == 0:
        m = 1
    # go through rows and cols in A
    rows, cols = A.shape
    score = 0
    for i in range(rows):
        for j in range(cols):
            Aij = A[i][j]
            ki = A[i].sum()
            kj = A[j].sum()
            # Calculate the value of each node pair
            score += Aij - (ki * kj) / (2 * m)
    return score


# Here we determine communities by betweeness of edges
def GN_BBC(oriG):
    # Here we create a empty list to hold our communites
    communities = list(nx.connected_components(oriG))
    # Here we copy the input graph so that we don't modify the original graph
    G = oriG.copy()
    score = 0
    for comm in communities:
        temp_subgraph = G.subgraph(comm)
        temp_score = ModularityScore(temp_subgraph)
        if temp_score == 'Empty_Graph':
            continue
        score += temp_score
    old_score = score - 1
    # We will continue finding communities until we reach the desired number
    # of communities

    while old_score < score:
        old_score = score
        old_communities = communities
        # Here we get the current number of communities, which is just the
        # number of connected components
        comm_size = nx.number_connected_components(G)
        # Here we calculate the edge betweeness of every edge in the graph
        edge_btw_dict = nx.edge_betweenness_centrality(G)
        # Here we convert the resulting ditconary to a list and then sort
        # by highest to lowest edge centrality
        edge_btw_list = [(edge, btw) for edge, btw in edge_btw_dict.items()]
        edge_btw_list.sort(key=btw_sorter, reverse=True)
        # We continue to delete edges (highest to lowest betweeness)
        # until we have created a new community
        while nx.number_connected_components(G) == comm_size:
            edge = edge_btw_list[0][0]
            G.remove_edge(edge[0], edge[1])
            del (edge_btw_list[0])
        communities = list(nx.connected_components(G))
    # Here we calculate the final modularity score
    score = 0
    for comm in communities:
        temp_subgraph = G.subgraph(comm)
        temp_score = ModularityScore(temp_subgraph)
        if temp_score == 'Empty_Graph':
            continue
        score += temp_score
    return old_communities, old_score


# Here we create a sorter function to sort a list
def eigs_sorter(elem):
    # Return the second element in the element tuple
    return elem[0]


# Here we find the eigenvalues of the modularity matrix of a graph
def B_EigvenValue(G):
    # Here we get the Adjaceny matrix
    A = nx.adjacency_matrix(G).todense()
    # Here we get teh degree matrix and compute dd^T/2m
    d = np.array([degree for (node, degree) in G.degree()])
    dd_2m = (d * d.transpose()) / (2 * len(list(G.edges())))
    # Here we compute the modularity matrix
    B = A - dd_2m
    # Here we get the eigenvalues and vectors of B
    eigenvalues, eigenvectors = np.linalg.eig(B)
    # Transpose the eigenvectors such that a row corresponds to an eigenvalue
    eigenvectors = eigenvectors.transpose()
    # Here we get the eigenvector corresponding to the max eigenvalue
    max_eigenvalue_index = np.argmax(eigenvalues)
    max_eigenvector = np.array(eigenvectors[max_eigenvalue_index])[0]

    return np.max(eigenvalues), max_eigenvector


# Here we find communities based on modularity maximization
def MM(oriG):
    # Here we copy the input graph so that we don't modify the original graph
    G = oriG.copy()
    current_communities = list(nx.connected_components(G))
    score = 0
    for comm in current_communities:
        temp_subgraph = G.subgraph(comm)
        temp_score = ModularityScore(temp_subgraph)
        if temp_score == 'Empty_Graph':
            continue
        score += temp_score
    old_score = score - 1
    # Continue you creating communites until we get k communities
    while old_score < score:
        old_score = score
        old_communities = current_communities
        # Here we create temp arrays to hold the proposed communities and their modularity score
        modularity_scores = []
        temp_communities = []
        # Go through all the current decided communites, split each one in two, and calculate the modularity
        # score of each possible divsion
        for node_set in current_communities:
            # Here we make the current community into a sperate subgraph so that it is properly isolated
            subgraph = G.subgraph(node_set)
            # Here we compute the eigenvalues/vectors of the modularity matrix of the subgraph
            temp_eig, temp_eigV = B_EigvenValue(subgraph)
            # Here we split the current subgraph into two new subgraphs
            temp_community1 = [node for (index, node) in enumerate(list(subgraph.nodes())) if temp_eigV[index] >= 0]
            temp_community2 = [node for (index, node) in enumerate(list(subgraph.nodes())) if temp_eigV[index] < 0]
            # Add all communites that are not being split to a temp list
            temp_total_community = [other_node_set for other_node_set in current_communities if
                                    other_node_set != node_set]
            # Add the two new communites to the list
            temp_total_community.append(temp_community1)
            temp_total_community.append(temp_community2)
            score = 0
            for comm in temp_total_community:
                temp_subgraph = G.subgraph(comm)
                score += ModularityScore(temp_subgraph)
                temp_score = ModularityScore(temp_subgraph)
                if temp_score == 'Empty_Graph':
                    continue
                score += temp_score

            modularity_scores.append(score)
            temp_communities.append(temp_total_community)
        # Find the highest modularity score, and make the corresponding community the current community
        max_modularity_index = np.argmax(modularity_scores)
        current_communities = temp_communities[max_modularity_index]

    # Here we calculate the final modularity score
    score = 0
    for comm in current_communities:
        temp_subgraph = G.subgraph(comm)
        score += ModularityScore(temp_subgraph)
        temp_score = ModularityScore(temp_subgraph)
        if temp_score == 'Empty_Graph':
            continue
        score += temp_score

    return old_communities, old_score


# Here we create a sorter function to sort a list
def eigs_sorter(elem):
    # Return the second element in the element tuple
    return elem[0]


# Here we find the eigenvalues of the modularity matrix of a graph
def L_EigvenValue(G):
    # Here we get the Laplacian matrix
    L = nx.laplacian_matrix(G).todense()
    # Here we get the eigenvalues and vectors of L
    eigenvalues, eigenvectors = np.linalg.eig(L)
    # Transpose the eigenvectors such that a row corresponds to an eigenvalue
    eigenvectors = np.array(eigenvectors.transpose())
    # Here we get the eigenvector corresponding to the smallest eigenvalue and toss it by setting it to infinity
    min_eigenvalue_index = np.argmin(eigenvalues)
    eigenvalues[min_eigenvalue_index] = np.inf
    eigenvectors[min_eigenvalue_index] = np.inf

    # Now we get the second smallest eigenvalue and correspoonding vector, which is now the smallest in the array
    min_eigenvalue_index = np.argmin(eigenvalues)
    min_eigenvector = eigenvectors[min_eigenvalue_index]
    return np.min(eigenvalues), min_eigenvector


# Here we find communities based on spectral clustering
def SpecClust(oriG):
    # Here we copy the input graph so that we don't modify the original graph
    G = oriG.copy()
    current_communities = list(nx.connected_components(G))
    score = 0;

    for comm in current_communities:
        temp_subgraph = G.subgraph(comm)
        temp_score = ModularityScore(temp_subgraph)
        if temp_score == 'Empty_Graph':
            continue
        score += temp_score
    old_score = score - 1
    # Continue you creating communites until we get k communities
    counter = 0
    while old_score < score:
        print("Counter:", counter)
        old_score = score
        old_communities = current_communities
        # Here we create temp arrays to hold the proposed communities and their modularity score
        modularity_scores = []
        temp_communities = []
        # Go through all the current decided communites, split each one in two, and calculate the modularity
        # score of each possible divsion
        for node_set in current_communities:
            # Here we make the current community into a sperate subgraph so that it is properly isolated
            subgraph = G.subgraph(node_set)
            # Here we compute the eigenvalues/vectors of the modularity matrix of the subgraph
            temp_eig, temp_eigV = L_EigvenValue(subgraph)
            # Here we split the current subgraph into two new subgraphs
            temp_community1 = [node for (index, node) in enumerate(list(subgraph.nodes())) if temp_eigV[index] >= 0]
            temp_community2 = [node for (index, node) in enumerate(list(subgraph.nodes())) if temp_eigV[index] < 0]
            # Add all communites that are not being split to a temp list
            temp_total_community = [other_node_set for other_node_set in current_communities if
                                    other_node_set != node_set]
            # Add the two new communites to the list
            temp_total_community.append(temp_community1)
            temp_total_community.append(temp_community2)
            score = 0
            for comm in temp_total_community:
                temp_subgraph = G.subgraph(comm)
                temp_score = ModularityScore(temp_subgraph)
                if temp_score == 'Empty_Graph':
                    continue
                score += temp_score

            modularity_scores.append(score)
            temp_communities.append(temp_total_community)
        # Find the highest modularity score, and make the corresponding community the current community
        max_modularity_index = np.argmax(modularity_scores)
        current_communities = temp_communities[max_modularity_index]
        counter = counter + 1
    # Here we calculate the modularity score of the final clustering
    score = 0
    for comm in current_communities:
        temp_subgraph = G.subgraph(comm)
        # score += ModularityScore(temp_subgraph)
        temp_score = ModularityScore(temp_subgraph)
        if temp_score == 'Empty_Graph':
            continue
        score += temp_score

    return old_communities, old_score


# Here we find the size of each community in a list of communities
def Size_of_Comms(Comm_list):
    Comm_Size_List = []
    for community in Comm_list:
        Comm_Size_List.append(len(community))
    return Comm_Size_List


# Here we convert a dictonary of nodes labled by community into a list of communities
def Get_Community(comms_dict):
    total_comms = max(list(comms_dict.values()))
    comm_list = [[] for i in range(total_comms + 1)]

    for (name, comm) in comms_dict.items():
        comm_list[comm].append(name)

    return comm_list


# Here we get the k most frequent properties in each community in a list of communities
def Get_K_Properties(comm_list, struct_dict, k, property):
    prop_freq_dict = {}
    comm_size = len(comm_list)
    k_prop_set = set()
    for (name, obj) in struct_dict.items():
        if name in comm_list:
            obj_prop = get_property(obj, property)
            for prop in obj_prop:
                if prop in prop_freq_dict:
                    prop_freq_dict[prop] += 1
                else:
                    prop_freq_dict[prop] = 1
    seen_set = set()
    while len(k_prop_set) < k and len(prop_freq_dict) != 0:
        prop = max(prop_freq_dict.items(), key=operator.itemgetter(1))[0]
        if prop not in seen_set:
            freq = round(prop_freq_dict[prop] / comm_size, 3)
            k_prop_set.add((prop, freq))
        del prop_freq_dict[prop]
    return k_prop_set


# Here we get the similarity score of a node
def similarity_score(node, k_prop_list, property):
    node_props = get_property(node, property)
    pure_k_prop = [prop for (prop, freq) in k_prop_list]
    similarity = 0
    for prop in node_props:
        if prop in pure_k_prop:
            similarity += 1

    try:
        score = round(similarity / len(k_prop_list), 3)
    except:
        score = 0

    return score


# Here we get the similarity score of a community
def community_score(comm_list, struct_dict, k_prop_list, property):
    comm_score = 0
    for node in comm_list:
        node_obj = struct_dict[node]
        sim_score = similarity_score(node_obj, k_prop_list, property)
        comm_score += sim_score
    try:
        comm_score = round(comm_score / len(comm_list), 3)
    except:
        comm_score = 0
    return comm_score


# Here we randomly assign nodes to communities. There is alread a set number of communities and community sizes
def create_random_communites(comms_dict):
    total_comms = max(list(comms_dict.values()))
    comm_list = Get_Community(comms_dict)
    allocation_limit_list = [len(comm) for comm in comm_list]
    allocation_list = [0 for comm in comm_list]
    rand_comm_dict = {}
    for (name, value) in comms_dict.items():
        allocated = False
        while allocated == False:
            rand_comm = random.randint(0, total_comms)
            if allocation_list[rand_comm] < allocation_limit_list[rand_comm]:
                rand_comm_dict[name] = rand_comm
                allocation_list[rand_comm] += 1
                allocated = True
    return rand_comm_dict


# Here we randomly assign nodes to communities. There is alread a set number of communities but sizes are random
def create_true_random_communites(comm_list):
    total_num_comms = len(comm_list)
    rand_comms = []
    for comm in comm_list:
        rand_comms.append([])
    for comm in comm_list:
        for node in comm:
            assignment = random.randrange(0, total_num_comms, 1)
            rand_comms[assignment].append(node)

    return rand_comms


# Here we get the similarity score of the graph
def score_graph(comms_dict, struc_dict, k, property, already_list=False):
    if already_list == True:
        comm_list = comms_dict
    else:
        comm_list = Get_Community(comms_dict)
    scores = []
    score_weight = []
    for comm in comm_list:
        k_prop_list = Get_K_Properties(comm, struc_dict, k, property)
        comm_score = community_score(comm, struc_dict, k_prop_list, property)
        scores.append(comm_score)
        score_weight.append(len(comm))
    graph_score = sum(scores[i] * score_weight[i] for i in range(len(scores))) / sum(score_weight)
    return round(graph_score, 3)


# Here we compare the score of a graph to the score of a random graph
def compared_to_random(comms_dict, struc_dict, k, property):
    rand_comm_dict = create_random_communites(comms_dict)
    comm_score = score_graph(comms_dict, struc_dict, k, property)
    rand_score = score_graph(rand_comm_dict, struc_dict, k, property)
    diff = round(comm_score - rand_score, 3)

    return (diff, (comm_score, rand_score)), rand_comm_dict


# A sorter
def tuple_sorter(item):
    return item[1][0]


# A sorter
def size_score_sorter(item):
    return item[0]


# Here we optimize the louvian resolution to get the best graph score
def optimize_louv(G, struct_dict, start_res, step_size, property, k):
    opt_list = []
    comm_size_score_list = []
    old_score = (-500, 0)
    diff_score = (-499, 0)
    res = start_res
    while old_score[0] < diff_score[0]:
        old_score = diff_score
        temp_louv = community.best_partition(G, resolution=res)
        num_coms = len(Get_Community(temp_louv))
        diff_score, rand_dict = compared_to_random(temp_louv, struct_dict, k, property)
        comm_size_score_list.append((num_coms, diff_score[1][0]))
        opt_list.append((temp_louv, diff_score, res, rand_dict, num_coms))
        # print("Resolution Score", res, diff_score)
        res += step_size

    opt_list.sort(reverse=True, key=tuple_sorter)
    comm_size_score_list.sort(reverse=True, key=size_score_sorter)
    louv = opt_list[0][0]
    diff_score = opt_list[0][1]
    res = opt_list[0][2]
    rand_dict = opt_list[0][3]
    num_coms = opt_list[0][4]
    return (louv, diff_score, num_coms, res, rand_dict, comm_size_score_list)


# Here we get the score of each individual community
def community_score_list(comm_dict, struct_dict, k_property, property, already_list=False):
    if already_list == False:
        comm_list = Get_Community(comm_dict)
    else:
        comm_list = comm_dict
    comm_score_list = []
    for comm in comm_list:
        size = len(comm)
        k_prop_list = Get_K_Properties(comm, struct_dict, k_property, property)
        score = community_score(comm, struct_dict, k_prop_list, property)
        comm_score_list.append((size, score))
    return comm_score_list


# Here we create a projected list, in which edges exist between protein nodes if the share some k number of
# specific traits in common
def create_projected_graph(Struct_Dict, k_common, property, percent=False):
    G = nx.Graph()
    edge_list = []
    seen_edges = set()
    seen_struct = set()
    for id_1, struct_1 in Struct_Dict.items():
        for id_2, struct_2 in Struct_Dict.items():
            if id_1 == id_2:
                continue
            prop_1 = get_property(struct_1, property)
            prop_2 = get_property(struct_2, property)
            prop_inter = prop_1.intersection(prop_2)

            if percent == True:
                prop_union = prop_1.union(prop_2)

                IOU = len(prop_inter) / len(prop_union)

                if IOU >= k_common:
                    if (id_1, id_2) in seen_edges or (id_2, id_1) in seen_edges:
                        continue
                    seen_struct.add(id_1)
                    seen_struct.add(id_2)
                    edge_list.append((id_1, id_2))
                    seen_edges.add((id_1, id_2))
                    seen_edges.add((id_2, id_1))
                    G.add_edge(id_1, id_2)

            else:
                if len(prop_inter) >= k_common:
                    if (id_1, id_2) in seen_edges or (id_2, id_1) in seen_edges:
                        continue
                    seen_struct.add(id_1)
                    seen_struct.add(id_2)
                    edge_list.append((id_1, id_2))
                    seen_edges.add((id_1, id_2))
                    seen_edges.add((id_2, id_1))
                    G.add_edge(id_1, id_2)
    for id, struct in Struct_Dict.items():
        if id not in seen_struct:
            G.add_node(id)

    return G


# Get the number of nodes in a community
def num_nodes(comm_list):
    counter = 0
    for comm in comm_list:
        for node in comm:
            counter += 1
    return counter


# Delete communitites based off their modularity score and its distance from the mean
def delete_comms(G, comm_list, allowed_std):
    modularity_score_list = []
    refined_comm_list = []
    for comm in comm_list:
        comm_graph = G.subgraph(comm)
        mod_score = ModularityScore(comm_graph)
        if mod_score == "Empty_Graph":
            mod_score = 0
        modularity_score_list.append(mod_score)

    mod_std = np.std(modularity_score_list)
    mod_mean = np.mean(modularity_score_list)
    for index, comm in enumerate(comm_list):
        if modularity_score_list[index] >= mod_mean + mod_std * allowed_std:
            refined_comm_list.append(comm)

    return refined_comm_list


# Plot community scores vs. the scores of random communities
def plot_vs_random(comm_list, struct_dict, k_property, property, plot_save_name, show_plots=True):
    try:
        os.remove(plot_save_name + ".png")
    except:
        skip = 1
    rand_comm = create_true_random_communites(comm_list)
    size_score_list_comm = community_score_list(comm_list, struct_dict, k_property, property, already_list=True)
    size_score_list_rand = community_score_list(rand_comm, struct_dict, k_property, property, already_list=True)
    comm_size_axis = [index for index, score in enumerate(size_score_list_comm)]
    comm_score_axis = [score for (size, score) in size_score_list_comm]
    rand_score_axis = [score for (size, score) in size_score_list_rand]

    plt.figure()

    plt.plot(comm_size_axis, comm_score_axis, 'bo', label='Communities')
    plt.plot(comm_size_axis, rand_score_axis, 'ro', label='Random Communities')

    plt.title("Community Similarity Scores")
    plt.xlabel("Community")
    plt.ylabel("Community Score")
    plt.legend()
    plt.savefig(plot_save_name + '.png')
    if show_plots == True:
        plt.show()


# Convert a list of communities in to a dictonary in which each node has a community number
def list_to_dict(comm_list):
    comm_dict = {}
    for comm_id, comm in enumerate(comm_list):
        for node in comm:
            comm_dict[node] = comm_id

    return comm_dict


# A sorter
def opt_sorter(item):
    return item[1]


# Find best k for the k_clique community detection
def opt_k_clique(G, start_k, end_k, num_trials):
    super_comm_list = []
    temp_k = []
    temp_score = []
    break_flag = False
    k = start_k
    while k <= end_k:
        temp_k.append(k)
        score = 0
        for i in range(num_trials):
            try:
                comm_list = nx.algorithms.community.k_clique_communities(G, k)
            except:
                break_flag = True
                break
            comm_list = list(list(comm_list))
            for comm in comm_list:
                comm_graph = G.subgraph(comm)
                mod_score = ModularityScore(comm_graph)
                if mod_score == "Empty_Graph":
                    mod_score = 0
                score += mod_score
        if break_flag == True:
            del (temp_k[-1])
            break
        score = score / num_trials
        temp_score.append(score)
        super_comm_list.append((comm_list, score, k))
        k += 1
    super_comm_list.sort(reverse=True, key=opt_sorter)
    comm_list = super_comm_list[0][0]
    score = super_comm_list[0][1]
    k = super_comm_list[0][2]
    plt.figure()
    plt.plot(temp_k, temp_score, 'bo')
    plt.title("Modularity Score Per K Clique Size")
    plt.xlabel("K")
    plt.ylabel("Modularity Score")
    plt.savefig('Optimize_K_Clique.png')
    plt.show()
    print(k, score)
    return comm_list


# Find best number of communities for Fluid community detection
def opt_fluid(G, start_comms, end_comms, step_size, num_trials):
    super_comm_list = []
    temp_comms = []
    temp_score = []
    break_flag = False
    num_comms = start_comms
    while num_comms <= end_comms:
        temp_comms.append(num_comms)
        score = 0
        for i in range(num_trials):
            try:
                comm_list = nx.algorithms.community.asyn_fluid.asyn_fluidc(G, num_comms)
            except:
                break_flag = True
                break
            comm_list = list(list(comm_list))
            for comm in comm_list:
                comm_graph = G.subgraph(comm)
                mod_score = ModularityScore(comm_graph)
                if mod_score == "Empty_Graph":
                    mod_score = 0
                score += mod_score
        if break_flag == True:
            del (temp_comms[-1])
            break
        score = score / num_trials
        temp_score.append(score)
        super_comm_list.append((comm_list, score, num_comms))
        num_comms += step_size

    super_comm_list.sort(reverse=True, key=opt_sorter)
    comm_list = super_comm_list[0][0]
    score = super_comm_list[0][1]
    num_comms = super_comm_list[0][2]
    plt.figure()
    plt.plot(temp_comms, temp_score, 'bo')
    plt.title("Modularity Score Per Number of Communities")
    plt.xlabel("Number of Communities")
    plt.ylabel("Modularity Score")
    plt.savefig('Optimize_Fluid.png')
    plt.show()
    print(num_comms, score)
    return comm_list

# %% importing things
import numpy as np
import networkx as nx
import helper_functions as hf
from helper_functions import Structure
import importlib
import matplotlib.pyplot as plt
import time
import community

#This make sure that we are always using the most updated version of helper_functions.py
hf = importlib.reload(hf)

#This function serves as an easy way to run various tests /  algorithims on the data
def master(struct_save_name="ProteinDict_ten_thousand", edge_type="ligands", edge_comm_num=3, property="processes",
           graph_filename="Protein-Protein_Graph_Default_Name", load_graph=False, print_dict_props=False,
           bipart_graph=False, bipartite_filename="Bipartite_Default_Name", show_plots=False, avg_clust=False,
           print_graph_props=False, degree_dist=False, k_clique=False, mod_max=False, fluid=False, louv=False,
           k_property=20, num_k_cliques=7, num_fluid_comms=100, std_val=-0.5, k_clique_opt = False, start_k_clique_opt = 3,
           end_k_clique_opt = 10, num_trials_k = 3, opt_fluid=False, start_fluid_comms= 100, end_fluid_comms=300,
           fluid_step_size=20, fluid_num_trials=3):
    Structure_Dict = {}
    Structure_Dict = hf.readDict(struct_save_name, Structure_Dict)

    #Here we print out some helpful information about the dataset we are using
    if print_dict_props == True:
        avg_ligands = hf.get_mean_property(Structure_Dict, "ligands")
        print("Average Number of Ligands:", avg_ligands)

        avg_subunits = hf.get_mean_property(Structure_Dict, "subunits")
        print("Average Number of Subunits:", avg_subunits)

        avg_functions = hf.get_mean_property(Structure_Dict, "functions")
        print("Average Number of Functions:", avg_functions)

        avg_processes = hf.get_mean_property(Structure_Dict, "processes")
        print("Average Number of Processes:", avg_processes)

        # Get Total Number of Ligands, Functions, Proccesses and Subunits
        num_ligands = len(hf.get_all_property(Structure_Dict, "ligands"))
        print("Number of Ligands:", num_ligands)

        num_subunits = len(hf.get_all_property(Structure_Dict, "subunits"))
        print("Number of Subunits:", num_subunits)

        num_functions = len(hf.get_all_property(Structure_Dict, "functions"))
        print("Number of Functions:", num_functions)

        num_processes = len(hf.get_all_property(Structure_Dict, "processes"))
        print("Number of Processes:", num_processes)

    #Here we create a bipartite graph of ligands and proteins, which can be analyzed on its own, or used to
    #to create a projected graph.
    if bipart_graph == True:
        Protein_Bipartite_Graph = nx.Graph()
        struct_name_set = set()
        # Create a bipartite graph in which there are structure nodes and ligand ndoes
        for (struct_name, struct) in Structure_Dict.items():
            struct_name_set.add(struct_name)
            hf.create_Edge(struct, Protein_Bipartite_Graph, property)

        print('Bipartite Nodes:', len(Protein_Bipartite_Graph.nodes()))
        print('Bipartite Edges:', len(Protein_Bipartite_Graph.edges()))
        nx.write_gml(Protein_Bipartite_Graph, bipartite_filename)

    #Here we create a new projected graph
    if load_graph == False:
        # Create a projected graph from the bipartite
        Protein_Graph = hf.create_projected_graph(Structure_Dict, edge_comm_num, edge_type)
        # Get the Giant Component of graph
        Protein_Graph_GC = Protein_Graph.subgraph(
            sorted(nx.connected_components(Protein_Graph), key=len, reverse=True)[0])
        nx.write_gml(Protein_Graph, graph_filename)

    #If the garph has already been created, load in the graph to save time
    if load_graph == True:
        Protein_Graph = nx.read_gml(graph_filename)
        Protein_Graph_GC = Protein_Graph.subgraph(
            sorted(nx.connected_components(Protein_Graph), key=len, reverse=True)[0])

    #Print out some useful informatoion about the graph
    if print_graph_props == True:
        print('Protein_Graph Nodes:', len(Protein_Graph.nodes()))
        print('Protein_Graph Edges:', len(Protein_Graph.edges()))
        print('Protein_Graph Num connected Components:', nx.number_connected_components(Protein_Graph))
        print('Protein_Graph Num edges in largest Components:', len(Protein_Graph_GC.edges()))
        print('Protein_Graph Num nodes in largest Components:', len(Protein_Graph_GC.nodes()))


    # K-Clique Implementation
    if k_clique == True:
        print('Begin K_Clique')
        #Create a copy of the graph, which will be used when we lable nodes by community
        k_clique_graph = Protein_Graph_GC.copy()

        #You can use a predetermined k, or optimize the k for the graph
        if k_clique_opt ==False:
            k_clique_comms_pre_del = nx.algorithms.community.k_clique_communities(Protein_Graph_GC, num_k_cliques)
            k_clique_comms_pre_del = list(list(k_clique_comms_pre_del))
        else:
            k_clique_comms_pre_del = hf.opt_k_clique(Protein_Graph_GC, start_k_clique_opt, end_k_clique_opt, num_trials_k)

        # Get the average size of found communities
        avg_comm_pre_del = sum([len(comm) for comm in k_clique_comms_pre_del]) / len(k_clique_comms_pre_del)

        #Get the graph similiarty score
        K_clique_score_pre_del = hf.score_graph(k_clique_comms_pre_del, Structure_Dict, k_property, property,
                                                already_list=True)
        print(K_clique_score_pre_del, len(k_clique_comms_pre_del), avg_comm_pre_del, hf.num_nodes(k_clique_comms_pre_del))

        #Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(k_clique_comms_pre_del, Structure_Dict, k_property, property,
                          "K_Clique_" + str(k_property) + "_" + property + "_Pre_Del_Comms_"+edge_type+"_edges", show_plots=show_plots)

        #Delete some communities based off there modularity score, and the standard deviation of community scores in the graph
        k_clique_comms = hf.delete_comms(Protein_Graph_GC, k_clique_comms_pre_del, std_val)

        #Get the graph similiarty score after deletion
        k_clique_score = hf.score_graph(k_clique_comms, Structure_Dict, k_property, property, already_list=True)

        # Get the average size of found communities after deleting 'bad' communities
        avg_comm= sum([len(comm) for comm in k_clique_comms]) / len(k_clique_comms)
        print(k_clique_score, len(k_clique_comms), avg_comm, hf.num_nodes(k_clique_comms))

        # Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(k_clique_comms, Structure_Dict, k_property, property,
                          "K_Clique_" + str(k_property) + "_" + property + "_Comms_"+edge_type+"_edges", show_plots=show_plots)

        # Label nodes by community
        nx.set_node_attributes(k_clique_graph, hf.list_to_dict(k_clique_comms_pre_del), "Community")

        #Save the graph with nodes labled by community
        nx.write_gml(k_clique_graph, "K_Clique_Protein_Protein_" + edge_type + "_edges_Network_" + str(
            k_property) + "_" + property + ".gml")
        print('End K_Clique')

    # Modularity Maximization Implementation
    if mod_max == True:
        print('Begin Modularity Maximization')
        # Create a copy of the graph, which will be used when we lable nodes by community
        mod_graph = Protein_Graph_GC.copy()

        #Find communities using modularity maximization
        mod_max_comms_pre_del = nx.algorithms.community.modularity_max.greedy_modularity_communities(Protein_Graph_GC)
        mod_max_comms_pre_del = list(list(mod_max_comms_pre_del))

        # Get the average size of found communities
        avg_comm_pre_del = sum([len(comm) for comm in mod_max_comms_pre_del]) / len(mod_max_comms_pre_del)

        # Get the graph similiarty score
        mod_max_score_pre_del = hf.score_graph(mod_max_comms_pre_del, Structure_Dict, k_property, property,
                                               already_list=True)
        print(mod_max_score_pre_del, len(mod_max_comms_pre_del), avg_comm_pre_del, hf.num_nodes(mod_max_comms_pre_del))

        # Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(mod_max_comms_pre_del, Structure_Dict, k_property, property,
                          "Mod_Max" + str(k_property) + "_" + property + "_Pre_Del_Comms_"+edge_type+"_edges", show_plots=show_plots)

        # Delete some communities based off there modularity score, and the standard deviation of community scores in the graph
        mod_max_comms = hf.delete_comms(Protein_Graph_GC, mod_max_comms_pre_del, std_val)

        # Get the graph similiarty score after deletion
        mod_max_score = hf.score_graph(mod_max_comms, Structure_Dict, k_property, property, already_list=True)

        # Get the average size of found communities after deleting 'bad' communities
        avg_comm = sum([len(comm) for comm in mod_max_comms]) / len(mod_max_comms)
        print(mod_max_score, len(mod_max_comms), avg_comm, hf.num_nodes(mod_max_comms))

        # Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(mod_max_comms, Structure_Dict, k_property, property,
                          "Mod_Max" + str(k_property) + "_" + property + "_Comms_"+edge_type+"_edges", show_plots=show_plots)

        # Label nodes by community
        nx.set_node_attributes(mod_graph, hf.list_to_dict(mod_max_comms_pre_del), "Community")

        # Save the graph with nodes labled by community
        nx.write_gml(mod_graph, "Mod_Max_Protein_Protein_" + edge_type + "_edges_Network_" + str(
            k_property) + "_" + property + ".gml")
        print('End Modularity Maximization')

    # Fluid Implementation
    if fluid == True:
        print('Begin Fluid')
        # Create a copy of the graph, which will be used when we lable nodes by community
        fluid_graph = Protein_Graph_GC.copy()

        # You can use a predetermined number of communities, or optimize the number of communieis for the graph
        if opt_fluid == False:
            fluid_comms_pre_del = nx.algorithms.community.asyn_fluid.asyn_fluidc(Protein_Graph_GC, num_fluid_comms)
            fluid_comms_pre_del = list(list(fluid_comms_pre_del))
        else:
            fluid_comms_pre_del = hf.opt_fluid(Protein_Graph_GC, start_fluid_comms, end_fluid_comms, fluid_step_size, fluid_num_trials)

        # Get the average size of found communities
        avg_comm_pre_del = sum([len(comm) for comm in fluid_comms_pre_del]) / len(fluid_comms_pre_del)

        # Get the graph similiarty score
        fluid_score_pre_del = hf.score_graph(fluid_comms_pre_del, Structure_Dict, k_property, property,
                                             already_list=True)
        print(fluid_score_pre_del, len(fluid_comms_pre_del), avg_comm_pre_del, hf.num_nodes(fluid_comms_pre_del))

        # Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(fluid_comms_pre_del, Structure_Dict, k_property, property,
                          "Fluid" + str(k_property) + "_" + property + "_Pre_Del_Comms_"+edge_type+"_edges", show_plots=show_plots)

        # Delete some communities based off there modularity score, and the standard deviation of community scores in
        # the graph
        fluid_comms = hf.delete_comms(Protein_Graph_GC, fluid_comms_pre_del, std_val)

        # Get the graph similiarty score after deletion
        fluid_score = hf.score_graph(fluid_comms, Structure_Dict, k_property, property, already_list=True)

        # Get the average size of found communities after deleting 'bad' communities
        avg_comm = sum([len(comm) for comm in fluid_comms]) / len(fluid_comms)
        print(fluid_score, len(fluid_comms), avg_comm, hf.num_nodes(fluid_comms))

        # Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(fluid_comms, Structure_Dict, k_property, property,
                          "Fluid" + str(k_property) + "_" + property + "_Comms_"+edge_type+"_edges", show_plots=show_plots)

        # Label nodes by community
        nx.set_node_attributes(fluid_graph, hf.list_to_dict(fluid_comms_pre_del), "Community")

        # Save the graph with nodes labled by community
        nx.write_gml(fluid_graph, "Fluid_Protein_Protein_" + edge_type + "_edges_Network_" + str(
            k_property) + "_" + property + ".gml")
        print('End Fluid')

    # louvian Implmentation
    if louv == True:
        print('Begin Louvain')
        # Create a copy of the graph, which will be used when we lable nodes by community
        louv_graph = Protein_Graph_GC.copy()

        #Create communities using the louvian
        opt_louv = hf.optimize_louv(Protein_Graph_GC, Structure_Dict, 100, 1, property, k_property)
        louv_comm_pre_del = hf.Get_Community(opt_louv[0])

        # Get the average size of found communities
        avg_comm_pre_del = sum([len(comm) for comm in louv_comm_pre_del]) / len(louv_comm_pre_del)

        # Get the graph similiarty score
        louv_score_pre_del = hf.score_graph(louv_comm_pre_del, Structure_Dict, k_property, property, already_list=True)
        print(louv_score_pre_del, len(louv_comm_pre_del), avg_comm_pre_del, hf.num_nodes(louv_comm_pre_del))

        # Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(louv_comm_pre_del, Structure_Dict, k_property, property,
                          "Louv" + str(k_property) + "_" + property + "_Pre_Del_Comms_"+edge_type+"_edges", show_plots=show_plots)

        # Delete some communities based off there modularity score, and the standard deviation of community scores in the graph
        louv_comms = hf.delete_comms(Protein_Graph_GC, louv_comm_pre_del, std_val)

        # Get the graph similiarty score after deletion
        louv_score = hf.score_graph(louv_comms, Structure_Dict, k_property, property, already_list=True)

        # Get the average size of found communities after deleting 'bad' communities
        avg_comm = sum([len(comm) for comm in louv_comms]) / len(louv_comms)
        print(louv_score, len(louv_comms), avg_comm, hf.num_nodes(louv_comms))

        # Create a plot of each community similiarity score vs. a random communities similarity score
        hf.plot_vs_random(louv_comms, Structure_Dict, k_property, property,
                          "Louv" + str(k_property) + "_" + property + "_Comms_"+edge_type+"_edges", show_plots=show_plots)

        # Label nodes by community
        nx.set_node_attributes(louv_graph, hf.list_to_dict(louv_comm_pre_del), "Community")

        # Save the graph with nodes labled by community
        nx.write_gml(louv_graph, "Louv_Protein_Protein_" + edge_type + "_edges_Network_" + str(
            k_property) + "_" + property + ".gml")
        print('End Louvain')

    # Create Degree Distribution Plot and print out the expexted degree of the node
    if degree_dist == True:
        x, y, expected_degree = hf.degree_dist(Protein_Graph_GC)
        print("Expected Degree:", expected_degree)
        plt.figure()
        plt.loglog(x, y, 'bo')
        plt.title("Degree distribution")
        plt.xlabel("log(degree values)")
        plt.ylabel("log(degree frequencies)")
        plt.savefig('degree_dist_'+edge_type+'.png')
        plt.show()

    #Find the average clustering coefficient of the graph
    if avg_clust == True:
        average_clustering = nx.average_clustering(Protein_Graph_GC)
        print("Average Clustering Coefficient:", average_clustering)


#An example call of the function
master(show_plots=False, load_graph=True, graph_filename="Protein_Protein_Multi_Ligand_Network_k_3.gml",
       edge_type="ligands", edge_comm_num=3, property="functions", print_dict_props=False, print_graph_props=True,
       k_clique=False, mod_max=False, fluid=True, louv=False, k_property=7, num_k_cliques=7, num_fluid_comms=100,
       degree_dist=False, avg_clust=False, k_clique_opt = True, start_k_clique_opt = 3, end_k_clique_opt = 30, num_trials_k = 3,
       start_fluid_comms= 100, end_fluid_comms=300, fluid_step_size=5, fluid_num_trials=5, opt_fluid=True)

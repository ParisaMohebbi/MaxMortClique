# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 17:35:30 2024

@author: pvaghfi
"""

########################################### Module for Finding Max Clique ###################################################

##import packages
import copy
import networkx as nx
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import *
import pandas as pd
import os.path
import csv 
import json
from networkx.readwrite import json_graph
import time

excluded_cliques = []

### Accessories ###
## Graph Preprocessing
def FindExp(data):      
    deadly_diseases = set()                           
    for patient_id, patient_data in data.items():
        status, diseases = patient_data
        if status == 1:
            deadly_diseases.update(diseases)
    return deadly_diseases

## Function for generating a graph from deadly diseases
def MakeGraph(nodes_set, edge_list):
    ## Create an empty graph
    G = nx.Graph()
    ## Add nodes to the graph
    G.add_nodes_from(nodes_set)    
    ## Read the edge list from the file and add edges to the graph
    with open(edge_list, 'r') as f:
        for line in f:
            source, target, attributes = line.split()
            if source in nodes_set and target in nodes_set:
                G.add_edge(source, target)
    return G

def MIP_model(data, graph, excluded_cliques=[]):
    try:
        ## Create a new model
        ## Get connected components in graph
        Comp = [c for c in sorted(nx.connected_components(graph), key=len, reverse=True)]
        numComp = nx.number_connected_components(graph)
        print("The graph has ", graph.number_of_nodes(), "nodes and ", graph.number_of_edges(), " edges, and ", numComp, "connected components.")
        model = gp.Model("clique_mortality")
        C = list(range(numComp))
        M = list(data.keys())
        N = graph.nodes()

        ## Create variables
        x = model.addVars(N, vtype=GRB.BINARY, name="x")
        f = model.addVars(C, vtype=GRB.BINARY, name="f")
        
        ## Set objective
        model.setObjective(quicksum(x[i] for i in N), GRB.MAXIMIZE)
        
        ## Add clique constraints for each component
        for H in Comp:
            for i in H:
                for j in H:
                    if j not in graph.neighbors(i) and int(i) < int(j):
                        model.addConstr(x[i] + x[j] <= 1, "clique_constr" + str(i) + "," + str(j))

        ## Add exclusion constraints for previously found cliques
        for clique in excluded_cliques:
            model.addConstr(quicksum(x[v] for v in clique) <= len(clique) - 1, name="exclude_clique")

        ## Add additional constraints
        model.addConstr(quicksum(f[i] for i in C) == 1)
        model.addConstrs(x[i] <= f[j] for j in C for i in Comp[j])

        # Optimize model
        model.Params.timeLimit = 10500
        model.optimize()

        status = model.Status

        if status == GRB.Status.TIME_LIMIT or status == GRB.Status.OPTIMAL:
            x_solution = model.getAttr('x', x)
            selected_vertices = [i for i in N if x_solution[i] > 0.9]
            print(f"Clique found: {selected_vertices}")
            model.write('formulation.lp') #write MIP constraints
            return selected_vertices
        else:
            print('Optimization stopped with status ' + str(status))
            return None

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        print('Encountered an attribute error')

def main():
    with open('random_sample_1k_1.json', 'r') as f:
        data = json.load(f)

    start_pp = time.time()
    deadly_diseases = FindExp(data)
    
    ## Making the graph from the edgelist file
    nodes_set = deadly_diseases
    edge_list = "Comorbidity_Network_SCI=0.05.edgelist"
    graph = MakeGraph(nodes_set, edge_list)
    print("Preprocessing Time:", time.time() - start_pp)
    
    start_solve = time.time()

    excluded_cliques = []
    clique_sizes = []

    for i in range(1,4):
        ## Find maximum clique
        clique = MIP_model(data, graph, excluded_cliques)
        if clique is None:
            break

        excluded_cliques.append(clique)
        clique_sizes.append(len(clique))

        ## THIS WORKS BUT IT IS WRONG ## Remove nodes and edges of the clique, 
# =============================================================================
#         graph.remove_nodes_from(clique)
# =============================================================================
#         graph.remove_edges_from([(u, v) for u in clique for v in clique if u != v])
# =============================================================================

    print("Solving Time:", time.time() - start_solve)

if __name__ == "__main__":
    main()

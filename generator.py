from scipy.spatial.distance import pdist, squareform
import pandas as pd
import numpy as np


# Setup config
amount_of_neighbors = 4
# List of LA-arcs P
P = []
R = []
r_dict = {}
infinity = 1000000
# Set up LAs with nodes representing the numbers on a clock
nearest_neighbors = []
nodes = [i for i in range(12)]
for time in nodes:
    LA = (time, [])
    for neighbor in range(1, int(amount_of_neighbors / 2) + 1):    
        positive_neighbor = ((time + neighbor) % 12)
        negative_neighbor = ((time - neighbor) % 12)
        LA[1].append(positive_neighbor)
        LA[1].append(negative_neighbor)
    nearest_neighbors.append(LA)
    

clock_coordinates = [
    (3, 6),
    (4, 5),
    (5, 4),
    (6, 3),
    (5, 2),
    (4, 1),
    (3, 0),
    (2, 1),
    (1, 2),
    (0, 3),
    (1, 4),
    (2, 5),
    (3, 3),
    ]

time_windows = [
    (3, 6),
    (4, 5),
    (5, 4),
    (6, 3),
    (5, 2),
    (4, 1),
    (3, 0),
    (2, 1),
    (1, 2),
    (0, 3),
    (1, 4),
    (2, 5),
    (0, 10),
    ]

distances = pdist(clock_coordinates)
distance_matrix = squareform(distances)
distance_df = pd.DataFrame(distance_matrix, index=range(len(clock_coordinates)), columns=range(len(clock_coordinates)))


def get_info(nodes):
    key = tuple(sorted(nodes))
    if key in r_dict:
        return r_dict[key]
    else:
        return None

    
def add_route(nodes, cost, start_node, end_node, **kwargs):
    # Convert the set of nodes to a tuple and add it to the dictionary for unique id
    key = tuple(sorted(nodes))
    r_dict[key] = {'cost': cost, 'start_node': start_node, 'end_node': end_node}


# Function to get the power set of a given set. 
def get_power_set(s):
  power_set=[[]]
  for elem in s:
    # iterate over the sub sets so far
    for sub_set in power_set:
      # add a new subset consisting of the subset at hand added elem to it
      # effectively doubling the sets resulting in the 2^n sets in the powerset of s.
      power_set=power_set+[list(sub_set)+[elem]]
  return power_set


# Transform a list of lists of lists to a list of lists with only unique lists
def get_unique_sets(power_set):
    unique_sets = []
    for subset in power_set:
        for subsubset in subset:
            if subsubset not in unique_sets and subsubset != []:
                unique_sets.append(subsubset)
    return unique_sets


def get_P_from_u_and_N_u(tuple_with_u_and_N_u, all_nodes):
    # Possible v_p (end customer) nodes, cannot be equal to u or in N_u
    nodes_not_in_N_u_nor_u = [x for x in all_nodes if x not in tuple_with_u_and_N_u[1]]
    # Add depot to possible v_p (end customers)
    nodes_not_in_N_u_nor_u.append(12)
    # Remove u_p (start customer) from possible v_p
    nodes_not_in_N_u_nor_u.remove(tuple_with_u_and_N_u[0])
    # Create power set of LA neighbors of u
    power_set_of_N_u = get_power_set(tuple_with_u_and_N_u[1])
    #power_set_of_N_u.remove([])
    # Generate each p in P by iterating through all possible v_p and LA neighbors
#    for v_p in nodes:
#        if tuple_with_u_and_N_u[0] != v_p:
#            P.append((tuple_with_u_and_N_u[0], v_p, []))
    for v_p in tuple_with_u_and_N_u[1]:
        P.append((tuple_with_u_and_N_u[0], v_p, []))
        # Increased 
#        for neighbor_node in all_nodes: # tuple_with_u_and_N_u[1]
#            if v_p != neighbor_node:
#                P.append((tuple_with_u_and_N_u[0], v_p, [neighbor_node]))
    for v_p in nodes_not_in_N_u_nor_u:
        if v_p != tuple_with_u_and_N_u[0]:
            for subset in power_set_of_N_u:
#                if subset != []:
                P.append((tuple_with_u_and_N_u[0], v_p, subset))
        
        
def calculate_cost_for_p(p: tuple):
    # If there are no intermediate customers, the cost of the edge from start customer to end customer is used
    invalid_arc_count = 0
    if len(p[2]) == 0:
        cost = distance_df[p[0]][p[1]]
        p_tuple = (p[0], p[1], tuple(p[2]))
        r_dict[p_tuple] = cost
    if len(p[2]) == 1:
        # If there is one intermediate customer, the cost of the edges from start customer to the intermediate
        # customer plus the from the intermediate customer to the end customer are used
        cost = distance_df[p[0]][p[2][0]] + distance_df[p[1]][p[2][0]]
        p_tuple = (p[0], p[1], tuple(p[2]))
        r_dict[p_tuple] = cost
    if len(p[2]) > 1:
        # If there are multiple intermediate customers, the previously calculated paths are used to concatenate
        # LA-arcs with the same end/start customer.
        r_minus_dict = {}
        for node in p[2]:
            first_arc_cost   = distance_df[p[0]][node]
            neighbors_without_new_start_node = p[2][:]
            neighbors_without_new_start_node.remove(node)
            try:
                second_arc_cost = r_dict[node, p[1], tuple(neighbors_without_new_start_node)]
            except KeyError:
#                try: 
                second_arc_cost = r_dict[p[1], node, tuple(neighbors_without_new_start_node)]
#                except KeyError:
#                    break
#                    second_arc_cost = infinity
            cost = first_arc_cost + second_arc_cost
            r_minus_dict[p[0], node, tuple(neighbors_without_new_start_node)] = cost
        r_minus_series = pd.Series(r_minus_dict)
        lowest_cost_in_r_minus_key = r_minus_series.idxmin()
        lowest_cost_in_r_minus_cost = r_minus_series.min()
        r_dict[lowest_cost_in_r_minus_key] = lowest_cost_in_r_minus_cost
               

print("Nearest Neighbors: ", nearest_neighbors)

# Generate every p for each node u
for time in nodes:
    get_P_from_u_and_N_u(nearest_neighbors[time], nodes)
    
# Sort P by the amount of neighbors
P = sorted(P, key=lambda x: len(x[2]))
print("P: ", P)
print("Cardinality of P: ", len(P))
# Sort neighbors (N_p) in P for saving and accessing in dict
P = [(x[0], x[1], np.sort(x[2]).tolist()) for x in P]

print("P: ", P)
print("Cardinality of P: ", len(P))

for p in P:
    calculate_cost_for_p(p)


invalid_count = 0
valid_count = 0
for key, value in r_dict.items():
    if value > 10000:
        invalid_count += 1
    if value < 10000:
        valid_count += 1
        
#print("Invalid count: ", invalid_count)
#print("Valid count: ", valid_count)

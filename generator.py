from scipy.spatial.distance import pdist, squareform
import pandas as pd


# Setup config
amount_of_neighbors = 4
# List of LA-arcs P
P = []
R = []
r_dict = {}
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


# Generate P by adding start and end customers (including depot)
def generate_P_from_N_p(N_p: list, all_nodes: list):
    # get elements from all_nodes that are not in N_p
    unique_elements_all_nodes = [x for x in all_nodes if x not in N_p]
    # Take node with index 12 as a depot
    unique_elements_all_nodes_and_depot = unique_elements_all_nodes
    unique_elements_all_nodes_and_depot.append(12)
    # start_customer represents u_p and end_customer represents v_p
    # Iterate over all possible nodes and the depot excluding of N_p for possible start_customers
    for start_customer in unique_elements_all_nodes_and_depot:
        # Iterate over all possible nodes and the depot excluding of N_p for possible start_customers
        for end_customer in unique_elements_all_nodes_and_depot:
            # u_p and v_p cannot be the same node
            if start_customer != end_customer:
                # Create LA-arc p by adding start and end customers for the set P
                p = (start_customer, end_customer, N_p)
                P.append(p)
          

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
    for v_p in nodes_not_in_N_u_nor_u:
        for subset in power_set_of_N_u:
            P.append((tuple_with_u_and_N_u[0], v_p, subset))
            
    
def get_efficient_frontier(p: tuple) -> tuple:
    # Denoted as R_p^* in the literature
    efficient_frontier = []
    # neighbor denoted as w in literature, p[2] represents N_p
    for neighbor in p[2]:
        N_p_without_neighbor = p[2]
        N_p_without_neighbor.remove(neighbor)
        
        p_hat = (neighbor, p[1], N_p_without_neighbor)
               

print(nearest_neighbors)

# Generate every p for each node u
for time in nodes:
    get_P_from_u_and_N_u(nearest_neighbors[time], nodes)
    
# Sort P by the amount of neighbors
P = sorted(P, key=lambda x: len(x[2]))

print(P)
print("Cardinality of P: ", len(P))

print("Distance df: ", distance_df)

print(distance_df[0][1])
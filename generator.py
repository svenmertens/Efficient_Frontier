# Setup config
amount_of_neighbors = 4
# List of LA-arcs P
P = []
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
    
    nodes_not_in_N_u_nor_u = [x for x in all_nodes if x not in tuple_with_u_and_N_u[1]]
    nodes_not_in_N_u_nor_u.append(12)
    nodes_not_in_N_u_nor_u.remove(tuple_with_u_and_N_u[0])
    power_set_of_N_u = get_power_set(tuple_with_u_and_N_u[1])
#    power_set_of_N_u.remove([])
    for v_p in nodes_not_in_N_u_nor_u:
        for subset in power_set_of_N_u:
            P.append((tuple_with_u_and_N_u[0], v_p, subset))
               

# Generate the power set from each set of nearest neighbors (with one for each node)
# N_p represents the set of neighborhoods
#N_p = []
#for base_set in nearest_neighbors:
#    power_set = get_power_set(base_set)
#    N_p.append(power_set)
      
#N_p = get_unique_sets(N_p)
#print("N_p with unique sets: ")
#print(N_p)
#print("Amount of outer sets: ", len(N_p))
#print(nodes)

#print("\n P: \n")
#for neighborhood in N_p:
#    generate_P_from_N_p(neighborhood, nodes)
#print(P)
#print("Cardinality of P: ", len(P))

print(nearest_neighbors)

print("\nTest: \n")
for time in nodes:
    get_P_from_u_and_N_u(nearest_neighbors[time], nodes)
print(P)
print("Cardinality of P: ", len(P))

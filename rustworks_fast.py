import rustworkx as rx
from tqdm import tqdm
import argparse
import numpy as np
import os


# def remove_low_degree_nodes(digraph, degree_threshold):
#     nodes_to_remove = []
#     for node in digraph.node_indexes():
#         degree = digraph.in_degree(node) + digraph.degree(node)
#         if degree < degree_threshold:
#             nodes_to_remove.append(node)
            
#     digraph.remove_nodes_from(nodes_to_remove)


def remove_low_degree_nodes(graph, degree_threshold):
    nodes_to_remove = []
    for node in graph.node_indexes():
        degree = graph.degree(node)
        if degree < degree_threshold:
            nodes_to_remove.append(node)
            
    graph.remove_nodes_from(nodes_to_remove)


def k_core_decomposition(graph, k):
    core_nodes = []
    for node in graph.node_indexes():
        degree = graph.degree(node)
        if degree >= k:
            core_nodes.append(node)
    return core_nodes

parser = argparse.ArgumentParser(description='Pairwise graph clustering')
parser.add_argument('-c', '--cutoff', type=int, required=True,
                    help="clustering threshold (0:100)%%")
parser.add_argument('-p', '--edges', type=str,
                    required=True, help="pairwise TSV file")
args = parser.parse_args()

input_edges_file = args.edges

edges_file = input_edges_file

id_to_name_file = "gtdb_nodeID_to_acc.tsv"

ANI_THRESHOLD = float(args.cutoff)

output = f"{input_edges_file}_kSpider_clusters_{ANI_THRESHOLD}%.tsv"


# parse kmer counts
id_to_kmer_count = {}
with open("gtdb_genomic_k31_s1k_groupID_to_kmerCount.tsv") as F:
    for line in F:
        line = line.strip().split('\t')
        id_to_kmer_count[int(line[0]) - 1] = int(line[1])


# no_edges = 2726667056 #gtdb
# no_edges = 6928010548 #nasa_no_filtering
no_edges = 738_181_634 #5371116604 #nasa

# loading id_to_group_name
id_to_name = {}
with open(id_to_name_file) as F:
    next(F)
    for line in F:
        line = line.strip().split('\t')
        id_to_name[int(line[0]) - 1] = line[1]


distance_col_idx = -1 # avg ani
# distance_col_idx = -2 # avg cont


graph = rx.PyGraph()
nodes_indeces = graph.add_nodes_from(list(id_to_name.keys()))

batch_size = 10000000
batch_counter = 0
edges_tuples = []

# EXPERIMENTAL load with numpy
print("[i] loading edges")
# save if not saved
if not os.path.exists("all_edges.npy"):
    np_edges = np.loadtxt(edges_file, delimiter='\t', skiprows=1, usecols=(0, 1), dtype=int)
    print("saving numpy edges")
    np.save("all_edges.npy", np_edges)
else:
    np_edges = np.load("all_edges.npy")
print("Loading complete!")



print("[i] constructing graph")

for row in tqdm(np_edges):
    seq1, seq2 = row
    seq1 -= 1
    seq2 -= 1

    if batch_counter < batch_size:
        batch_counter += 1
        if id_to_kmer_count[seq1] < id_to_kmer_count[seq2]:
            edges_tuples.append((seq1, seq2))
        else:
            edges_tuples.append((seq2, seq1))
    else:
        graph.add_edges_from_no_data(edges_tuples)
        batch_counter = 0
        edges_tuples.clear()
if len(edges_tuples):
    graph.add_edges_from_no_data(edges_tuples)



# print("[i] constructing graph")
# with open(edges_file, 'r') as pairwise_tsv:
#     next(pairwise_tsv)  # skip header
#     for row in tqdm(pairwise_tsv, total=no_edges):
#         row = row.strip().split('\t')
#         seq1 = int(row[0]) - 1
#         seq2 = int(row[1]) - 1
#         # distance = float(row[distance_col_idx]) * 100

#         # don't make graph edge
#         # if distance < ANI_THRESHOLD:
#         #     continue

#         if batch_counter < batch_size:
#             batch_counter += 1
#             if id_to_kmer_count[seq1] < id_to_kmer_count[seq2]:
#                 edges_tuples.append((seq1, seq2))
#             else:
#                 edges_tuples.append((seq2, seq1))
#         else:
#             graph.add_edges_from_no_data(edges_tuples)
#             batch_counter = 0
#             edges_tuples.clear()

#     if len(edges_tuples):
#         graph.add_edges_from_no_data(edges_tuples)


remove_low_degree_nodes(graph, 5)


print("clustering...")
connected_components = rx.connected_components(graph)
print(f"number of clusters: {len(connected_components)}")
print("printing results")
single_components = 0
with open(output, 'w') as CLUSTERS:
    for component in connected_components:
        # uncomment to exclude single genome clusters from exporting
        # if len(component) == 1:
        #     single_components += 1
        #     continue
        named_component = [id_to_name[node] for node in component]
        CLUSTERS.write(','.join(named_component) + '\n')

# print(f"skipped clusters with single node: {single_components}")
import math
from random import random

from DirectMethod import Bender
from GraphProcessor import GraphProcessor

args = input("").split()
if len(args) < 2:
    raise Exception("Invalid input")

print("1. Process Input Parameters and Initialize")
directed = True
SMPLS = 100000
G_N = 3
fileName = args[0]

for i in range(1, len(args)):
    s = args[i]
    if s == "-s":
        i += 1
        G_N = int(args[i])
    elif s == "-t":
        i += 1
        SMPLS = int(args[i])
    elif s == "-d" or s == "-directed":
        directed = True
    elif s == "-ud" or s == "-undirected":
        directed = False

# Generate Graph whose node name starts with 0
myGP = GraphProcessor()
graph = myGP.loadGraph(fileName, directed=directed)
file = open("tempFiles\\newGraph.txt", "w")
s = set()
for u, v in graph.getEdgeList():
    if (u, v) not in s and (v, u) not in s:
        toWrite = str(u) + " " + str(v) + "\n"
        file.write(toWrite)
    s.add((u, v))
file.close()

print("2. Build Graph")
# read graph
inputGraph = myGP.loadGraph("tempFiles\\newGraph.txt", directed=False)

n = len(inputGraph.getVertexList().keys())
m = len(inputGraph.getEdgeList())
edges = inputGraph.getEdgeList()
num_neighbours = []

print(" - Building degree sequence for graph with ")
print(" {} vertices and {} connections".format(n, m))

# info for Direct Method
degree_r = [0 for i in range(n)]
degree_c = [0 for i in range(n)]
print(degree_r == degree_c)
for u, v in edges:
    if directed:
        degree_r[u] += 1
        degree_c[v] += 1
    else:
        degree_r[u] += 1
        degree_r[v] += 1
        degree_c[u] += 1
        degree_c[v] += 1

print("3. Sampling")
total_num_subgraphs_undir = [1, 1, 1, 2, 6, 21, 112, 853, 11117, 261080]

subgraph_library3_undir = [78, 238]
subgraph_library4_undir = [4382, 4958, 8598, 13278, 27030, 31710]
subgraph_library3_dir = [6, 12, 14, 36, 38, 46, 78, 102, 140, 164, 166, 174, 238]
subgraph_library4_dir = [14, 28, 30, 74, 76, 78, 90, 92, 94, 204, 206, 222, 280, 282, 286,
                         328, 330, 332, 334, 344, 346, 348, 350, 390, 392, 394, 396, 398,
                         404, 406, 408, 410, 412, 414, 454, 456, 458, 460, 462, 468, 470,
                         472, 474, 476, 478, 856, 858, 862, 904, 906, 908, 910, 922, 924,
                         926, 972, 974, 990, 2184, 2186, 2190, 2202, 2204, 2206, 2252, 2254,
                         2270, 2458, 2462, 2506, 2510, 2524, 2526, 3038, 4370, 4374, 4382,
                         4418, 4420, 4422, 4424, 4426, 4428, 4430, 4434, 4436, 4438, 4440,
                         4442, 4444, 4446, 4546, 4548, 4550, 4556, 4558, 4562, 4564, 4566,
                         4572, 4574, 4678, 4682, 4686, 4692, 4694, 4698, 4700, 4702, 4740,
                         4742, 4748, 4750, 4758, 4764, 4766, 4812, 4814, 4830, 4946, 4950,
                         4952, 4954, 4958, 4994, 4998, 5002, 5004, 5006, 5010, 5012, 5014,
                         5016, 5018, 5020, 5022, 5058, 5062, 5064, 5066, 5068, 5070, 5074,
                         5076, 5078, 5080, 5082, 5084, 5086, 6342, 6348, 6350, 6356, 6358,
                         6364, 6366, 6550, 6552, 6554, 6558, 6598, 6602, 6604, 6606, 6614,
                         6616, 6618, 6620, 6622, 6854, 6858, 6862, 6870, 6874, 6876, 6878,
                         7126, 7128, 7130, 7134, 13142, 13146, 13148, 13150, 13260, 13262,
                         13278, 14678, 14686, 14790, 14798, 14810, 14812, 14814, 15258,
                         15262, 15310, 15326, 31710]
num_samples = SMPLS
bender = Bender(n, degree_r, None, random)

target_arr_str = "subgraph_library" + str(G_N) + "_"
if directed:
    target_arr_str += "dir"
else:
    target_arr_str += "undir"

target_arr = eval(target_arr_str)

res = [0] * len(target_arr)
num_of_subgraphs = len(res)
for i in range(num_of_subgraphs):
    print("   Motif {} of {}".format(i + 1, num_of_subgraphs))
    res[i] = bender.motif_sample_log(target_arr[i], G_N, num_samples)

max_value = -50000000.0
# normalization the results
for i in range(num_of_subgraphs):
    if res[i] > max_value and not (res[i] != res[i]):
        max_value = res[i]
summation = 0.0
for i in range(num_of_subgraphs):
    if not (res[i] != res[i]):
        res[i] = math.exp(res[i] - max_value)
        summation += res[i]
print("4. Results")
for i in range(num_of_subgraphs):
    print("   Motif {} : ".format(target_arr[i]),
          end="")
    if not (res[i] != res[i]):
        print(res[i] / summation)
        # print(res[i])
    else:
        print("Not found")



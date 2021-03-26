# DirectMethodPy
This is a prototype of Direct Method implementation by Python.This program supports the Direct Method calculation for directed and undirected subgraphs in size 3 and 4.

By calling this program, the input parameter format is as follows:

**GRAPH_NAME -s SUBGRAPH_SIZE -t NUM_SAMPLE -GRAPH_TYPE**

GRAPH_NAME is the original graph applied this program.
SUBGRAPH_SIZE is the size of subgraph, either 3 or 4 available.
NUM_SAMPLE is the number of random iterations, at least 100,000 recommended.
GRAPH_TYPE means the graph is directed or undirected, -d or -ud.

Here is an example of the input format:
**Y2k.txt -s 4 -t 10000 -ud**

This means I am looking for the estimated concentration of subgraph in size 4 for undirected graph Y2k.txt.

For estimated Z-score, it will need to be compared with the concentration of subgraphs in the original graphs. Such techniques could be found on : https://bioresearch.css.uwb.edu/biores/nemo/

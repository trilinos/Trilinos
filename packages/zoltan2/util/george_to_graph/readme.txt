A tool to convert George Slota's graph files into Chaco/METIS format for XtraPuLP. By default, this tool appends up to 3 vertex weights to each vertices: the first is just integer 1, the second is the vertex index, the third is a random integer between 1 and 100. More weights can be added or changed through source code. 

To run by itself and to generate a graph file with up to 3 vertex weights, the command line format is:

george_to_graph.exe <name of graph file>.graph <# of vertex components>

Example: george_to_graph.exe kkt_power.graph 3

Running the above command will generate a new graph file with the format <original graph name>_<# of vertex weights>vw.graph.

Example: kkt_power_3vw.graph

To modify and build: Import the files george_to_graph.c, mmio.c, and mmio.h into a C project. Modify the function generate_weights in george_to_graph.c at the end of the file to modify weights. 
INFO Records
Inline mesh specification requested: 
	5 Elements 
	10 Nodes and
 	14 Edges.
Using BISECTION LAYOUT decomposition.
Number of elements/segment in directions I/X/R 		2
Number of elements/segment in directions J/Y/THETA 	2
Number of elements/segment in directions K/Z/PHI 	1
Number of mesh segments in directions I/X/R 		1
Number of mesh segments in directions J/Y/THETA 	1
Number of mesh segments in directions K/Z/PHI 	1
Using  BISECTION decomposition.
 No cuts will be made in radial direction.

Exodus header info:
Title: PAMGEN Inline Mesh
Dimension 2 
Number of Nodes 10 
Number of Elements 5 
Number of Element Blocks 2 
Number of Node Sets 1 
Number of Side Sets 2 

num node set nodes 0
num node set dfs 0
num side set elements 0
num side set nodes 8
num side set dfs 0
num block properties 0
num node set properties 0
num side set properties 0
A taste of coords
X 0.000000 Y 0.000000 
X 0.000707 Y 0.000000 
X 0.000000 Y 0.000707 
X 0.000604 Y 0.000604 
X 0.002000 Y 0.000000 
X 0.005000 Y 0.000000 
X 0.001414 Y 0.001414 
X 0.003536 Y 0.003536 
X 0.000000 Y 0.002000 
X 0.000000 Y 0.005000 
coord name 0 X
coord name 1 Y
A tast of map
map i=0, val=1
map i=1, val=2
map i=2, val=3
map i=3, val=4
map i=4, val=5
A tast of global elem numbers
global el num  i=0, val=1
global el num  i=1, val=2
global el num  i=2, val=3
global el num  i=3, val=4
global el num  i=4, val=5
A tast of global elem numbers
global node num  i=0, val=1
global node num  i=1, val=2
global node num  i=2, val=3
global node num  i=3, val=4
global node num  i=4, val=5
block i = 0 has id 1 
block i = 0
block id 1
element_type QUAD
num elements 3
nodes per element 4
element attributes 0
block i = 1 has id 2 
block i = 1
block id 2
element_type QUAD
num elements 2
nodes per element 4
element attributes 0
block 1 element 0 connectivty 1 2 4 3 
block 1 element 1 connectivty 2 5 7 4 
block 1 element 2 connectivty 4 7 9 3 
block 2 element 0 connectivty 5 6 8 7 
block 2 element 1 connectivty 7 8 10 9 
Nodeset i = 0 id = 1 has 1 nodes
nodeset node i=0 = 1
Side set index 0 id 2 has 2 elements
element 4 and face 2
element 5 and face 2
Side set index 1 id 3 has 2 elements
element 2 and face 2
element 3 and face 2
num qa records 1

QA Record 0
 PAMGEN
PArallel Mesh GENerator
Num Info Records 0
Nemesis data
Num nodes global 10
Num elems global 5
Num elm_blks global 2
Num node sets global 1
Num side sets global 2
Num total proc 1
Num proc in file 1
element block index 0 has id 1 and 3 elements
element block index 1 has id 2 and 2 elements
global ns info for ns index 0 id 1 num_nodes = 0 num_ns_df = 0
global ss info for ss index 0 id 2 num_elements = 2 num_ss_df = 0
global ss info for ss index 1 id 3 num_elements = 2 num_ss_df = 0
Loadbal params:
num_internal_nodes 10
num_border_nodes0
num_external_nodes0
num_internal_elems5
num_border_elems0
num_node_comm_maps0
num_elem_comm_maps0
internal node i=0 = 1
internal node i=1 = 2
internal node i=2 = 3
internal node i=3 = 4
internal node i=4 = 5
internal node i=5 = 6
internal node i=6 = 7
internal node i=7 = 8
internal node i=8 = 9
internal node i=9 = 10
internal elem i=0 = 1
internal elem i=1 = 2
internal elem i=2 = 3
internal elem i=3 = 4
internal elem i=4 = 5

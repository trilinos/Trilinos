INFO Records
Inline mesh specification requested: 
	448 Elements 
	465 Nodes and
 	912 Edges.
Using BISECTION LAYOUT decomposition.
Number of elements/segment in directions I/X/R 		12
Number of elements/segment in directions J/Y/THETA 	16
Number of elements/segment in directions K/Z/PHI 	1
Number of mesh segments in directions I/X/R 		1
Number of mesh segments in directions J/Y/THETA 	2
Number of mesh segments in directions K/Z/PHI 	1
Using  BISECTION decomposition.
 No cuts will be made in radial direction.

Exodus header info:
Title: PAMGEN Inline Mesh
Dimension 2 
Number of Nodes 249 
Number of Elements 224 
Number of Element Blocks 3 
Number of Node Sets 1 
Number of Side Sets 1 

num node set nodes 0
num node set dfs 0
num side set elements 0
num side set nodes 32
num side set dfs 0
num block properties 0
num node set properties 0
num side set properties 0
A taste of coords
X 0.000000 Y 0.000000 
X 0.282843 Y 0.000000 
X 0.565685 Y 0.000000 
X 0.848528 Y 0.000000 
X 1.131371 Y 0.000000 
X -0.282843 Y 0.000000 
X -0.565685 Y 0.000000 
X -0.848528 Y 0.000000 
X -1.131371 Y 0.000000 
X -0.000000 Y -0.282843 
coord name 0 X
coord name 1 Y
A tast of map
map i=0, val=33
map i=1, val=34
map i=2, val=35
map i=3, val=36
map i=4, val=37
map i=5, val=38
map i=6, val=39
map i=7, val=40
map i=8, val=41
map i=9, val=42
A tast of global elem numbers
global el num  i=0, val=33
global el num  i=1, val=34
global el num  i=2, val=35
global el num  i=3, val=36
global el num  i=4, val=37
global el num  i=5, val=38
global el num  i=6, val=39
global el num  i=7, val=40
global el num  i=8, val=41
global el num  i=9, val=42
A tast of global elem numbers
global node num  i=0, val=1
global node num  i=1, val=2
global node num  i=2, val=3
global node num  i=3, val=4
global node num  i=4, val=5
global node num  i=5, val=26
global node num  i=6, val=31
global node num  i=7, val=36
global node num  i=8, val=41
global node num  i=9, val=46
block i = 0 has id 1 
block i = 0
block id 1
element_type QUAD
num elements 96
nodes per element 4
element attributes 0
block i = 1 has id 2 
block i = 1
block id 2
element_type QUAD
num elements 64
nodes per element 4
element attributes 0
block i = 2 has id 3 
block i = 2
block id 3
element_type QUAD
num elements 64
nodes per element 4
element attributes 0
block 1 element 0 connectivty 1 6 11 10 
block 1 element 1 connectivty 6 7 12 11 
block 1 element 2 connectivty 7 8 13 12 
block 1 element 3 connectivty 8 9 14 13 
block 1 element 4 connectivty 10 11 16 15 
block 1 element 5 connectivty 11 12 17 16 
block 1 element 6 connectivty 12 13 18 17 
block 1 element 7 connectivty 13 14 19 18 
block 1 element 8 connectivty 15 16 21 20 
block 1 element 9 connectivty 16 17 22 21 
block 2 element 0 connectivty 61 62 74 73 
block 2 element 1 connectivty 62 63 75 74 
block 2 element 2 connectivty 63 64 76 75 
block 2 element 3 connectivty 64 65 77 76 
block 2 element 4 connectivty 73 74 86 85 
block 2 element 5 connectivty 74 75 87 86 
block 2 element 6 connectivty 75 76 88 87 
block 2 element 7 connectivty 76 77 89 88 
block 2 element 8 connectivty 85 86 98 97 
block 2 element 9 connectivty 86 87 99 98 
block 3 element 0 connectivty 65 66 78 77 
block 3 element 1 connectivty 66 67 79 78 
block 3 element 2 connectivty 67 68 80 79 
block 3 element 3 connectivty 68 69 81 80 
block 3 element 4 connectivty 77 78 90 89 
block 3 element 5 connectivty 78 79 91 90 
block 3 element 6 connectivty 79 80 92 91 
block 3 element 7 connectivty 80 81 93 92 
block 3 element 8 connectivty 89 90 102 101 
block 3 element 9 connectivty 90 91 103 102 
Nodeset i = 0 id = 40 has 17 nodes
nodeset node i=0 = 57
nodeset node i=1 = 69
nodeset node i=2 = 81
nodeset node i=3 = 93
nodeset node i=4 = 105
nodeset node i=5 = 117
nodeset node i=6 = 129
nodeset node i=7 = 141
nodeset node i=8 = 153
nodeset node i=9 = 165
Side set index 0 id 45 has 16 elements
element 100 and face 2
element 104 and face 2
element 108 and face 2
element 112 and face 2
element 116 and face 2
element 120 and face 2
element 124 and face 2
element 128 and face 2
element 132 and face 2
element 136 and face 2
num qa records 1

QA Record 0
 PAMGEN
PArallel Mesh GENerator
Num Info Records 0
Nemesis data
Num nodes global 429
Num elems global 448
Num elm_blks global 3
Num node sets global 1
Num side sets global 1
Num total proc 2
Num proc in file 1
element block index 0 has id 1 and 192 elements
element block index 1 has id 2 and 128 elements
element block index 2 has id 3 and 128 elements
global ns info for ns index 0 id 40 num_nodes = 33 num_ns_df = 0
global ss info for ss index 0 id 45 num_elements = 32 num_ss_df = 0
Loadbal params:
num_internal_nodes 216
num_border_nodes33
num_external_nodes0
num_internal_elems192
num_border_elems32
num_node_comm_maps1
num_elem_comm_maps1
internal node i=0 = 10
internal node i=1 = 11
internal node i=2 = 12
internal node i=3 = 13
internal node i=4 = 14
internal node i=5 = 15
internal node i=6 = 16
internal node i=7 = 17
internal node i=8 = 18
internal node i=9 = 19
border node i=0 = 1
border node i=1 = 2
border node i=2 = 3
border node i=3 = 4
border node i=4 = 5
border node i=5 = 6
border node i=6 = 7
border node i=7 = 8
border node i=8 = 9
border node i=9 = 46
internal elem i=0 = 5
internal elem i=1 = 6
internal elem i=2 = 7
internal elem i=3 = 8
internal elem i=4 = 9
internal elem i=5 = 10
internal elem i=6 = 11
internal elem i=7 = 12
internal elem i=8 = 13
internal elem i=9 = 14
border elem i=0 = 1
border elem i=1 = 2
border elem i=2 = 3
border elem i=3 = 4
border elem i=4 = 17
border elem i=5 = 21
border elem i=6 = 25
border elem i=7 = 29
border elem i=8 = 33
border elem i=9 = 34
node_cmap_id i = 0 node_cmap_id = 1 node_cmap_node_cnts = 33
elem_cmap_id i = 0 elem_cmap_id = 1 elem_cmap_elem_cnts = 32
node_cmap_id i=0 = 1 comm_node_ids = 1 comm_node_proc_ids = 1
node_cmap_id i=1 = 1 comm_node_ids = 2 comm_node_proc_ids = 1
node_cmap_id i=2 = 1 comm_node_ids = 3 comm_node_proc_ids = 1
node_cmap_id i=3 = 1 comm_node_ids = 4 comm_node_proc_ids = 1
node_cmap_id i=4 = 1 comm_node_ids = 5 comm_node_proc_ids = 1
node_cmap_id i=5 = 1 comm_node_ids = 6 comm_node_proc_ids = 1
node_cmap_id i=6 = 1 comm_node_ids = 7 comm_node_proc_ids = 1
node_cmap_id i=7 = 1 comm_node_ids = 8 comm_node_proc_ids = 1
node_cmap_id i=8 = 1 comm_node_ids = 9 comm_node_proc_ids = 1
node_cmap_id i=9 = 1 comm_node_ids = 46 comm_node_proc_ids = 1
elem_cmap_id i=0 = 1 comm_elem_ids = 1 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=1 = 1 comm_elem_ids = 2 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=2 = 1 comm_elem_ids = 3 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=3 = 1 comm_elem_ids = 4 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=4 = 1 comm_elem_ids = 17 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=5 = 1 comm_elem_ids = 21 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=6 = 1 comm_elem_ids = 25 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=7 = 1 comm_elem_ids = 29 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=8 = 1 comm_elem_ids = 33 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=9 = 1 comm_elem_ids = 34 comm_side_ids = 1 comm_elem_proc_ids = 1

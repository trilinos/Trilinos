INFO Records
Inline mesh specification requested: 
	448 Elements 
	465 Nodes and
 	912 Edges.
Using BISECTION LAYOUT decomposition.
Number of elements/segment in directions I/X/R 		12
Number of elements/segment in directions J/Y/THETA 	8
Number of elements/segment in directions K/Z/PHI 	1
Number of mesh segments in directions I/X/R 		1
Number of mesh segments in directions J/Y/THETA 	4
Number of mesh segments in directions K/Z/PHI 	1
Using  BISECTION decomposition.
 No cuts will be made in radial direction.

Exodus header info:
Title: PAMGEN Inline Mesh
Dimension 2 
Number of Nodes 133 
Number of Elements 112 
Number of Element Blocks 3 
Number of Node Sets 1 
Number of Side Sets 1 

num node set nodes 0
num node set dfs 0
num side set elements 0
num side set nodes 16
num side set dfs 0
num block properties 0
num node set properties 0
num side set properties 0
A taste of coords
X 0.000000 Y 0.000000 
X -0.282843 Y 0.000000 
X -0.565685 Y 0.000000 
X -0.848528 Y 0.000000 
X -1.131371 Y 0.000000 
X -0.000000 Y -0.282843 
X -0.266274 Y -0.266274 
X -0.532548 Y -0.249706 
X -0.798823 Y -0.233137 
X -1.065097 Y -0.216569 
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
global node num  i=1, val=26
global node num  i=2, val=31
global node num  i=3, val=36
global node num  i=4, val=41
global node num  i=5, val=46
global node num  i=6, val=47
global node num  i=7, val=48
global node num  i=8, val=49
global node num  i=9, val=50
block i = 0 has id 1 
block i = 0
block id 1
element_type QUAD
num elements 48
nodes per element 4
element attributes 0
block i = 1 has id 2 
block i = 1
block id 2
element_type QUAD
num elements 32
nodes per element 4
element attributes 0
block i = 2 has id 3 
block i = 2
block id 3
element_type QUAD
num elements 32
nodes per element 4
element attributes 0
block 1 element 0 connectivty 1 2 7 6 
block 1 element 1 connectivty 2 3 8 7 
block 1 element 2 connectivty 3 4 9 8 
block 1 element 3 connectivty 4 5 10 9 
block 1 element 4 connectivty 6 7 12 11 
block 1 element 5 connectivty 7 8 13 12 
block 1 element 6 connectivty 8 9 14 13 
block 1 element 7 connectivty 9 10 15 14 
block 1 element 8 connectivty 11 12 17 16 
block 1 element 9 connectivty 12 13 18 17 
block 2 element 0 connectivty 29 30 42 41 
block 2 element 1 connectivty 30 31 43 42 
block 2 element 2 connectivty 31 32 44 43 
block 2 element 3 connectivty 32 33 45 44 
block 2 element 4 connectivty 41 42 54 53 
block 2 element 5 connectivty 42 43 55 54 
block 2 element 6 connectivty 43 44 56 55 
block 2 element 7 connectivty 44 45 57 56 
block 2 element 8 connectivty 53 54 66 65 
block 2 element 9 connectivty 54 55 67 66 
block 3 element 0 connectivty 33 34 46 45 
block 3 element 1 connectivty 34 35 47 46 
block 3 element 2 connectivty 35 36 48 47 
block 3 element 3 connectivty 36 37 49 48 
block 3 element 4 connectivty 45 46 58 57 
block 3 element 5 connectivty 46 47 59 58 
block 3 element 6 connectivty 47 48 60 59 
block 3 element 7 connectivty 48 49 61 60 
block 3 element 8 connectivty 57 58 70 69 
block 3 element 9 connectivty 58 59 71 70 
Nodeset i = 0 id = 40 has 9 nodes
nodeset node i=0 = 37
nodeset node i=1 = 49
nodeset node i=2 = 61
nodeset node i=3 = 73
nodeset node i=4 = 85
nodeset node i=5 = 97
nodeset node i=6 = 109
nodeset node i=7 = 121
nodeset node i=8 = 133
Side set index 0 id 45 has 8 elements
element 52 and face 2
element 56 and face 2
element 60 and face 2
element 64 and face 2
element 68 and face 2
element 72 and face 2
element 76 and face 2
element 80 and face 2
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
Num total proc 4
Num proc in file 1
element block index 0 has id 1 and 192 elements
element block index 1 has id 2 and 128 elements
element block index 2 has id 3 and 128 elements
global ns info for ns index 0 id 40 num_nodes = 33 num_ns_df = 0
global ss info for ss index 0 id 45 num_elements = 32 num_ss_df = 0
Loadbal params:
num_internal_nodes 100
num_border_nodes33
num_external_nodes0
num_internal_elems81
num_border_elems31
num_node_comm_maps3
num_elem_comm_maps2
internal node i=0 = 7
internal node i=1 = 8
internal node i=2 = 9
internal node i=3 = 10
internal node i=4 = 12
internal node i=5 = 13
internal node i=6 = 14
internal node i=7 = 15
internal node i=8 = 17
internal node i=9 = 18
border node i=0 = 1
border node i=1 = 2
border node i=2 = 3
border node i=3 = 4
border node i=4 = 5
border node i=5 = 6
border node i=6 = 11
border node i=7 = 16
border node i=8 = 21
border node i=9 = 26
internal elem i=0 = 6
internal elem i=1 = 7
internal elem i=2 = 8
internal elem i=3 = 10
internal elem i=4 = 11
internal elem i=5 = 12
internal elem i=6 = 14
internal elem i=7 = 15
internal elem i=8 = 16
internal elem i=9 = 21
border elem i=0 = 1
border elem i=1 = 2
border elem i=2 = 3
border elem i=3 = 4
border elem i=4 = 5
border elem i=5 = 9
border elem i=6 = 13
border elem i=7 = 17
border elem i=8 = 18
border elem i=9 = 19
node_cmap_id i = 0 node_cmap_id = 0 node_cmap_node_cnts = 1
node_cmap_id i = 1 node_cmap_id = 1 node_cmap_node_cnts = 17
node_cmap_id i = 2 node_cmap_id = 2 node_cmap_node_cnts = 17
elem_cmap_id i = 0 elem_cmap_id = 1 elem_cmap_elem_cnts = 16
elem_cmap_id i = 1 elem_cmap_id = 2 elem_cmap_elem_cnts = 16
node_cmap_id i=0 = 0 comm_node_ids = 1 comm_node_proc_ids = 0
node_cmap_id i=0 = 1 comm_node_ids = 1 comm_node_proc_ids = 1
node_cmap_id i=1 = 1 comm_node_ids = 2 comm_node_proc_ids = 1
node_cmap_id i=2 = 1 comm_node_ids = 3 comm_node_proc_ids = 1
node_cmap_id i=3 = 1 comm_node_ids = 4 comm_node_proc_ids = 1
node_cmap_id i=4 = 1 comm_node_ids = 5 comm_node_proc_ids = 1
node_cmap_id i=5 = 1 comm_node_ids = 26 comm_node_proc_ids = 1
node_cmap_id i=6 = 1 comm_node_ids = 27 comm_node_proc_ids = 1
node_cmap_id i=7 = 1 comm_node_ids = 28 comm_node_proc_ids = 1
node_cmap_id i=8 = 1 comm_node_ids = 29 comm_node_proc_ids = 1
node_cmap_id i=9 = 1 comm_node_ids = 30 comm_node_proc_ids = 1
node_cmap_id i=0 = 2 comm_node_ids = 1 comm_node_proc_ids = 2
node_cmap_id i=1 = 2 comm_node_ids = 6 comm_node_proc_ids = 2
node_cmap_id i=2 = 2 comm_node_ids = 11 comm_node_proc_ids = 2
node_cmap_id i=3 = 2 comm_node_ids = 16 comm_node_proc_ids = 2
node_cmap_id i=4 = 2 comm_node_ids = 21 comm_node_proc_ids = 2
node_cmap_id i=5 = 2 comm_node_ids = 122 comm_node_proc_ids = 2
node_cmap_id i=6 = 2 comm_node_ids = 123 comm_node_proc_ids = 2
node_cmap_id i=7 = 2 comm_node_ids = 124 comm_node_proc_ids = 2
node_cmap_id i=8 = 2 comm_node_ids = 125 comm_node_proc_ids = 2
node_cmap_id i=9 = 2 comm_node_ids = 126 comm_node_proc_ids = 2
elem_cmap_id i=0 = 1 comm_elem_ids = 1 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=1 = 1 comm_elem_ids = 2 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=2 = 1 comm_elem_ids = 3 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=3 = 1 comm_elem_ids = 4 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=4 = 1 comm_elem_ids = 17 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=5 = 1 comm_elem_ids = 18 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=6 = 1 comm_elem_ids = 19 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=7 = 1 comm_elem_ids = 20 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=8 = 1 comm_elem_ids = 49 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=9 = 1 comm_elem_ids = 50 comm_side_ids = 1 comm_elem_proc_ids = 1
elem_cmap_id i=0 = 2 comm_elem_ids = 1 comm_side_ids = 4 comm_elem_proc_ids = 2
elem_cmap_id i=1 = 2 comm_elem_ids = 5 comm_side_ids = 4 comm_elem_proc_ids = 2
elem_cmap_id i=2 = 2 comm_elem_ids = 9 comm_side_ids = 4 comm_elem_proc_ids = 2
elem_cmap_id i=3 = 2 comm_elem_ids = 13 comm_side_ids = 4 comm_elem_proc_ids = 2
elem_cmap_id i=4 = 2 comm_elem_ids = 45 comm_side_ids = 3 comm_elem_proc_ids = 2
elem_cmap_id i=5 = 2 comm_elem_ids = 46 comm_side_ids = 3 comm_elem_proc_ids = 2
elem_cmap_id i=6 = 2 comm_elem_ids = 47 comm_side_ids = 3 comm_elem_proc_ids = 2
elem_cmap_id i=7 = 2 comm_elem_ids = 48 comm_side_ids = 3 comm_elem_proc_ids = 2
elem_cmap_id i=8 = 2 comm_elem_ids = 77 comm_side_ids = 3 comm_elem_proc_ids = 2
elem_cmap_id i=9 = 2 comm_elem_ids = 78 comm_side_ids = 3 comm_elem_proc_ids = 2

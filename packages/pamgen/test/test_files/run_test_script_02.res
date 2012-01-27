INFO Records
Inline mesh specification requested: 
	32000 Elements 
	35640 Nodes and
 	103240 Edges.
Using BISECTION LAYOUT decomposition.
Number of elements/segment in directions I/X/R 		26
Number of elements/segment in directions J/Y/THETA 	40
Number of elements/segment in directions K/Z/PHI 	10
Number of mesh segments in directions I/X/R 		3
Number of mesh segments in directions J/Y/THETA 	1
Number of mesh segments in directions K/Z/PHI 	1

Exodus header info:
Title: PAMGEN Inline Mesh
Dimension 3 
Number of Nodes 12320 
Number of Elements 10800 
Number of Element Blocks 1 
Number of Node Sets 0 
Number of Side Sets 0 

num node set nodes 0
num node set dfs 0
num side set elements 0
num side set nodes 0
num side set dfs 0
num block properties 0
num node set properties 0
num side set properties 0
A taste of coords
X 0.002071 Y 0.000000 Z -0.000750 
X 0.002109 Y 0.000000 Z -0.000750 
X 0.002147 Y 0.000000 Z -0.000750 
X 0.002185 Y 0.000000 Z -0.000750 
X 0.002223 Y 0.000000 Z -0.000750 
X 0.002261 Y 0.000000 Z -0.000750 
X 0.002299 Y 0.000000 Z -0.000750 
X 0.002338 Y 0.000000 Z -0.000750 
X 0.002376 Y 0.000000 Z -0.000750 
X 0.002414 Y 0.000000 Z -0.000750 
coord name 0 X
coord name 1 Y
coord name 2 Z
A tast of map
map i=0, val=54
map i=1, val=55
map i=2, val=56
map i=3, val=57
map i=4, val=58
map i=5, val=59
map i=6, val=60
map i=7, val=61
map i=8, val=62
map i=9, val=63
A tast of global elem numbers
global el num  i=0, val=54
global el num  i=1, val=55
global el num  i=2, val=56
global el num  i=3, val=57
global el num  i=4, val=58
global el num  i=5, val=59
global el num  i=6, val=60
global el num  i=7, val=61
global el num  i=8, val=62
global el num  i=9, val=63
A tast of global elem numbers
global node num  i=0, val=54
global node num  i=1, val=55
global node num  i=2, val=56
global node num  i=3, val=57
global node num  i=4, val=58
global node num  i=5, val=59
global node num  i=6, val=60
global node num  i=7, val=61
global node num  i=8, val=62
global node num  i=9, val=63
block i = 0 has id 1 
block i = 0
block id 1
element_type HEX
num elements 10800
nodes per element 8
element attributes 0
block 1 element 0 connectivty 1 2 30 29 1121 1122 1150 1149 
block 1 element 1 connectivty 2 3 31 30 1122 1123 1151 1150 
block 1 element 2 connectivty 3 4 32 31 1123 1124 1152 1151 
block 1 element 3 connectivty 4 5 33 32 1124 1125 1153 1152 
block 1 element 4 connectivty 5 6 34 33 1125 1126 1154 1153 
block 1 element 5 connectivty 6 7 35 34 1126 1127 1155 1154 
block 1 element 6 connectivty 7 8 36 35 1127 1128 1156 1155 
block 1 element 7 connectivty 8 9 37 36 1128 1129 1157 1156 
block 1 element 8 connectivty 9 10 38 37 1129 1130 1158 1157 
block 1 element 9 connectivty 10 11 39 38 1130 1131 1159 1158 
num qa records 1

QA Record 0
 PAMGEN
PArallel Mesh GENerator
Num Info Records 0
Nemesis data
Num nodes global 35640
Num elems global 32000
Num elm_blks global 1
Num node sets global 0
Num side sets global 0
Num total proc 3
Num proc in file 1
element block index 0 has id 1 and 32000 elements
Loadbal params:
num_internal_nodes 11880
num_border_nodes440
num_external_nodes0
num_internal_elems10400
num_border_elems400
num_node_comm_maps1
num_elem_comm_maps1
internal node i=0 = 2
internal node i=1 = 3
internal node i=2 = 4
internal node i=3 = 5
internal node i=4 = 6
internal node i=5 = 7
internal node i=6 = 8
internal node i=7 = 9
internal node i=8 = 10
internal node i=9 = 11
border node i=0 = 1
border node i=1 = 29
border node i=2 = 57
border node i=3 = 85
border node i=4 = 113
border node i=5 = 141
border node i=6 = 169
border node i=7 = 197
border node i=8 = 225
border node i=9 = 253
internal elem i=0 = 2
internal elem i=1 = 3
internal elem i=2 = 4
internal elem i=3 = 5
internal elem i=4 = 6
internal elem i=5 = 7
internal elem i=6 = 8
internal elem i=7 = 9
internal elem i=8 = 10
internal elem i=9 = 11
border elem i=0 = 1
border elem i=1 = 28
border elem i=2 = 55
border elem i=3 = 82
border elem i=4 = 109
border elem i=5 = 136
border elem i=6 = 163
border elem i=7 = 190
border elem i=8 = 217
border elem i=9 = 244
node_cmap_id i = 0 node_cmap_id = 1 node_cmap_node_cnts = 440
elem_cmap_id i = 0 elem_cmap_id = 1 elem_cmap_elem_cnts = 400
node_cmap_id i=0 = 1 comm_node_ids = 1 comm_node_proc_ids = 1
node_cmap_id i=1 = 1 comm_node_ids = 29 comm_node_proc_ids = 1
node_cmap_id i=2 = 1 comm_node_ids = 57 comm_node_proc_ids = 1
node_cmap_id i=3 = 1 comm_node_ids = 85 comm_node_proc_ids = 1
node_cmap_id i=4 = 1 comm_node_ids = 113 comm_node_proc_ids = 1
node_cmap_id i=5 = 1 comm_node_ids = 141 comm_node_proc_ids = 1
node_cmap_id i=6 = 1 comm_node_ids = 169 comm_node_proc_ids = 1
node_cmap_id i=7 = 1 comm_node_ids = 197 comm_node_proc_ids = 1
node_cmap_id i=8 = 1 comm_node_ids = 225 comm_node_proc_ids = 1
node_cmap_id i=9 = 1 comm_node_ids = 253 comm_node_proc_ids = 1
elem_cmap_id i=0 = 1 comm_elem_ids = 1 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=1 = 1 comm_elem_ids = 28 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=2 = 1 comm_elem_ids = 55 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=3 = 1 comm_elem_ids = 82 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=4 = 1 comm_elem_ids = 109 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=5 = 1 comm_elem_ids = 136 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=6 = 1 comm_elem_ids = 163 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=7 = 1 comm_elem_ids = 190 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=8 = 1 comm_elem_ids = 217 comm_side_ids = 4 comm_elem_proc_ids = 1
elem_cmap_id i=9 = 1 comm_elem_ids = 244 comm_side_ids = 4 comm_elem_proc_ids = 1

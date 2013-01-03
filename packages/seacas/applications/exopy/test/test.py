#!/opt/local/bin/python

import exopy

path = 'Cyl_Hex.exo'
#mode = int(0)
#comp_ws = int(0)
#io_ws = int(0)
(exoid,comp_ws,io_ws,version) = exopy.ex_open(path) #exopy.ex_open(path,mode,comp_ws,io_ws)
print "exoid = ", exoid
print "comp_ws =", comp_ws
print "io_ws =", io_ws
print "version =", version

(title,num_dim,num_node,num_elem,num_elem_blk,num_node_sets,num_side_sets) = exopy.ex_get_init(exoid)
print "title = ", title
print "num_dim = ", num_dim
print "num_node = ", num_node
print "num_elem = ", num_elem
print "num_elem_blk = ", num_elem_blk
print "num_node_sets = ", num_node_sets
print "num_side_sets = ", num_side_sets

error = exopy.ex_close(exoid)
print "ex_close() error = ", error


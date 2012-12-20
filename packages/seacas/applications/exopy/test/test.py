#!/opt/local/bin/python

import exopy

path = 'Cyl_Hex.exo'
mode = int(0)

(exoid,comp_ws,io_ws,version) = exopy.ex_open_test(path,mode)

print "exoid = ", exoid
print "mode = ", mode
print "comp_ws =", comp_ws
print "io_ws =", io_ws
print "version =", version

title = exopy.ref()
num_dim = exopy.ref()
num_node = exopy.ref()
num_elem = exopy.ref()
num_elem_blk = exopy.ref()
num_node_sets = exopy.ref()
num_side_sets = exopy.ref()
error = exopy.ex_get_init(exoid,title,num_dim,num_node,num_elem,num_elem_blk,num_node_sets,num_side_sets)
print "error = ", error
print "title = ", title.value
print "num_dim = ", num_dim.value
print "num_node = ", num_node.value
print "num_elem = ", num_elem.value
print "num_elem_blk = ", num_elem_blk.value
print "num_node_sets = ", num_node_sets.value
print "num_side_sets = ", num_side_sets.value

error = exopy.ex_close(exoid)

print "ex_close() error = ", error


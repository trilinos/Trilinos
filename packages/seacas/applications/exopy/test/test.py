#!/opt/local/bin/python

import exopy

path = 'Cyl_Hex.exo'
mode = int(0)
comp_ws = exopy.ref(0)   ## int
io_ws = exopy.ref(0)     ## int
version = exopy.ref(0.0) ## float or double

exoid = exopy.ex_open(path,mode,comp_ws,io_ws,version)

print "exoid = ", exoid
print "mode = ", mode
print "comp_ws =", comp_ws.value
print "io_ws =", io_ws.value
print "version =", version.value

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


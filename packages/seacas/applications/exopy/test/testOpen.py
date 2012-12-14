#!/opt/local/bin/python

import exopy

path = 'Cyl_Hex.exo'
mode = int(0)
comp_ws = exopy.ref(0)
io_ws = exopy.ref(0)
version = exopy.ref(0.0)

exoid = exopy.ex_open(path,mode,comp_ws,io_ws,version)

print "exoid = ", exoid
print "mode = ", mode
print "comp_ws =", comp_ws.value
print "io_ws =", io_ws.value
print "version =", version.value




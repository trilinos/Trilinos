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

error = exopy.ex_close(exoid)

print "ex_close() error = ", error


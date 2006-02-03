#! /usr/bin/env python

import setpath
import Epetra   # This will initialize MPI if we are running the MPI test harness
import Thyra

vs = Thyra.SerialVectorSpaceStd_double(4)
print "vs.dim() =", vs.dim()
print "dir(vs) =", dir(vs)
v = Thyra.createMember_double(vs)
print "type(v) =", type(v)
print "v.this =", v.this
print "dir(v) =", dir(v)
print "v.space().dim() =", v.space().dim()
v2 = Thyra.createMember_double(v.space())
print "v2.space().dim() =", v2.space().dim()
Thyra.V_S_double(v.__deref__(),1.0)
#Thyra.V_S_double(v,1.0)
print "sum_double(v) =", Thyra.sum_double(v.__deref__())
#print "sum_double(v) =", Thyra.sum_double(v)

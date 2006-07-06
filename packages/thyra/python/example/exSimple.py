#! /usr/bin/env python

from optparse import *

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
    import setpath
    import Epetra, Thyra
else:
    try:
        import setpath
        import Epetra, Thyra
    except ImportError:
        from PyTrilinos import Epetra, Thyra
        print >>sys.stderr, "Using system-installed Epetra, Thyra"

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

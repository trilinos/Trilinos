#! /usr/bin/env python

import setpath
import Thyra

vs = Thyra.SerialVectorSpaceStdDouble(4)
print vs.dim()
print dir(vs)

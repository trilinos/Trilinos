# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Chan

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


import LAPACK
import Abstract
import Continuation
import NOX.Abstract
import NOX.StatusTest
import MultiContinuation
import Homotopy
import TimeDependent
import Bifurcation
import NOX.LAPACK
class ProblemInterface(LAPACK.Interface):
    __swig_setmethods__ = {}
    for _s in [LAPACK.Interface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ProblemInterface, name, value)
    __swig_getmethods__ = {}
    for _s in [LAPACK.Interface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ProblemInterface, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ChanProblemInterface instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ProblemInterface, 'this', _Chan.new_ProblemInterface(*args))
        _swig_setattr(self, ProblemInterface, 'thisown', 1)
    def __del__(self, destroy=_Chan.delete_ProblemInterface):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getInitialGuess(*args): return _Chan.ProblemInterface_getInitialGuess(*args)
    def computeF(*args): return _Chan.ProblemInterface_computeF(*args)
    def computeJacobian(*args): return _Chan.ProblemInterface_computeJacobian(*args)
    def setParams(*args): return _Chan.ProblemInterface_setParams(*args)
    def printSolution(*args): return _Chan.ProblemInterface_printSolution(*args)

class ProblemInterfacePtr(ProblemInterface):
    def __init__(self, this):
        _swig_setattr(self, ProblemInterface, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ProblemInterface, 'thisown', 0)
        _swig_setattr(self, ProblemInterface,self.__class__,ProblemInterface)
_Chan.ProblemInterface_swigregister(ProblemInterfacePtr)



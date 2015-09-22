# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _LAPACK

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


import Abstract
import Continuation
import NOX.Abstract
import NOX.StatusTest
import MultiContinuation
import Homotopy
import TimeDependent
import Bifurcation
import NOX.LAPACK
class Interface(NOX.LAPACK.Interface):
    __swig_setmethods__ = {}
    for _s in [NOX.LAPACK.Interface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Interface, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.LAPACK.Interface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Interface, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::LAPACK::Interface instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_LAPACK.delete_Interface):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setParams(*args): return _LAPACK.Interface_setParams(*args)
    def printSolution(*args): return _LAPACK.Interface_printSolution(*args)
    def computeMass(*args): return _LAPACK.Interface_computeMass(*args)
    def projectToDraw(*args): return _LAPACK.Interface_projectToDraw(*args)
    def projectToDrawDimension(*args): return _LAPACK.Interface_projectToDrawDimension(*args)

class InterfacePtr(Interface):
    def __init__(self, this):
        _swig_setattr(self, Interface, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Interface, 'thisown', 0)
        _swig_setattr(self, Interface,self.__class__,Interface)
_LAPACK.Interface_swigregister(InterfacePtr)

class Group(NOX.LAPACK.Group,Abstract.Group):
    __swig_setmethods__ = {}
    for _s in [NOX.LAPACK.Group,Abstract.Group]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Group, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.LAPACK.Group,Abstract.Group]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Group, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::LAPACK::Group instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Group, 'this', _LAPACK.new_Group(*args))
        _swig_setattr(self, Group, 'thisown', 1)
    def __del__(self, destroy=_LAPACK.delete_Group):
        try:
            if self.thisown: destroy(self)
        except: pass

    def clone(*args): return _LAPACK.Group_clone(*args)
    def computeF(*args): return _LAPACK.Group_computeF(*args)
    def computeJacobian(*args): return _LAPACK.Group_computeJacobian(*args)
    def setParams(*args): return _LAPACK.Group_setParams(*args)
    def setParam(*args): return _LAPACK.Group_setParam(*args)
    def getParams(*args): return _LAPACK.Group_getParams(*args)
    def getParam(*args): return _LAPACK.Group_getParam(*args)
    def computeScaledDotProduct(*args): return _LAPACK.Group_computeScaledDotProduct(*args)
    def scaleVector(*args): return _LAPACK.Group_scaleVector(*args)
    def applyJacobianInverseMulti(*args): return _LAPACK.Group_applyJacobianInverseMulti(*args)
    def projectToDraw(*args): return _LAPACK.Group_projectToDraw(*args)
    def projectToDrawDimension(*args): return _LAPACK.Group_projectToDrawDimension(*args)
    def applyBorderedJacobianInverse(*args): return _LAPACK.Group_applyBorderedJacobianInverse(*args)
    def computeMassMatrix(*args): return _LAPACK.Group_computeMassMatrix(*args)
    def isMassMatrix(*args): return _LAPACK.Group_isMassMatrix(*args)
    def applyComplexInverse(*args): return _LAPACK.Group_applyComplexInverse(*args)
    def applyComplexInverseMulti(*args): return _LAPACK.Group_applyComplexInverseMulti(*args)
    def augmentJacobianForHomotopy(*args): return _LAPACK.Group_augmentJacobianForHomotopy(*args)
    def printSolution(*args): return _LAPACK.Group_printSolution(*args)
    def applyMassMatrix(*args): return _LAPACK.Group_applyMassMatrix(*args)
    def hasMass(*args): return _LAPACK.Group_hasMass(*args)
    def getJacobianMatrix(*args): return _LAPACK.Group_getJacobianMatrix(*args)
    def getMassMatrix(*args): return _LAPACK.Group_getMassMatrix(*args)

class GroupPtr(Group):
    def __init__(self, this):
        _swig_setattr(self, Group, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Group, 'thisown', 0)
        _swig_setattr(self, Group,self.__class__,Group)
_LAPACK.Group_swigregister(GroupPtr)



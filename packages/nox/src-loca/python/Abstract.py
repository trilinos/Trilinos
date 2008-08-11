# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Abstract

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


import Continuation
import NOX.Abstract
import NOX.StatusTest
import MultiContinuation
import Homotopy
import TimeDependent
import Bifurcation
class Group(Bifurcation.HopfBordFiniteDifferenceGroup,Bifurcation.TPBordSingularSolveGroup,Homotopy.AbstractGroup,MultiContinuation.FiniteDifferenceGroup):
    __swig_setmethods__ = {}
    for _s in [Bifurcation.HopfBordFiniteDifferenceGroup,Bifurcation.TPBordSingularSolveGroup,Homotopy.AbstractGroup,MultiContinuation.FiniteDifferenceGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Group, name, value)
    __swig_getmethods__ = {}
    for _s in [Bifurcation.HopfBordFiniteDifferenceGroup,Bifurcation.TPBordSingularSolveGroup,Homotopy.AbstractGroup,MultiContinuation.FiniteDifferenceGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Group, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Abstract::Group instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Abstract.delete_Group):
        try:
            if self.thisown: destroy(self)
        except: pass

    def applyBorderedJacobianInverse(*args): return _Abstract.Group_applyBorderedJacobianInverse(*args)
    def applyShiftedMatrixInverse(*args): return _Abstract.Group_applyShiftedMatrixInverse(*args)
    def applyComplexInverse(*args): return _Abstract.Group_applyComplexInverse(*args)
    def augmentJacobianForHomotopy(*args): return _Abstract.Group_augmentJacobianForHomotopy(*args)
    def setParamsMulti(*args): return _Abstract.Group_setParamsMulti(*args)

class GroupPtr(Group):
    def __init__(self, this):
        _swig_setattr(self, Group, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Group, 'thisown', 0)
        _swig_setattr(self, Group,self.__class__,Group)
_Abstract.Group_swigregister(GroupPtr)



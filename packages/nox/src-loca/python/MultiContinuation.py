# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _MultiContinuation

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
class AbstractGroup(Continuation.AbstractGroup):
    __swig_setmethods__ = {}
    for _s in [Continuation.AbstractGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, AbstractGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [Continuation.AbstractGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, AbstractGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::MultiContinuation::AbstractGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_MultiContinuation.delete_AbstractGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setParamsMulti(*args): return _MultiContinuation.AbstractGroup_setParamsMulti(*args)
    def computeDfDpMulti(*args): return _MultiContinuation.AbstractGroup_computeDfDpMulti(*args)
    def projectToDraw(*args): return _MultiContinuation.AbstractGroup_projectToDraw(*args)
    def projectToDrawDimension(*args): return _MultiContinuation.AbstractGroup_projectToDrawDimension(*args)

class AbstractGroupPtr(AbstractGroup):
    def __init__(self, this):
        _swig_setattr(self, AbstractGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, AbstractGroup, 'thisown', 0)
        _swig_setattr(self, AbstractGroup,self.__class__,AbstractGroup)
_MultiContinuation.AbstractGroup_swigregister(AbstractGroupPtr)

class FiniteDifferenceGroup(AbstractGroup):
    __swig_setmethods__ = {}
    for _s in [AbstractGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FiniteDifferenceGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [AbstractGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FiniteDifferenceGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::MultiContinuation::FiniteDifferenceGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_MultiContinuation.delete_FiniteDifferenceGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def computeDfDpMulti(*args): return _MultiContinuation.FiniteDifferenceGroup_computeDfDpMulti(*args)

class FiniteDifferenceGroupPtr(FiniteDifferenceGroup):
    def __init__(self, this):
        _swig_setattr(self, FiniteDifferenceGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FiniteDifferenceGroup, 'thisown', 0)
        _swig_setattr(self, FiniteDifferenceGroup,self.__class__,FiniteDifferenceGroup)
_MultiContinuation.FiniteDifferenceGroup_swigregister(FiniteDifferenceGroupPtr)



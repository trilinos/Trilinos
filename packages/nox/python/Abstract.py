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


class Group(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Group, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Group, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Abstract::Group instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    Ok = _Abstract.Group_Ok
    NotDefined = _Abstract.Group_NotDefined
    BadDependency = _Abstract.Group_BadDependency
    NotConverged = _Abstract.Group_NotConverged
    Failed = _Abstract.Group_Failed
    def __del__(self, destroy=_Abstract.delete_Group):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setX(*args): return _Abstract.Group_setX(*args)
    def computeX(*args): return _Abstract.Group_computeX(*args)
    def computeF(*args): return _Abstract.Group_computeF(*args)
    def computeJacobian(*args): return _Abstract.Group_computeJacobian(*args)
    def computeGradient(*args): return _Abstract.Group_computeGradient(*args)
    def computeNewton(*args): return _Abstract.Group_computeNewton(*args)
    def applyJacobian(*args): return _Abstract.Group_applyJacobian(*args)
    def applyJacobianTranspose(*args): return _Abstract.Group_applyJacobianTranspose(*args)
    def applyJacobianInverse(*args): return _Abstract.Group_applyJacobianInverse(*args)
    def applyRightPreconditioning(*args): return _Abstract.Group_applyRightPreconditioning(*args)
    def isF(*args): return _Abstract.Group_isF(*args)
    def isJacobian(*args): return _Abstract.Group_isJacobian(*args)
    def isGradient(*args): return _Abstract.Group_isGradient(*args)
    def isNewton(*args): return _Abstract.Group_isNewton(*args)
    def getX(*args): return _Abstract.Group_getX(*args)
    def getF(*args): return _Abstract.Group_getF(*args)
    def getNormF(*args): return _Abstract.Group_getNormF(*args)
    def getGradient(*args): return _Abstract.Group_getGradient(*args)
    def getNewton(*args): return _Abstract.Group_getNewton(*args)
    def getNormLastLinearSolveResidual(*args): return _Abstract.Group_getNormLastLinearSolveResidual(*args)
    def clone(*args): return _Abstract.Group_clone(*args)

class GroupPtr(Group):
    def __init__(self, this):
        _swig_setattr(self, Group, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Group, 'thisown', 0)
        _swig_setattr(self, Group,self.__class__,Group)
_Abstract.Group_swigregister(GroupPtr)

DeepCopy = _Abstract.DeepCopy
ShapeCopy = _Abstract.ShapeCopy
class Vector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vector, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Abstract::Vector instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    TwoNorm = _Abstract.Vector_TwoNorm
    OneNorm = _Abstract.Vector_OneNorm
    MaxNorm = _Abstract.Vector_MaxNorm
    def __del__(self, destroy=_Abstract.delete_Vector):
        try:
            if self.thisown: destroy(self)
        except: pass

    def init(*args): return _Abstract.Vector_init(*args)
    def random(*args): return _Abstract.Vector_random(*args)
    def abs(*args): return _Abstract.Vector_abs(*args)
    def reciprocal(*args): return _Abstract.Vector_reciprocal(*args)
    def scale(*args): return _Abstract.Vector_scale(*args)
    def update(*args): return _Abstract.Vector_update(*args)
    def clone(*args): return _Abstract.Vector_clone(*args)
    def norm(*args): return _Abstract.Vector_norm(*args)
    def dot(*args): return _Abstract.Vector_dot(*args)
    def length(*args): return _Abstract.Vector_length(*args)

class VectorPtr(Vector):
    def __init__(self, this):
        _swig_setattr(self, Vector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vector, 'thisown', 0)
        _swig_setattr(self, Vector,self.__class__,Vector)
_Abstract.Vector_swigregister(VectorPtr)



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
class Vector(Abstract.Vector):
    __swig_setmethods__ = {}
    for _s in [Abstract.Vector]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vector, name, value)
    __swig_getmethods__ = {}
    for _s in [Abstract.Vector]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Vector, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::LAPACK::Vector instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Vector, 'this', _LAPACK.new_Vector(*args))
        _swig_setattr(self, Vector, 'thisown', 1)
    def __del__(self, destroy=_LAPACK.delete_Vector):
        try:
            if self.thisown: destroy(self)
        except: pass

    def init(*args): return _LAPACK.Vector_init(*args)
    def random(*args): return _LAPACK.Vector_random(*args)
    def abs(*args): return _LAPACK.Vector_abs(*args)
    def reciprocal(*args): return _LAPACK.Vector_reciprocal(*args)
    def scale(*args): return _LAPACK.Vector_scale(*args)
    def update(*args): return _LAPACK.Vector_update(*args)
    def clone(*args): return _LAPACK.Vector_clone(*args)
    def norm(*args): return _LAPACK.Vector_norm(*args)
    def dot(*args): return _LAPACK.Vector_dot(*args)
    def length(*args): return _LAPACK.Vector_length(*args)
    def __call__(*args): return _LAPACK.Vector___call__(*args)
    def leftshift(*args): return _LAPACK.Vector_leftshift(*args)

class VectorPtr(Vector):
    def __init__(self, this):
        _swig_setattr(self, Vector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vector, 'thisown', 0)
        _swig_setattr(self, Vector,self.__class__,Vector)
_LAPACK.Vector_swigregister(VectorPtr)
cvar = _LAPACK.cvar
d_one = cvar.d_one
d_mone = cvar.d_mone
d_zero = cvar.d_zero
i_one = cvar.i_one
i_zero = cvar.i_zero

class Interface(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Interface, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Interface, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::LAPACK::Interface instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_LAPACK.delete_Interface):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getInitialGuess(*args): return _LAPACK.Interface_getInitialGuess(*args)
    def computeF(*args): return _LAPACK.Interface_computeF(*args)
    def computeJacobian(*args): return _LAPACK.Interface_computeJacobian(*args)

class InterfacePtr(Interface):
    def __init__(self, this):
        _swig_setattr(self, Interface, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Interface, 'thisown', 0)
        _swig_setattr(self, Interface,self.__class__,Interface)
_LAPACK.Interface_swigregister(InterfacePtr)

class Group(Abstract.Group):
    __swig_setmethods__ = {}
    for _s in [Abstract.Group]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Group, name, value)
    __swig_getmethods__ = {}
    for _s in [Abstract.Group]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Group, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::LAPACK::Group instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Group, 'this', _LAPACK.new_Group(*args))
        _swig_setattr(self, Group, 'thisown', 1)
    def __del__(self, destroy=_LAPACK.delete_Group):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setX(*args): return _LAPACK.Group_setX(*args)
    def computeX(*args): return _LAPACK.Group_computeX(*args)
    def computeF(*args): return _LAPACK.Group_computeF(*args)
    def computeJacobian(*args): return _LAPACK.Group_computeJacobian(*args)
    def computeGradient(*args): return _LAPACK.Group_computeGradient(*args)
    def computeNewton(*args): return _LAPACK.Group_computeNewton(*args)
    def applyJacobian(*args): return _LAPACK.Group_applyJacobian(*args)
    def applyJacobianTranspose(*args): return _LAPACK.Group_applyJacobianTranspose(*args)
    def applyJacobianInverse(*args): return _LAPACK.Group_applyJacobianInverse(*args)
    def isF(*args): return _LAPACK.Group_isF(*args)
    def isJacobian(*args): return _LAPACK.Group_isJacobian(*args)
    def isGradient(*args): return _LAPACK.Group_isGradient(*args)
    def isNewton(*args): return _LAPACK.Group_isNewton(*args)
    def getX(*args): return _LAPACK.Group_getX(*args)
    def getF(*args): return _LAPACK.Group_getF(*args)
    def getNormF(*args): return _LAPACK.Group_getNormF(*args)
    def getGradient(*args): return _LAPACK.Group_getGradient(*args)
    def getNewton(*args): return _LAPACK.Group_getNewton(*args)
    def clone(*args): return _LAPACK.Group_clone(*args)

class GroupPtr(Group):
    def __init__(self, this):
        _swig_setattr(self, Group, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Group, 'thisown', 0)
        _swig_setattr(self, Group,self.__class__,Group)
_LAPACK.Group_swigregister(GroupPtr)



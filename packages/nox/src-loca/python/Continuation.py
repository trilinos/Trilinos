# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Continuation

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


import NOX.Abstract
import NOX.StatusTest
class AbstractGroup(NOX.Abstract.Group):
    __swig_setmethods__ = {}
    for _s in [NOX.Abstract.Group]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, AbstractGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.Abstract.Group]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, AbstractGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Continuation::AbstractGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Continuation.delete_AbstractGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setParams(*args): return _Continuation.AbstractGroup_setParams(*args)
    def setParam(*args): return _Continuation.AbstractGroup_setParam(*args)
    def getParams(*args): return _Continuation.AbstractGroup_getParams(*args)
    def getParam(*args): return _Continuation.AbstractGroup_getParam(*args)
    def computeDfDp(*args): return _Continuation.AbstractGroup_computeDfDp(*args)
    def applyJacobianInverseMulti(*args): return _Continuation.AbstractGroup_applyJacobianInverseMulti(*args)
    def computeScaledDotProduct(*args): return _Continuation.AbstractGroup_computeScaledDotProduct(*args)
    def printSolution(*args): return _Continuation.AbstractGroup_printSolution(*args)
    def applyHouseholderJacobianInverse(*args): return _Continuation.AbstractGroup_applyHouseholderJacobianInverse(*args)
    def scaleVector(*args): return _Continuation.AbstractGroup_scaleVector(*args)

class AbstractGroupPtr(AbstractGroup):
    def __init__(self, this):
        _swig_setattr(self, AbstractGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, AbstractGroup, 'thisown', 0)
        _swig_setattr(self, AbstractGroup,self.__class__,AbstractGroup)
_Continuation.AbstractGroup_swigregister(AbstractGroupPtr)

class FiniteDifferenceGroup(AbstractGroup):
    __swig_setmethods__ = {}
    for _s in [AbstractGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FiniteDifferenceGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [AbstractGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FiniteDifferenceGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Continuation::FiniteDifferenceGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Continuation.delete_FiniteDifferenceGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def computeDfDp(*args): return _Continuation.FiniteDifferenceGroup_computeDfDp(*args)

class FiniteDifferenceGroupPtr(FiniteDifferenceGroup):
    def __init__(self, this):
        _swig_setattr(self, FiniteDifferenceGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FiniteDifferenceGroup, 'thisown', 0)
        _swig_setattr(self, FiniteDifferenceGroup,self.__class__,FiniteDifferenceGroup)
_Continuation.FiniteDifferenceGroup_swigregister(FiniteDifferenceGroupPtr)

class ParameterResidualNorm(NOX.StatusTest.Generic):
    __swig_setmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ParameterResidualNorm, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ParameterResidualNorm, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Continuation::StatusTest::ParameterResidualNorm instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ParameterResidualNorm, 'this', _Continuation.new_ParameterResidualNorm(*args))
        _swig_setattr(self, ParameterResidualNorm, 'thisown', 1)
    def __del__(self, destroy=_Continuation.delete_ParameterResidualNorm):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Continuation.ParameterResidualNorm_checkStatus(*args)
    def getStatus(*args): return _Continuation.ParameterResidualNorm_getStatus(*args)
    def getParameterResidualNorm(*args): return _Continuation.ParameterResidualNorm_getParameterResidualNorm(*args)
    def getRTOL(*args): return _Continuation.ParameterResidualNorm_getRTOL(*args)
    def getATOL(*args): return _Continuation.ParameterResidualNorm_getATOL(*args)
    def getTOL(*args): return _Continuation.ParameterResidualNorm_getTOL(*args)

class ParameterResidualNormPtr(ParameterResidualNorm):
    def __init__(self, this):
        _swig_setattr(self, ParameterResidualNorm, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ParameterResidualNorm, 'thisown', 0)
        _swig_setattr(self, ParameterResidualNorm,self.__class__,ParameterResidualNorm)
_Continuation.ParameterResidualNorm_swigregister(ParameterResidualNormPtr)

class ParameterUpdateNorm(NOX.StatusTest.Generic):
    __swig_setmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ParameterUpdateNorm, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ParameterUpdateNorm, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Continuation::StatusTest::ParameterUpdateNorm instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ParameterUpdateNorm, 'this', _Continuation.new_ParameterUpdateNorm(*args))
        _swig_setattr(self, ParameterUpdateNorm, 'thisown', 1)
    def __del__(self, destroy=_Continuation.delete_ParameterUpdateNorm):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Continuation.ParameterUpdateNorm_checkStatus(*args)
    def getStatus(*args): return _Continuation.ParameterUpdateNorm_getStatus(*args)
    def getParameterUpdateNorm(*args): return _Continuation.ParameterUpdateNorm_getParameterUpdateNorm(*args)
    def getRTOL(*args): return _Continuation.ParameterUpdateNorm_getRTOL(*args)
    def getATOL(*args): return _Continuation.ParameterUpdateNorm_getATOL(*args)
    def getTOL(*args): return _Continuation.ParameterUpdateNorm_getTOL(*args)

class ParameterUpdateNormPtr(ParameterUpdateNorm):
    def __init__(self, this):
        _swig_setattr(self, ParameterUpdateNorm, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ParameterUpdateNorm, 'thisown', 0)
        _swig_setattr(self, ParameterUpdateNorm,self.__class__,ParameterUpdateNorm)
_Continuation.ParameterUpdateNorm_swigregister(ParameterUpdateNormPtr)



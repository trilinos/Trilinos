# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Solver

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


Unevaluated = _Solver.Unevaluated
Unconverged = _Solver.Unconverged
Converged = _Solver.Converged
Failed = _Solver.Failed
Complete = _Solver.Complete
Minimal = _Solver.Minimal
StatusTest_None = _Solver.StatusTest_None
class StatusTest_Generic(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, StatusTest_Generic, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, StatusTest_Generic, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::Generic instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Solver.delete_StatusTest_Generic):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Solver.StatusTest_Generic_checkStatus(*args)
    def checkStatusEfficiently(*args): return _Solver.StatusTest_Generic_checkStatusEfficiently(*args)
    def getStatus(*args): return _Solver.StatusTest_Generic_getStatus(*args)

class StatusTest_GenericPtr(StatusTest_Generic):
    def __init__(self, this):
        _swig_setattr(self, StatusTest_Generic, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, StatusTest_Generic, 'thisown', 0)
        _swig_setattr(self, StatusTest_Generic,self.__class__,StatusTest_Generic)
_Solver.StatusTest_Generic_swigregister(StatusTest_GenericPtr)

class Generic(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Generic, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Generic, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Solver::Generic instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Solver.delete_Generic):
        try:
            if self.thisown: destroy(self)
        except: pass

    def reset(*args): return _Solver.Generic_reset(*args)
    def getStatus(*args): return _Solver.Generic_getStatus(*args)
    def iterate(*args): return _Solver.Generic_iterate(*args)
    def solve(*args): return _Solver.Generic_solve(*args)
    def getSolutionGroup(*args): return _Solver.Generic_getSolutionGroup(*args)
    def getPreviousSolutionGroup(*args): return _Solver.Generic_getPreviousSolutionGroup(*args)
    def getNumIterations(*args): return _Solver.Generic_getNumIterations(*args)
    def getParameterList(*args): return _Solver.Generic_getParameterList(*args)

class GenericPtr(Generic):
    def __init__(self, this):
        _swig_setattr(self, Generic, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Generic, 'thisown', 0)
        _swig_setattr(self, Generic,self.__class__,Generic)
_Solver.Generic_swigregister(GenericPtr)

class Manager(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Manager, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Manager, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Solver::Manager instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Manager, 'this', _Solver.new_Manager(*args))
        _swig_setattr(self, Manager, 'thisown', 1)
    def __del__(self, destroy=_Solver.delete_Manager):
        try:
            if self.thisown: destroy(self)
        except: pass

    def reset(*args): return _Solver.Manager_reset(*args)
    def getStatus(*args): return _Solver.Manager_getStatus(*args)
    def iterate(*args): return _Solver.Manager_iterate(*args)
    def solve(*args): return _Solver.Manager_solve(*args)
    def getSolutionGroup(*args): return _Solver.Manager_getSolutionGroup(*args)
    def getPreviousSolutionGroup(*args): return _Solver.Manager_getPreviousSolutionGroup(*args)
    def getNumIterations(*args): return _Solver.Manager_getNumIterations(*args)
    def getParameterList(*args): return _Solver.Manager_getParameterList(*args)

class ManagerPtr(Manager):
    def __init__(self, this):
        _swig_setattr(self, Manager, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Manager, 'thisown', 0)
        _swig_setattr(self, Manager,self.__class__,Manager)
_Solver.Manager_swigregister(ManagerPtr)



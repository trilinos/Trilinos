# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _TopLevel

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


class Iterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Iterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Iterator, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Abstract::Iterator instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    LastIteration = _TopLevel.Iterator_LastIteration
    Finished = _TopLevel.Iterator_Finished
    Failed = _TopLevel.Iterator_Failed
    NotFinished = _TopLevel.Iterator_NotFinished
    Successful = _TopLevel.Iterator_Successful
    Unsuccessful = _TopLevel.Iterator_Unsuccessful
    Provisional = _TopLevel.Iterator_Provisional
    def __del__(self, destroy=_TopLevel.delete_Iterator):
        try:
            if self.thisown: destroy(self)
        except: pass

    def resetIterator(*args): return _TopLevel.Iterator_resetIterator(*args)
    def getIteratorStatus(*args): return _TopLevel.Iterator_getIteratorStatus(*args)
    def getStepNumber(*args): return _TopLevel.Iterator_getStepNumber(*args)
    def getNumFailedSteps(*args): return _TopLevel.Iterator_getNumFailedSteps(*args)
    def getNumTotalSteps(*args): return _TopLevel.Iterator_getNumTotalSteps(*args)
    def run(*args): return _TopLevel.Iterator_run(*args)

class IteratorPtr(Iterator):
    def __init__(self, this):
        _swig_setattr(self, Iterator, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Iterator, 'thisown', 0)
        _swig_setattr(self, Iterator,self.__class__,Iterator)
_TopLevel.Iterator_swigregister(IteratorPtr)

class Stepper(Iterator):
    __swig_setmethods__ = {}
    for _s in [Iterator]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Stepper, name, value)
    __swig_getmethods__ = {}
    for _s in [Iterator]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Stepper, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Stepper instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Stepper, 'this', _TopLevel.new_Stepper(*args))
        _swig_setattr(self, Stepper, 'thisown', 1)
    def __del__(self, destroy=_TopLevel.delete_Stepper):
        try:
            if self.thisown: destroy(self)
        except: pass

    def reset(*args): return _TopLevel.Stepper_reset(*args)
    def getSolutionGroup(*args): return _TopLevel.Stepper_getSolutionGroup(*args)
    def getBifurcationGroup(*args): return _TopLevel.Stepper_getBifurcationGroup(*args)
    def getParameterList(*args): return _TopLevel.Stepper_getParameterList(*args)
    def getSolver(*args): return _TopLevel.Stepper_getSolver(*args)

class StepperPtr(Stepper):
    def __init__(self, this):
        _swig_setattr(self, Stepper, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Stepper, 'thisown', 0)
        _swig_setattr(self, Stepper,self.__class__,Stepper)
_TopLevel.Stepper_swigregister(StepperPtr)

class ParameterVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ParameterVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ParameterVector, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::ParameterVector instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ParameterVector, 'this', _TopLevel.new_ParameterVector(*args))
        _swig_setattr(self, ParameterVector, 'thisown', 1)
    def clone(*args): return _TopLevel.ParameterVector_clone(*args)
    def __del__(self, destroy=_TopLevel.delete_ParameterVector):
        try:
            if self.thisown: destroy(self)
        except: pass

    def addParameter(*args): return _TopLevel.ParameterVector_addParameter(*args)
    def init(*args): return _TopLevel.ParameterVector_init(*args)
    def scale(*args): return _TopLevel.ParameterVector_scale(*args)
    def update(*args): return _TopLevel.ParameterVector_update(*args)
    def setValue(*args): return _TopLevel.ParameterVector_setValue(*args)
    def getValue(*args): return _TopLevel.ParameterVector_getValue(*args)
    def getIndex(*args): return _TopLevel.ParameterVector_getIndex(*args)
    def getDoubleArrayPointer(*args): return _TopLevel.ParameterVector_getDoubleArrayPointer(*args)
    def isParameter(*args): return _TopLevel.ParameterVector_isParameter(*args)
    def getLabel(*args): return _TopLevel.ParameterVector_getLabel(*args)
    def length(*args): return _TopLevel.ParameterVector_length(*args)
    def Print(*args): return _TopLevel.ParameterVector_Print(*args)
    def getValuesVector(*args): return _TopLevel.ParameterVector_getValuesVector(*args)
    def getNamesVector(*args): return _TopLevel.ParameterVector_getNamesVector(*args)

class ParameterVectorPtr(ParameterVector):
    def __init__(self, this):
        _swig_setattr(self, ParameterVector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ParameterVector, 'thisown', 0)
        _swig_setattr(self, ParameterVector,self.__class__,ParameterVector)
_TopLevel.ParameterVector_swigregister(ParameterVectorPtr)

class Utils(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Utils, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Utils, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ Utils instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    Error = _TopLevel.Utils_Error
    Warning = _TopLevel.Utils_Warning
    StepperIteration = _TopLevel.Utils_StepperIteration
    StepperDetails = _TopLevel.Utils_StepperDetails
    Solver = _TopLevel.Utils_Solver
    SolverDetails = _TopLevel.Utils_SolverDetails
    Direction = _TopLevel.Utils_Direction
    Parameters = _TopLevel.Utils_Parameters
    def __init__(self, *args):
        _swig_setattr(self, Utils, 'this', _TopLevel.new_Utils(*args))
        _swig_setattr(self, Utils, 'thisown', 1)
    def __del__(self, destroy=_TopLevel.delete_Utils):
        try:
            if self.thisown: destroy(self)
        except: pass


class UtilsPtr(Utils):
    def __init__(self, this):
        _swig_setattr(self, Utils, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Utils, 'thisown', 0)
        _swig_setattr(self, Utils,self.__class__,Utils)
_TopLevel.Utils_swigregister(UtilsPtr)



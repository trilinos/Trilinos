# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _StatusTest

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


Unevaluated = _StatusTest.Unevaluated
Unconverged = _StatusTest.Unconverged
Converged = _StatusTest.Converged
Failed = _StatusTest.Failed
Complete = _StatusTest.Complete
Minimal = _StatusTest.Minimal
StatusTest_None = _StatusTest.StatusTest_None
class Generic(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Generic, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Generic, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::Generic instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_StatusTest.delete_Generic):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.Generic_checkStatus(*args)
    def checkStatusEfficiently(*args): return _StatusTest.Generic_checkStatusEfficiently(*args)
    def getStatus(*args): return _StatusTest.Generic_getStatus(*args)
    def __str__(*args): return _StatusTest.Generic___str__(*args)

class GenericPtr(Generic):
    def __init__(self, this):
        _swig_setattr(self, Generic, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Generic, 'thisown', 0)
        _swig_setattr(self, Generic,self.__class__,Generic)
_StatusTest.Generic_swigregister(GenericPtr)

class Combo(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Combo, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Combo, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::Combo instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    AND = _StatusTest.Combo_AND
    OR = _StatusTest.Combo_OR
    def __init__(self, *args):
        _swig_setattr(self, Combo, 'this', _StatusTest.new_Combo(*args))
        _swig_setattr(self, Combo, 'thisown', 1)
    def addStatusTest(*args): return _StatusTest.Combo_addStatusTest(*args)
    def __del__(self, destroy=_StatusTest.delete_Combo):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.Combo_checkStatus(*args)
    def checkStatusEfficiently(*args): return _StatusTest.Combo_checkStatusEfficiently(*args)
    def getStatus(*args): return _StatusTest.Combo_getStatus(*args)

class ComboPtr(Combo):
    def __init__(self, this):
        _swig_setattr(self, Combo, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Combo, 'thisown', 0)
        _swig_setattr(self, Combo,self.__class__,Combo)
_StatusTest.Combo_swigregister(ComboPtr)

class NormF(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, NormF, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, NormF, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::NormF instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    Unscaled = _StatusTest.NormF_Unscaled
    Scaled = _StatusTest.NormF_Scaled
    Relative = _StatusTest.NormF_Relative
    Absolute = _StatusTest.NormF_Absolute
    def __init__(self, *args):
        _swig_setattr(self, NormF, 'this', _StatusTest.new_NormF(*args))
        _swig_setattr(self, NormF, 'thisown', 1)
    def __del__(self, destroy=_StatusTest.delete_NormF):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.NormF_checkStatus(*args)
    def checkStatusEfficiently(*args): return _StatusTest.NormF_checkStatusEfficiently(*args)
    def getStatus(*args): return _StatusTest.NormF_getStatus(*args)
    def getNormF(*args): return _StatusTest.NormF_getNormF(*args)
    def getTrueTolerance(*args): return _StatusTest.NormF_getTrueTolerance(*args)
    def getSpecifiedTolerance(*args): return _StatusTest.NormF_getSpecifiedTolerance(*args)
    def getInitialTolerance(*args): return _StatusTest.NormF_getInitialTolerance(*args)

class NormFPtr(NormF):
    def __init__(self, this):
        _swig_setattr(self, NormF, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, NormF, 'thisown', 0)
        _swig_setattr(self, NormF,self.__class__,NormF)
_StatusTest.NormF_swigregister(NormFPtr)

class NormUpdate(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, NormUpdate, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, NormUpdate, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::NormUpdate instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    Unscaled = _StatusTest.NormUpdate_Unscaled
    Scaled = _StatusTest.NormUpdate_Scaled
    def __init__(self, *args):
        _swig_setattr(self, NormUpdate, 'this', _StatusTest.new_NormUpdate(*args))
        _swig_setattr(self, NormUpdate, 'thisown', 1)
    def __del__(self, destroy=_StatusTest.delete_NormUpdate):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.NormUpdate_checkStatus(*args)
    def checkStatusEfficiently(*args): return _StatusTest.NormUpdate_checkStatusEfficiently(*args)
    def getStatus(*args): return _StatusTest.NormUpdate_getStatus(*args)
    def getNormUpdate(*args): return _StatusTest.NormUpdate_getNormUpdate(*args)
    def getTolerance(*args): return _StatusTest.NormUpdate_getTolerance(*args)

class NormUpdatePtr(NormUpdate):
    def __init__(self, this):
        _swig_setattr(self, NormUpdate, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, NormUpdate, 'thisown', 0)
        _swig_setattr(self, NormUpdate,self.__class__,NormUpdate)
_StatusTest.NormUpdate_swigregister(NormUpdatePtr)

class NormWRMS(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, NormWRMS, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, NormWRMS, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::NormWRMS instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, NormWRMS, 'this', _StatusTest.new_NormWRMS(*args))
        _swig_setattr(self, NormWRMS, 'thisown', 1)
    def __del__(self, destroy=_StatusTest.delete_NormWRMS):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.NormWRMS_checkStatus(*args)
    def getStatus(*args): return _StatusTest.NormWRMS_getStatus(*args)
    def getNormWRMS(*args): return _StatusTest.NormWRMS_getNormWRMS(*args)
    def getTolerance(*args): return _StatusTest.NormWRMS_getTolerance(*args)
    def getRTOL(*args): return _StatusTest.NormWRMS_getRTOL(*args)
    def getATOL(*args): return _StatusTest.NormWRMS_getATOL(*args)
    def getBDFMultiplier(*args): return _StatusTest.NormWRMS_getBDFMultiplier(*args)
    def getAlpha(*args): return _StatusTest.NormWRMS_getAlpha(*args)
    def getBeta(*args): return _StatusTest.NormWRMS_getBeta(*args)

class NormWRMSPtr(NormWRMS):
    def __init__(self, this):
        _swig_setattr(self, NormWRMS, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, NormWRMS, 'thisown', 0)
        _swig_setattr(self, NormWRMS,self.__class__,NormWRMS)
_StatusTest.NormWRMS_swigregister(NormWRMSPtr)

class MaxIters(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, MaxIters, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, MaxIters, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::MaxIters instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, MaxIters, 'this', _StatusTest.new_MaxIters(*args))
        _swig_setattr(self, MaxIters, 'thisown', 1)
    def __del__(self, destroy=_StatusTest.delete_MaxIters):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.MaxIters_checkStatus(*args)
    def checkStatusEfficiently(*args): return _StatusTest.MaxIters_checkStatusEfficiently(*args)
    def getStatus(*args): return _StatusTest.MaxIters_getStatus(*args)
    def getMaxIters(*args): return _StatusTest.MaxIters_getMaxIters(*args)
    def getNumIters(*args): return _StatusTest.MaxIters_getNumIters(*args)

class MaxItersPtr(MaxIters):
    def __init__(self, this):
        _swig_setattr(self, MaxIters, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, MaxIters, 'thisown', 0)
        _swig_setattr(self, MaxIters,self.__class__,MaxIters)
_StatusTest.MaxIters_swigregister(MaxItersPtr)

class Stagnation(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Stagnation, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Stagnation, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::Stagnation instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Stagnation, 'this', _StatusTest.new_Stagnation(*args))
        _swig_setattr(self, Stagnation, 'thisown', 1)
    def __del__(self, destroy=_StatusTest.delete_Stagnation):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.Stagnation_checkStatus(*args)
    def getStatus(*args): return _StatusTest.Stagnation_getStatus(*args)
    def getMaxNumSteps(*args): return _StatusTest.Stagnation_getMaxNumSteps(*args)
    def getCurrentNumSteps(*args): return _StatusTest.Stagnation_getCurrentNumSteps(*args)
    def getTolerance(*args): return _StatusTest.Stagnation_getTolerance(*args)
    def getConvRate(*args): return _StatusTest.Stagnation_getConvRate(*args)

class StagnationPtr(Stagnation):
    def __init__(self, this):
        _swig_setattr(self, Stagnation, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Stagnation, 'thisown', 0)
        _swig_setattr(self, Stagnation,self.__class__,Stagnation)
_StatusTest.Stagnation_swigregister(StagnationPtr)

class FiniteValue(Generic):
    __swig_setmethods__ = {}
    for _s in [Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FiniteValue, name, value)
    __swig_getmethods__ = {}
    for _s in [Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FiniteValue, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::StatusTest::FiniteValue instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    FVector = _StatusTest.FiniteValue_FVector
    SolutionVector = _StatusTest.FiniteValue_SolutionVector
    def __init__(self, *args):
        _swig_setattr(self, FiniteValue, 'this', _StatusTest.new_FiniteValue(*args))
        _swig_setattr(self, FiniteValue, 'thisown', 1)
    def __del__(self, destroy=_StatusTest.delete_FiniteValue):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _StatusTest.FiniteValue_checkStatus(*args)
    def checkStatusEfficiently(*args): return _StatusTest.FiniteValue_checkStatusEfficiently(*args)
    def getStatus(*args): return _StatusTest.FiniteValue_getStatus(*args)
    def finiteNumberTest(*args): return _StatusTest.FiniteValue_finiteNumberTest(*args)

class FiniteValuePtr(FiniteValue):
    def __init__(self, this):
        _swig_setattr(self, FiniteValue, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FiniteValue, 'thisown', 0)
        _swig_setattr(self, FiniteValue,self.__class__,FiniteValue)
_StatusTest.FiniteValue_swigregister(FiniteValuePtr)



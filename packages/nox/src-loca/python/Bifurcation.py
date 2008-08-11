# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Bifurcation

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
import TimeDependent
class TPBordAbstractGroup(Continuation.AbstractGroup):
    __swig_setmethods__ = {}
    for _s in [Continuation.AbstractGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, TPBordAbstractGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [Continuation.AbstractGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, TPBordAbstractGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::TPBord::AbstractGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Bifurcation.delete_TPBordAbstractGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def applySingularJacobianInverse(*args): return _Bifurcation.TPBordAbstractGroup_applySingularJacobianInverse(*args)
    def applySingularJacobianInverseMulti(*args): return _Bifurcation.TPBordAbstractGroup_applySingularJacobianInverseMulti(*args)
    def computeDJnDp(*args): return _Bifurcation.TPBordAbstractGroup_computeDJnDp(*args)
    def computeDJnDxa(*args): return _Bifurcation.TPBordAbstractGroup_computeDJnDxa(*args)
    def innerProduct(*args): return _Bifurcation.TPBordAbstractGroup_innerProduct(*args)
    def applyBorderedJacobianInverse(*args): return _Bifurcation.TPBordAbstractGroup_applyBorderedJacobianInverse(*args)
    def applyBorderedJacobianInverseMulti(*args): return _Bifurcation.TPBordAbstractGroup_applyBorderedJacobianInverseMulti(*args)

class TPBordAbstractGroupPtr(TPBordAbstractGroup):
    def __init__(self, this):
        _swig_setattr(self, TPBordAbstractGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, TPBordAbstractGroup, 'thisown', 0)
        _swig_setattr(self, TPBordAbstractGroup,self.__class__,TPBordAbstractGroup)
_Bifurcation.TPBordAbstractGroup_swigregister(TPBordAbstractGroupPtr)

class TPBordFiniteDifferenceGroup(TPBordAbstractGroup,Continuation.FiniteDifferenceGroup):
    __swig_setmethods__ = {}
    for _s in [TPBordAbstractGroup,Continuation.FiniteDifferenceGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, TPBordFiniteDifferenceGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [TPBordAbstractGroup,Continuation.FiniteDifferenceGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, TPBordFiniteDifferenceGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::TPBord::FiniteDifferenceGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Bifurcation.delete_TPBordFiniteDifferenceGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def computeDJnDp(*args): return _Bifurcation.TPBordFiniteDifferenceGroup_computeDJnDp(*args)
    def computeDJnDxa(*args): return _Bifurcation.TPBordFiniteDifferenceGroup_computeDJnDxa(*args)

class TPBordFiniteDifferenceGroupPtr(TPBordFiniteDifferenceGroup):
    def __init__(self, this):
        _swig_setattr(self, TPBordFiniteDifferenceGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, TPBordFiniteDifferenceGroup, 'thisown', 0)
        _swig_setattr(self, TPBordFiniteDifferenceGroup,self.__class__,TPBordFiniteDifferenceGroup)
_Bifurcation.TPBordFiniteDifferenceGroup_swigregister(TPBordFiniteDifferenceGroupPtr)

class TPBordSingularSolveGroup(TPBordAbstractGroup):
    __swig_setmethods__ = {}
    for _s in [TPBordAbstractGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, TPBordSingularSolveGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [TPBordAbstractGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, TPBordSingularSolveGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::TPBord::SingularSolveGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Bifurcation.delete_TPBordSingularSolveGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def applySingularJacobianInverse(*args): return _Bifurcation.TPBordSingularSolveGroup_applySingularJacobianInverse(*args)
    def applySingularJacobianInverseMulti(*args): return _Bifurcation.TPBordSingularSolveGroup_applySingularJacobianInverseMulti(*args)

class TPBordSingularSolveGroupPtr(TPBordSingularSolveGroup):
    def __init__(self, this):
        _swig_setattr(self, TPBordSingularSolveGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, TPBordSingularSolveGroup, 'thisown', 0)
        _swig_setattr(self, TPBordSingularSolveGroup,self.__class__,TPBordSingularSolveGroup)
_Bifurcation.TPBordSingularSolveGroup_swigregister(TPBordSingularSolveGroupPtr)

class HopfBordAbstractGroup(TPBordAbstractGroup,TimeDependent.AbstractGroup):
    __swig_setmethods__ = {}
    for _s in [TPBordAbstractGroup,TimeDependent.AbstractGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, HopfBordAbstractGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [TPBordAbstractGroup,TimeDependent.AbstractGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, HopfBordAbstractGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::HopfBord::AbstractGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Bifurcation.delete_HopfBordAbstractGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def applyComplexInverse(*args): return _Bifurcation.HopfBordAbstractGroup_applyComplexInverse(*args)
    def computeDCeDp(*args): return _Bifurcation.HopfBordAbstractGroup_computeDCeDp(*args)
    def computeDCeDxa(*args): return _Bifurcation.HopfBordAbstractGroup_computeDCeDxa(*args)
    def applyComplex(*args): return _Bifurcation.HopfBordAbstractGroup_applyComplex(*args)
    def applyComplexInverseMulti(*args): return _Bifurcation.HopfBordAbstractGroup_applyComplexInverseMulti(*args)

class HopfBordAbstractGroupPtr(HopfBordAbstractGroup):
    def __init__(self, this):
        _swig_setattr(self, HopfBordAbstractGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, HopfBordAbstractGroup, 'thisown', 0)
        _swig_setattr(self, HopfBordAbstractGroup,self.__class__,HopfBordAbstractGroup)
_Bifurcation.HopfBordAbstractGroup_swigregister(HopfBordAbstractGroupPtr)

class HopfBordFiniteDifferenceGroup(HopfBordAbstractGroup,TPBordFiniteDifferenceGroup):
    __swig_setmethods__ = {}
    for _s in [HopfBordAbstractGroup,TPBordFiniteDifferenceGroup]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, HopfBordFiniteDifferenceGroup, name, value)
    __swig_getmethods__ = {}
    for _s in [HopfBordAbstractGroup,TPBordFiniteDifferenceGroup]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, HopfBordFiniteDifferenceGroup, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_Bifurcation.delete_HopfBordFiniteDifferenceGroup):
        try:
            if self.thisown: destroy(self)
        except: pass

    def computeDCeDp(*args): return _Bifurcation.HopfBordFiniteDifferenceGroup_computeDCeDp(*args)
    def computeDCeDxa(*args): return _Bifurcation.HopfBordFiniteDifferenceGroup_computeDCeDxa(*args)

class HopfBordFiniteDifferenceGroupPtr(HopfBordFiniteDifferenceGroup):
    def __init__(self, this):
        _swig_setattr(self, HopfBordFiniteDifferenceGroup, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, HopfBordFiniteDifferenceGroup, 'thisown', 0)
        _swig_setattr(self, HopfBordFiniteDifferenceGroup,self.__class__,HopfBordFiniteDifferenceGroup)
_Bifurcation.HopfBordFiniteDifferenceGroup_swigregister(HopfBordFiniteDifferenceGroupPtr)

class TPBordNullVectorNormWRMS(NOX.StatusTest.Generic):
    __swig_setmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, TPBordNullVectorNormWRMS, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, TPBordNullVectorNormWRMS, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, TPBordNullVectorNormWRMS, 'this', _Bifurcation.new_TPBordNullVectorNormWRMS(*args))
        _swig_setattr(self, TPBordNullVectorNormWRMS, 'thisown', 1)
    def __del__(self, destroy=_Bifurcation.delete_TPBordNullVectorNormWRMS):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Bifurcation.TPBordNullVectorNormWRMS_checkStatus(*args)
    def getStatus(*args): return _Bifurcation.TPBordNullVectorNormWRMS_getStatus(*args)
    def getNullVectorNormWRMS(*args): return _Bifurcation.TPBordNullVectorNormWRMS_getNullVectorNormWRMS(*args)
    def getRTOL(*args): return _Bifurcation.TPBordNullVectorNormWRMS_getRTOL(*args)
    def getATOL(*args): return _Bifurcation.TPBordNullVectorNormWRMS_getATOL(*args)
    def getTOL(*args): return _Bifurcation.TPBordNullVectorNormWRMS_getTOL(*args)

class TPBordNullVectorNormWRMSPtr(TPBordNullVectorNormWRMS):
    def __init__(self, this):
        _swig_setattr(self, TPBordNullVectorNormWRMS, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, TPBordNullVectorNormWRMS, 'thisown', 0)
        _swig_setattr(self, TPBordNullVectorNormWRMS,self.__class__,TPBordNullVectorNormWRMS)
_Bifurcation.TPBordNullVectorNormWRMS_swigregister(TPBordNullVectorNormWRMSPtr)

class TPBordParameterUpdateNorm(NOX.StatusTest.Generic):
    __swig_setmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, TPBordParameterUpdateNorm, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, TPBordParameterUpdateNorm, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, TPBordParameterUpdateNorm, 'this', _Bifurcation.new_TPBordParameterUpdateNorm(*args))
        _swig_setattr(self, TPBordParameterUpdateNorm, 'thisown', 1)
    def __del__(self, destroy=_Bifurcation.delete_TPBordParameterUpdateNorm):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Bifurcation.TPBordParameterUpdateNorm_checkStatus(*args)
    def getStatus(*args): return _Bifurcation.TPBordParameterUpdateNorm_getStatus(*args)
    def getParameterUpdateNorm(*args): return _Bifurcation.TPBordParameterUpdateNorm_getParameterUpdateNorm(*args)
    def getRTOL(*args): return _Bifurcation.TPBordParameterUpdateNorm_getRTOL(*args)
    def getATOL(*args): return _Bifurcation.TPBordParameterUpdateNorm_getATOL(*args)
    def getTOL(*args): return _Bifurcation.TPBordParameterUpdateNorm_getTOL(*args)

class TPBordParameterUpdateNormPtr(TPBordParameterUpdateNorm):
    def __init__(self, this):
        _swig_setattr(self, TPBordParameterUpdateNorm, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, TPBordParameterUpdateNorm, 'thisown', 0)
        _swig_setattr(self, TPBordParameterUpdateNorm,self.__class__,TPBordParameterUpdateNorm)
_Bifurcation.TPBordParameterUpdateNorm_swigregister(TPBordParameterUpdateNormPtr)

class PitchforkBordNullVectorNormWRMS(NOX.StatusTest.Generic):
    __swig_setmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, PitchforkBordNullVectorNormWRMS, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, PitchforkBordNullVectorNormWRMS, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, PitchforkBordNullVectorNormWRMS, 'this', _Bifurcation.new_PitchforkBordNullVectorNormWRMS(*args))
        _swig_setattr(self, PitchforkBordNullVectorNormWRMS, 'thisown', 1)
    def __del__(self, destroy=_Bifurcation.delete_PitchforkBordNullVectorNormWRMS):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Bifurcation.PitchforkBordNullVectorNormWRMS_checkStatus(*args)
    def getStatus(*args): return _Bifurcation.PitchforkBordNullVectorNormWRMS_getStatus(*args)
    def getNullVectorNormWRMS(*args): return _Bifurcation.PitchforkBordNullVectorNormWRMS_getNullVectorNormWRMS(*args)
    def getRTOL(*args): return _Bifurcation.PitchforkBordNullVectorNormWRMS_getRTOL(*args)
    def getATOL(*args): return _Bifurcation.PitchforkBordNullVectorNormWRMS_getATOL(*args)
    def getTOL(*args): return _Bifurcation.PitchforkBordNullVectorNormWRMS_getTOL(*args)

class PitchforkBordNullVectorNormWRMSPtr(PitchforkBordNullVectorNormWRMS):
    def __init__(self, this):
        _swig_setattr(self, PitchforkBordNullVectorNormWRMS, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, PitchforkBordNullVectorNormWRMS, 'thisown', 0)
        _swig_setattr(self, PitchforkBordNullVectorNormWRMS,self.__class__,PitchforkBordNullVectorNormWRMS)
_Bifurcation.PitchforkBordNullVectorNormWRMS_swigregister(PitchforkBordNullVectorNormWRMSPtr)

class PitchforkBordParameterUpdateNorm(NOX.StatusTest.Generic):
    __swig_setmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, PitchforkBordParameterUpdateNorm, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, PitchforkBordParameterUpdateNorm, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, PitchforkBordParameterUpdateNorm, 'this', _Bifurcation.new_PitchforkBordParameterUpdateNorm(*args))
        _swig_setattr(self, PitchforkBordParameterUpdateNorm, 'thisown', 1)
    def __del__(self, destroy=_Bifurcation.delete_PitchforkBordParameterUpdateNorm):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Bifurcation.PitchforkBordParameterUpdateNorm_checkStatus(*args)
    def getStatus(*args): return _Bifurcation.PitchforkBordParameterUpdateNorm_getStatus(*args)
    def getParameterUpdateNorm(*args): return _Bifurcation.PitchforkBordParameterUpdateNorm_getParameterUpdateNorm(*args)
    def getRTOL(*args): return _Bifurcation.PitchforkBordParameterUpdateNorm_getRTOL(*args)
    def getATOL(*args): return _Bifurcation.PitchforkBordParameterUpdateNorm_getATOL(*args)
    def getTOL(*args): return _Bifurcation.PitchforkBordParameterUpdateNorm_getTOL(*args)

class PitchforkBordParameterUpdateNormPtr(PitchforkBordParameterUpdateNorm):
    def __init__(self, this):
        _swig_setattr(self, PitchforkBordParameterUpdateNorm, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, PitchforkBordParameterUpdateNorm, 'thisown', 0)
        _swig_setattr(self, PitchforkBordParameterUpdateNorm,self.__class__,PitchforkBordParameterUpdateNorm)
_Bifurcation.PitchforkBordParameterUpdateNorm_swigregister(PitchforkBordParameterUpdateNormPtr)

class PitchforkBordSlackUpdateNorm(NOX.StatusTest.Generic):
    __swig_setmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, PitchforkBordSlackUpdateNorm, name, value)
    __swig_getmethods__ = {}
    for _s in [NOX.StatusTest.Generic]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, PitchforkBordSlackUpdateNorm, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, PitchforkBordSlackUpdateNorm, 'this', _Bifurcation.new_PitchforkBordSlackUpdateNorm(*args))
        _swig_setattr(self, PitchforkBordSlackUpdateNorm, 'thisown', 1)
    def __del__(self, destroy=_Bifurcation.delete_PitchforkBordSlackUpdateNorm):
        try:
            if self.thisown: destroy(self)
        except: pass

    def checkStatus(*args): return _Bifurcation.PitchforkBordSlackUpdateNorm_checkStatus(*args)
    def getStatus(*args): return _Bifurcation.PitchforkBordSlackUpdateNorm_getStatus(*args)
    def getSlackUpdateNorm(*args): return _Bifurcation.PitchforkBordSlackUpdateNorm_getSlackUpdateNorm(*args)
    def getRTOL(*args): return _Bifurcation.PitchforkBordSlackUpdateNorm_getRTOL(*args)
    def getATOL(*args): return _Bifurcation.PitchforkBordSlackUpdateNorm_getATOL(*args)
    def getTOL(*args): return _Bifurcation.PitchforkBordSlackUpdateNorm_getTOL(*args)

class PitchforkBordSlackUpdateNormPtr(PitchforkBordSlackUpdateNorm):
    def __init__(self, this):
        _swig_setattr(self, PitchforkBordSlackUpdateNorm, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, PitchforkBordSlackUpdateNorm, 'thisown', 0)
        _swig_setattr(self, PitchforkBordSlackUpdateNorm,self.__class__,PitchforkBordSlackUpdateNorm)
_Bifurcation.PitchforkBordSlackUpdateNorm_swigregister(PitchforkBordSlackUpdateNormPtr)



# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Parameter

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


class List(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, List, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, List, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Parameter::List instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, List, 'this', _Parameter.new_List(*args))
        _swig_setattr(self, List, 'thisown', 1)
    def __del__(self, destroy=_Parameter.delete_List):
        try:
            if self.thisown: destroy(self)
        except: pass

    def unused(*args): return _Parameter.List_unused(*args)
    def sublist(*args): return _Parameter.List_sublist(*args)
    def setParameter(*args): return _Parameter.List_setParameter(*args)
    def getArbitraryParameter(*args): return _Parameter.List_getArbitraryParameter(*args)
    def isParameter(*args): return _Parameter.List_isParameter(*args)
    def isParameterBool(*args): return _Parameter.List_isParameterBool(*args)
    def isParameterInt(*args): return _Parameter.List_isParameterInt(*args)
    def isParameterDouble(*args): return _Parameter.List_isParameterDouble(*args)
    def isParameterString(*args): return _Parameter.List_isParameterString(*args)
    def isParameterSublist(*args): return _Parameter.List_isParameterSublist(*args)
    def isParameterArbitrary(*args): return _Parameter.List_isParameterArbitrary(*args)
    def isParameterEqual(*args): return _Parameter.List_isParameterEqual(*args)
    def getParameter(*args): return _Parameter.List_getParameter(*args)
    def __str__(*args): return _Parameter.List___str__(*args)

class ListPtr(List):
    def __init__(self, this):
        _swig_setattr(self, List, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, List, 'thisown', 0)
        _swig_setattr(self, List,self.__class__,List)
_Parameter.List_swigregister(ListPtr)

class Utils(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Utils, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Utils, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ Utils instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    Error = _Parameter.Utils_Error
    Warning = _Parameter.Utils_Warning
    OuterIteration = _Parameter.Utils_OuterIteration
    InnerIteration = _Parameter.Utils_InnerIteration
    Parameters = _Parameter.Utils_Parameters
    Details = _Parameter.Utils_Details
    OuterIterationStatusTest = _Parameter.Utils_OuterIterationStatusTest
    LinearSolverDetails = _Parameter.Utils_LinearSolverDetails
    TestDetails = _Parameter.Utils_TestDetails
    def __init__(self, *args):
        _swig_setattr(self, Utils, 'this', _Parameter.new_Utils(*args))
        _swig_setattr(self, Utils, 'thisown', 1)
    def __del__(self, destroy=_Parameter.delete_Utils):
        try:
            if self.thisown: destroy(self)
        except: pass


class UtilsPtr(Utils):
    def __init__(self, this):
        _swig_setattr(self, Utils, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Utils, 'thisown', 0)
        _swig_setattr(self, Utils,self.__class__,Utils)
_Parameter.Utils_swigregister(UtilsPtr)



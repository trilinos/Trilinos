# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _StringFunction
import new
new_instancemethod = new.instancemethod
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class StringFunction(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, StringFunction, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, StringFunction, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _StringFunction.new_StringFunction(*args)
        try: self.this.append(this)
        except: self.this = this
    def resolve(*args): return _StringFunction.StringFunction_resolve(*args)
    def getFunctionString(*args): return _StringFunction.StringFunction_getFunctionString(*args)
    def derivative(*args): return _StringFunction.StringFunction_derivative(*args)
    def __call__(*args): return _StringFunction.StringFunction___call__(*args)
    __swig_destroy__ = _StringFunction.delete_StringFunction
    __del__ = lambda self : None;
StringFunction_swigregister = _StringFunction.StringFunction_swigregister
StringFunction_swigregister(StringFunction)

__add__ = _StringFunction.__add__
__div__ = _StringFunction.__div__
__mul__ = _StringFunction.__mul__

__sub__ = _StringFunction.__sub__


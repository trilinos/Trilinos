# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Epetra

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


import PyTrilinos.RawEpetra
import PyTrilinos.EpetraExt
import Abstract
class Interface(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Interface, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Interface, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Epetra::Interface instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    F = _Epetra.Interface_F
    Jacobian = _Epetra.Interface_Jacobian
    Preconditioner = _Epetra.Interface_Preconditioner
    FiniteDifferenceF = _Epetra.Interface_FiniteDifferenceF
    MatrixFreeF = _Epetra.Interface_MatrixFreeF
    def __del__(self, destroy=_Epetra.delete_Interface):
        try:
            if self.thisown: destroy(self)
        except: pass

    def computeF(*args): return _Epetra.Interface_computeF(*args)
    def computeJacobian(*args): return _Epetra.Interface_computeJacobian(*args)
    def computePrecMatrix(*args): return _Epetra.Interface_computePrecMatrix(*args)
    def computePreconditioner(*args): return _Epetra.Interface_computePreconditioner(*args)

class InterfacePtr(Interface):
    def __init__(self, this):
        _swig_setattr(self, Interface, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Interface, 'thisown', 0)
        _swig_setattr(self, Interface,self.__class__,Interface)
_Epetra.Interface_swigregister(InterfacePtr)

class Group(Abstract.Group):
    __swig_setmethods__ = {}
    for _s in [Abstract.Group]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Group, name, value)
    __swig_getmethods__ = {}
    for _s in [Abstract.Group]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Group, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Epetra::Group instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    Group_None = _Epetra.Group_Group_None
    EpetraOperator = _Epetra.Group_EpetraOperator
    EpetraRowMatrix = _Epetra.Group_EpetraRowMatrix
    NoxOperator = _Epetra.Group_NoxOperator
    NoxFiniteDifferenceRowMatrix = _Epetra.Group_NoxFiniteDifferenceRowMatrix
    NoxMatrixFreeOperator = _Epetra.Group_NoxMatrixFreeOperator
    def __init__(self, *args):
        _swig_setattr(self, Group, 'this', _Epetra.new_Group(*args))
        _swig_setattr(self, Group, 'thisown', 1)
    def __del__(self, destroy=_Epetra.delete_Group):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setX(*args): return _Epetra.Group_setX(*args)
    def computeX(*args): return _Epetra.Group_computeX(*args)
    def computeF(*args): return _Epetra.Group_computeF(*args)
    def computeJacobian(*args): return _Epetra.Group_computeJacobian(*args)
    def computeGradient(*args): return _Epetra.Group_computeGradient(*args)
    def computeNewton(*args): return _Epetra.Group_computeNewton(*args)
    def applyJacobian(*args): return _Epetra.Group_applyJacobian(*args)
    def applyJacobianTranspose(*args): return _Epetra.Group_applyJacobianTranspose(*args)
    def applyJacobianInverse(*args): return _Epetra.Group_applyJacobianInverse(*args)
    def applyRightPreconditioning(*args): return _Epetra.Group_applyRightPreconditioning(*args)
    def isF(*args): return _Epetra.Group_isF(*args)
    def isJacobian(*args): return _Epetra.Group_isJacobian(*args)
    def isGradient(*args): return _Epetra.Group_isGradient(*args)
    def isNewton(*args): return _Epetra.Group_isNewton(*args)
    def isNormNewtonSolveResidual(*args): return _Epetra.Group_isNormNewtonSolveResidual(*args)
    def isPreconditioner(*args): return _Epetra.Group_isPreconditioner(*args)
    def getX(*args): return _Epetra.Group_getX(*args)
    def getF(*args): return _Epetra.Group_getF(*args)
    def getNormF(*args): return _Epetra.Group_getNormF(*args)
    def getGradient(*args): return _Epetra.Group_getGradient(*args)
    def getNewton(*args): return _Epetra.Group_getNewton(*args)
    def getNormLastLinearSolveResidual(*args): return _Epetra.Group_getNormLastLinearSolveResidual(*args)
    def clone(*args): return _Epetra.Group_clone(*args)
    def getSharedJacobian(*args): return _Epetra.Group_getSharedJacobian(*args)
    def getSharedPreconditioner(*args): return _Epetra.Group_getSharedPreconditioner(*args)
    def getUserInterface(*args): return _Epetra.Group_getUserInterface(*args)
    def setLinearSolveScaling(*args): return _Epetra.Group_setLinearSolveScaling(*args)
    def getOperatorType(*args): return _Epetra.Group_getOperatorType(*args)
    def getSoln(*args): return _Epetra.Group_getSoln(*args)

class GroupPtr(Group):
    def __init__(self, this):
        _swig_setattr(self, Group, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Group, 'thisown', 0)
        _swig_setattr(self, Group,self.__class__,Group)
_Epetra.Group_swigregister(GroupPtr)

class FiniteDifference(PyTrilinos.RawEpetra.RowMatrix):
    __swig_setmethods__ = {}
    for _s in [PyTrilinos.RawEpetra.RowMatrix]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FiniteDifference, name, value)
    __swig_getmethods__ = {}
    for _s in [PyTrilinos.RawEpetra.RowMatrix]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FiniteDifference, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Epetra::FiniteDifference instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    Forward = _Epetra.FiniteDifference_Forward
    Backward = _Epetra.FiniteDifference_Backward
    Centered = _Epetra.FiniteDifference_Centered
    def __init__(self, *args):
        _swig_setattr(self, FiniteDifference, 'this', _Epetra.new_FiniteDifference(*args))
        _swig_setattr(self, FiniteDifference, 'thisown', 1)
    def __del__(self, destroy=_Epetra.delete_FiniteDifference):
        try:
            if self.thisown: destroy(self)
        except: pass

    def Label(*args): return _Epetra.FiniteDifference_Label(*args)
    def SetUseTranspose(*args): return _Epetra.FiniteDifference_SetUseTranspose(*args)
    def Apply(*args): return _Epetra.FiniteDifference_Apply(*args)
    def ApplyInverse(*args): return _Epetra.FiniteDifference_ApplyInverse(*args)
    def UseTranspose(*args): return _Epetra.FiniteDifference_UseTranspose(*args)
    def HasNormInf(*args): return _Epetra.FiniteDifference_HasNormInf(*args)
    def OperatorDomainMap(*args): return _Epetra.FiniteDifference_OperatorDomainMap(*args)
    def OperatorRangeMap(*args): return _Epetra.FiniteDifference_OperatorRangeMap(*args)
    def Filled(*args): return _Epetra.FiniteDifference_Filled(*args)
    def NumMyRowEntries(*args): return _Epetra.FiniteDifference_NumMyRowEntries(*args)
    def MaxNumEntries(*args): return _Epetra.FiniteDifference_MaxNumEntries(*args)
    def ExtractMyRowCopy(*args): return _Epetra.FiniteDifference_ExtractMyRowCopy(*args)
    def ExtractDiagonalCopy(*args): return _Epetra.FiniteDifference_ExtractDiagonalCopy(*args)
    def Multiply(*args): return _Epetra.FiniteDifference_Multiply(*args)
    def Solve(*args): return _Epetra.FiniteDifference_Solve(*args)
    def InvRowSums(*args): return _Epetra.FiniteDifference_InvRowSums(*args)
    def LeftScale(*args): return _Epetra.FiniteDifference_LeftScale(*args)
    def InvColSums(*args): return _Epetra.FiniteDifference_InvColSums(*args)
    def RightScale(*args): return _Epetra.FiniteDifference_RightScale(*args)
    def NormInf(*args): return _Epetra.FiniteDifference_NormInf(*args)
    def NormOne(*args): return _Epetra.FiniteDifference_NormOne(*args)
    def NumGlobalNonzeros(*args): return _Epetra.FiniteDifference_NumGlobalNonzeros(*args)
    def NumGlobalRows(*args): return _Epetra.FiniteDifference_NumGlobalRows(*args)
    def NumGlobalCols(*args): return _Epetra.FiniteDifference_NumGlobalCols(*args)
    def NumGlobalDiagonals(*args): return _Epetra.FiniteDifference_NumGlobalDiagonals(*args)
    def NumMyNonzeros(*args): return _Epetra.FiniteDifference_NumMyNonzeros(*args)
    def NumMyRows(*args): return _Epetra.FiniteDifference_NumMyRows(*args)
    def NumMyCols(*args): return _Epetra.FiniteDifference_NumMyCols(*args)
    def NumMyDiagonals(*args): return _Epetra.FiniteDifference_NumMyDiagonals(*args)
    def LowerTriangular(*args): return _Epetra.FiniteDifference_LowerTriangular(*args)
    def UpperTriangular(*args): return _Epetra.FiniteDifference_UpperTriangular(*args)
    def Comm(*args): return _Epetra.FiniteDifference_Comm(*args)
    def RowMatrixRowMap(*args): return _Epetra.FiniteDifference_RowMatrixRowMap(*args)
    def RowMatrixColMap(*args): return _Epetra.FiniteDifference_RowMatrixColMap(*args)
    def RowMatrixImporter(*args): return _Epetra.FiniteDifference_RowMatrixImporter(*args)
    def Map(*args): return _Epetra.FiniteDifference_Map(*args)
    def computeJacobian(*args): return _Epetra.FiniteDifference_computeJacobian(*args)
    def computePreconditioner(*args): return _Epetra.FiniteDifference_computePreconditioner(*args)
    def setDifferenceMethod(*args): return _Epetra.FiniteDifference_setDifferenceMethod(*args)
    def getUnderlyingMatrix(*args): return _Epetra.FiniteDifference_getUnderlyingMatrix(*args)
    def Print(*args): return _Epetra.FiniteDifference_Print(*args)

class FiniteDifferencePtr(FiniteDifference):
    def __init__(self, this):
        _swig_setattr(self, FiniteDifference, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FiniteDifference, 'thisown', 0)
        _swig_setattr(self, FiniteDifference,self.__class__,FiniteDifference)
_Epetra.FiniteDifference_swigregister(FiniteDifferencePtr)

class FiniteDifferenceColoring(FiniteDifference):
    __swig_setmethods__ = {}
    for _s in [FiniteDifference]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, FiniteDifferenceColoring, name, value)
    __swig_getmethods__ = {}
    for _s in [FiniteDifference]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, FiniteDifferenceColoring, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Epetra::FiniteDifferenceColoring instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, FiniteDifferenceColoring, 'this', _Epetra.new_FiniteDifferenceColoring(*args))
        _swig_setattr(self, FiniteDifferenceColoring, 'thisown', 1)
    def __del__(self, destroy=_Epetra.delete_FiniteDifferenceColoring):
        try:
            if self.thisown: destroy(self)
        except: pass

    def computeJacobian(*args): return _Epetra.FiniteDifferenceColoring_computeJacobian(*args)

class FiniteDifferenceColoringPtr(FiniteDifferenceColoring):
    def __init__(self, this):
        _swig_setattr(self, FiniteDifferenceColoring, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FiniteDifferenceColoring, 'thisown', 0)
        _swig_setattr(self, FiniteDifferenceColoring,self.__class__,FiniteDifferenceColoring)
_Epetra.FiniteDifferenceColoring_swigregister(FiniteDifferenceColoringPtr)

class Vector(Abstract.Vector):
    __swig_setmethods__ = {}
    for _s in [Abstract.Vector]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vector, name, value)
    __swig_getmethods__ = {}
    for _s in [Abstract.Vector]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, Vector, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ NOX::Epetra::Vector instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Vector, 'this', _Epetra.new_Vector(*args))
        _swig_setattr(self, Vector, 'thisown', 1)
    def __del__(self, destroy=_Epetra.delete_Vector):
        try:
            if self.thisown: destroy(self)
        except: pass

    def getEpetraVector(*args): return _Epetra.Vector_getEpetraVector(*args)
    def init(*args): return _Epetra.Vector_init(*args)
    def random(*args): return _Epetra.Vector_random(*args)
    def abs(*args): return _Epetra.Vector_abs(*args)
    def reciprocal(*args): return _Epetra.Vector_reciprocal(*args)
    def scale(*args): return _Epetra.Vector_scale(*args)
    def update(*args): return _Epetra.Vector_update(*args)
    def clone(*args): return _Epetra.Vector_clone(*args)
    def norm(*args): return _Epetra.Vector_norm(*args)
    def dot(*args): return _Epetra.Vector_dot(*args)
    def length(*args): return _Epetra.Vector_length(*args)

class VectorPtr(Vector):
    def __init__(self, this):
        _swig_setattr(self, Vector, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Vector, 'thisown', 0)
        _swig_setattr(self, Vector,self.__class__,Vector)
_Epetra.Vector_swigregister(VectorPtr)

class Callback(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Callback, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Callback, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ Callback instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, Callback, 'this', _Epetra.new_Callback(*args))
        _swig_setattr(self, Callback, 'thisown', 1)
    def __del__(self, destroy=_Epetra.delete_Callback):
        try:
            if self.thisown: destroy(self)
        except: pass

    def setFunction(*args): return _Epetra.Callback_setFunction(*args)
    def getFunction(*args): return _Epetra.Callback_getFunction(*args)

class CallbackPtr(Callback):
    def __init__(self, this):
        _swig_setattr(self, Callback, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Callback, 'thisown', 0)
        _swig_setattr(self, Callback,self.__class__,Callback)
_Epetra.Callback_swigregister(CallbackPtr)

class PyInterface(Interface):
    __swig_setmethods__ = {}
    for _s in [Interface]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, PyInterface, name, value)
    __swig_getmethods__ = {}
    for _s in [Interface]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, PyInterface, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ PyInterface instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, PyInterface, 'this', _Epetra.new_PyInterface(*args))
        _swig_setattr(self, PyInterface, 'thisown', 1)
    def __del__(self, destroy=_Epetra.delete_PyInterface):
        try:
            if self.thisown: destroy(self)
        except: pass

    def computeF(*args): return _Epetra.PyInterface_computeF(*args)
    def computeJacobian(*args): return _Epetra.PyInterface_computeJacobian(*args)
    def computePrecMatrix(*args): return _Epetra.PyInterface_computePrecMatrix(*args)
    def computePreconditioner(*args): return _Epetra.PyInterface_computePreconditioner(*args)
    def setComputeF(*args): return _Epetra.PyInterface_setComputeF(*args)
    def setComputeJacobian(*args): return _Epetra.PyInterface_setComputeJacobian(*args)
    def setComputePrecMatrix(*args): return _Epetra.PyInterface_setComputePrecMatrix(*args)
    def setComputePreconditioner(*args): return _Epetra.PyInterface_setComputePreconditioner(*args)
    def unloadX(*args): return _Epetra.PyInterface_unloadX(*args)
    def loadRHS(*args): return _Epetra.PyInterface_loadRHS(*args)
    def unloadRHS(*args): return _Epetra.PyInterface_unloadRHS(*args)

class PyInterfacePtr(PyInterface):
    def __init__(self, this):
        _swig_setattr(self, PyInterface, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, PyInterface, 'thisown', 0)
        _swig_setattr(self, PyInterface,self.__class__,PyInterface)
_Epetra.PyInterface_swigregister(PyInterfacePtr)



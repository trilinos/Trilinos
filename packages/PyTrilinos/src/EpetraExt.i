// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %epetraext_docstring
"
PyTrilinos.EpetraExt is the python interface to the Trilinos linear
algebra services extension package EpetraExt:

    http://trilinos.sandia.gov/packages/epetraext

The purpose of EpetraExt is to provide various extensions to Epetra
that were not considered appropriate for the Epetra package.  These
extensions include I/O, matrix-matrix operations and graph coloring.

Currently, only a subset of EpetraExt classes and functions have
python interfaces, including the following user-level classes:

    * XMLReader                      - Read Epetra data from an XML file
    * XMLWriter                      - Write Epetra data as an XML file
    * CrsGraph_MapColoring           - Compute a graph coloring
    * CrsGraph_MapColoringIndex      - Compute indexes for a graph coloring

and functions:

    * MatrixMarketFileToBlockMap     - Read a BlockMap from an MM file
    * MatrixMarketFileToBlockMaps    - Read BlockMaps from an MM file
    * MatrixMarketFileToMap          - Read a Map from an MM file
    * MatrixMarketFileToMultiVector  - Read a MultiVector from an MM file
    * MatrixMarketFileToCrsMatrix    - Read a CrsMatrix from an MM file
    * MatlabFileToCrsMatrix          - Read a CrsMatrix from an ML file
    * BlockMapToMatrixMarketFile     - Write a BlockMap to an MM file
    * MultiVectorToMatrixMarketFile  - Write a MultiVector to an MM file
    * MultiVectorToMatlabFile        - Write a MultiVector to an ML file
    * RowMatrixToMatrixMarketFile    - Write a RowMatrix to an MM file
    * RowMatrixToMatlabFile          - Write a RowMatrix to an ML file
    * Add                            - Add two CrsMatrix objects
    * Multiply                       - Multiply two CrsMatrix objects

For examples of usage, please consult the following scripts in the
example subdirectory of the PyTrilinos package:

    * exEpetraExt_IO_MatrixMarket.py
    * exEpetraExt_IO_XML.py
    * exEpetraExt_MatrixMatrix.py
"
%enddef

%module(package   = "PyTrilinos",
	directors = "1",
	docstring = %epetraext_docstring) EpetraExt

%{
// System includes
#include <vector>

// Configuration includes
#include "PyTrilinos_config.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "Teuchos_XMLObject.hpp"
#include "PyTrilinos_Teuchos_Util.h"

// Epetra includes
#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_LocalMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_InvOperator.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_VbrRowMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_OffsetIndex.h"
#include "Epetra_Time.h"

// Epetra python includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
#include "PyTrilinos_Epetra_Util.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialSymDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"

// EpetraExt includes
#include "EpetraExt_config.h"
#include "EpetraExt_Version.h"
#include "EpetraExt_Exception.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_HDF5.h"
#include "EpetraExt_XMLReader.h"
#include "EpetraExt_XMLWriter.h"
#include "EpetraExt_ModelEvaluator.h"

// EpetraExt python includes
#include "PyTrilinos_EpetraExt_Util.h"
%}

// PyTrilinos configuration
%include "PyTrilinos_config.h"
%include "EpetraExt_config.h"

// Standard exception handling
%include "exception.i"

// Auto-documentation feature
%feature("autodoc", "1");

// Include EpetraExt documentation
%include "EpetraExt_dox.i"

// C++ STL support
%include "stl.i"

// Trilinos interface support
%import "Teuchos.i"
%import "Epetra.i"

// General exception handling
%feature("director:except")
{
  if ($error != NULL)
  {
    throw Swig::DirectorMethodException();
  }
}

%exception
{
  try
  {
    $action
  }
  catch(EpetraExt::Exception & eee)
  {
    eee.Print();
    SWIG_exception(SWIG_UnknownError, "EpetraExt exception");
  }
  catch(Teuchos::EmptyXMLError & e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(Swig::DirectorException &e)
  {
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch(...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown C++ exception");
  }
}

////////////
// Macros //
////////////

// The overloaded HDF5 and XMLReader Read() methods cannot be
// type-disambiguated in python.  We therefore replace selected
// overloaded Read() methods with a python version that has the type
// in the method name.  For example,
//
//   void HDF5::Read(const std::string &, Epetra_Map *&)
//
// is translated from C++ to python as
//
//   HDF5.ReadMap(str) -> Epetra.Map
//
// These translations are made possible by the following macros:
%define %epetraext_primitive_read_method(readtype,TypeName,initVal)
readtype Read ## TypeName(const std::string & group, const std::string & name)
{
  readtype obj = initVal;
  self->Read(group, name, obj);
  return obj;
}
%enddef
%define %epetraext_epetra_read_method(ClassName)
Epetra_ ## ClassName * Read ## ClassName(const std::string & name)
{
  Epetra_ ## ClassName * obj = NULL;
  self->Read(name, obj);
  return obj;
}
%enddef

///////////////////////////////
// EpetraExt_Version support //
///////////////////////////////
%include "EpetraExt_Version.h"
%pythoncode
%{
__version__ = EpetraExt_Version().split()[2]
%}

////////////////////////////
// EpetraExt_HDF5 support //
////////////////////////////
%feature("autodoc",
"ReadBlockMapProperties(str groupName) -> (numGlobalElements, numGlobalPoints,
                                           indexBase, numProc)

Return the properties of a BlockMap from an HDF5 file.")
EpetraExt::HDF5::ReadBlockMapProperties;
%feature("docstring")
EpetraExt::HDF5::ReadBlockMap
"
Return a BlockMap read from an open HDF5 file with group name 'name'.
"
%feature("autodoc",
"ReadMapProperties(str groupName) -> (numGlobalElements, indexBase, numProc)

Return the properties of a Map from an HDF5 file.")
EpetraExt::HDF5::ReadMapProperties;
%feature("docstring")
EpetraExt::HDF5::ReadMap
"
Return a Map read from an open HDF5 file with group name 'name'.
"
%feature("autodoc",
"ReadIntVectorProperties(str groupName) -> int globalLength

Return the global length of an IntVector from an HDF5 file.")
EpetraExt::HDF5::ReadIntVectorProperties;
%feature("docstring")
EpetraExt::HDF5::ReadIntVector
"
Return an IntVector read from an open HDF5 file with group name 'name'.
"
%feature("autodoc",
"ReadMultiVectorProperties(str groupName) -> (globalLength, numVectors)

Return the properties of a MultiVector from an HDF5 file.")
EpetraExt::HDF5::ReadMultiVectorProperties;
%feature("docstring")
EpetraExt::HDF5::ReadMultiVector
"
Return a MultiVector read from an open HDF5 file with group name 'name'.
"
%feature("autodoc",
"ReadCrsGraphProperties(str groupName) -> (numGlobalRows, numGlobalCols,
                                           numGlobalNonzeros, numGlobalDiagonals,
                                           maxNumIndices)

Return the properties of a CrsGraph from an HDF5 file.")
EpetraExt::HDF5::ReadCrsGraphProperties;
%feature("docstring")
EpetraExt::HDF5::ReadCrsGraph
"
Return a CrsGraph read from an open HDF5 file with group name 'name'.
"
%feature("autodoc",
"ReadCrsMatrixProperties(str groupName) -> (numGlobalRows, numGlobalCols,
                                           numNonzeros, numGlobalDiagonals,
                                           maxNumEntries, normOne, normInf)

Return the properties of a CrsMatrix from an HDF5 file.")
EpetraExt::HDF5::ReadCrsMatrixProperties;
%feature("docstring")
EpetraExt::HDF5::ReadCrsMatrix
"
Return a CrsMatrix read from an open HDF5 file with group name 'name'.
"
// There are a series of referenced arguments for the ReadComment()
// and Read*Properties() methods that need to be changed to python
// output arguments
%typemap(in,numinputs=0) std::string& Comment (std::string temp = "") { $1 = &temp; }
%typemap(argout) std::string& Comment
{
  $result = SWIG_Python_AppendOutput($result, PyString_FromString($1->c_str()));
}
%typemap(in,numinputs=0) int& NUM (int temp = 0) { $1 = &temp; }
%typemap(argout) int& NUM
{
  $result = SWIG_Python_AppendOutput($result, PyInt_FromLong((long)(*$1)));
}
%apply int& NUM {int& GlobalLength,
                 int& IndexBase,
                 int& MaxNumEntries,
                 int& MaxNumIndices,
                 int& NumGlobalCols,
                 int& NumGlobalDiagonals,
                 int& NumGlobalElements,
                 int& NumGlobalNonzeros,
                 int& NumGlobalPoints,
                 int& NumGlobalRows,
                 int& NumNonzeros,
                 int& NumProc,
                 int& NumVectors,
                 int& RowSize};
%typemap(in,numinputs=0) double& NORM (double temp = 0.0) { $1 = &temp; }
%typemap(argout) double& NORM
{
  $result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble(*$1));
}
%apply double& NORM {double& NormOne,
                     double& NormInf};

#ifdef HAVE_EPETRAEXT_HDF5
namespace EpetraExt
{
%extend HDF5
{
  %epetraext_primitive_read_method(int        , Int   , 0  )
  %epetraext_primitive_read_method(double     , Float , 0.0)
  %epetraext_primitive_read_method(std::string, String, "" )
  %epetraext_epetra_read_method(BlockMap   )
  %epetraext_epetra_read_method(Map        )
  %epetraext_epetra_read_method(MultiVector)
  %epetraext_epetra_read_method(CrsGraph   )
  %epetraext_epetra_read_method(CrsMatrix  )
  %epetraext_epetra_read_method(IntVector  )

  Teuchos::ParameterList ReadParameterList(const std::string & name)
  {
    Teuchos::ParameterList tpl;
    self->Read(name, tpl);
    return tpl;
  }
} // HDF5
}
%ignore EpetraExt::HDF5::Write(const string &, const string &, const int, const int, void *);
%ignore EpetraExt::HDF5::Write(const string &, const string &, int, int, int, const void*);
%ignore EpetraExt::HDF5::Write(const string &, const DistArray<int> &);
%ignore EpetraExt::HDF5::Write(const string &, const DistArray<double> &);
%ignore EpetraExt::HDF5::Write(const string &, const Handle &);
%ignore EpetraExt::HDF5::Read;
%ignore EpetraExt::HDF5::ReadIntDistArrayProperties;
%ignore EpetraExt::HDF5::ReadDoubleDistArrayProperties;
%ignore EpetraExt::HDF5::ReadHandleProperties;
%include "EpetraExt_HDF5.h"
#endif

/////////////////////////////////
// EpetraExt_XMLReader support //
/////////////////////////////////
%feature("docstring")
EpetraExt::XMLReader::ReadMap
"
Return a Map read from an XML file specified by filename 'name'.
"
%feature("docstring")
EpetraExt::XMLReader::ReadMultiVector
"
Return a MultiVector read from an XML file specified by filename
'name'.
"
%feature("docstring")
EpetraExt::XMLReader::ReadCrsGraph
"
Return a CrsGraph read from an XML file specified by filename 'name'.
"
%feature("docstring")
EpetraExt::XMLReader::ReadCrsMatrix
"
Return a CrsMatrix read from an XML file specified by filename 'name'.
"
%ignore EpetraExt::XMLReader::Read;
%include "EpetraExt_XMLReader.h"
namespace EpetraExt
{
  %extend XMLReader
  {
    %epetraext_epetra_read_method(Map        )
    %epetraext_epetra_read_method(MultiVector)
    %epetraext_epetra_read_method(CrsGraph   )
    %epetraext_epetra_read_method(CrsMatrix  )
  } // XMLReader
}

/////////////////////////////////
// EpetraExt_XMLWriter support //
/////////////////////////////////
%include "EpetraExt_XMLWriter.h"

/////////////////////////////////
// EpetraExt_Transform support //
/////////////////////////////////
%include "EpetraExt_Transform.h"
%template () std::vector<Epetra_IntVector>;
%template () EpetraExt::Transform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::Transform<Epetra_CrsGraph,
				  std::vector<Epetra_IntVector,
					      std::allocator<Epetra_IntVector> > >;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph,
					    std::vector<Epetra_IntVector> >;

///////////////////////////////////
// EpetraExt_MapColoring support //
///////////////////////////////////
%include "EpetraExt_MapColoring.h"

////////////////////////////////////////
// EpetraExt_MapColoringIndex support //
////////////////////////////////////////
%include "EpetraExt_MapColoringIndex.h"

/////////////////////////////////////
// EpetraExt_MultiVectorIn support //
/////////////////////////////////////
%feature("autodoc",
"MatrixMarketFileToMultiVector(str filename, Epetra.BlockMap) ->
    Epetra.MultiVector

Return a MultiVector read from a matix market file.")
EpetraExt::MatrixMarketFileToMultiVector;
%include "EpetraExt_MultiVectorIn.h"

//////////////////////////////////////
// EpetraExt_MultiVectorOut support //
//////////////////////////////////////
%include "EpetraExt_MultiVectorOut.h"

///////////////////////////////////
// EpetraExt_CrsMatrixIn support //
///////////////////////////////////
%feature("autodoc",
"MatlabFileToCrsMatrix(str filename, Epetra.Comm) -> Epetra.CrsMatrix

Return a CrsMatrix read from a matlab file.")
EpetraExt::MatlabFileToCrsMatrix;
%feature("autodoc",
"MatrixMarketFileToCrsMatrix(str filename, Epetra.Map rowMap, Epetra.Map
    colMap=None, Epetra.Map rangeMap=None, Epetra.Map domainMap=None) ->
    Epetra.CrsMatrix

Return a CrsMatrix read from a matrix market file.")
EpetraExt::MatrixMarketFileToCrsMatrix;
%include "EpetraExt_CrsMatrixIn.h"

////////////////////////////////////
// EpetraExt_RowMatrixOut support //
////////////////////////////////////
%include "EpetraExt_RowMatrixOut.h"

//////////////////////////////////
// EpetraExt_BlockMapIn support //
//////////////////////////////////
%feature("autodoc",
"MatrixMarketFileToBlockMap(str filename, Epetra.Comm) -> Epetra.BlockMap

Return a BlockMap read from a matrix market file.")
EpetraExt::MatrixMarketFileToBlockMap;
%feature("autodoc",
"MatrixMarketFileToBlockMaps(str filename, Epetra.Comm) ->
    (Epetra.BlockMap rowMap, Epetra.BlockMap colMap, Epetra.BlockMap rangeMap,
     Epetra.BlockMap domainMap)

Return a tuple of BlockMaps read from a matrix market file.  The
BlockMaps, listed in order, are the row map, the column map, the range
map and the domain map.")
EpetraExt::MatrixMarketFileToBlockMaps;
%feature("autodoc",
"MatrixMarketFileToMap(str filename, Epetra.Comm) -> Epetra.Map

Return a Map read from a matrix market file.")
EpetraExt::MatrixMarketFileToMap;
%include "EpetraExt_BlockMapIn.h"

///////////////////////////////////
// EpetraExt_BlockMapOut support //
///////////////////////////////////
%include "EpetraExt_BlockMapOut.h"

////////////////////////////////////////////
// EpetraExt.Add() and Multiply() support //
////////////////////////////////////////////
%feature("autodoc",
"Add(Epetra.CrsMatrix A, bool flag, float valA, Epetra.CrsMatrix B,
    float valB) -> int

Compute B <- valA * A + valB * B.  If flag is True, use the transpose
of A.  B must either have the structure of A+B or not yet have
FillComplete() called on it.")
EpetraExt::Add;
%feature("autodoc",
"Multiply(Epetra.CrsMatrix A, bool transposeA, Epetra.CrsMatrix B, bool
    transposeB, Epetra.CrsMatrix C) -> int

Compute C <- A * B, where transposeA and transposeB control the
transposition of A and B respectively.  C must have the structure of A
* B, or not yet have FillComplete() called on it.")
EpetraExt::Multiply;
%inline
%{
  namespace EpetraExt
  {
    int Add(Epetra_CrsMatrix& A, const bool flag, const double ValA,
            Epetra_CrsMatrix& B, const double ValB)
    {
      EpetraExt::MatrixMatrix M;
      return(M.Add(A, flag, ValA, B, ValB));
    }

    int Multiply(const Epetra_CrsMatrix& A, bool transposeA, 
                 const Epetra_CrsMatrix& B, bool transposeB, 
                 Epetra_CrsMatrix& C)
    {
      EpetraExt::MatrixMatrix M;
      return(M.Multiply(A, transposeA, B, transposeB, C));
    }
  }
%}

///////////////////////////////////////////
// EpetraExt ModelEvaluator support code //
///////////////////////////////////////////
//
// The EpetraExt::ModelEvaluator class is difficult to wrap because it
// has nested classes, and SWIG does not support those.  These classes
// are also resistent to the work-around presented in the SWIG manual
// because of the use of nested enumerations.  So the fix is to write
// python-only proxy classes that are not nested and converter
// typemaps so that the ModelEvaluator class can use them.

// We start by importing utility base classes
%pythoncode
%{
from PyTrilinos import Epetra
from PyTrilinos import PropertyBase
from PyTrilinos import typed_tuple
from PyTrilinos import tuple_of_int
from PyTrilinos import tuple_of_str
tuple_of_Vector = typed_tuple(Epetra.Vector)

%}

//////////////////
// InArgs class //
//////////////////
%pythoncode
%{
class InArgs(PropertyBase):
    """
    InArgs proxy class.

    This is a 'Property' class restricted to specific attributes that are
    type-checked. These properties are:

    description  - string description of associated ModelEvaluation
    x            - bool or Epetra.Vector: solution vector support.  If True, the
                   solver should allocate the vector.  If False, x is not
                   supported.  If a vector, the solver should use the
                   user-provided data.  (default False)
    x_dot        - bool or Epetra.Vector: time derivative of solution vector
                   support.  If True, the solver should allocate the vector.  If
                   False, x_dot is not supported.  If a vector, the solver
                   should use the user- provided data.  (default False)
    p            - int or tuple_of_Vector: VARIABLE P support.  If an int, the
                   solver should allocate an array of the given number of
                   vectors.  If 0, p is not supported.  If a tuple_of_Vector,
                   the solver should use the user-provided data.  (default 0)
    t            - float: time (default None)
    alpha        - float: VARIABLE ALPHA (default None)
    beta         - float: VARIABLE BETA (default None)
    """
    props = {'description' : str,
             'x'           : (bool, Epetra.Vector),
             #'x_poly'      : Teuchos.Polynomial,
             'x_dot'       : (bool, Epetra.Vector),
             #'x_dot_poly'  : Teuchos.Polynomial,
             'p'           : (int, tuple_of_Vector),
             't'           : float,
             'alpha'       : float,
             'beta'        : float
             }
    defaults = {'description' : 'WARNING!  THIS INARGS OBJECT IS UNINITALIZED!',
                'x'           : False,
                'x_dot'       : False,
                'p'           : 0
                }
    def __init__(self, **kwargs):
        PropertyBase.__init__(self, **kwargs)

%}

%typemap(out) EpetraExt::ModelEvaluator::InArgs
{
  $result = PyTrilinos::convertInArgsToPython($1);
  if (PyErr_Occurred()) SWIG_fail;
}

%typemap(directorin) const EpetraExt::ModelEvaluator::InArgs &
{
  $input = PyTrilinos::convertInArgsToPython($1_name);
  if (PyErr_Occurred()) throw PyTrilinos::PythonException();
}

%typemap(directorout) EpetraExt::ModelEvaluator::InArgs
{
  $result = EpetraExt::convertInArgsFromPython($1);
  if (PyErr_Occurred()) throw PyTrilinos::PythonException();
}

%typemap(in) const EpetraExt::ModelEvaluator::InArgs &
{
  *$1 = EpetraExt::convertInArgsFromPython($input);
  if (PyErr_Occurred()) SWIG_fail;
}

//////////////////////
// Evaluation class //
//////////////////////
%pythoncode
%{
class Evaluation(PropertyBase):
    """
    Evaluation< RCP<Epetra_Vector> > proxy class

    This is a 'Property' class restricted to specific attributes that are
    type-checked. These properties are:

    vector  - Epetra.Vector: (default None)
    type    - str: an enumeration limited to 'exact', 'approx_deriv', and
              'very_approx_deriv' (default None)
    """
    props = {'vector' : Epetra.Vector,
             'type'   : ('exact', 'approx_deriv', 'very_approx_deriv')
             }
    def __init__(self, **kwargs):
        PropertyBase.__init__(self, **kwargs)

tuple_of_Evaluation = typed_tuple(Evaluation)

%}

/////////////////////////////
// DerivativeSupport class //
/////////////////////////////
%pythoncode
%{
class DerivativeSupport(PropertyBase):
    """
    DerivativeSupport proxy class

    This is a 'Property' class restricted to specific attributes that are
    type-checked. These properties are:

    linearOp      - bool: True indicates that derivative is a linear operator
                    (default False)
    mVByCol       - bool: True indicates that derivative is a MultiVector stored
                    by column (defualt False)
    transMVByRow  - bool: True indicates that derivative is a transpose
                    MultiVector stored by row (default False)
    """
    props = {'linearOp'     : bool,
             'mVByCol'      : bool,
             'transMVByRow' : bool
             }
    defaults = {'linearOp'     : False,
                'mVByCol'      : False,
                'transMVByRow' : False
                }
    def __init__(self, **kwargs):
        PropertyBase.__init__(self, **kwargs)
    def none(self):
        noTrue = [True] * len(self.props)
        return (props.values() == noTrue)

tuple_of_DerivativeSupport = typed_tuple(DerivativeSupport)

%}

////////////////////////////////
// DerivativeProperties class //
////////////////////////////////
%pythoncode
%{
class DerivativeProperties(PropertyBase):
    """
    DerivativeProperties proxy class

    This is a 'Property' class restricted to specific attributes that are
    type-checked. These properties are:

    linearity        - str: an enumeration limited to 'unknown', 'const' and
                       'nonconst' (default 'unknown')
    rank             - str: an enumeration limited to 'unknown', 'full' and
                       'deficient' (default 'unknown')
    supportsAdjoint  - bool: True indicates that the adjoint is supported
                       (default False)
    """
    props = {'linearity'       : ('unknown', 'const', 'nonconst'),
             'rank'            : ('unknown', 'full', 'deficient'),
             'supportsAdjoint' : bool}
    defaults = {'linearity'       : 'unknown',
                'rank'            : 'unknown',
                'supportsAdjoint' : False}
    def __init__(self, **kwargs):
        PropertyBase.__init__(self, **kwargs)

tuple_of_DerivativeProperties = typed_tuple(DerivativeProperties)

%}

/////////////////////////////////
// DerivativeMultiVector class //
/////////////////////////////////
%pythoncode
%{
class DerivativeMultiVector(PropertyBase):
    """
    DerivativeMultiVector proxy class

    This is a 'Property' class restricted to specific attributes that are
    type-checked. These properties are:

    multiVector   - Epetra.MultiVector: (default None)
    orientation   - str: an enumeration limited to 'mv_by_col', and
                    'trans_mv_by_row' (default None)
    paramIndexes  - tuple_of_int: (default None)
    """
    props = {'multiVector'  : Epetra.MultiVector,
             'orientation'  : ('mv_by_col', 'trans_mv_by_row'),
             'paramIndexes' : tuple_of_int
             }
    def __init__(self, **kwargs):
        PropertyBase.__init__(self, **kwargs)

%}

//////////////////////
// Derivative class //
//////////////////////
%pythoncode
%{
class Derivative(PropertyBase):
    """
    Derivative proxy class

    This is a 'Property' class restricted to specific attributes that are
    type-checked. These properties are:

    operator               - Epetra.Operator (default None)
    derivativeMultiVector  - DerivativeMultiVector (default None)

    Only one or the other of these two attributes should be set, to indicate the
    nature of the derivitive evaluation.
    """
    props = {'operator'              : Epetra.Operator,
             'derivativeMultiVector' : DerivativeMultiVector
             }
    def __init__(self, **kwargs):
        PropertyBase.__init__(self, **kwargs)
    def isEmpty(self):
        return (self.operator is None) and (self.derivativeMultiVector is None)

tuple_of_Derivative = typed_tuple(Derivative)

%}

///////////////////
// OutArgs class //
///////////////////
%pythoncode
%{
class OutArgs(PropertyBase):
    """
    OutArgs proxy class

    This is a 'Property' class restricted to specific attributes that are
    type-checked. These properties are:

    description          - string description of associated ModelEvaluation
    g                    - (int, tuple_of_Evaluation): VARIABLE G support.  If
                           an int, the solver should allocate an array of the
                           given number of Evaluations.  If 0, g is not
                           supported.  If a tuple_of_Evaluation, the solver
                           should use the user-provided data.  (default 0)
    f                    - (bool, Evaluation): VARIABLE F support.  If True, the
                           solver should allocate the Evaluation.  If False, f
                           is not supported.  If an Evaluation, the solver
                           should use the user-supplied data.  (default False)
    W                    - (bool, Epetra.Operator): VARIABLE W support.  If True, the
                           solver should allocate the operator.  If False, W
                           is not supported.  If an operator, the solver
                           should use the user-supplied data.  (default False)
    W_properties         - DerivativeProperties: derivative properties for
                           VARIABLE W.  (default None)
    DfDp                 - (int, tuple_of_Derivative): VARIABLE DFDP support.  If
                           an int, the solver should allocate an array of the
                           given number of Derivatives.  If 0, DfDp is not
                           supported.  If a tuple_of_Derivative, the solver
                           should use the user-provided data.  (default 0)
    DfDp_properties      - tuple_of_DerivativeProperties: derivative properties
                           for VARIABLE DFDP.  (default None)
    DgDx                 - (int, tuple_of_Derivative): VARIABLE DGDX support.  If
                           an int, the solver should allocate an array of the
                           given number of Derivatives.  If 0, DgDx is not
                           supported.  If a tuple_of_Derivative, the solver
                           should use the user-provided data.  (default 0)
    DgDx_properties      - tuple_of_DerivativeProperties: derivative properties
                           for VARIABLE DGDX.  (default None)
    DgDx_dot             - (int, tuple_of_Derivative): VARIABLE DGDX_DOT support.
                           If an int, the solver should allocate an array of the
                           given number of Derivatives.  If 0, DgDx_dot is not
                           supported.  If a tuple_of_Derivative, the solver
                           should use the user-provided data.  (default 0)
    DgDx_dot_properties  - tuple_of_DerivativeProperties: derivative properties
                           for VARIABLE DGDX_DOT.  (default None)
    DgDp                 - (int, tuple_of_Derivative): VARIABLE DGDP support.  If
                           an int, the solver should allocate an array of the
                           given number of Derivatives.  If 0, DgDp is not
                           supported.  If a tuple_of_Derivative, the solver
                           should use the user-provided data.  (default 0)
    DgDp_properties      - tuple_of_DerivativeProperties: derivative properties
                           for VARIABLE DGDP.  (default None)
    """
    props = {'description'         : str,
             'g'                   : (int, tuple_of_Evaluation),
             'f'                   : (bool, Evaluation),
             'W'                   : (bool, Epetra.Operator),
             'W_properties'        : DerivativeProperties,
             #'f_poly'              : Teuchos.Polynomial,
             'DfDp'                : (int, tuple_of_Derivative),
             'DfDp_properties'     : tuple_of_DerivativeProperties,
             'DgDx'                : (int, tuple_of_Derivative),
             'DgDx_properties'     : tuple_of_DerivativeProperties,
             'DgDx_dot'            : (int, tuple_of_Derivative),
             'DgDx_dot_properties' : tuple_of_DerivativeProperties,
             'DgDp'                : (int, tuple_of_Derivative),
             'DgDp_properties'     : tuple_of_DerivativeProperties
             }
    defaults = {'description' : 'WARNING!  THIS OUTARGS OBJECT IS UNINITALIZED!',
                'g'           : 0,
                'f'           : False,
                'W'           : False,
                'DfDp'        : 0,
                'DgDx'        : 0,
                'DgDx_dot'    : 0,
                'DgDp'        : 0
                }
    def __init__(self, **kwargs):
        PropertyBase.__init__(self, **kwargs)

%}

%typemap(out) EpetraExt::ModelEvaluator::OutArgs
{
  $result = PyTrilinos::convertOutArgsToPython($1);
  if (PyErr_Occurred()) SWIG_fail;
}

%typemap(directorin) const EpetraExt::ModelEvaluator::OutArgs &
{
  $input = PyTrilinos::convertOutArgsToPython($1_name);
  if (PyErr_Occurred()) throw PyTrilinos::PythonException();
}

%typemap(directorout) EpetraExt::ModelEvaluator::OutArgs
{
  $result = EpetraExt::convertOutArgsFromPython($1);
  if (PyErr_Occurred()) throw PyTrilinos::PythonException();
}

%typemap(in) const EpetraExt::ModelEvaluator::OutArgs &
{
  *$1 = EpetraExt::convertOutArgsFromPython($input);
  if (PyErr_Occurred()) SWIG_fail;
}

//////////////////////////////////////
// EpetraExt ModelEvaluator support //
//////////////////////////////////////
//
// The EpetraExt::ModelEvaluator class is sufficiently complex,
// including nested classes, that it confuses the SWIG code parser (in
// part because SWIG does not support nested classes, it is a use-case
// SWIG developers do not test against).  Therefore, I provide here a
// stripped-down declaration of the class for wrapping purposes only,
// with enumerations and all but two nested classes removed.
//
%feature("director") EpetraExt::ModelEvaluator;

namespace EpetraExt {
class ModelEvaluator : virtual public Teuchos::Describable
{
public:
  class InArgs;
  class OutArgs;

  virtual ~ModelEvaluator();
  virtual Teuchos::RCP<	const Epetra_Map > get_x_map() const = 0;
  virtual Teuchos::RCP<	const Epetra_Map > get_f_map() const = 0;
  virtual Teuchos::RCP<	const Epetra_Map > get_p_map(int l) const;
  virtual Teuchos::RCP<	const Teuchos::Array<std::string> > get_p_names(int l) const;
  virtual Teuchos::RCP<	const Epetra_Map > get_g_map(int j) const;
  virtual Teuchos::RCP<	const Epetra_Vector > get_x_init() const;
  virtual Teuchos::RCP<	const Epetra_Vector > get_x_dot_init() const;
  virtual Teuchos::RCP<	const Epetra_Vector > get_p_init(int l) const;
  virtual double get_t_init() const;
  virtual double getInfBound() const;
  virtual Teuchos::RCP<	const Epetra_Vector > get_x_lower_bounds() const;
  virtual Teuchos::RCP<	const Epetra_Vector > get_x_upper_bounds() const;
  virtual Teuchos::RCP<	const Epetra_Vector > get_p_lower_bounds(int l) const;
  virtual Teuchos::RCP<	const Epetra_Vector > get_p_upper_bounds(int l) const;
  virtual double get_t_lower_bound() const;
  virtual double get_t_upper_bound() const;
  virtual Teuchos::RCP<	Epetra_Operator	> create_W() const;
  virtual Teuchos::RCP<	Epetra_Operator	> create_DfDp_op(int l) const;
  virtual Teuchos::RCP<	Epetra_Operator	> create_DgDx_dot_op(int j) const;
  virtual Teuchos::RCP<	Epetra_Operator	> create_DgDx_op(int j) const;
  virtual Teuchos::RCP<	Epetra_Operator	> create_DgDp_op( int j, int l ) const;
  virtual InArgs createInArgs() const = 0;
  virtual OutArgs createOutArgs() const = 0;
  virtual void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const = 0;
};
}

// Notes:
//
// * Teuchos::Polynomial is not yet wrapped, so the following have
//   been ignored:
//   + InArgs::x_poly
//   + InArgs::x_dot_poly
//   + OutArgs::f_poly
// * Evaluation should derive from Epetra.Vector, rather than simply
//   containing one.
// * Should DerivativeMultiVector derive from Epetra.MultiVector?
// * Allow derivative properties to be sparse by utilizing
//   dictionaries
// * Add documentation to Epetra_ModelEvaluator.h (based upon Thyra
//   ModelEvaluator documentation)

// Turn off the exception handling
%exception;

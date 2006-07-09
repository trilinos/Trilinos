// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//         PyTrilinos.EpetraExt: Python Interface to EpetraExt
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define EPETRAEXT_DOCSTRING
"The EpetraExt module allows access to The Trilinos package EpetraExt.  Note
that the 'EpetraExt_' prefix has been stripped from all EpetraExt objects,
but that if imported with 'from PyTrilinos import EpetraExt', these
objects exist in the 'EpetraExt' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.

The most important classes of the EpetraExt module are:
*) Graph coloring classes:
   - MapColoring
   - MapColoringIndex
*) Input functions:
   - BlockMapToMatrixMarketFile()
   - MatrixMarketFileToMultiVector()
   - MatrixMarketFileToBlockMap()
   - MatrixMarketFileToMap()
   - MatrixMarketFileToCrsMatrix()
*) Output functions:
   - BlockMapToMatrixMarketFile()
   - RowMatrixToMatrixMarketFile()
   - MultiVectorToMatrixMarketFile()
*) Input/Output classes:
   - HDF5
   - XML
*) Matrix-Matrix functions:
   - Add()
   - Multiply()
"
%enddef

%module(package="PyTrilinos", directors="1", docstring=EPETRAEXT_DOCSTRING) EpetraExt

%{
// System includes
#include <vector>

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_LocalMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"

// Epetra python includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_JadMatrix.h"

// EpetraExt includes
#include "EpetraExt_Version.h"
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

#include "EpetraExt_HDF5.h" // FIXME: memory management still scary...
#include "EpetraExt_XMLReader.h" // FIXME: memory management still scary...
#include "EpetraExt_XMLWriter.h" // FIXME: memory management still scary...

namespace EpetraExt {
  int MatrixMarketFileToCrsMatrix(const char *filename, Epetra_CrsMatrix* A)
  {
    return(MatrixMarketFileToCrsMatrixHandle(filename, A));
  }
  int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_BlockMap & map, 
                               Epetra_CrsMatrix * & OutCrsMatrix)
  {
    // MS // this works out only for serial runs
    if (map.Comm().NumProc() != 1)
      EPETRA_CHK_ERR(-2);
    if (map.NumGlobalElements() != map.NumGlobalPoints())
      EPETRA_CHK_ERR(-1);
    Epetra_Map map2(map.NumGlobalElements(), map.IndexBase(), map.Comm());
    return(MatrixMarketFileToCrsMatrix(filename, map2, OutCrsMatrix));
  }
}
%}

// Ignore directives
// %ignore Epetra_CrsGraph::operator[](int);
// %ignore Epetra_CrsGraph::operator[](int) const;
// %ignore Epetra_CrsGraph::operator=(const Epetra_CrsGraph &);
// %ignore Epetra_IntVector::operator=(const Epetra_IntVector &);
// %ignore Epetra_IntVector::operator[](int);
// %ignore Epetra_IntVector::operator[](int) const;
// %ignore Epetra_MapColoring::operator[](int);
// %ignore Epetra_MapColoring::operator[](int) const;

using namespace std;
// C++ STL support
%include "std_string.i"
%include "std_vector.i"

// Epetra interface import
%import "Epetra.i"

// Typemaps
%typemap(argout) Epetra_BlockMap*& OutBlockMap {
  PyObject *o1, *oBlockMap;
  oBlockMap = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_BlockMap, 1);
  if (!PyTuple_Check($result)) $result = Py_BuildValue("(O)", $result);
  if (result >= 0)
    o1 = Py_BuildValue("(O)", oBlockMap);
  else
    o1 = Py_BuildValue("(i)", 0);
  $result = PySequence_Concat($result,o1);
}

%typemap(in,numinputs=0) Epetra_BlockMap *&OutBlockMap(Epetra_BlockMap* _BlockMap) {
  $1 = &_BlockMap;
}

%typemap(argout) Epetra_Map*& OutMap {
  PyObject *o1, *oMap;
  oMap = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_Map, 1);
  if (!PyTuple_Check($result)) $result = Py_BuildValue("(O)", $result);
  if (result >= 0)
    o1 = Py_BuildValue("(O)", oMap);
  else
    o1 = Py_BuildValue("(i)", 0);
  $result = PySequence_Concat($result,o1);
}

%typemap(in,numinputs=0) Epetra_Map *&OutMap(Epetra_Map* _Map) {
  $1 = &_Map;
}

%typemap(argout) Epetra_MultiVector*& OutVector {
  PyObject *o1, *oVector;
  Epetra_NumPyMultiVector * npmv;
  static swig_type_info *ty = SWIG_TypeQuery("Epetra_NumPyMultiVector *");
  npmv = new Epetra_NumPyMultiVector(**$1);
  oVector = SWIG_NewPointerObj(npmv, ty, 1);
  if (!PyTuple_Check($result)) $result = Py_BuildValue("(O)", $result);
  if (!result) {
    o1 = Py_BuildValue("(O)", oVector);}
  else {
    o1 = Py_BuildValue("(i)", 0);
  }
  $result = PySequence_Concat($result, o1);
}

%typemap(in,numinputs=0) Epetra_MultiVector *&OutVector(Epetra_MultiVector* _Vector) {
  $1 = &_Vector;
}

%typemap(argout) Epetra_Vector*& OutVector {
  PyObject *o1, *oVector;
  oVector = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_Vector, 1);
  if (!PyTuple_Check($result)) $result = Py_BuildValue("(O)", $result);
  if (!result) 
    o1 = Py_BuildValue("(O)",oVector);
  else
    o1 = Py_BuildValue("(i)", 0);
  $result = PySequence_Concat($result, o1);
}

%typemap(in,numinputs=0) Epetra_Vector *&OutVector(Epetra_Vector* _Vector) {
  $1 = &_Vector;
}

%typemap(argout) Epetra_CrsMatrix*& OutCrsMatrix {
  PyObject *o1, *oCrsMatrix;
  oCrsMatrix = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_CrsMatrix, 1);
  if (!PyTuple_Check($result)) $result = Py_BuildValue("(O)", $result);
  if (!result) 
    o1 = Py_BuildValue("(O)", oCrsMatrix);
  else
    o1 = Py_BuildValue("(i)", 0);
  $result = PySequence_Concat($result, o1);
}

%typemap(in,numinputs=0) Epetra_CrsMatrix *&OutCrsMatrix(Epetra_CrsMatrix* _CrsMatrix) {
  $1 = &_CrsMatrix;
}

// EpetraExt interface includes
%include "EpetraExt_Version.h"
%include "EpetraExt_Transform.h"
%template () std::vector<Epetra_IntVector>;
%template () EpetraExt::Transform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::Transform<Epetra_CrsGraph, std::vector<Epetra_IntVector,
							       std::allocator<Epetra_IntVector> > >;
//%template () EpetraExt::Transform<Epetra_CrsGraph, std::vector<Epetra_IntVector> >;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, std::vector<Epetra_IntVector> >;

%include "EpetraExt_MapColoring.h"
%include "EpetraExt_MapColoringIndex.h"
%include "EpetraExt_HDF5.h"
%include "EpetraExt_XMLReader.h"
%include "EpetraExt_XMLWriter.h"

namespace EpetraExt 
{
  int VectorToMatrixMarketFile(const char *filename, const Epetra_Vector & A,
                               const char * matrixName=0,
                               const char *matrixDescription=0,
                               bool writeHeader=true);
  int VectorToMatlabFile(const char *filename, const Epetra_Vector & A);

  int MultiVectorToMatrixMarketFile(const char *filename, const Epetra_MultiVector & A,
                                    const char * matrixName=0,
                                    const char *matrixDescription=0,
                                    bool writeHeader=true);
  int MultiVectorToMatlabFile(const char *filename, const Epetra_MultiVector & A);

  int RowMatrixToMatrixMarketFile(const char *filename, const Epetra_RowMatrix & A,
                                  const char * matrixName=0,
                                  const char *matrixDescription=0,
                                  bool writeHeader=true);
  int RowMatrixToMatlabFile(const char *filename, const Epetra_RowMatrix & A);

  int BlockMapToMatrixMarketFile(const char *filename, const Epetra_BlockMap & blockMap,
                                 const char * mapName=0,
                                 const char *mapDescription=0,
                                 bool writeHeader=true);

  int MatrixMarketFileToBlockMap(const char *filename, const Epetra_Comm &
                                 comm, Epetra_BlockMap * & OutBlockMap);

  int MatrixMarketFileToMap(const char *filename, const Epetra_Comm &
                                 comm, Epetra_Map * & OutMap);

  int MatrixMarketFileToMultiVector(const char *filename, const Epetra_BlockMap & map, 
                               Epetra_MultiVector * & OutVector);

  int MatrixMarketFileToVector(const char *filename, const Epetra_BlockMap & map, 
                               Epetra_Vector * & OutVector);

  int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_BlockMap & map, 
                               Epetra_CrsMatrix * & OutCrsMatrix);
  int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_Map & map, 
                               Epetra_CrsMatrix * & OutCrsMatrix);

}

%extend EpetraExt::HDF5 
{
  Epetra_BlockMap* ReadBlockMap(string Name)
  {
    Epetra_BlockMap* obj = 0;
    self->Read(Name, obj);
    return(obj);
  }

  Epetra_Map* ReadMap(string Name)
  {
    Epetra_Map* obj = 0;
    self->Read(Name, obj);
    return(obj);
  }

  Epetra_MultiVector* ReadMultiVector(string Name)
  {
    Epetra_MultiVector* obj = 0;
    self->Read(Name, obj);
    return(obj);
  }

  Epetra_CrsGraph* ReadCrsGraph(string Name)
  {
    Epetra_CrsGraph* obj = 0;
    self->Read(Name, obj);
    return(obj);
  }

  Epetra_CrsMatrix* ReadCrsMatrix(string Name)
  {
    Epetra_CrsMatrix* obj = 0;
    self->Read(Name, obj);
    return(obj);
  }

  Epetra_IntVector* ReadIntVector(string Name)
  {
    Epetra_IntVector* obj = 0;
    self->Read(Name, obj);
    return(obj);
  }
}

%extend EpetraExt::XMLReader 
{
  Epetra_Map* ReadMap(string Name)
  {
    Epetra_Map* obj = 0;
    try {
      self->Read(Name, obj);
    } catch(...) {
      cout << "Caught generic exception, maybe the specified class does not exist" << endl;
    }
    return(obj);
  }

  Epetra_MultiVector* ReadMultiVector(string Name)
  {
    Epetra_MultiVector* obj = 0;
    try {
      self->Read(Name, obj);
    } catch(...) {
      cout << "Caught generic exception, maybe the specified class does not exist" << endl;
    }
    return(obj);
  }

  Epetra_CrsGraph* ReadCrsGraph(string Name)
  {
    Epetra_CrsGraph* obj = 0;
    try {
      self->Read(Name, obj);
    } catch(...) {
      cout << "Caught generic exception, maybe the specified class does not exist" << endl;
    }
    return(obj);
  }

  Epetra_CrsMatrix* ReadCrsMatrix(string Name)
  {
    Epetra_CrsMatrix* obj = 0;
    try {
      self->Read(Name, obj);
    } catch(...) {
      cout << "Caught generic exception, maybe the specified class does not exist" << endl;
    }
    return(obj);
  }
}

%inline %{
  namespace EpetraExt
  {
    int Add(Epetra_CrsMatrix& A, const bool flag, const double ValA,
            Epetra_CrsMatrix& B, const double ValB)
    {
      EpetraExt::MatrixMatrix M;
      return(M.Add(A, flag, ValA, B, ValB));
    }
  }
%}

// Python code.
%pythoncode %{
  __version__ = EpetraExt_Version().split()[2]
%}

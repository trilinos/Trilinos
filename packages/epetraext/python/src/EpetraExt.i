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
*) Graph coloring:
   - MapColoring
   - MapColoringIndex
*) Input functions:
   - BlockMapToMatrixMarketFile()
   - MatrixMarketFileToMultiVector()
   - MatrixMarketFileToBlockMap()
   - MatrixMarketFileToCrsMatrix()
*) Output functions:
   - BlockMapToMatrixMarketFile()
   - RowMatrixToMatrixMarketFile()
   - MultiVectorToMatrixMarketFile()
*) Matrix-Matrix operations:
   - Add()
   - Multiply()
"
%enddef

%module(package="PyTrilinos", docstring=EPETRAEXT_DOCSTRING) EpetraExt

%{
// System includes
#include <vector>

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"
#include "Epetra_NumPyVector.h"

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
%ignore Epetra_CrsGraph::operator[](int);
%ignore Epetra_CrsGraph::operator[](int) const;
%ignore Epetra_CrsGraph::operator=(const Epetra_CrsGraph &);
%ignore Epetra_IntVector::operator=(const Epetra_IntVector &);
%ignore Epetra_IntVector::operator[](int);
%ignore Epetra_IntVector::operator[](int) const;
%ignore Epetra_MapColoring::operator[](int);
%ignore Epetra_MapColoring::operator[](int) const;

using namespace std;
// C++ STL support
%include "std_string.i"
%include "std_vector.i"

// Epetra interface import
%import "Epetra.i"

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

  %typemap(argout) (Epetra_BlockMap*& OutBlockMap)
  {
    PyObject *o3, *oBlockMap;
    oBlockMap = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_BlockMap, 1);

    if (!PyTuple_Check($result)) {
      PyObject *o2 = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($target,0,o2);
    }

    o3 = PyTuple_New(1);
    // add matrix
    if (result >= 0) 
      PyTuple_SetItem(o3,0,oBlockMap);
    else
      PyTuple_SetItem(o3,0,PyInt_FromLong(0));

    $result = PySequence_Concat($result,o3);
  }

  %typemap(in,numinputs=0) Epetra_BlockMap *&OutBlockMap(Epetra_BlockMap* _BlockMap) {
    $1 = &_BlockMap;
  }
  int MatrixMarketFileToBlockMap(const char *filename, const Epetra_Comm &
                                 comm, Epetra_BlockMap * & OutBlockMap);

  %typemap(argout) (Epetra_MultiVector*& OutVector)
  {
    PyObject *o3, *oVector;
    oVector = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_MultiVector, 1);

    if (!PyTuple_Check($result)) {
      PyObject *o2 = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($target,0,o2);
    }

    o3 = PyTuple_New(1);
    // add vector to the tuple
    if (!result) 
      PyTuple_SetItem(o3,0,oVector);
    else
      PyTuple_SetItem(o3,0,PyInt_FromLong(0));
    $result = PySequence_Concat($result,o3);
  }

  %typemap(in,numinputs=0) Epetra_MultiVector *&OutVector(Epetra_MultiVector* _Vector) {
    $1 = &_Vector;
  }
  int MatrixMarketFileToMultiVector(const char *filename, const Epetra_BlockMap & map, 
                               Epetra_MultiVector * & OutVector);

  %typemap(argout) (Epetra_Vector*& OutVector)
  {
    PyObject *o3, *oVector;
    oVector = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_Vector, 1);

    if (!PyTuple_Check($result)) {
      PyObject *o2 = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($target,0,o2);
    }

    o3 = PyTuple_New(1);
    // add vector to the tuple
    if (!result) 
      PyTuple_SetItem(o3,0,oVector);
    else
      PyTuple_SetItem(o3,0,PyInt_FromLong(0));
    $result = PySequence_Concat($result,o3);
  }

  %typemap(in,numinputs=0) Epetra_Vector *&OutVector(Epetra_Vector* _Vector) {
    $1 = &_Vector;
  }
  int MatrixMarketFileToVector(const char *filename, const Epetra_BlockMap & map, 
                               Epetra_Vector * & OutVector);

  %typemap(argout) (Epetra_CrsMatrix*& OutCrsMatrix)
  {
    PyObject *o3, *oCrsMatrix;
    oCrsMatrix = SWIG_NewPointerObj((void*)(*$1), SWIGTYPE_p_Epetra_CrsMatrix, 1);

    if (!PyTuple_Check($result)) {
      PyObject *o2 = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($target,0,o2);
    }

    o3 = PyTuple_New(1);
    // add vector to the tuple
    if (!result) 
      PyTuple_SetItem(o3,0,oCrsMatrix);
    else
      PyTuple_SetItem(o3,0,PyInt_FromLong(0));
    $result = PySequence_Concat($result,o3);
  }

  %typemap(in,numinputs=0) Epetra_CrsMatrix *&OutCrsMatrix(Epetra_CrsMatrix* _CrsMatrix) {
    $1 = &_CrsMatrix;
  }
  int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_BlockMap & map, 
                               Epetra_CrsMatrix * & OutCrsMatrix);
  int MatrixMarketFileToCrsMatrix(const char *filename, const Epetra_Map & map, 
                               Epetra_CrsMatrix * & OutCrsMatrix);

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

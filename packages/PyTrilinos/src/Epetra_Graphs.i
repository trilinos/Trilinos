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

%{
#include "Epetra_CrsGraph.h"
#include "Epetra_OffsetIndex.h"
%}

////////////////////////
// Typemap directives //
////////////////////////

// Apply the numpy.i (int DIM1, int* IN_ARRAY1) typemap to
// the (int NumIndices, int* Indices) argument signature
%apply (int DIM1, int* IN_ARRAY1) {(int NumIndices, int* Indices)};

// Begin argout typemap collection for (int & NumIndices, int *& Indices)
%typecheck(SWIG_TYPECHECK_INT32_ARRAY) (int & NumIndices, int *& Indices) {
  $1 = ($input != 0);
}
%typemap(in,numinputs=0) (int & NumIndices, int *& Indices) { }
%typemap(argout)         (int & NumIndices, int *& Indices)
{
  Py_XDECREF($result);
  if (result == -1) SWIG_exception(SWIG_ValueError,   "Invalid row index"  );
  if (result == -2) SWIG_exception(SWIG_RuntimeError, "Graph not completed");
  intp dims[ ] = { *$1 };
  $result = PyArray_SimpleNewFromData(1,dims,PyArray_INT,(void*)(*$2));
  if ($result == NULL) SWIG_exception(SWIG_RuntimeError, "Error creating integer array");
}
// End argout typemap collection for (int & NumIndices, int * Indices)

/////////////////////////////
// Epetra_CrsGraph support //
/////////////////////////////
%ignore Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess,
					 const Epetra_BlockMap &,
					 const int *,
					 bool);
%ignore Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess,
					 const Epetra_BlockMap &,
					 const Epetra_BlockMap &,
					 const int *,
					 bool);
%ignore Epetra_CrsGraph::ExtractGlobalRowCopy(int, int, int&, int*) const;
%ignore Epetra_CrsGraph::ExtractMyRowCopy(int, int, int&, int*) const;
%ignore Epetra_CrsGraph::ExtractGlobalRowView(int, int&, int*&) const;
%ignore Epetra_CrsGraph::ExtractMyRowView(int, int&, int*&) const;
%rename(CrsGraph) Epetra_CrsGraph;
%epetra_exception(Epetra_CrsGraph,Epetra_CrsGraph)
%epetra_exception(Epetra_CrsGraph,OptimizeStorage)
%include "Epetra_CrsGraph.h"
%extend Epetra_CrsGraph {

  Epetra_CrsGraph(Epetra_DataAccess       CV,
		  const Epetra_BlockMap & rowMap,
		  PyObject              * numIndicesList,
		  bool                    staticProfile=false) {
    // Declarations
    PyArrayObject   * numIndicesArray    = NULL;
    int             * numIndicesPerRow   = NULL;
    Epetra_CrsGraph * returnCrsGraph     = NULL;
    int               constIndicesPerRow = 0;
    int               listSize           = 0;
    int               is_new             = 0;

    if (PyInt_Check(numIndicesList)) {
      constIndicesPerRow = (int) PyInt_AsLong(numIndicesList);
      returnCrsGraph     = new Epetra_CrsGraph(CV,rowMap,constIndicesPerRow,staticProfile);
    } else {
      numIndicesArray = obj_to_array_contiguous_allow_conversion(numIndicesList,
								 PyArray_INT, &is_new);
      if (!numIndicesArray || !require_dimensions(numIndicesArray,1)) goto fail;
      numIndicesPerRow = (int*) (numIndicesArray->data);
      listSize = (int) array_size(numIndicesArray,0);
      if (listSize != rowMap.NumMyElements()) {
	PyErr_Format(PyExc_ValueError,
		     "Row map has %d elements, list of number of indices has %d",
		     rowMap.NumMyElements(), listSize);
	goto fail;
      }
      returnCrsGraph = new Epetra_CrsGraph(CV,rowMap,numIndicesPerRow,staticProfile);
      if (is_new) Py_DECREF(numIndicesArray);
    }
    return returnCrsGraph;

  fail:
    if (is_new) Py_XDECREF(numIndicesArray);
    return NULL;
  }

  Epetra_CrsGraph(Epetra_DataAccess       CV,
		  const Epetra_BlockMap & rowMap,
		  const Epetra_BlockMap & colMap,
		  PyObject              * numIndicesList,
		  bool                    staticProfile=false) {
    // Declarations
    PyArrayObject   * numIndicesArray    = NULL;
    int             * numIndicesPerRow   = NULL;
    Epetra_CrsGraph * returnCrsGraph     = NULL;
    int               constIndicesPerRow = 0;
    int               listSize           = 0;
    int               is_new             = 0;

    if (PyInt_Check(numIndicesList)) {
      constIndicesPerRow = (int) PyInt_AsLong(numIndicesList);
      returnCrsGraph     = new Epetra_CrsGraph(CV,rowMap,colMap,constIndicesPerRow,staticProfile);
    } else {
      numIndicesArray = obj_to_array_contiguous_allow_conversion(numIndicesList,
								 PyArray_INT, &is_new);
      if (!numIndicesArray || !require_dimensions(numIndicesArray,1)) goto fail;
      numIndicesPerRow = (int*) (numIndicesArray->data);
      listSize = (int) array_size(numIndicesArray,0);
      if (listSize != rowMap.NumMyElements()) {
	PyErr_Format(PyExc_ValueError,
		     "Row map has %d elements, list of number of indices has %d",
		     rowMap.NumMyElements(), listSize);
	goto fail;
      }
      returnCrsGraph = new Epetra_CrsGraph(CV,rowMap,colMap,numIndicesPerRow,staticProfile);
      if (is_new) Py_DECREF(numIndicesArray);
    }
    return returnCrsGraph;

  fail:
    if (is_new) Py_XDECREF(numIndicesArray);
    return NULL;
  }

  PyObject * ExtractGlobalRowCopy(int globalRow) const {
    int        lrid          = 0;
    int        numIndices    = 0;
    int        result        = 0;
    intp       dimensions[ ] = { 0 };
    int      * indices       = NULL;
    PyObject * indicesArray  = NULL;

    lrid = self->LRID(globalRow);
    if (lrid == -1) {
      PyErr_Format(PyExc_ValueError, "Invalid global row index = %d", globalRow);
      goto fail;
    }
    dimensions[0] = self->NumMyIndices(lrid);
    indicesArray  = PyArray_SimpleNew(1,dimensions,PyArray_INT);
    indices       = (int *) ((PyArrayObject *)indicesArray)->data;
    result        = self->ExtractGlobalRowCopy(globalRow, dimensions[0], numIndices, indices);
    if (result == -2) {
      PyErr_SetString(PyExc_RuntimeError, "Graph not completed");
      goto fail;
    }
    return PyArray_Return((PyArrayObject*)indicesArray);
  fail:
    Py_XDECREF(indicesArray);
    return NULL;
  }

  int * __getitem__(int i) {
    return self->operator[](i);
  }

  PyObject * ExtractMyRowCopy(int localRow) const {
    int        numIndices    = 0;
    int        result        = 0;
    intp       dimensions[ ] = { 0 };
    int      * indices       = NULL;
    PyObject * indicesArray  = NULL;

    if (localRow < 0 || localRow >= self->NumMyRows()) {
      PyErr_Format(PyExc_ValueError, "Invalid local row index = %d", localRow);
      goto fail;
    }
    dimensions[0] = self->NumMyIndices(localRow);
    indicesArray  = PyArray_SimpleNew(1,dimensions,PyArray_INT);
    indices       = (int *) ((PyArrayObject *)indicesArray)->data;
    result        = self->ExtractMyRowCopy(localRow, dimensions[0], numIndices, indices);
    if (result == -2) {
      PyErr_SetString(PyExc_RuntimeError, "Graph not completed");
      goto fail;
    }
    return PyArray_Return((PyArrayObject*)indicesArray);
  fail:
    Py_XDECREF(indicesArray);
    return NULL;
  }

  int * __getitem__(int i) {
    return self->operator[](i);
  }
}

////////////////////////////////
// Epetra_OffsetIndex support //
////////////////////////////////
%rename(OffsetIndex) Epetra_OffsetIndex;
%include "Epetra_OffsetIndex.h"

///////////////////////////////
// Epetra_FECrsGraph support //
///////////////////////////////
// ** Epetra_FECrsGraph is apparently not built **
//#include "Epetra_FECrsGraph.h"
//%rename(FECrsGraph) Epetra_FECrsGraph;
//%include "Epetra_FECrsGraph.h"

// Clear the typemaps
%clear (int   NumIndices, int *  Indices);
%clear (int & NumIndices, int *& Indices);

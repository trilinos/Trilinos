// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//            PyTrilinos.Epetra: Python Interface to Epetra
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

%{
#include "Epetra_CrsGraph.h"
%}

// Rename directives
%rename(CrsGraph) Epetra_CrsGraph;

// Local interface includes
%include "NumPyArray.h"

// Typemap directives
%typemap(in) (int * Indices)
{
  PyArrayObject * array = NumPyWrapper::contiguous_typed_array($input, PyArray_INT,
							       0, 0);
  assert (array && "Function should have already checked this.");
  $1 = (int *) array->data;
}

%typemap(argout) (int & NumIndices, int *& Indices)
{
  // Decrement the current result object
  Py_XDECREF($result);
  // Check for bad result
  if (result == -1)
    PyErr_SetString(PyExc_RuntimeError, "Invalid row index");
  if (result == -2)
    PyErr_SetString(PyExc_RuntimeError, "Graph not completed");
  // Copy the indices into a python tuple
  $result = PyTuple_New(*$1);
  for (int i=0; i<*$1; ++i)
    PyTuple_SetItem($result, i, PyInt_FromLong((long) *$2[i])); 
}

// Include directives
%include "Epetra_CrsGraph.h"

// Extend directives
%extend Epetra_CrsGraph {
  int * __getitem__(int i) {
    return self->operator[](i);
  }
}

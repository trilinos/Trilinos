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

// Begin argout typemap collection for (int & NumIndices, int *& Indices)
%typecheck(SWIG_TYPECHECK_INT32_ARRAY) (int & NumIndices, int *& Indices)
{
  $1 = ($input != 0);
}
%typemap(in,numinputs=0) (int & NumIndices, int *& Indices)
{ }
%typemap(argout)         (int & NumIndices, int *& Indices)
{
  Py_XDECREF($result);
  if (result == -1) SWIG_exception(SWIG_ValueError,   "Invalid row index"  );
  if (result == -2) SWIG_exception(SWIG_RuntimeError, "Graph not completed");
  npy_intp dims[ ] = { *$1 };
  $result = PyArray_SimpleNewFromData(1,dims,NPY_INT,(void*)(*$2));
  if ($result == NULL) SWIG_exception(SWIG_RuntimeError, "Error creating integer array");
}
// End argout typemap collection for (int & NumIndices, int * Indices)

/////////////////////////////
// Epetra_CrsGraph support //
/////////////////////////////
%teuchos_rcp(Epetra_CrsGraph)
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, int numIndicesPerRow,
    bool staticProfile=False) -> CrsGraph

  Constructor with implicit column map and constant indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    numIndicesPerRow  - Integer number of indices per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 int, bool);
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, PySequence
    numIndicesPerRow, bool staticProfile=False) -> CrsGraph

  Constructor with implicit column map and variable indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    numIndicesPerRow  - Sequence of integers representing the number of indices
                        per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 int, const int*, bool);
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    int numIndicesPerRow, bool staticProfile=False) -> CrsGraph

  Constructor with specified column map and constant indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    colMap            - Map describing distribution of columns across processors
    numIndicesPerRow  - Integer number of indices per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 const Epetra_BlockMap&, int, bool);
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    PySequence numIndicesPerRow, bool staticProfile=False) -> CrsGraph

  Constructor with specified column map and variable indices per row.
  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - Map describing distribution of rows across processors
    colMap            - Map describing distribution of columns across processors
    numIndicesPerRow  - Sequence of integers representing the number of indices
                        per row
    staticProfile     - Static profile flag")
Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess, const Epetra_BlockMap&,
				 const Epetra_BlockMap&, int, const int*, bool);
%feature("autodoc",
"
__init__(self, CrsGraph graph) -> CrsGraph

  Copy constructor.  Arguments:

    graph - Source graph for copy constructor")
Epetra_CrsGraph::Epetra_CrsGraph(const Epetra_CrsGraph&);
%feature("autodoc",
"InsertGlobalIndices(self, int globalRow, PySequence indices) -> int

Insert a sequence of global indices into the set of nonzero columns
for the specified global row.  Argument indices can be a numpy array
of integers or any python sequence that can be converted to a numpy
array of integers.  The integers represent global IDs that are to be
inserted into the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::InsertGlobalIndices;
%feature("autodoc",
"RemoveGlobalIndices(self, int globalRow, PySequence indices) -> int

Remove a sequence of global indices from the set of nonzero columns
for the specified global row.  Argument indices can be a numpy array
of integers or any python sequence that can be converted to a numpy
array of integers.  The integers represent global IDs that are to be
removed from the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::RemoveGlobalIndices(int, int, int*);
%feature("autodoc",
"InsertMyIndices(self, int localRow, PySequence indices) -> int

Insert a sequence of local indices into the set of nonzero columns for
the specified local row.  Argument indices can be a numpy array of
integers or any python sequence that can be converted to a numpy array
of integers.  The integers represent local IDs that are to be inserted
into the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::InsertMyIndices;
%feature("autodoc",
"RemoveMyIndices(self, int localRow, PySequence indices) -> int

Remove a sequence of local indices from the set of nonzero columns for
the specified local row.  Argument indices can be a numpy array of
integers or any python sequence that can be converted to a numpy array
of integers.  The integers represent local IDs that are to be removed
from the graph.  An integer error/warning code is returned.
")
Epetra_CrsGraph::RemoveMyIndices(int, int, int*);
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
%apply (int DIM1, int* IN_ARRAY1) {(int NumRows,    const int * NumIndicesPerRow)};
%apply (int DIM1, int* IN_ARRAY1) {(int NumIndices, int       * Indices         )};
%extend Epetra_CrsGraph
{
  Epetra_CrsGraph(Epetra_DataAccess       CV,
		  const Epetra_BlockMap & RowMap,
		  int                     NumRows,
		  const int             * NumIndicesPerRow,
		  bool                    StaticProfile=false)
  {
    if (NumRows != RowMap.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
		   "Row map has %d elements, NumIndicesPerRow has %d",
		   RowMap.NumMyElements(), NumRows);
      return NULL;
    }
    return new Epetra_CrsGraph(CV, RowMap, NumIndicesPerRow, StaticProfile);
  }

  Epetra_CrsGraph(Epetra_DataAccess       CV,
		  const Epetra_BlockMap & RowMap,
		  const Epetra_BlockMap & ColMap,
		  int                     NumRows,
		  const int             * NumIndicesPerRow,
		  bool                    StaticProfile=false)
  {
    if (NumRows != RowMap.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
		   "Row map has %d elements, NumIndicesPerRow has %d",
		   RowMap.NumMyElements(), NumRows);
      return NULL;
    }
    return new Epetra_CrsGraph(CV, RowMap, ColMap, NumIndicesPerRow, StaticProfile);
  }

  PyObject * ExtractGlobalRowCopy(int globalRow) const
  {
    int        lrid          = 0;
    int        numIndices    = 0;
    int        result        = 0;
    npy_intp   dimensions[ ] = { 0 };
    int      * indices       = NULL;
    PyObject * indicesArray  = NULL;

    lrid = self->LRID(globalRow);
    if (lrid == -1)
    {
      PyErr_Format(PyExc_ValueError, "Invalid global row index = %d", globalRow);
      goto fail;
    }
    dimensions[0] = self->NumMyIndices(lrid);
    indicesArray  = PyArray_SimpleNew(1,dimensions,NPY_INT);
    if (indicesArray == NULL) goto fail;
    indices = (int*) array_data(indicesArray);
    result  = self->ExtractGlobalRowCopy(globalRow, dimensions[0], numIndices, indices);
    if (result == -2)
    {
      PyErr_SetString(PyExc_RuntimeError, "Graph not completed");
      goto fail;
    }
    return PyArray_Return((PyArrayObject*)indicesArray);

  fail:
    Py_XDECREF(indicesArray);
    return NULL;
  }

  int * __getitem__(int i)
  {
    return self->operator[](i);
  }

  PyObject * ExtractMyRowCopy(int localRow) const
  {
    int        numIndices    = 0;
    int        result        = 0;
    npy_intp   dimensions[ ] = { 0 };
    int      * indices       = NULL;
    PyObject * indicesArray  = NULL;

    if (localRow < 0 || localRow >= self->NumMyRows())
    {
      PyErr_Format(PyExc_ValueError, "Invalid local row index = %d", localRow);
      goto fail;
    }
    dimensions[0] = self->NumMyIndices(localRow);
    indicesArray  = PyArray_SimpleNew(1,dimensions,NPY_INT);
    if (indicesArray == NULL) goto fail;
    indices = (int*) array_data(indicesArray);
    result  = self->ExtractMyRowCopy(localRow, dimensions[0], numIndices, indices);
    if (result == -2)
    {
      PyErr_SetString(PyExc_RuntimeError, "Graph not completed");
      goto fail;
    }
    return PyArray_Return((PyArrayObject*)indicesArray);

  fail:
    Py_XDECREF(indicesArray);
    return NULL;
  }

  int * __getitem__(int i)
  {
    return self->operator[](i);
  }
}
%include "Epetra_CrsGraph.h"
%clear (const int * NumIndicesPerRow, int    NumRows);
%clear (int         NumIndices,       int *  Indices);
%clear (int       & NumIndices,       int *& Indices);

////////////////////////////////
// Epetra_OffsetIndex support //
////////////////////////////////
%teuchos_rcp(Epetra_OffsetIndex)
%rename(OffsetIndex) Epetra_OffsetIndex;
%include "Epetra_OffsetIndex.h"

///////////////////////////////
// Epetra_FECrsGraph support //
///////////////////////////////
// ** Epetra_FECrsGraph is apparently not built **
//#include "Epetra_FECrsGraph.h"
//%teuchos_rcp(Epetra_FECrsGraph)
//%rename(FECrsGraph) Epetra_FECrsGraph;
//%include "Epetra_FECrsGraph.h"

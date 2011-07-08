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
// Epetra includes
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_MapColoring.h"
%}

//////////////////////////////////
// Epetra_SrcDistObject support //
//////////////////////////////////
%teuchos_rcp(Epetra_SrcDistObject)
%rename(SrcDistObject) Epetra_SrcDistObject;
%include "Epetra_SrcDistObject.h"

///////////////////////////////
// Epetra_DistObject support //
///////////////////////////////
%teuchos_rcp(Epetra_DistObject)
%rename(DistObject) Epetra_DistObject;
%include "Epetra_DistObject.h"

////////////////////////////////
// Epetra_MapColoring support //
////////////////////////////////
%teuchos_rcp(Epetra_MapColoring)
%rename(MapColoring) Epetra_MapColoring;
%ignore Epetra_MapColoring::Epetra_MapColoring(const Epetra_BlockMap &, int*, const int);
%ignore Epetra_MapColoring::operator()(int) const;
%ignore Epetra_MapColoring::ListOfColors() const;
%ignore Epetra_MapColoring::ColorLIDList(int) const;
%ignore Epetra_MapColoring::ElementColors() const;
%apply (int DIM1, int* IN_ARRAY1) {(int numColors, int* elementColors)};
%extend Epetra_MapColoring
{
  Epetra_MapColoring(const Epetra_BlockMap & map,
		     int numColors, int* elementColors,
		     const int defaultColor=0)
  {
    Epetra_MapColoring * mapColoring;
    if (numColors != map.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
		   "Epetra.BlockMap has %d elements, while elementColors has %d",
		   map.NumMyElements(), numColors);
      goto fail;
    }
    mapColoring = new Epetra_MapColoring(map, elementColors, defaultColor);
    return mapColoring;
  fail:
    return NULL;
  }

  int __getitem__(int i)
  {
    return self->operator[](i);
  }

  void __setitem__(int i, int color)
  {
    self->operator[](i) = color;
  }

  PyObject * ListOfColors()
  {
    int      * list    = self->ListOfColors();
    npy_intp   dims[ ] = { self->NumColors() };
    int      * data;
    PyObject * retObj  = PyArray_SimpleNew(1,dims,NPY_INT);
    if (retObj == NULL) goto fail;
    data = (int*) array_data(retObj);
    for (int i = 0; i<dims[0]; i++) data[i] = list[i];
    return PyArray_Return((PyArrayObject*)retObj);
  fail:
    Py_XDECREF(retObj);
    return NULL;
  }

  PyObject * ColorLIDList(int color)
  {
    int      * list    = self->ColorLIDList(color);
    npy_intp   dims[ ] = { self->NumElementsWithColor(color) };
    int      * data;
    PyObject * retObj  = PyArray_SimpleNew(1,dims,NPY_INT);
    if (retObj == NULL) goto fail;
    data = (int*) array_data(retObj);
    for (int i = 0; i<dims[0]; i++) data[i] = list[i];
    return PyArray_Return((PyArrayObject*)retObj);
  fail:
    Py_XDECREF(retObj);
    return NULL;
  }

  PyObject * ElementColors()
  {
    int      * list    = self->ElementColors();
    npy_intp   dims[ ] = { self->Map().NumMyElements() };
    int      * data;
    PyObject * retObj  = PyArray_SimpleNew(1,dims,NPY_INT);
    if (retObj == NULL) goto fail;
    data = (int*) array_data(retObj);
    for (int i = 0; i<dims[0]; i++) data[i] = list[i];
    return PyArray_Return((PyArrayObject*)retObj);
  fail:
    Py_XDECREF(retObj);
    return NULL;
  }
}
%include "Epetra_MapColoring.h"
%clear (int numColors, int* elementColors);

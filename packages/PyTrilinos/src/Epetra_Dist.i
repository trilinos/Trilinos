// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%{
// Epetra include files
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_MapColoring.h"
%}

/////////////////////////////////////////////////////////
// Teuchos::RCP<> support for all classes in this file //
/////////////////////////////////////////////////////////
%teuchos_rcp(Epetra_SrcDistObject)
%teuchos_rcp(Epetra_DistObject)
%teuchos_rcp(Epetra_MapColoring)

//////////////////////////////////
// Epetra_SrcDistObject support //
//////////////////////////////////
%rename(SrcDistObject) Epetra_SrcDistObject;
%include "Epetra_SrcDistObject.h"

///////////////////////////////
// Epetra_DistObject support //
///////////////////////////////
%rename(DistObject) Epetra_DistObject;
%include "Epetra_DistObject.h"

////////////////////////////////
// Epetra_MapColoring support //
////////////////////////////////
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

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
// Epetra includes
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
%}

// Ignore directives
// These are replaced with extensions
%ignore Epetra_BlockMap::Epetra_BlockMap(int,int,int*,int, int,const Epetra_Comm&);
%ignore Epetra_BlockMap::Epetra_BlockMap(int,int,int*,int*,int,const Epetra_Comm&);

// Rename directives
%rename(BlockMap) Epetra_BlockMap;
%rename(Map     ) Epetra_Map;
%rename(LocalMap) Epetra_LocalMap;

// Include directives
%include "Epetra_BlockMap.h"
%include "Epetra_Map.h"
%include "Epetra_LocalMap.h"

// Extend directives
%extend Epetra_BlockMap {

  Epetra_BlockMap(int                 numGlobalElements,
		  PyObject *          myGlobalElementArray,
		  int                 elementSize,
		  int                 indexBase,
		  const Epetra_Comm & comm                ) {
    // The wrapper for this constructor cannot be called because
    // (under the SWIG-generated logic) it is masked by the next
    // constructor; However, it does serve a documentation purpose
    // [via %feature("autodoc")], so I keep it as a shell.
    return NULL;
  }

  Epetra_BlockMap(int                 numGlobalElements,
		  PyObject *          myGlobalElementArray,
		  PyObject *          myElementSizes,
		  int                 indexBase,
		  const Epetra_Comm & comm                ) {
    // Declarations
    int              numMyElements;
    int              numMyElementSizes;
    int             *myGlobalElements = NULL;
    int             *elementSizeList  = NULL;
    PyObject        *elementArray     = NULL;
    PyObject        *elementSizeArray = NULL;
    Epetra_BlockMap *returnBlockMap   = NULL;
    // Check for integer PyObjects in the argument list
    if (PyInt_Check(myElementSizes)) {
      int elementSize = (int) PyInt_AsLong(myElementSizes);
      if (PyInt_Check(myGlobalElementArray)) {
	numMyElements   = (int) PyInt_AsLong(myGlobalElementArray);
	// Constructor for user-defined, linear distribution of constant-size elements
	returnBlockMap = new Epetra_BlockMap(numGlobalElements,numMyElements,elementSize,
					     indexBase,comm);
      } else {
	elementArray = PyArray_ContiguousFromObject(myGlobalElementArray,'i',1,1);
	if (elementArray == NULL) goto fail;
	numMyElements    = ((PyArrayObject*)elementArray)->dimensions[0];
	myGlobalElements = (int *) (((PyArrayObject*)elementArray)->data);
	// Constructor for user-defined, arbitrary distribution of constant-size elements
	returnBlockMap = new Epetra_BlockMap(numGlobalElements,numMyElements,myGlobalElements,
					     elementSize,indexBase,comm);
	Py_DECREF(elementArray);
      }
    } else {
      // Obtain a Numeric element array and check
      elementArray = PyArray_ContiguousFromObject(myGlobalElementArray,'i',1,1);
      if (elementArray == NULL) goto fail;
      numMyElements    = ((PyArrayObject*)elementArray)->dimensions[0];
      myGlobalElements = (int *) (((PyArrayObject*)elementArray)->data);
      // Obtain a Numric element size array and check
      elementSizeArray = PyArray_ContiguousFromObject(myElementSizes,'i',1,1);
      if (elementArray == NULL) goto fail;
      numMyElementSizes = ((PyArrayObject*)elementSizeArray)->dimensions[0];
      if (numMyElements != numMyElementSizes) {
	PyErr_Format(PyExc_ValueError,
		     "Element and element size arrays must have same lengths\n"
		     "Lengths = %d, %d", numMyElements, numMyElementSizes);
	goto fail;
      }
      elementSizeList = (int *) (((PyArrayObject*)elementSizeArray)->data);
      // Obtain a new Epetra_BlockMap
      returnBlockMap = new Epetra_BlockMap(numGlobalElements,numMyElements,myGlobalElements,
					   elementSizeList,indexBase,comm);
      Py_DECREF(elementArray);
      Py_DECREF(elementSizeArray);
    }
    return returnBlockMap;
  fail:
    Py_XDECREF(elementArray);
    Py_XDECREF(elementSizeArray);
    return NULL;
  }

}


%extend Epetra_Map {
  Epetra_Map(const int NumGlobalElements,
             const Epetra_IntSerialDenseVector& MyGlobalElements,
             const int IndexBase, const Epetra_Comm& Comm)
  {
    return(new Epetra_Map(NumGlobalElements, MyGlobalElements.Length(),
                         (int*)MyGlobalElements.Values(), IndexBase, Comm));
  }

  Epetra_Map(const int NumGlobalElements,
             PyObject* MyGlobalElements, const int IndexBase,
             const Epetra_Comm& Comm)
  {
    if (PyList_Check(MyGlobalElements) == 0)
    {
      cerr << "Input object is not a list" << endl;
      return NULL;
    }

    int len = PyList_Size(MyGlobalElements);

    vector<int> list(len);

    for (int i = 0 ; i < len ; ++i)
    {
      PyObject* Index;
      Index = PyList_GetItem(MyGlobalElements, i);

      if (PyInt_Check(Index) == 0)
      {
        cerr << "Indices must be integers" << endl;
        return NULL;
      }

      list[i] = PyLong_AsLong(Index);
    }
    return(new Epetra_Map(NumGlobalElements, len, &list[0], IndexBase, Comm));
  }

  PyObject*  MyGlobalElements()
  {
    int* MyGlobalElements_Epetra = self->MyGlobalElements();
    PyObject* MyGlobalElements_Python,* item;
    int size = self->NumMyElements();
    if (size <= 0)
      goto fail;

    MyGlobalElements_Python = PyList_New(size);

    for (int i = 0 ; i < size ; ++i)
    {
      item = PyInt_FromLong(MyGlobalElements_Epetra[i]);
      PyList_SetItem(MyGlobalElements_Python, i, item);
    }

    return(MyGlobalElements_Python);
fail:
    Py_INCREF(Py_None);
    return Py_None;
  }
}

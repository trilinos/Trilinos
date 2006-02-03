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
%ignore Epetra_BlockMap::RemoteIDList(int,const int *,int*,int*) const;
%ignore Epetra_BlockMap::RemoteIDList(int,const int *,int*,int*,int*) const;
%ignore Epetra_BlockMap::FindLocalElementID(int,int&,int&) const;
%ignore Epetra_BlockMap::MyGlobalElements(int*) const;
%ignore Epetra_BlockMap::MyGlobalElements() const;
%ignore Epetra_BlockMap::FirstPointInElementList(int*) const;
%ignore Epetra_BlockMap::FirstPointInElementList() const;
%ignore Epetra_BlockMap::ElementSizeList(int*) const;
%ignore Epetra_BlockMap::ElementSizeList() const;
%ignore Epetra_BlockMap::PointToElementList(int*) const;
%ignore Epetra_BlockMap::PointToElementList() const;
%ignore Epetra_Map::Epetra_Map(int,int,int*,int,const Epetra_Comm&);
// These are expert methods not wrapped
%ignore Epetra_BlockMap::ReferenceCount() const;
%ignore Epetra_BlockMap::DataPtr() const;

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
    // The wrapper for this constructor cannot be called from python
    // because (under the SWIG-generated logic) it is masked by the
    // next constructor; However, it does serve a documentation
    // purpose [via %feature("autodoc")], so I keep it as a shell.
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

  PyObject * RemoteIDList(PyObject * GIDList) {
    int        numIDs[1];
    int        result;
    int      * GIDData 	 = NULL;
    int      * PIDData 	 = NULL;
    int      * LIDData 	 = NULL;
    int      * sizeData  = NULL;
    PyObject * GIDArray  = NULL;
    PyObject * PIDArray  = NULL;
    PyObject * LIDArray  = NULL;
    PyObject * sizeArray = NULL;
    PyObject * returnObj = NULL;
    GIDArray = PyArray_ContiguousFromObject(GIDList,'i',1,1);
    if (GIDArray == NULL) goto fail;
    numIDs[0] = ((PyArrayObject*)GIDArray)->dimensions[0];
    PIDArray  = PyArray_FromDims(1,numIDs,PyArray_INT);
    if (PIDArray == NULL) goto fail;
    LIDArray  = PyArray_FromDims(1,numIDs,PyArray_INT);
    if (LIDArray == NULL) goto fail;
    sizeArray = PyArray_FromDims(1,numIDs,PyArray_INT);
    if (sizeArray == NULL) goto fail;
    GIDData  = (int*) (((PyArrayObject*)GIDArray)->data);
    PIDData  = (int*) (((PyArrayObject*)PIDArray)->data);
    LIDData  = (int*) (((PyArrayObject*)LIDArray)->data);
    sizeData = (int*) (((PyArrayObject*)sizeArray)->data);
    result   = self->RemoteIDList(numIDs[0],GIDData,PIDData,LIDData,sizeData);
    if (result != 0) {
      PyErr_Format(PyExc_RuntimeError,"Bad RemoteIDList return code = %d", result);
      goto fail;
    }
    returnObj = Py_BuildValue("(OOO)",PIDArray,LIDArray,sizeArray);
    Py_DECREF(GIDArray );
    return returnObj;

  fail:
    Py_XDECREF(GIDArray );
    Py_XDECREF(PIDArray );
    Py_XDECREF(LIDArray );
    Py_XDECREF(sizeArray);
    return NULL;
  }

  PyObject * FindLocalElementID(int pointID) {
    int result;
    int elementID;
    int elementOffset;
    result   = self->FindLocalElementID(pointID,elementID,elementOffset);
    if (result != 0) {
      PyErr_Format(PyExc_RuntimeError,"Bad FindLocalElementID return code = %d", result);
      goto fail;
    }
    return Py_BuildValue("(ii)",elementID,elementOffset);

  fail:
    return NULL;
  }

  PyObject * MyGlobalElements() {
    int        numEls[1];
    int        result;
    int      * geData  = NULL;
    PyObject * geArray = NULL;
    numEls[0] = self->NumMyElements();
    geArray   = PyArray_FromDims(1,numEls,PyArray_INT);
    if (geArray == NULL) goto fail;
    geData = (int*) (((PyArrayObject*)geArray)->data);
    result = self->MyGlobalElements(geData);
    if (result != 0) {
      PyErr_Format(PyExc_RuntimeError,"Bad MyGlobalElements return code = %d", result);
      goto fail;
    }
    return geArray;

  fail:
    Py_XDECREF(geArray );
    return NULL;
  }

  PyObject * FirstPointInElementList() {
    int        numEls[1];
    int        result;
    int      * fpeData  = NULL;
    PyObject * fpeArray = NULL;
    numEls[0] = self->NumMyElements();
    fpeArray  = PyArray_FromDims(1,numEls,PyArray_INT);
    if (fpeArray == NULL) goto fail;
    fpeData = (int*) (((PyArrayObject*)fpeArray)->data);
    result = self->FirstPointInElementList(fpeData);
    if (result != 0) {
      PyErr_Format(PyExc_RuntimeError,"Bad FirstPointInElementList return code = %d", result);
      goto fail;
    }
    return fpeArray;

  fail:
    Py_XDECREF(fpeArray );
    return NULL;
  }

  PyObject * ElementSizeList() {
    int        numEls[1];
    int        result;
    int      * eslData  = NULL;
    PyObject * eslArray = NULL;
    numEls[0] = self->NumMyElements();
    eslArray  = PyArray_FromDims(1,numEls,PyArray_INT);
    if (eslArray == NULL) goto fail;
    eslData = (int*) (((PyArrayObject*)eslArray)->data);
    result = self->ElementSizeList(eslData);
    if (result != 0) {
      PyErr_Format(PyExc_RuntimeError,"Bad ElementSizeList return code = %d", result);
      goto fail;
    }
    return eslArray;

  fail:
    Py_XDECREF(eslArray );
    return NULL;
  }

  PyObject * PointToElementList() {
    int        numPts[1];
    int        result;
    int      * pteData  = NULL;
    PyObject * pteArray = NULL;
    numPts[0] = self->NumMyPoints();
    pteArray  = PyArray_FromDims(1,numPts,PyArray_INT);
    if (pteArray == NULL) goto fail;
    pteData = (int*) (((PyArrayObject*)pteArray)->data);
    result = self->PointToElementList(pteData);
    if (result != 0) {
      PyErr_Format(PyExc_RuntimeError,"Bad PointToElementList return code = %d", result);
      goto fail;
    }
    return pteArray;

  fail:
    Py_XDECREF(pteArray );
    return NULL;
  }

}


%extend Epetra_Map {

  Epetra_Map(int                 numGlobalElements,
	     PyObject          * myGlobalElementArray,
	     int                 indexBase,
	     const Epetra_Comm & comm) {
    // Declarations
    int          numMyElements;
    int        * myGlobalElements = NULL;
    PyObject   * elementArray     = NULL;
    Epetra_Map * returnMap        = NULL;
    // Check for integer PyObject in the argument list
    if (PyInt_Check(myGlobalElementArray)) {
      numMyElements = (int) PyInt_AsLong(myGlobalElementArray);
      returnMap = new Epetra_Map(numGlobalElements,numMyElements,indexBase,comm);
    } else {
      // Obtain a Numeric element array and check
      elementArray = PyArray_ContiguousFromObject(myGlobalElementArray,'i',1,1);
      if (elementArray == NULL) goto fail;
      numMyElements    = ((PyArrayObject*)elementArray)->dimensions[0];
      myGlobalElements = (int *) (((PyArrayObject*)elementArray)->data);
      returnMap = new Epetra_Map(numGlobalElements,numMyElements,myGlobalElements,
				 indexBase,comm);
      Py_DECREF(elementArray);
    }
    return returnMap;
  fail:
    Py_XDECREF(elementArray);
    return NULL;
  }

}

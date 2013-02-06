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
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Directory.h"
#include "Epetra_BasicDirectory.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
%}

// General ignore directives
%ignore *::PermuteFromLIDs() const;
%ignore *::PermuteToLIDs() const;
%ignore *::RemoteLIDs() const;
%ignore *::ExportLIDs() const;
%ignore *::ExportPIDs() const;

/////////////////////////////
// Epetra_BlockMap support //
/////////////////////////////
%teuchos_rcp(Epetra_BlockMap)
%teuchos_rcp_epetra_argout(Epetra_BlockMap)
%feature("autodoc",
"
__init__(self, int numGlobalElements, int elementSize, int indexBase,
     Comm comm) -> BlockMap

BlockMap constructor with implicit local elements and constant element
size.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     elementSize        - The number of degrees of freedom associated
                          with every element.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, int, const Epetra_Comm &);
%feature("autodoc",
"
__init__(self, int numGlobalElements, int numMyElements, int elementSize,
     int indexBase, Comm comm) -> BlockMap

BlockMap constructor with specified number of local elements and
constant element size.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     numMyElements      - Number of local elements on this processor.
     elementSize        - The number of degrees of freedom associated
                          with every element.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, int, int, const Epetra_Comm &);
%feature("autodoc",
"
__init__(self, int numGlobalElements, PySequence myGlobalElements,
     int elementSize, int indexBase, Comm comm) -> BlockMap

BlockMap constructor with specified list of local elements and
constant element size.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     myGlobalElements   - A sequence of integers specifying the global
                          element indexes on this processor.
     elementSize        - The number of degrees of freedom associated
                          with every element.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, const int*, int, int,
				 const Epetra_Comm &);
%feature("autodoc",
"
__init__(self, int numGlobalElements, PySequence myGlobalElements,
     PySequence elementsSizes, int indexBase, Comm comm) -> BlockMap

BlockMap constructor with specified list of local elements and
specified list of element sizes.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     myGlobalElements   - A sequence of integers specifying the global
                          element indexes on this processor.
     elementSizes       - A sequence of integers specifying the number of
                          degrees of freedom associated with each element
                          on this processor.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_BlockMap::Epetra_BlockMap(int, int, const int*, int, const int*, int,
				 const Epetra_Comm &);
%feature("autodoc",
"
__init__(self, BlockMap map) -> BlockMap

BlockMap copy constructor.")
Epetra_BlockMap::Epetra_BlockMap(const Epetra_BlockMap &);
%feature("docstring")
Epetra_BlockMap::RemoteIDList
"``GIDList`` is a sequence of integer global IDs, and the return
argument is the three-tuple ``(PIDList, LIDList, sizeList)``, which
are ``numpy.ndarray`` objects of integers representing the processor
IDs, local IDs and element sizes, respectively."
%feature("docstring")
Epetra_BlockMap::FindLocalElementID
"Returns a tuple containing the local ID of the element that contains
the given local pointID, and the offset of the point in that element."
%feature("docstring")
Epetra_BlockMap::MyGlobalElements
"Returns a numpy array of integers specifying the list of global IDs on
the processor."
%feature("docstring")
Epetra_BlockMap::FirstPointInElementList
"Returns a numpy array of integer first local point numbers for all of
the local elements."
%feature("docstring")
Epetra_BlockMap::ElementSizeList
"Returns a numpy array of integer sizes for each local element."
%feature("docstring")
Epetra_BlockMap::PointToElementList
"Returns a numpy array of integers such that for each local point, it
indicates the local element ID that the point belongs to."
%ignore Epetra_BlockMap::Epetra_BlockMap(int,int,const int*,const int*,int,
					 const Epetra_Comm&);
%ignore Epetra_BlockMap::RemoteIDList(int,const int*,int*,int*) const;
%ignore Epetra_BlockMap::RemoteIDList(int,const int*,int*,int*,int*) const;
%ignore Epetra_BlockMap::FindLocalElementID(int,int&,int&) const;
%ignore Epetra_BlockMap::MyGlobalElements() const;
%ignore Epetra_BlockMap::MyGlobalElements(int*) const;
%ignore Epetra_BlockMap::MyGlobalElements(int const *&,long long const *&) const;
%ignore Epetra_BlockMap::FirstPointInElementList(int*) const;
%ignore Epetra_BlockMap::FirstPointInElementList() const;
%ignore Epetra_BlockMap::ElementSizeList(int*) const;
%ignore Epetra_BlockMap::ElementSizeList() const;
%ignore Epetra_BlockMap::PointToElementList(int*) const;
%ignore Epetra_BlockMap::PointToElementList() const;
%ignore Epetra_BlockMap::ReferenceCount() const;
%ignore Epetra_BlockMap::DataPtr() const;
%rename(BlockMap) Epetra_BlockMap;
%apply (int DIM1, int * IN_ARRAY1) {(int NumMyElements, const int * MyGlobalElements)};
%apply (int DIM1, int * IN_ARRAY1) {(int NumElements,   const int * ElementSizeList )};
%extend Epetra_BlockMap
{
  Epetra_BlockMap(int                 NumGlobalElements,
		  int                 NumMyElements,
		  const int         * MyGlobalElements,
		  int                 NumElements,
		  const int         * ElementSizeList,
		  int                 IndexBase,
		  const Epetra_Comm & Comm              )
  {
    if (NumMyElements != NumElements)
    {
      PyErr_Format(PyExc_ValueError,
		   "MyGlobalElements and ElementSizeList must have same lengths\n"
		   "Lengths = %d, %d", NumMyElements, NumElements);
      return NULL;
    }
    return new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,
			       ElementSizeList, IndexBase, Comm);
  }

  PyObject * RemoteIDList(PyObject * GIDList)
  {
    npy_intp        numIDs[1];
    int             result;
    int           * GIDData   = NULL;
    int           * PIDData   = NULL;
    int           * LIDData   = NULL;
    int           * sizeData  = NULL;
    PyArrayObject * GIDArray  = NULL;
    PyArrayObject * PIDArray  = NULL;
    PyArrayObject * LIDArray  = NULL;
    PyArrayObject * sizeArray = NULL;
    PyObject      * returnObj = NULL;
    int             is_new    = 0;

    GIDArray = obj_to_array_contiguous_allow_conversion(GIDList, NPY_INT, &is_new);
    if (!GIDArray || !require_dimensions(GIDArray,1)) goto fail;
    numIDs[0] = (int) array_size(GIDArray,0);
    PIDArray  = (PyArrayObject*) PyArray_SimpleNew(1,numIDs,NPY_INT);
    if (PIDArray == NULL) goto fail;
    LIDArray  = (PyArrayObject*) PyArray_SimpleNew(1,numIDs,NPY_INT);
    if (LIDArray == NULL) goto fail;
    sizeArray = (PyArrayObject*) PyArray_SimpleNew(1,numIDs,NPY_INT);
    if (sizeArray == NULL) goto fail;
    GIDData  = (int*) array_data(GIDArray);
    PIDData  = (int*) array_data(PIDArray);
    LIDData  = (int*) array_data(LIDArray);
    sizeData = (int*) array_data(sizeArray);
    result   = self->RemoteIDList(numIDs[0],GIDData,PIDData,LIDData,sizeData);
    if (result != 0)
    {
      PyErr_Format(PyExc_RuntimeError,"Bad RemoteIDList return code = %d", result);
      goto fail;
    }
    returnObj = Py_BuildValue("(OOO)",PIDArray,LIDArray,sizeArray);
    if (is_new) { Py_DECREF(GIDArray); }
    return returnObj;

  fail:
    if (is_new) { Py_XDECREF(GIDArray); }
    Py_XDECREF(PIDArray );
    Py_XDECREF(LIDArray );
    Py_XDECREF(sizeArray);
    return NULL;
  }

  PyObject * FindLocalElementID(int pointID)
  {
    int result;
    int elementID;
    int elementOffset;
    result = self->FindLocalElementID(pointID,elementID,elementOffset);
    if (result != 0)
    {
      PyErr_Format(PyExc_RuntimeError,"Bad FindLocalElementID return code = %d", result);
      goto fail;
    }
    return Py_BuildValue("(ii)",elementID,elementOffset);

  fail:
    return NULL;
  }

  PyObject * MyGlobalElements()
  {
    npy_intp   numEls[1];
    int        result;
    int      * geData  = NULL;
    PyObject * geArray = NULL;
    numEls[0] = self->NumMyElements();
    geArray   = PyArray_SimpleNew(1,numEls,NPY_INT);
    if (geArray == NULL) goto fail;
    geData = (int*) array_data(geArray);
    result = self->MyGlobalElements(geData);
    if (result != 0)
    {
      PyErr_Format(PyExc_RuntimeError,"Bad MyGlobalElements return code = %d", result);
      goto fail;
    }
    return PyArray_Return((PyArrayObject*)geArray);

  fail:
    Py_XDECREF(geArray);
    return NULL;
  }

  PyObject * FirstPointInElementList()
  {
    npy_intp   numEls[1];
    int      * result;
    int      * fpeData  = NULL;
    PyObject * fpeArray = NULL;
    numEls[0] = self->NumMyElements();
    fpeArray  = PyArray_SimpleNew(1,numEls,NPY_INT);
    if (fpeArray == NULL) goto fail;
    fpeData = (int*) array_data(fpeArray);
    result  = self->FirstPointInElementList();
    for (int i=0; i<numEls[0]; ++i) fpeData[i] = result[i];
    return PyArray_Return((PyArrayObject*)fpeArray);
  fail:
    Py_XDECREF(fpeArray);
    return NULL;
  }

  PyObject * ElementSizeList()
  {
    npy_intp   numEls[1];
    int        result;
    int      * eslData  = NULL;
    PyObject * eslArray = NULL;
    numEls[0] = self->NumMyElements();
    eslArray  = PyArray_SimpleNew(1,numEls,NPY_INT);
    if (eslArray == NULL) goto fail;
    eslData = (int*) array_data(eslArray);
    result = self->ElementSizeList(eslData);
    if (result != 0)
    {
      PyErr_Format(PyExc_RuntimeError,"Bad ElementSizeList return code = %d", result);
      goto fail;
    }
    return PyArray_Return((PyArrayObject*)eslArray);

  fail:
    Py_XDECREF(eslArray );
    return NULL;
  }

  PyObject * PointToElementList()
  {
    npy_intp   numPts[1];
    int        result;
    int      * pteData  = NULL;
    PyObject * pteArray = NULL;
    numPts[0] = self->NumMyPoints();
    pteArray  = PyArray_SimpleNew(1,numPts,NPY_INT);
    if (pteArray == NULL) goto fail;
    pteData = (int*) array_data(pteArray);
    result = self->PointToElementList(pteData);
    if (result != 0)
    {
      PyErr_Format(PyExc_RuntimeError,"Bad PointToElementList return code = %d", result);
      goto fail;
    }
    return PyArray_Return((PyArrayObject*)pteArray);

  fail:
    Py_XDECREF(pteArray );
    return NULL;
  }
}
%include "Epetra_BlockMap.h"
%clear (int NumElements, const int * ElementSizeList);

////////////////////////
// Epetra_Map support //
////////////////////////
%teuchos_rcp(Epetra_Map)
%teuchos_rcp_epetra_argout(Epetra_Map)
%feature("autodoc",
"
__init__(self, int numGlobalElements, int indexBase, Comm comm) -> Map

Map constructor with implicit number of elements per processor.
Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_Map::Epetra_Map(int, int, const Epetra_Comm &);
%feature("autodoc",
"
__init__(self, int numGlobalElements, int numMyElements, int indexBase,
     Comm comm) -> Map

Map constructor with specified number of elements per processor.
Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     numMyElements      - Number of local elements on this processor.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_Map::Epetra_Map(int, int, int, const Epetra_Comm &);
%feature("autodoc",
"
__init__(self, int numGlobalElements, PySequence myGlobalElements,
     int indexBase, Comm comm) -> Map

Map constructor with specified list of global element IDs for each
processor.  Arguments are:

     numGlobalElements  - Total number of elements over all processors.
                          Specify -1 to have the constructor compute
                          the number of global elements
     myGlobalElements   - A sequence of integers specifying the global
                          element indexes on this processor.
     indexBase          - The base integer value for indexed array
                          references.  Typically this is 0 for C/C++ and 1
                          for Fortran, but it can be set to any integer
                          value.
     comm               - The Epetra.Comm communicator. This communicator
                          can in turn be queried for processor rank and
                          size information.")
Epetra_Map::Epetra_Map(int, int, const int*, int, const Epetra_Comm &);
%feature("autodoc",
"
__init__(self, Map map) -> Map

Map copy constructor.")
Epetra_Map::Epetra_Map(const Epetra_Map &);
%rename(Map) Epetra_Map;
%include "Epetra_Map.h"
%clear (int NumMyElements, const int * MyGlobalElements);

/////////////////////////////
// Epetra_LocalMap support //
/////////////////////////////
%teuchos_rcp(Epetra_LocalMap)
%teuchos_rcp_epetra_argout(Epetra_LocalMap)
%rename(LocalMap) Epetra_LocalMap;
%include "Epetra_LocalMap.h"

//////////////////////////////
// Epetra_Directory support //
//////////////////////////////
%rename(Directory) Epetra_Directory;
%include "Epetra_Directory.h"

///////////////////////////////////
// Epetra_BasicDirectory support //
///////////////////////////////////
%rename(BasicDirectory) Epetra_BasicDirectory;
%include "Epetra_BasicDirectory.h"

///////////////////////////
// Import/Export support //
///////////////////////////

// Import/Export extensions are done with a couple of nested macros
%define %epetra_mover_method(methodName, numMethod)
PyObject * methodName()
{
  npy_intp   numIDs[ ]   = { self->numMethod() };
  int      * ids         = NULL;
  int      * returnData  = NULL;
  PyObject * returnArray = PyArray_SimpleNew(1,numIDs,NPY_INT);
  if (returnArray == NULL) goto fail;
  ids        = self->methodName();
  returnData = (int*) array_data(returnArray);
  for (int i=0; i<numIDs[0]; i++) returnData[i] = ids[i];
  return PyArray_Return((PyArrayObject*)returnArray);
  fail:
  return NULL;
}
%enddef

%define %epetra_mover_class(CLASS)
%extend Epetra_##CLASS
{
  %epetra_mover_method(PermuteFromLIDs,	NumPermuteIDs)
  %epetra_mover_method(PermuteToLIDs,   NumPermuteIDs)
  %epetra_mover_method(RemoteLIDs,      NumRemoteIDs )
  %epetra_mover_method(ExportLIDs,     	NumExportIDs )
  %epetra_mover_method(ExportPIDs,     	NumExportIDs )
}
%enddef

///////////////////////////
// Epetra_Import support //
///////////////////////////
%teuchos_rcp(Epetra_Import)
%rename(Import) Epetra_Import;
%include "Epetra_Import.h"
%epetra_mover_class(Import)

///////////////////////////
// Epetra_Export support //
///////////////////////////
%teuchos_rcp(Epetra_Export)
%rename(Export) Epetra_Export;
%include "Epetra_Export.h"
%epetra_mover_class(Export)

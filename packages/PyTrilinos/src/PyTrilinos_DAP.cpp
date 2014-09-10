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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// Include the PyTrilinos DistArray Protocol header
#include "PyTrilinos_DAP.hpp"

// Python exception handling
#include "PyTrilinos_PythonException.h"

// Verbosity
// #define PYTRILINOS_DAP_VERBOSE

////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

DimensionDictionary::DimensionDictionary(PyObject * dim_dict)
{
  // Check that dim_dict is a dictionary
  if (!PyDict_Check(dim_dict))
  {
    PyErr_SetString(PyExc_TypeError,
                    "'dim_dict' element is not a dictionary");
    throw PythonException();
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "dim_dict = " << PyString_AsString(PyObject_Str(dim_dict))
            << std::endl;
#endif

  ////////////////////////////////
  // Get the 'dist_type' attribute
  ////////////////////////////////
  PyObject * distTypeObj = PyDict_GetItemString(dim_dict, "dist_type");
  if (distTypeObj == NULL)
  {
    PyErr_SetString(PyExc_KeyError,
                    "'dim_dict' has no 'dist_type' key");
    throw PythonException();
  }
  if (!PyString_Check(distTypeObj))
  {
    PyErr_SetString(PyExc_TypeError,
                    "'dist_type' is not a string");
    throw PythonException();
  }
  char * distTypeStr = PyString_AsString(distTypeObj);
  if (strcmp(distTypeStr, "b") == 0)
    dist_type = BLOCK;
  else if (strcmp(distTypeStr, "c") == 0)
    dist_type = CYCLIC;
  else if (strcmp(distTypeStr, "u") == 0)
    dist_type = UNSTRUCTURED;
  else
  {
    PyErr_Format(PyExc_ValueError,
                 "'dist_type' value '%s' is unrecognized",
                 distTypeStr);
    throw PythonException();
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  dist_type = '" << distTypeStr << "'" << std::endl;
#endif

  ///////////////////////////
  // Get the 'size' attribute
  ///////////////////////////
  PyObject * sizeObj = PyDict_GetItemString(dim_dict, "size");
  if (sizeObj == NULL)
  {
    PyErr_SetString(PyExc_KeyError,
                    "'dim_dict' required key 'size' not present");
    throw PythonException();
  }
  if (!PyInt_Check(sizeObj))
  {
    PyErr_Format(PyExc_TypeError, "'size' is not an int");
    throw PythonException();
  }
  size = PyInt_AsLong(sizeObj);
  if (size < 0)
  {
    PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'size' = %d is less "
                 "than zero", size);
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  size = " << size << std::endl;
#endif

  /////////////////////////////////////
  // Get the 'proc_grid_size' attribute
  /////////////////////////////////////
  PyObject * commDimObj = PyDict_GetItemString(dim_dict, "proc_grid_size");
  if (commDimObj == NULL)
  {
    PyErr_SetString(PyExc_KeyError,
                    "'dim_dict' required key 'proc_grid_size' not present");
    throw PythonException();
  }
  if (!PyInt_Check(commDimObj))
  {
    PyErr_Format(PyExc_TypeError,
                 "'proc_grid_size' is not an int");
    throw PythonException();
  }
  proc_grid_size = PyInt_AsLong(commDimObj);
  if (proc_grid_size < 1)
  {
    PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'proc_grid_size' = "
                 "%d is less than one", proc_grid_size);
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  proc_grid_size = " << proc_grid_size << std::endl;
#endif

  /////////////////////////////////////
  // Get the 'proc_grid_rank' attribute
  /////////////////////////////////////
  PyObject * commRankObj = PyDict_GetItemString(dim_dict, "proc_grid_rank");
  if (commRankObj == NULL)
  {
    PyErr_SetString(PyExc_KeyError,
                    "'dim_dict' required key 'proc_grid_rank' not present");
    throw PythonException();
  }
  if (!PyInt_Check(commRankObj))
  {
    PyErr_Format(PyExc_TypeError,
                 "'proc_grid_rank' is not an int");
    throw PythonException();
  }
  proc_grid_rank = PyInt_AsLong(commRankObj);
  if ((proc_grid_rank < 0) || (proc_grid_rank >= proc_grid_size))
  {
    PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'proc_grid_rank' = "
                 "%d is out of range [0, %d)", proc_grid_rank, proc_grid_size);
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  proc_grid_rank = " << proc_grid_rank << std::endl;
#endif

  ////////////////////////////
  // Get the 'start' attribute
  ////////////////////////////
  PyObject * startObj = PyDict_GetItemString(dim_dict, "start");
  if (startObj == NULL)
  {
    if ((dist_type == BLOCK) || (dist_type == CYCLIC))
    {
      PyErr_SetString(PyExc_KeyError,
                      "'dim_dict' is of type BLOCK or CYCLIC but has no "
                      "'start' key");
      throw PythonException();
    }
    else
    {
      start = -1;
    }
  }
  else
  {
    if (!PyInt_Check(startObj))
    {
      PyErr_SetString(PyExc_TypeError,
                      "'start' is not an int");
      throw PythonException();
    }
    start = PyInt_AsLong(startObj);
    if ((start < 0) || (start > size))
    {
      PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'start' = "
                   "%d is out of range [0, %d]", start, size);
    }
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  start = " << start << std::endl;
#endif

  ///////////////////////////
  // Get the 'stop' attribute
  ///////////////////////////
  PyObject * stopObj = PyDict_GetItemString(dim_dict, "stop");
  if (stopObj == NULL)
  {
    if (dist_type == BLOCK)
    {
      PyErr_SetString(PyExc_KeyError,
                      "'dim_dict' is of type BLOCK but has no 'stop' key");
      throw PythonException();
    }
    else
    {
      stop = -1;
    }
  }
  else
  {
    if (!PyInt_Check(stopObj))
    {
      PyErr_SetString(PyExc_TypeError,
                      "'stop' for axis %d is not an int");
      throw PythonException();
    }
    stop = PyInt_AsLong(stopObj);
    if ((stop < 0) || (stop > size))
    {
      PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'stop' = "
                   "%d is out of range [0, %d]", stop, size);
    }
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  stop = " << stop << std::endl;
#endif

  /////////////////////////////////
  // Get the 'block_size' attribute
  /////////////////////////////////
  PyObject * blockSizeObj = PyDict_GetItemString(dim_dict, "block_size");
  if (blockSizeObj == NULL)
  {
    block_size = 1;
  }
  else
  {
    if (!PyInt_Check(blockSizeObj))
    {
      PyErr_SetString(PyExc_TypeError,
                      "'block_size' is not an int");
      throw PythonException();
    }
    block_size = PyInt_AsLong(stopObj);
    if (block_size < 1)
    {
      PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'block_size' = "
                   "%d is less than one", block_size);
    }
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  block_size = " << block_size << std::endl;
#endif

  //////////////////////////////
  // Get the 'padding' attribute
  //////////////////////////////
  PyObject * paddingObj = PyDict_GetItemString(dim_dict, "padding");
  if (paddingObj == NULL)
  {
    padding[0] = 0;
    padding[1] = 0;
  }
  else
  {
    if (!PySequence_Check(paddingObj))
    {
      PyErr_SetString(PyExc_ValueError, "'padding' attribute must be a "
                      "sequence of two integers");
      throw PythonException();
    }
    if (PySequence_Size(paddingObj) != 2)
    {
      PyErr_SetString(PyExc_ValueError, "'padding' attribute must be a "
                      "sequence of two integers");
      throw PythonException();
    }
    PyObject * padObj = PySequence_GetItem(padObj, 0);
    if (!PyInt_Check(padObj))
    {
      PyErr_SetString(PyExc_TypeError, "Lower 'padding' object is not an "
                      "integer");
      throw PythonException();
    }
    padding[0] = PyInt_AsLong(padObj);
    if (padding[0] < 0)
    {
      PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'padding'[0] = "
                   "%d is less than zero", padding[0]);
    }
    padObj = PySequence_GetItem(padObj, 1);
    if (!PyInt_Check(padObj))
    {
      PyErr_SetString(PyExc_TypeError, "Upper 'padding' object is not an "
                      "integer");
      throw PythonException();
    }
    padding[1] = PyInt_AsLong(padObj);
    if (padding[1] < 0)
    {
      PyErr_Format(PyExc_ValueError, "'dim_dict' attribute 'padding'[1] = "
                   "%d is less than zero", padding[1]);
    }
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  padding = " << padding << std::endl;
#endif

  ///////////////////////////////
  // Get the 'periodic' attribute
  ///////////////////////////////
  PyObject * periodicObj = PyDict_GetItemString(dim_dict, "periodic");
  if (periodicObj == NULL)
  {
    periodic = false;
  }
  else
  {
    if (!PyBool_Check(periodicObj))
    {
      PyErr_SetString(PyExc_TypeError,
                      "'periodic' is not a bool");
      throw PythonException();
    }
    if (periodicObj == Py_True) periodic = true;
    else                        periodic = false;
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  periodic = " << (periodic ? "true" : "false") << std::endl;
#endif

  /////////////////////////////////
  // Get the 'one_to_one' attribute
  /////////////////////////////////
  PyObject * oneToOneObj = PyDict_GetItemString(dim_dict, "one_to_one");
  if (oneToOneObj == NULL)
  {
    one_to_one = false;
  }
  else
  {
    if (!PyBool_Check(oneToOneObj))
    {
      PyErr_SetString(PyExc_TypeError,
                      "'one_to_one' is not a bool");
      throw PythonException();
    }
    if (oneToOneObj == Py_True) one_to_one = true;
    else                        one_to_one = false;
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  one_to_one = " << (one_to_one ? "true" : "false")
            << std::endl;
#endif

  //////////////////////////////
  // Get the 'indices' attribute
  //////////////////////////////
  PyObject * indicesObj = PyDict_GetItemString(dim_dict, "indices");
  if (indicesObj == NULL)
  {
    if (dist_type == UNSTRUCTURED)
    {
      PyErr_SetString(PyExc_KeyError, "'dim_dict' with type UNSTRUCTURED "
                      "requires 'indices' key");
      throw PythonException();
    }
  }
  else
  {
    PyObject * indicesArray =
      PyArray_ContiguousFromAny(indicesObj, NPY_INT, 1, 1);
    if (indicesArray == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "'dim_dict' attribute 'indices' "
                      "specified illegally");
      throw PythonException();
    }
    npy_intp numIndices = PyArray_Size(indicesArray);
    indices.resize(numIndices);
    int * source = (int*) PyArray_DATA((PyArrayObject*)indicesArray);
    for (IndicesType::size_type i = 0; i < numIndices; ++i)
      indices[i] = *(source++);
    Py_DECREF(indicesArray);
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "  indices = " << indices() << std::endl;
#endif
}

////////////////////////////////////////////////////////////////////////

DimensionDictionary::~DimensionDictionary()
{
}

////////////////////////////////////////////////////////////////////////

DistArrayProtocol::DistArrayProtocol(PyObject * distarray)
{
  // Check for a dictionary
  if (!PyDict_Check(distarray))
  {
    PyObject * pyType = PyObject_Type(distarray);
    char * typeStr = PyString_AsString(PyObject_Str(pyType));
    //PyErr_SetString(PyExc_ValueError, "distarray object must be a dictionary");
    PyErr_Format(PyExc_ValueError, "distarray object must be a dictionary, "
                 "given '%s'", typeStr);
    throw PythonException();
  }
  __distarray__ = distarray;
  Py_INCREF(__distarray__);

  // Check for an instance method
  // if (!PyMethod_Check(distarray))
  // {
  //   PyObject * pyType = PyObject_Type(distarray);
  //   char * typeStr = PyString_AsString(PyObject_Str(pyType));
  //   PyErr_Format(PyExc_ValueError, "distarray attribute must be an instance "
  //                "method, given '%s'", typeStr);
  //   throw PythonException();    
  // }
  // std::cout << "a" << std::endl;
  // //PyObject * self = PyMethod_Self(distarray);
  // //std::cout << "b" << std::endl;
  // PyObject * function = PyMethod_Function(distarray);
  // std::cout << "b" << std::endl;
  // if (function == NULL) throw PythonException();
  // std::cout << "d" << std::endl;
  // __distarray__ = PyObject_CallObject(function, NULL);
  // std::cout << "__distarray__ = "
  //           << PyString_AsString(PyObject_Str(__distarray__)) << std::endl;

  //////////////////
  // Get the version
  //////////////////
  PyObject * versionObj = PyDict_GetItemString(distarray, "__version__");
  if (versionObj == NULL)
  {
    PyErr_SetString(PyExc_KeyError, "distarray required key '__version__' not "
                    "present");
    throw PythonException();
  }
  __version__ = PyString_AsString(versionObj);
  if (PyErr_Occurred())
  {
    Py_XDECREF(__distarray__);
    throw PythonException();
  }
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "__version__ = " << __version__ << std::endl;
#endif

  /////////////////
  // Get the buffer
  /////////////////
  __buffer__ = PyDict_GetItemString(distarray, "buffer");
  if (__buffer__ == NULL)
  {
    PyErr_SetString(PyExc_KeyError, "distarray required key 'buffer' not "
                    "present");
    Py_XDECREF(__distarray__);
    throw PythonException();
  }
  // What is the appropriate type check for buffer?
  Py_INCREF(__buffer__);
#ifdef PYTRILINOS_DAP_VERBOSE
  std::cout << "buffer = " << PyString_AsString(PyObject_Str(__buffer__))
            << std::endl;
#endif

  /////////////////////////////////
  // Get the dim_data Python object
  /////////////////////////////////
  PyObject * dimDataObj = PyDict_GetItemString(distarray, "dim_data");
  if (dimDataObj == NULL)
  {
    PyErr_SetString(PyExc_KeyError, "distarray required key 'dim_data' not "
                    "present");
    Py_XDECREF(__distarray__);
    Py_XDECREF(__buffer__);
    throw PythonException();
  }
  if (!PySequence_Check(dimDataObj))
  {
    PyErr_SetString(PyExc_TypeError, "distarray key 'dim_data' must have a "
                    "value that is a sequence");
    Py_XDECREF(__distarray__);
    Py_XDECREF(__buffer__);
    throw PythonException();
  }

  ////////////////////////////////////////////////////////////
  // Loop over dim_data and extract the dimension dictionaries
  ////////////////////////////////////////////////////////////
  dim_data.clear();
  int numDims = PySequence_Size(dimDataObj);
  for (int axis = 0; axis < numDims; ++axis)
  {
    try
    {
      PyObject * dim_dict = PySequence_GetItem(dimDataObj, axis);
#ifdef PYTRILINOS_DAP_VERBOSE
      std::cout << "axis " << axis << ": ";
#endif
      dim_data.push_back(DimensionDictionary(dim_dict));
    }
    catch (PythonException e)
    {
      // One of the dim_dicts was improperly formed.  That messes up
      // this constructor, so de-reference the stored Python objects
      // and re-throw.
      Py_XDECREF(__distarray__);
      Py_XDECREF(__buffer__);
      throw e;
    }
  }
}

////////////////////////////////////////////////////////////////////////

DistArrayProtocol::~DistArrayProtocol()
{
  Py_XDECREF(__distarray__);
  Py_XDECREF(__buffer__);
}

////////////////////////////////////////////////////////////////////////

int DistArrayProtocol::num_dims() const
{
  return dim_data.size();
}

////////////////////////////////////////////////////////////////////////

DistributionType DistArrayProtocol::dist_type(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].dist_type;
}

////////////////////////////////////////////////////////////////////////

int DistArrayProtocol::size(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].size;
}

////////////////////////////////////////////////////////////////////////

int DistArrayProtocol::proc_grid_size(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].proc_grid_size;
}

////////////////////////////////////////////////////////////////////////

int DistArrayProtocol::proc_grid_rank(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].proc_grid_rank;
}

////////////////////////////////////////////////////////////////////////

int DistArrayProtocol::start(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].start;
}

////////////////////////////////////////////////////////////////////////

int DistArrayProtocol::stop(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].stop;
}

////////////////////////////////////////////////////////////////////////

IndicesType DistArrayProtocol::indices(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].indices;
}

////////////////////////////////////////////////////////////////////////

PaddingType DistArrayProtocol::padding(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].padding;
}

////////////////////////////////////////////////////////////////////////

int DistArrayProtocol::block_size(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].block_size;
}

////////////////////////////////////////////////////////////////////////

bool DistArrayProtocol::periodic(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].periodic;
}

////////////////////////////////////////////////////////////////////////

bool DistArrayProtocol::one_to_one(int axis) const
{
  checkAxis(axis);
  return dim_data[axis].one_to_one;
}

////////////////////////////////////////////////////////////////////////

PyObject * DistArrayProtocol::buffer() const
{
  return __buffer__;
}

////////////////////////////////////////////////////////////////////////

std::string DistArrayProtocol::version() const
{
  return __version__;
}

////////////////////////////////////////////////////////////////////////

void DistArrayProtocol::checkAxis(int axis) const
{
  if ((axis < 0) || (axis >= dim_data.size()))
  {
    PyErr_Format(PyExc_IndexError, "axis = %d out of range [0, %ld)",
                 axis, dim_data.size());
    throw PythonException();
  }
}

////////////////////////////////////////////////////////////////////////

}

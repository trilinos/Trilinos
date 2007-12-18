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
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_CrsSingletonFilter.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_LinearProblem.h"

%}

//////////////
// Typemaps //
//////////////
%epetra_argout_typemaps(Epetra_CrsMatrix)
%epetra_argout_typemaps(Epetra_VbrMatrix)

////////////////
// Macro code //
////////////////
%define %epetra_global_row_method(className,methodName)
%apply (double * IN_ARRAY1, int DIM1) {(double * Values,  int NumValues )};
%apply (int    * IN_ARRAY1, int DIM1) {(int    * Indices, int NumIndices)};
%extend className
{
  int methodName(int Row, double* Values, int NumValues, int* Indices,
		 int NumIndices)
  {
    if (NumValues != NumIndices)
    {
      PyErr_Format(PyExc_ValueError,
		   "Values length %d not equal to Indices length %d", 
		   NumValues, NumIndices);
      return -1;
    }
    return self->methodName(Row, NumValues, Values, Indices);
  }

  int methodName(PyObject* Rows, PyObject* Cols, PyObject* Values)
  {
    int numValEntries;
    int numRowEntries;
    int numColEntries;
    int result=0;
    PyArrayObject * rowArray = NULL;
    PyArrayObject * colArray = NULL;
    PyArrayObject * valArray = NULL;

    // Create the array of values
    valArray = (PyArrayObject*) PyArray_ContiguousFromObject(Values,'d',0,0);
    if (valArray == NULL) goto fail;
    numValEntries = (int) PyArray_MultiplyList(valArray->dimensions,valArray->nd);

    // Create the array of rows
    rowArray = (PyArrayObject*) PyArray_ContiguousFromObject(Rows,'i',0,0);
    if (rowArray == NULL) goto fail;
    numRowEntries = (int) PyArray_MultiplyList(rowArray->dimensions,rowArray->nd);

    // Create the array of cols
    colArray = (PyArrayObject*) PyArray_ContiguousFromObject(Cols,'i',0,0);
    if (colArray == NULL) goto fail;
    numColEntries = (int) PyArray_MultiplyList(colArray->dimensions,colArray->nd);

    if(numValEntries != numColEntries || numValEntries != numRowEntries ||
       numRowEntries != numColEntries)
    {
      PyErr_Format(PyExc_ValueError, 
		   "lengths of Rows, Cols, Values not equal: %d, %d, %d", 
		   numRowEntries, numColEntries, numValEntries);
      goto fail;
    }

    for(int i = 0 ; i < numValEntries ; ++i)
    {
      double Value = ((double*)valArray->data)[i];
      int Row = ((int*)rowArray->data)[i];
      int Col = ((int*)colArray->data)[i];

      result = self->methodName(Row, 1, &Value, &Col);
      if(result < 0) goto fail;
    }
    Py_DECREF(valArray);
    Py_DECREF(rowArray);
    Py_DECREF(colArray);
    return result;
  fail:
    Py_XDECREF(valArray);
    Py_XDECREF(rowArray);
    Py_XDECREF(colArray);
    return -1;
  }
}
%ignore className::methodName;
%enddef

%define %epetra_my_row_method(className,methodName)
%ignore className::methodName(int,int,double*,int*);
%apply (double * IN_ARRAY1, int DIM1) {(double * Values,  int NumValues )};
%apply (int    * IN_ARRAY1, int DIM1) {(int    * Indices, int NumIndices)};
%extend className
{
  int methodName(int Row, double* Values, int NumValues, int* Indices,
		 int NumIndices)
  {
    // Check for column map
    if (!self->HaveColMap())
    {
      PyErr_SetString(PyExc_RuntimeError,
		      "methodName" " cannot be called on " "className"
		      " that does not have a column map");
      return -2;
    }
    if (NumValues != NumIndices)
    {
      PyErr_Format(PyExc_ValueError,
		   "Values length %d not equal to Indices length %d", 
		   NumValues, NumIndices);
      return -1;
    }
    return self->methodName(Row, NumValues, Values, Indices);
  }

  int methodName(PyObject* Rows, PyObject* Cols, PyObject* Values)
  {
    int numValEntries;
    int numRowEntries;
    int numColEntries;
    int result=0;
    PyArrayObject * rowArray = NULL;
    PyArrayObject * colArray = NULL;
    PyArrayObject * valArray = NULL;

    if (!self->HaveColMap())
    {
      PyErr_SetString(PyExc_RuntimeError, "methodName" " cannot be called on"
		      "className" " that does not have a column map");
      goto fail;
    }

    // Create the array of values
    valArray = (PyArrayObject*) PyArray_ContiguousFromObject(Values,'d',0,0);
    if (valArray == NULL) goto fail;
    numValEntries = (int) PyArray_MultiplyList(valArray->dimensions,valArray->nd);

    // Create the array of rows
    rowArray = (PyArrayObject*) PyArray_ContiguousFromObject(Rows,'i',0,0);
    if (rowArray == NULL) goto fail;
    numRowEntries = (int) PyArray_MultiplyList(rowArray->dimensions,rowArray->nd);

    // Create the array of cols
    colArray = (PyArrayObject*) PyArray_ContiguousFromObject(Cols,'i',0,0);
    if (colArray == NULL) goto fail;
    numColEntries = (int) PyArray_MultiplyList(colArray->dimensions,colArray->nd);

    if(numValEntries != numColEntries || numValEntries != numRowEntries || 
       numRowEntries != numColEntries)
    {
      PyErr_Format(PyExc_ValueError, 
		   "lengths of Rows, Cols, Values not equal: %d, %d, %d", 
		   numRowEntries, numColEntries, numValEntries);
      goto fail;
    }

    for(int i = 0 ; i < numValEntries ; ++i)
    {
      double Value = ((double*)valArray->data)[i];
      int Row = ((int*)rowArray->data)[i];
      int Col = ((int*)colArray->data)[i];

      result = self->methodName(Row, 1, &Value, &Col);
      if(result < 0) goto fail;
    }
    Py_DECREF(valArray);
    Py_DECREF(rowArray);
    Py_DECREF(colArray);
    return result;
  fail:
    Py_XDECREF(valArray);
    Py_XDECREF(rowArray);
    Py_XDECREF(colArray);
    return -1;
  }
}
%enddef

/////////////////////////////
// Epetra_Operator support //
/////////////////////////////
%warnfilter(473)     Epetra_Operator;
%feature("director") Epetra_Operator;
%rename(Operator)    Epetra_Operator;
%include "Epetra_Operator.h"

////////////////////////////////
// Epetra_InvOperator support //
////////////////////////////////
%warnfilter(473)     Epetra_InvOperator;
%feature("director") Epetra_InvOperator;
%rename(InvOperator) Epetra_InvOperator;
%include "Epetra_InvOperator.h"

//////////////////////////////
// Epetra_RowMatrix support //
//////////////////////////////
%include "Epetra_RowMatrix_Utils.i"
%warnfilter(473)     Epetra_RowMatrix;
%feature("director") Epetra_RowMatrix;
%rename(RowMatrix)   Epetra_RowMatrix;
// These typemaps are specifically for the NumMyRowEntries() and
// ExtractMyRowCopy() methods.  They turn output arguments NumEntries,
// Values and Indices into numpy arrays with appropriate data buffers.
// The user therefore has access to these buffers via the [ ]
// (bracket) operators (__setitem__ and __getitem__).  NumEntries,
// even though it is a scalar in C++, must be accessed as
// NumEntries[0].
%typemap(directorin) int &NumEntries
%{
  intp dims$argnum[ ] = { (intp) 1 };
  $input = PyArray_SimpleNewFromData(1, dims$argnum, NPY_INT, (void*)&$1_name);
%}
%typemap(directorin) double *Values
%{
  intp dims$argnum[ ] = { (intp) Length };
  $input = PyArray_SimpleNewFromData(1, dims$argnum, NPY_DOUBLE, (void*)$1_name);
%}
%typemap(directorin) int *Indices
%{
  intp dims$argnum[ ] = { (intp) Length };
  $input = PyArray_SimpleNewFromData(1, dims$argnum, NPY_INT, (void*)$1_name);
%}
%include "Epetra_RowMatrix.h"

///////////////////////////////////
// Epetra_BasicRowMatrix support //
///////////////////////////////////
%warnfilter(473)        Epetra_BasicRowMatrix;
%feature("director")    Epetra_BasicRowMatrix;
%rename(BasicRowMatrix) Epetra_BasicRowMatrix;
%ignore Epetra_BasicRowMatrix::ExtractMyEntryView(int,const double*&,int&,int&) const;
// Typemap for double * & Value
%typemap(directorin) double *&Value
%{
  intp dims$argnum[ ] = { (intp) 1 };
  $input = PyArray_SimpleNewFromData(1, dims$argnum, NPY_DOUBLE, (void*)$1_name);
%}
// Apply the referenced int typemap for NumEntries to RowIndex and ColIndex
%apply int &NumEntries {int &RowIndex, int &ColIndex};
%include "Epetra_BasicRowMatrix.h"
// Make sure the following typemaps specific to Epetra_BasicRowMap
// are not applied after this point.
%clear double *& Value;
%clear int    &  RowIndex;
%clear int    &  ColIndex;
// The typemaps below were defined for the Epetra_RowMatrix support
// just above, but they apply to the Epetra_BasicRowMap as well.
// However, we also want to make sure they do not get applied after
// this point.
%clear int    & NumEntries;
%clear double * Values;
%clear int    * Indices;

//////////////////////////////
// Epetra_CrsMatrix support //
//////////////////////////////
%rename(CrsMatrix) Epetra_CrsMatrix;
%epetra_global_row_method(Epetra_CrsMatrix, InsertGlobalValues )
%epetra_global_row_method(Epetra_CrsMatrix, ReplaceGlobalValues)
%epetra_global_row_method(Epetra_CrsMatrix, SumIntoGlobalValues)
%epetra_my_row_method(Epetra_CrsMatrix, InsertMyValues )
%epetra_my_row_method(Epetra_CrsMatrix, ReplaceMyValues)
%epetra_my_row_method(Epetra_CrsMatrix, SumIntoMyValues)
%apply (int * IN_ARRAY1, int DIM1) {(const int* NumEntriesPerRow, int NumRows)};
%extend Epetra_CrsMatrix
{
  // Add NumRows to the constructor argument list, so that I can use
  // an appropriate numpy.i typemap, and also check the length of
  // NumEntriesPerRow
  Epetra_CrsMatrix(Epetra_DataAccess   CV,
		   const Epetra_Map  & RowMap,
		   const int         * NumEntriesPerRow,
		   int                 NumRows,
		   bool                StaticProfile=false)
  {
    if (NumRows != RowMap.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
		   "RowMap has %d rows and NumEntriesPerRow has %d elements",
		   RowMap.NumMyElements(), NumRows);
      return NULL;
    }
    return new Epetra_CrsMatrix(CV, RowMap, NumEntriesPerRow, StaticProfile);
  }

  // Add NumRows to the constructor argument list, so that I can use
  // an appropriate numpy.i typemap, and also check the length of
  // NumEntriesPerRow
  Epetra_CrsMatrix(Epetra_DataAccess   CV,
		   const Epetra_Map  & RowMap,
		   const Epetra_Map  & ColMap,
		   const int         * NumEntriesPerRow,
		   int                 NumRows,
		   bool                StaticProfile=false)
  {
    if (NumRows != RowMap.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
		   "RowMap has %d rows and NumEntriesPerRow has %d elements",
		   RowMap.NumMyElements(), NumRows);
      return NULL;
    }
    return new Epetra_CrsMatrix(CV, RowMap, ColMap, NumEntriesPerRow, StaticProfile);
  }

  PyObject * ExtractGlobalRowCopy(int globalRow) const
  {
    int        lrid          = 0;
    int        numEntries    = 0;
    int        result        = 0;
    intp       dimensions[ ] = { 0 };
    double   * values        = NULL;
    int      * indices       = NULL;
    PyObject * valuesArray   = NULL;
    PyObject * indicesArray  = NULL;

    lrid = self->LRID(globalRow);
    if (lrid == -1)
    {
      PyErr_Format(PyExc_ValueError, "Invalid global row index = %d", globalRow);
      goto fail;
    }
    dimensions[0] = self->NumMyEntries(lrid);
    valuesArray   = PyArray_SimpleNew(1,dimensions,NPY_DOUBLE);
    indicesArray  = PyArray_SimpleNew(1,dimensions,NPY_INT   );
    values        = (double*) array_data(valuesArray );
    indices       = (int   *) array_data(indicesArray);
    result        = self->ExtractGlobalRowCopy(globalRow, dimensions[0], numEntries,
					       values, indices);
    if (result == -2)
    {
      PyErr_SetString(PyExc_RuntimeError, "Matrix not completed");
      goto fail;
    }
    return Py_BuildValue("(OO)",valuesArray,indicesArray);
  fail:
    Py_XDECREF(valuesArray );
    Py_XDECREF(indicesArray);
    return NULL;
  }

  PyObject * ExtractMyRowCopy(int localRow) const
  {
    int        numEntries    = 0;
    int        result        = 0;
    intp       dimensions[ ] = { 0 };
    double   * values        = NULL;
    int      * indices       = NULL;
    PyObject * valuesArray   = NULL;
    PyObject * indicesArray  = NULL;

    if (localRow < 0 || localRow >= self->NumMyRows())
    {
      PyErr_Format(PyExc_ValueError, "Invalid local row index = %d", localRow);
      goto fail;
    }
    dimensions[0] = self->NumMyEntries(localRow);
    valuesArray   = PyArray_SimpleNew(1,dimensions,NPY_DOUBLE);
    indicesArray  = PyArray_SimpleNew(1,dimensions,NPY_INT   );
    values        = (double*) array_data(valuesArray );
    indices       = (int   *) array_data(indicesArray);
    result        = self->ExtractMyRowCopy(localRow, dimensions[0], numEntries,
					   values, indices);
    if (result == -2)
    {
      PyErr_SetString(PyExc_RuntimeError, "Matrix not completed");
      goto fail;
    }
    return Py_BuildValue("(OO)",valuesArray,indicesArray);
  fail:
    Py_XDECREF(valuesArray );
    Py_XDECREF(indicesArray);
    return NULL;
  }

  PyObject * __setitem__(PyObject* args, double val)
  {
    PyObject * rowObj = NULL;
    PyObject * colObj = NULL;
    int row           = 0;
    int col           = 0;
    if (!(PyArg_ParseTuple(args, "OO:Epetra_CrsMatrix___setitem__",
			   &rowObj, &colObj)      &&
	  SWIG_IsOK(SWIG_AsVal_int(rowObj, &row)) &&
	  SWIG_IsOK(SWIG_AsVal_int(colObj, &col))    ))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    if (self->ReplaceGlobalValues(row, 1, &val, &col))
      self->InsertGlobalValues(row, 1, &val, &col);
    return Py_BuildValue("");
  }

  PyObject* __getitem__(PyObject* args) const
  {
    int        grid          = 0;
    int        lrid          = 0;
    int        gcid          = 0;
    int        lcid          = 0;
    int        error         = 0;
    int        numEntries    = 0;
    intp       dimensions[ ] = { 0 };
    int      * indices       = NULL;
    double     result        = 0.0;
    double   * values        = NULL;
    double   * data          = NULL;
    PyObject * rowObj        = NULL;
    PyObject * colObj        = NULL;
    PyObject * returnObj     = NULL;

    // If the argument is an integer, get the global row ID, construct
    // a return PyArray, and obtain the data pointer
    if (SWIG_IsOK(SWIG_AsVal_int(args, &grid)))
    {
      dimensions[0] = self->NumMyCols();
      returnObj = PyArray_SimpleNew(1,dimensions,NPY_DOUBLE);
      if (returnObj == NULL) goto fail;
      data = (double*) array_data(returnObj);
      for (int i=0; i<dimensions[0]; ++i) data[i] = 0.0;
      returnObj = PyArray_Return((PyArrayObject*)returnObj);

      // If the matrix is FillComplete()-ed, obtain the local row data
      // and copy it into the data buffer
      if (self->Filled())
      {
	lrid = self->LRID(grid);
	if (lrid == -1)
	{
	  PyErr_Format(PyExc_IndexError, "Global row index %d not on processor",
		       grid);
	  goto fail;
	}
	error = self->ExtractMyRowView(lrid, numEntries, values, indices);
	if (error)
	{
	  PyErr_Format(PyExc_RuntimeError, "ExtractMyRowView error code %d", error);
	  goto fail;
	}
	for (int i=0; i<numEntries; ++i)
	{
	  lcid = indices[i];
	  data[lcid] = values[i];
	}

      // If the matrix is not FillComplete()-ed, raise an exception
      }
      else
      {
	if (self->Comm().NumProc() > 1)
	{
	  PyErr_SetString(PyExc_IndexError,
			  "__getitem__ cannot be called with single "
			  "index unless CrsMatrix has been filled");
	  goto fail;
	}
	else
	{
	  error = self->ExtractGlobalRowView(grid, numEntries, values, indices);
	  if (error)
	  {
	    if (error == -1)
	    {
	      PyErr_Format(PyExc_IndexError, "Global row %d not on processor", grid);
	    }
	    else
	    {
	      PyErr_Format(PyExc_RuntimeError, "ExtractGlobalRowView error code %d",
			   error);
	    }
	    goto fail;
	  }
	  for (int i=0; i<numEntries; ++i)
	  {
	    gcid = indices[i];
	    lcid = self->LCID(gcid);
	    if (lcid == -1)
	    {
	      PyErr_Format(PyExc_IndexError,
			   "Global column index %d not on processor", gcid);
	      goto fail;
	    }
	    data[lcid] = values[i];
	  }
	}
      }

    // If the arguments are two integers, obtain a single result value
    }
    else if (PyArg_ParseTuple(args, "OO:Epetra_CrsMatrix___getitem__",
				&rowObj, &colObj)       &&
	       SWIG_IsOK(SWIG_AsVal_int(rowObj, &grid)) &&
	       SWIG_IsOK(SWIG_AsVal_int(colObj, &gcid))   )
    {
      lrid = self->LRID(grid);
      if (lrid == -1)
      {
	PyErr_Format(PyExc_IndexError, "Global row %d not on processor", grid);
	goto fail;
      }
      lcid = self->LCID(gcid);
      if (lcid == -1)
      {
	PyErr_Format(PyExc_IndexError, "Global column %d not on processor", gcid);
	goto fail;
      }

      // If the matrix is FillComplete()-ed, obtain the local row data
      // and column data
      if (self->Filled())
      {
	error = self->ExtractMyRowView(lrid, numEntries, values, indices);
	if (error)
	{
	  PyErr_Format(PyExc_RuntimeError, "ExtractMyRowView error code %d", error);
	  goto fail;
	}
	for (int i=0; i<numEntries; ++i)
	{
	  if (indices[i] == lcid)
	  {
	    result = values[i];
	    break;
	  }
	}

      // If the matrix is not FillComplete()-ed, obtain the local row data
      // and column data
      }
      else
      {
	error = self->ExtractGlobalRowView(grid, numEntries, values, indices);
	if (error)
	{
	  PyErr_Format(PyExc_RuntimeError, "ExtractGlobalRowView error code %d",
		       error);
	  goto fail;
	}
	for (int i=0; i<numEntries; ++i)
	{
	  if (indices[i] == gcid)
	  {
	    result = values[i];
	    break;
	  }
	}
      }
      returnObj = PyFloat_FromDouble(result);
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      goto fail;
    }
    return returnObj;
  fail:
    return NULL;
  }
}
%ignore Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess,
					   const Epetra_Map &,
					   const int *,
					   bool);
%ignore Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess,
					   const Epetra_Map &,
					   const Epetra_Map &,
					   const int *,
					   bool);
%ignore Epetra_CrsMatrix::ExtractGlobalRowCopy;
%ignore Epetra_CrsMatrix::ExtractGlobalRowView;
%ignore Epetra_CrsMatrix::ExtractMyRowCopy;
%ignore Epetra_CrsMatrix::ExtractMyRowView;
%ignore Epetra_CrsMatrix::ExtractCrsDataPointers;
%include "Epetra_CrsMatrix.h"
%clear (const int* NumEntriesPerRow, int NumRows   );
%clear (double*    Values,           int NumValues );
%clear (int*       Indices,          int NumIndices);

////////////////////////////////
// Epetra_FECrsMatrix support //
////////////////////////////////
%rename(FECrsMatrix) Epetra_FECrsMatrix;
%extend Epetra_FECrsMatrix
{
  void __setitem__(PyObject* args, double val) 
  {
    int row = 0;
    int col = 0;
    PyObject *rowObj = NULL;
    PyObject *colObj = NULL;

    if (!(PyArg_ParseTuple(args, "OO:Epetra_FECrsMatrix___setitem__",
			   &rowObj, &colObj)      &&
	  SWIG_IsOK(SWIG_AsVal_int(rowObj, &row)) &&
	  SWIG_IsOK(SWIG_AsVal_int(colObj, &col))    ))
    {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }

    if (self->ReplaceGlobalValues(1, &row, 1, &col, &val))
      self->InsertGlobalValues(1, &row, 1, &col, &val);
  }

  PyObject* __getitem__(PyObject* args) 
  {
    int row = 0;
    int col = 0;
    PyObject *rowObj = NULL;
    PyObject *colObj = NULL;

    if (SWIG_IsOK(SWIG_AsVal_int(args, &row)))
    {
      return(Epetra_RowMatrix_GetEntries(*self, row));
    }
    else if (PyArg_ParseTuple(args, "OO:Epetra_FECrsMatrix___getitem__",
			      &rowObj, &colObj)      && 
	     SWIG_IsOK(SWIG_AsVal_int(rowObj, &row)) &&
	     SWIG_IsOK(SWIG_AsVal_int(colObj, &col))    )
    {
      return(Epetra_RowMatrix_GetEntry(*self, row, col));
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      return NULL;
    }
  }

  int InsertGlobalValues(const int Row, const int Size, 
                         const Epetra_SerialDenseVector& Values,
                         const Epetra_IntSerialDenseVector& Entries)
  {
    return self->InsertGlobalValues(1, &Row,
                                    Size, (int*)Entries.Values(),
                                    Values.Values());
  }

  int InsertGlobalValue(int i, int j, double val)
  {
    double val2 = val;
    int j2 = j;
    return self->InsertGlobalValues(1, &i, 1, &j2, &val2);
  }
}
%include "Epetra_FECrsMatrix.h"

///////////////////////////////////////
// Epetra_CrsSingletonFilter support //
///////////////////////////////////////
%rename(CrsSingletonFilter) Epetra_CrsSingletonFilter;
%include "Epetra_CrsSingletonFilter.h"

//////////////////////////////
// Epetra_VbrMatrix support //
//////////////////////////////
%rename(VbrMatrix) Epetra_VbrMatrix;
%apply (int * IN_ARRAY1, int DIM1) {(int * NumBlockEntriesPerRow, int NumRows)};
%apply (int DIM1, int * IN_ARRAY1) {(int NumBlockEntries, int * BlockIndices)};
%extend Epetra_VbrMatrix
{
  // Add NumRows to the constructor argument list, so that I can use
  // an appropriate numpy.i typemap, and also check the length of
  // NumBlockEntriesPerRow
  Epetra_VbrMatrix(Epetra_DataAccess 	   CV,
		   const Epetra_BlockMap & RowMap,
		   int                   * NumBlockEntriesPerRow,
		   int                     NumRows)
  {
    if (NumRows != RowMap.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
		   "RowMap has %d rows and NumBlockEntriesPerRow has %d elements",
		   RowMap.NumMyElements(), NumRows);
      return NULL;
    }
    return new Epetra_VbrMatrix(CV, RowMap, NumBlockEntriesPerRow);
  }

  // Add NumRows to the constructor argument list, so that I can use
  // an appropriate numpy.i typemap, and also check the length of
  // NumBlockEntriesPerRow
  Epetra_VbrMatrix(Epetra_DataAccess 	   CV,
		   const Epetra_BlockMap & RowMap,
		   const Epetra_BlockMap & ColMap,
		   int                   * NumBlockEntriesPerRow,
		   int                     NumRows)
  {
    if (NumRows != RowMap.NumMyElements())
    {
      PyErr_Format(PyExc_ValueError,
		   "RowMap has %d rows and NumBlockEntriesPerRow has %d elements",
		   RowMap.NumMyElements(), NumRows);
      return NULL;
    }
    return new Epetra_VbrMatrix(CV, RowMap, ColMap, NumBlockEntriesPerRow);
  }

}
%ignore Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess 	   CV,
					   const Epetra_BlockMap & RowMap,
					   int                   * NumBlockEntriesPerRow);
%ignore Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess 	   CV,
					   const Epetra_BlockMap & RowMap,
					   const Epetra_BlockMap & ColMap,
					   int                   * NumBlockEntriesPerRow);
%ignore Epetra_VbrMatrix::Solve(bool, bool, bool,
				Epetra_Vector const&, Epetra_Vector&) const;
%include "Epetra_VbrMatrix.h"
%clear (int * NumBlockEntriesPerRow, int NumRows);
%clear (int NumBlockEntries, int * BlockIndices);

////////////////////////////////
// Epetra_FEVbrMatrix support //
////////////////////////////////
%rename(FEVbrMatrix) Epetra_FEVbrMatrix;
%include "Epetra_FEVbrMatrix.h"

//////////////////////////////
// Epetra_JadMatrix support //
//////////////////////////////
%ignore Epetra_JadMatrix::ExtractMyEntryView(int,double*&,int&,int&);
%rename(JadMatrix) Epetra_JadMatrix;
%include "Epetra_JadMatrix.h"

//////////////////////////////////
// Epetra_LinearProblem support //
//////////////////////////////////
%rename(LinearProblem) Epetra_LinearProblem;
%include "Epetra_LinearProblem.h"

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

////////////////
// Macro code //
////////////////
%define %epetra_global_row_method(method)
int method(int row, PyObject * values, PyObject * indices) {
  int numValEntries;
  int numIndEntries;
  int result;
  PyArrayObject * valArray = NULL;
  PyArrayObject * indArray = NULL;
  // Create the array of values
  valArray = (PyArrayObject*) PyArray_ContiguousFromObject(values,'d',0,0);
  if (valArray == NULL) goto fail;
  numValEntries = (int) PyArray_MultiplyList(valArray->dimensions,valArray->nd);
  // Create the array of indeces
  indArray = (PyArrayObject*) PyArray_ContiguousFromObject(indices,'i',0,0);
  if (indArray == NULL) goto fail;
  numIndEntries = (int) PyArray_MultiplyList(indArray->dimensions,indArray->nd);
  if (numIndEntries != numValEntries) {
    PyErr_Format(PyExc_ValueError, "values length of %d not equal to indices length %d", 
		 numValEntries, numIndEntries);
    goto fail;
  }
  // Manipulate the row values
  result = self->method(row, numValEntries, (double*)valArray->data,
			(int*)indArray->data);
  Py_DECREF(valArray);
  Py_DECREF(indArray);
  return result;
 fail:
  Py_XDECREF(valArray);
  Py_XDECREF(indArray);
  return -1;
}
%enddef

%define %epetra_my_row_method(method)
int method(int row, PyObject * values, PyObject * indices) {
  int numValEntries;
  int numIndEntries;
  int result;
  PyArrayObject * valArray = NULL;
  PyArrayObject * indArray = NULL;
  // Check for column map
  if (!self->HaveColMap()) {
    PyErr_SetString(PyExc_RuntimeError, "method" " cannot be called on a CrsMatrix"
		    " that does not have a column map");
    goto fail;
  }
  // Create the array of values
  valArray = (PyArrayObject*) PyArray_ContiguousFromObject(values,'d',0,0);
  if (valArray == NULL) goto fail;
  numValEntries = (int) PyArray_MultiplyList(valArray->dimensions,valArray->nd);
  // Create the array of indeces
  indArray = (PyArrayObject*) PyArray_ContiguousFromObject(indices,'i',0,0);
  if (indArray == NULL) goto fail;
  numIndEntries = (int) PyArray_MultiplyList(indArray->dimensions,indArray->nd);
  if (numIndEntries != numValEntries) {
    PyErr_Format(PyExc_ValueError, "values length of %d not equal to indices length %d", 
		 numValEntries, numIndEntries);
    goto fail;
  }
  // Manipulate the row values
  result = self->method(row, numValEntries, (double*)valArray->data,
			(int*)indArray->data);
  Py_DECREF(valArray);
  Py_DECREF(indArray);
  return result;
 fail:
  Py_XDECREF(valArray);
  Py_XDECREF(indArray);
  return -1;
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
%warnfilter(473)     Epetra_RowMatrix;
%feature("director") Epetra_RowMatrix;
%rename(RowMatrix)   Epetra_RowMatrix;
%include "Epetra_RowMatrix.h"

///////////////////////////////////
// Epetra_BasicRowMatrix support //
///////////////////////////////////
%ignore Epetra_BasicRowMatrix::ExtractMyEntryView(int,double*&,int&,int&);
%rename(BasicRowMatrix) Epetra_BasicRowMatrix;
%include "Epetra_BasicRowMatrix.h"

//////////////////////////////
// Epetra_CrsMatrix support //
//////////////////////////////
%ignore Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess,
					   const Epetra_Map &,
					   const int *,
					   bool);
%ignore Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess,
					   const Epetra_Map &,
					   const Epetra_Map &,
					   const int *,
					   bool);
%ignore Epetra_CrsMatrix::InsertGlobalValues( int,int,double*,int);
%ignore Epetra_CrsMatrix::ReplaceGlobalValues(int,int,double*,int);
%ignore Epetra_CrsMatrix::SumIntoGlobalValues(int,int,double*,int);
%ignore Epetra_CrsMatrix::InsertMyValues(     int,int,double*,int);
%ignore Epetra_CrsMatrix::ReplaceMyValues(    int,int,double*,int);
%ignore Epetra_CrsMatrix::SumIntoMyValues(    int,int,double*,int);
%ignore Epetra_CrsMatrix::ExtractGlobalRowCopy(int,int,int&,double* ,int*);
%ignore Epetra_CrsMatrix::ExtractMyRowCopy(    int,int,int&,double* ,int*);
%ignore Epetra_CrsMatrix::ExtractGlobalRowCopy(int,int,int&,double*);
%ignore Epetra_CrsMatrix::ExtractMyRowCopy(    int,int,int&,double*);
%ignore Epetra_CrsMatrix::ExtractGlobalRowView(int,int,int&,double*&,int*&);
%ignore Epetra_CrsMatrix::ExtractMyRowView(    int,int,int&,double*&,int*&);
%ignore Epetra_CrsMatrix::ExtractGlobalRowView(int,int,int&,double*&);
%ignore Epetra_CrsMatrix::ExtractMyRowView(    int,int,int&,double*&);
%ignore Epetra_CrsMatrix::ExtractCrsDataPointers(int*&,int*&,double*&);
%rename(CrsMatrix) Epetra_CrsMatrix;
%epetra_exception(Epetra_CrsMatrix, Epetra_CrsMatrix   )
%epetra_exception(Epetra_CrsMatrix, InsertGlobalValues )
%epetra_exception(Epetra_CrsMatrix, ReplaceGlobalValues)
%epetra_exception(Epetra_CrsMatrix, SumIntoGlobalValues)
%epetra_exception(Epetra_CrsMatrix, InsertMyValues     )
%epetra_exception(Epetra_CrsMatrix, ReplaceMyValues    )
%epetra_exception(Epetra_CrsMatrix, SumIntoMyValues    )
%epetra_exception(Epetra_CrsMatrix, OptimizeStorage    )
%epetra_exception(Epetra_CrsMatrix, __setitem__        )
%include "Epetra_CrsMatrix.h"
%extend Epetra_CrsMatrix {

  Epetra_CrsMatrix(Epetra_DataAccess   CV,
		   const Epetra_Map  & rowMap,
		   PyObject          * numEntriesList,
		   bool                staticProfile=false) {
    // Declarations
    PyArrayObject    * numEntriesArray    = NULL;
    int              * numEntriesPerRow   = NULL;
    Epetra_CrsMatrix * returnCrsMatrix    = NULL;
    int                constEntriesPerRow = 0;
    int                listSize           = 0;

    if (PyInt_Check(numEntriesList)) {
      constEntriesPerRow = (int) PyInt_AsLong(numEntriesList);
      returnCrsMatrix    = new Epetra_CrsMatrix(CV,rowMap,constEntriesPerRow,staticProfile);
    } else {
      numEntriesArray = (PyArrayObject*) PyArray_ContiguousFromObject(numEntriesList,'i',0,0);
      if (numEntriesArray == NULL) goto fail;
      numEntriesPerRow = (int*) (numEntriesArray->data);
      listSize = (int) PyArray_MultiplyList(numEntriesArray->dimensions,numEntriesArray->nd);
      if (listSize != rowMap.NumMyElements()) {
	PyErr_Format(PyExc_ValueError,
		     "Row map has %d elements, list of number of entries has %d",
		     rowMap.NumMyElements(), listSize);
	goto fail;
      }
      returnCrsMatrix = new Epetra_CrsMatrix(CV,rowMap,numEntriesPerRow,staticProfile);
      Py_DECREF(numEntriesArray);
    }
    return returnCrsMatrix;

  fail:
    Py_XDECREF(numEntriesArray);
    return NULL;
  }

  Epetra_CrsMatrix(Epetra_DataAccess   CV,
		   const Epetra_Map  & rowMap,
		   const Epetra_Map  & colMap,
		   PyObject          * numEntriesList,
		   bool                staticProfile=false) {
    // Declarations
    PyArrayObject    * numEntriesArray    = NULL;
    int              * numEntriesPerRow   = NULL;
    Epetra_CrsMatrix * returnCrsMatrix    = NULL;
    int                constEntriesPerRow = 0;
    int                listSize           = 0;

    if (PyInt_Check(numEntriesList)) {
      constEntriesPerRow = (int) PyInt_AsLong(numEntriesList);
      returnCrsMatrix    = new Epetra_CrsMatrix(CV,rowMap,colMap,constEntriesPerRow,staticProfile);
    } else {
      numEntriesArray = (PyArrayObject*) PyArray_ContiguousFromObject(numEntriesList,'i',0,0);
      if (numEntriesArray == NULL) goto fail;
      numEntriesPerRow = (int*) numEntriesArray->data;
      listSize = (int) PyArray_MultiplyList(numEntriesArray->dimensions,numEntriesArray->nd);
      if (listSize != rowMap.NumMyElements()) {
	PyErr_Format(PyExc_ValueError,
		     "Row map has %d elements, list of number of entries has %d",
		     rowMap.NumMyElements(), listSize);
	goto fail;
      }
      returnCrsMatrix = new Epetra_CrsMatrix(CV,rowMap,colMap,numEntriesPerRow,staticProfile);
      Py_DECREF(numEntriesArray);
    }
    return returnCrsMatrix;

  fail:
    Py_XDECREF(numEntriesArray);
    return NULL;
  }

  // These macros expand into code for the various methods
  %epetra_global_row_method(InsertGlobalValues)
  %epetra_global_row_method(ReplaceGlobalValues)
  %epetra_global_row_method(SumIntoGlobalValues)
  %epetra_my_row_method(InsertMyValues)
  %epetra_my_row_method(ReplaceMyValues)
  %epetra_my_row_method(SumIntoMyValues)

  PyObject * ExtractGlobalRowCopy(int globalRow) const {
    int        lrid          = 0;
    int        numEntries    = 0;
    int        result        = 0;
    intp       dimensions[ ] = { 0 };
    double   * values        = NULL;
    int      * indices       = NULL;
    PyObject * valuesArray   = NULL;
    PyObject * indicesArray  = NULL;

    lrid = self->LRID(globalRow);
    if (lrid == -1) {
      PyErr_Format(PyExc_ValueError, "Invalid global row index = %d", globalRow);
      goto fail;
    }
    dimensions[0] = self->NumMyEntries(lrid);
    valuesArray   = PyArray_SimpleNew(1,dimensions,'d');
    indicesArray  = PyArray_SimpleNew(1,dimensions,'i');
    values        = (double *) ((PyArrayObject *)valuesArray )->data;
    indices       = (int    *) ((PyArrayObject *)indicesArray)->data;
    result        = self->ExtractGlobalRowCopy(globalRow, dimensions[0], numEntries,
					       values, indices);
    if (result == -2) {
      PyErr_SetString(PyExc_RuntimeError, "Matrix not completed");
      goto fail;
    }
    return Py_BuildValue("(OO)",valuesArray,indicesArray);
  fail:
    Py_XDECREF(valuesArray );
    Py_XDECREF(indicesArray);
    return NULL;
  }

  PyObject * ExtractMyRowCopy(int localRow) const {
    int        numEntries    = 0;
    int        result        = 0;
    intp       dimensions[ ] = { 0 };
    double   * values        = NULL;
    int      * indices       = NULL;
    PyObject * valuesArray   = NULL;
    PyObject * indicesArray  = NULL;

    if (localRow < 0 || localRow >= self->NumMyRows()) {
      PyErr_Format(PyExc_ValueError, "Invalid local row index = %d", localRow);
      goto fail;
    }
    dimensions[0] = self->NumMyEntries(localRow);
    valuesArray   = PyArray_SimpleNew(1,dimensions,'d');
    indicesArray  = PyArray_SimpleNew(1,dimensions,'i');
    values        = (double *) ((PyArrayObject *)valuesArray )->data;
    indices       = (int    *) ((PyArrayObject *)indicesArray)->data;
    result        = self->ExtractMyRowCopy(localRow, dimensions[0], numEntries,
					   values, indices);
    if (result == -2) {
      PyErr_SetString(PyExc_RuntimeError, "Matrix not completed");
      goto fail;
    }
    return Py_BuildValue("(OO)",valuesArray,indicesArray);
  fail:
    Py_XDECREF(valuesArray );
    Py_XDECREF(indicesArray);
    return NULL;
  }

  void __setitem__(PyObject* args, double val) {
    int row;
    int col;
    if (!PyArg_ParseTuple(args, "ii", &row, &col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }
    if (self->ReplaceGlobalValues(row, 1, &val, &col))
      self->InsertGlobalValues(row, 1, &val, &col);
  }

  PyObject* __getitem__(PyObject* args) const {
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
    PyObject * returnObj     = NULL;

    // If the argument is an integer, get the global row ID, construct
    // a return PyArray, and obtain the data pointer
    if (PyInt_Check(args)) {
      grid = (int) PyInt_AsLong(args);
      dimensions[0] = self->NumMyCols();
      returnObj = PyArray_SimpleNew(1,dimensions,'d');
      if (returnObj == NULL) goto fail;
      data = (double*) ((PyArrayObject*) returnObj)->data;
      for (int i=0; i<dimensions[0]; ++i) data[i] = 0.0;
      returnObj = PyArray_Return((PyArrayObject*)returnObj);

      // If the matrix is FillComplete()-ed, obtain the local row data
      // and copy it into the data buffer
      if (self->Filled()) {
	lrid = self->LRID(grid);
	if (lrid == -1) {
	  PyErr_Format(PyExc_IndexError, "Global row index %d not on processor", grid);
	  goto fail;
	}
	error = self->ExtractMyRowView(lrid, numEntries, values, indices);
	if (error) {
	  PyErr_Format(PyExc_RuntimeError, "ExtractMyRowView error code %d", error);
	  goto fail;
	}
	for (int i=0; i<numEntries; ++i) {
	  lcid = indices[i];
	  data[lcid] = values[i];
	}

      // If the matrix is not FillComplete()-ed, raise an exception
      } else {
	if (self->Comm().NumProc() > 1) {
	  PyErr_SetString(PyExc_IndexError, "__getitem__ cannot be called with single "
			  "index unless CrsMatrix has been filled");
	  goto fail;
	} else {
	  error = self->ExtractGlobalRowView(grid, numEntries, values, indices);
	  if (error) {
	    if (error == -1) {
	      PyErr_Format(PyExc_IndexError, "Global row %d not on processor", grid);
	    } else {
	      PyErr_Format(PyExc_RuntimeError, "ExtractGlobalRowView error code %d", error);
	    }
	    goto fail;
	  }
	  for (int i=0; i<numEntries; ++i) {
	    gcid = indices[i];
	    lcid = self->LCID(gcid);
	    if (lcid == -1) {
	      PyErr_Format(PyExc_IndexError, "Global column index %d not on processor", gcid);
	      goto fail;
	    }
	    data[lcid] = values[i];
	  }
	}
      }

    // If the arguments are two integers, obtain a single result value
    } else if (PyArg_ParseTuple(args, "ii", &grid, &gcid)) {
      lrid = self->LRID(grid);
      if (lrid == -1) {
	PyErr_Format(PyExc_IndexError, "Global row %d not on processor", grid);
	goto fail;
      }
      lcid = self->LCID(gcid);
      if (lcid == -1) {
	PyErr_Format(PyExc_IndexError, "Global column %d not on processor", gcid);
	goto fail;
      }

      // If the matrix is FillComplete()-ed, obtain the local row data
      // and column data
      if (self->Filled()) {
	error = self->ExtractMyRowView(lrid, numEntries, values, indices);
	if (error) {
	  PyErr_Format(PyExc_RuntimeError, "ExtractMyRowView error code %d", error);
	  goto fail;
	}
	for (int i=0; i<numEntries; ++i) {
	  if (indices[i] == lcid) {
	    result = values[i];
	    break;
	  }
	}

      // If the matrix is FillComplete()-ed, obtain the local row data
      // and column data
      } else {
	error = self->ExtractGlobalRowView(grid, numEntries, values, indices);
	if (error) {
	  PyErr_Format(PyExc_RuntimeError, "ExtractGlobalRowView error code %d", error);
	  goto fail;
	}
	for (int i=0; i<numEntries; ++i) {
	  if (indices[i] == gcid) {
	    result = values[i];
	    break;
	  }
	}
      }
      returnObj = PyFloat_FromDouble(result);
    } else {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      goto fail;
    }
    return returnObj;
  fail:
    return NULL;
  }
}

////////////////////////////////
// Epetra_FECrsMatrix support //
////////////////////////////////
%rename(FECrsMatrix) Epetra_FECrsMatrix;
%include "Epetra_FECrsMatrix.h"
%extend Epetra_FECrsMatrix {
  void __setitem__(PyObject* args, double val) 
  {
    int Row, Col;
    if (!PyArg_ParseTuple(args, "ii", &Row, &Col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }

    if (self->ReplaceGlobalValues(1, &Row, 1, &Col, &val))
      self->InsertGlobalValues(1, &Row, 1, &Col, &val);
  }

  PyObject* __getitem__(PyObject* args) 
  {
    int Row, Col;
    if (PyInt_Check(args))
    {
      return(Epetra_RowMatrix_GetEntries(*self, PyLong_AsLong(args)));
    }
    else if (PyArg_ParseTuple(args, "ii", &Row, &Col))
    {
      return(Epetra_RowMatrix_GetEntry(*self, Row, Col));
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "Input argument not supported");
      return Py_BuildValue("");
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

  int InsertGlobalValue(int i, int j, double val) {
    double val2 = val;
    int j2 = j;
    return self->InsertGlobalValues(1, &i, 1, &j2, &val2);
  }
}

///////////////////////////////////////
// Epetra_CrsSingletonFilter support //
///////////////////////////////////////
%rename(CrsSingletonFilter) Epetra_CrsSingletonFilter;
%include "Epetra_CrsSingletonFilter.h"

//////////////////////////////
// Epetra_VbrMatrix support //
//////////////////////////////
%ignore Epetra_VbrMatrix::Solve(bool, bool, bool,
				Epetra_Vector const&, Epetra_Vector&) const;
%rename(VbrMatrix) Epetra_VbrMatrix;
%include "Epetra_VbrMatrix.h"

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
%epetra_exception(Epetra_JadMatrix, Epetra_JadMatrix)
%include "Epetra_JadMatrix.h"

//////////////////////////////////
// Epetra_LinearProblem support //
//////////////////////////////////
%rename(LinearProblem) Epetra_LinearProblem;
%include "Epetra_LinearProblem.h"

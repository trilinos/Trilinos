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
#include "Epetra_Util.h"

#include "PyTrilinos_Epetra_Util.h"

%}

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
    int numRowEntries;
    int numColEntries;
    int numValEntries;
    int result = 0;
    PyArrayObject * rowArray = NULL;
    PyArrayObject * colArray = NULL;
    PyArrayObject * valArray = NULL;
    int newRows = 0;
    int newCols = 0;
    int newVals = 0;

    // Create the array of rows
    rowArray = obj_to_array_allow_conversion(Rows, NPY_INT, &newRows);
    if (rowArray == NULL) goto fail;
    numRowEntries = (int) PyArray_MultiplyList(PyArray_DIMS(rowArray),
                                               PyArray_NDIM(rowArray));

    // Create the array of cols
    colArray = obj_to_array_allow_conversion(Cols, NPY_INT, &newCols);
    if (colArray == NULL) goto fail;
    numColEntries = (int) PyArray_MultiplyList(PyArray_DIMS(colArray),
                                               PyArray_NDIM(colArray));

    // Create the array of values
    valArray = obj_to_array_allow_conversion(Values, NPY_DOUBLE, &newVals);
    if (valArray == NULL) goto fail;
    numValEntries = (int) PyArray_MultiplyList(PyArray_DIMS(valArray),
                                               PyArray_NDIM(valArray));

    // Error checking
    if(numValEntries != numColEntries || numValEntries != numRowEntries ||
       numRowEntries != numColEntries)
    {
      PyErr_Format(PyExc_ValueError, 
		   "lengths of Rows, Cols, Values not equal: %d, %d, %d", 
		   numRowEntries, numColEntries, numValEntries);
      goto fail;
    }

    // Loop over the (row,col) pairs, assigning the corresponding
    // value to the matrix
    for(int i = 0 ; i < numValEntries ; ++i)
    {
      double Value = ((double*)PyArray_DATA(valArray))[i];
      int Row = ((int*)PyArray_DATA(rowArray))[i];
      int Col = ((int*)PyArray_DATA(colArray))[i];

      result = self->methodName(Row, 1, &Value, &Col);
      if(result < 0) goto fail;
    }

    // Object cleanup
    if (newRows) { Py_DECREF(rowArray); }
    if (newCols) { Py_DECREF(colArray); }
    if (newVals) { Py_DECREF(valArray); }
    return result;
  fail:
    if (newRows) { Py_XDECREF(rowArray); }
    if (newCols) { Py_XDECREF(colArray); }
    if (newVals) { Py_XDECREF(valArray); }
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
    int numRowEntries;
    int numColEntries;
    int numValEntries;
    int result=0;
    PyArrayObject * rowArray = NULL;
    PyArrayObject * colArray = NULL;
    PyArrayObject * valArray = NULL;
    int newRows = 0;
    int newCols = 0;
    int newVals = 0;

    if (!self->HaveColMap())
    {
      PyErr_SetString(PyExc_RuntimeError, "methodName" " cannot be called on"
		      "className" " that does not have a column map");
      goto fail;
    }

    // Create the array of rows
    rowArray = obj_to_array_allow_conversion(Rows,NPY_INT,&newRows);
    if (rowArray == NULL) goto fail;
    numRowEntries = (int) PyArray_MultiplyList(PyArray_DIMS(rowArray),
                                               PyArray_NDIM(rowArray));

    // Create the array of cols
    colArray = obj_to_array_allow_conversion(Cols,NPY_INT,&newCols);
    if (colArray == NULL) goto fail;
    numColEntries = (int) PyArray_MultiplyList(PyArray_DIMS(colArray),
                                               PyArray_NDIM(colArray));

    // Create the array of values
    valArray = obj_to_array_allow_conversion(Values,NPY_DOUBLE,&newVals);
    if (valArray == NULL) goto fail;
    numValEntries = (int) PyArray_MultiplyList(PyArray_DIMS(valArray),
                                               PyArray_NDIM(valArray));

    // Error checking
    if(numValEntries != numColEntries || numValEntries != numRowEntries || 
       numRowEntries != numColEntries)
    {
      PyErr_Format(PyExc_ValueError, 
		   "lengths of Rows, Cols, Values not equal: %d, %d, %d", 
		   numRowEntries, numColEntries, numValEntries);
      goto fail;
    }

    // Loop over the (row,col) pairs, assigning the corresponding
    // value to the matrix
    for(int i = 0 ; i < numValEntries ; ++i)
    {
      double Value = ((double*)PyArray_DATA(valArray))[i];
      int Row = ((int*)PyArray_DATA(rowArray))[i];
      int Col = ((int*)PyArray_DATA(colArray))[i];

      result = self->methodName(Row, 1, &Value, &Col);
      if(result < 0) goto fail;
    }

    // Object cleanup
    if (newRows) { Py_DECREF(rowArray); }
    if (newCols) { Py_DECREF(colArray); }
    if (newVals) { Py_DECREF(valArray); }
    return result;
  fail:
    if (newRows) { Py_XDECREF(rowArray); }
    if (newCols) { Py_XDECREF(colArray); }
    if (newVals) { Py_XDECREF(valArray); }
    return -1;
  }
}
%enddef

/////////////////////////////
// Epetra_Operator support //
/////////////////////////////
%teuchos_rcp(Epetra_Operator)
%feature("docstring")
Epetra_Operator
"
For cross-language polymorphism to work in python, you must call this
constructor::

    from PyTrilinos import Epetra
    class MyOperator(Epetra.Operator):
        def __init__(self):
            Epetra.Operator.__init__(self)

Other than that, the Epetra.Operator class is much more forgiving than
its C++ counterpart.  Often, you can override just the Label() and
Apply() methods.
"
%feature("autodoc",
"Apply(self, MultiVector x, MultiVector y) -> int

In C++, the Apply() method is pure virtual, thus intended to be
overridden by derived classes.  In python, cross-language polymorphism
is supported, and you are expected to derive classes from this base
class and redefine the Apply() method.  C++ code (e.g., AztecOO
solvers) can call back to your Apply() method as needed.  You must
support two arguments, labeled here MultiVector x and MultiVector y.
These will be converted from Epetra_MultiVector C++ objects to
numpy-hybrid Epetra.MultiVector objects before they are passed to you.
Thus, it is legal to use slice indexing and other numpy features to
compute y from x.

If application of your operator is successful, return 0; else return
some non-zero error code.

It is strongly suggested that you prevent Apply() from raising any
exceptions.  Accidental errors can be prevented by wrapping your code
in a try block:

    try:
        # Your code goes here...
    except Exception, e:
        print 'A python exception was raised by method Apply:'
        print e
        return -1

By returning a -1, you inform the calling routine that Apply() was
unsuccessful.
")
Epetra_Operator::Apply;
%feature("autodoc",
"ApplyInverse(self, MultiVector x, MultiVector y) -> int

In C++, the ApplyInverse() method is pure virtual, thus intended to be
overridden by derived classes.  In python, cross-language polymorphism
is supported, and you are expected to derive classes from this base
class and redefine the ApplyInverse() method.  C++ code (e.g., AztecOO
solvers) can call back to your ApplyInverse() method as needed.  You
must support two arguments, labeled here MultiVector x and MultiVector
y.  These will be converted from Epetra_MultiVector C++ objects to
numpy-hybrid Epetra.MultiVector objects before they are passed to you.
Thus, it is legal to use slice indexing and other numpy features to
compute y from x.

If application of your operator is successful, return 0; else return
some non-zero error code.

It is strongly suggested that you prevent ApplyInverse() from raising
any exceptions.  Accidental errors can be prevented by wrapping your
code in a try block:

    try:
        # Your code goes here...
    except Exception, e:
        print 'A python exception was raised by method ApplyInverse:'
        print e
        return -1

By returning a -1, you inform the calling routine that ApplyInverse()
was unsuccessful.
")
Epetra_Operator::ApplyInverse;
%warnfilter(473)     Epetra_Operator;
%feature("director") Epetra_Operator;
%rename(Operator)    Epetra_Operator;
%include "Epetra_Operator.h"
//////////////////////////////////
// Typemaps for Epetra_Operator //
//////////////////////////////////
#ifdef HAVE_TEUCHOS
%typemap(out) Teuchos::RCP< Epetra_Operator >
{
  if ($1 == Teuchos::null)
  {
    $result = Py_BuildValue("");
  }
  else
  {
    $result = PyTrilinos::convertEpetraOperatorToPython(&$1);
  }
}
%typemap(out) Teuchos::RCP< const Epetra_Operator >
{
  if ($1 == Teuchos::null)
  {
    $result = Py_BuildValue("");
  }
  else
  {
    $result = PyTrilinos::convertEpetraOperatorToPython(&$1);
  }
}
%typemap(directorin) Epetra_Operator &
{
  Teuchos::RCP< Epetra_Operator > *smartinput = new
    Teuchos::RCP< Epetra_Operator >(&$1_name, false);
  $input = PyTrilinos::convertEpetraOperatorToPython(smartinput);
  delete smartinput;
}
#else
%typemap(directorin) Epetra_Operator &
{
  $input = PyTrilinos::convertEpetraOperatorToPython(&$1_name, SWIG_POINTER_OWN);
}
#endif

////////////////////////////////
// Epetra_InvOperator support //
////////////////////////////////
%teuchos_rcp(Epetra_InvOperator)
%warnfilter(473)     Epetra_InvOperator;
%feature("director") Epetra_InvOperator;
%rename(InvOperator) Epetra_InvOperator;
%include "Epetra_InvOperator.h"

//////////////////////////////
// Epetra_RowMatrix support //
//////////////////////////////
%teuchos_rcp(Epetra_RowMatrix)
%feature("autodoc",
"NumMyRowEntries(int myRow, numpy.ndarray numEntries) -> int

In C++, numEntries in an int&.  In python, it is provided to you as a
numpy array of length one so that you can set its value in-place using
numEntries[0] = ....")
Epetra_RowMatrix::NumMyRowEntries;
%feature("autodoc",
"ExtractMyRowCopy(int myRow, int length, numpy.ndarray numEntries,
    numpy.ndarray values, numpy.ndarray indices) -> int

In C++, numEntries in an int&.  In python, it is provided to you as a
numpy array of length one so that you can set its value in-place using
numEntries[0] = ....

Arguments values and indices are double* and int*, respectively, in
C++.  In python, these are provided to you as numpy arrays of the
given length, so that you may alter their entries in-place.")
Epetra_RowMatrix::ExtractMyRowCopy;
%feature("autodoc",
"ExtractDiagonalCopy(Vector diagonal) -> int

Argument diagonal is provided to you as a numpy-hybrid Epetra.Vector,
giving you access to the numpy interface in addition to the
Epetra_Vector C++ interface.")
Epetra_RowMatrix::ExtractDiagonalCopy;
%feature("autodoc",
"Multiply(bool useTranspose, MultiVector x, MultiVector y) -> int

In C++, arguments x and y are Epetra_MultiVectors.  In python, they
are provided to you as numpy-hybrid Epetra.MultiVectors, giving you
access to the numpy interface in addition to the Epetra_MultiVector
C++ interface.")
Epetra_RowMatrix::Multiply;
%feature("autodoc",
"Solve((bool upper, bool trans, bool unitDiagonal, MultiVector x,
    MultiVector y) -> int

In C++, arguments x and y are Epetra_MultiVectors.  In python, they
are provided to you as numpy-hybrid Epetra.MultiVectors, giving you
access to the numpy interface in addition to the Epetra_MultiVector
C++ interface.")
Epetra_RowMatrix::Solve;
%feature("autodoc",
"InvRowSum(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::InvRowSum;
%feature("autodoc",
"LeftScale(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::LeftScale;
%feature("autodoc",
"InvColSums(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::InvColSums;
%feature("autodoc",
"InvRowSums(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::InvRowSums;
%feature("autodoc",
"RightScale(Vector x) -> int

Argument x is provided to you as a numpy-hybrid Epetra.Vector, giving
you access to the numpy interface in addition to the Epetra_Vector C++
interface.")
Epetra_RowMatrix::RightScale;
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
  npy_intp dims$argnum[ ] = { (npy_intp) 1 };
  $input = PyArray_SimpleNewFromData(1, dims$argnum, NPY_INT, (void*)&$1_name);
%}
%typemap(directorin) double *Values
%{
  npy_intp dims$argnum[ ] = { (npy_intp) Length };
  $input = PyArray_SimpleNewFromData(1, dims$argnum, NPY_DOUBLE, (void*)$1_name);
%}
%typemap(directorin) int *Indices
%{
  npy_intp dims$argnum[ ] = { (npy_intp) Length };
  $input = PyArray_SimpleNewFromData(1, dims$argnum, NPY_INT, (void*)$1_name);
%}
%include "Epetra_RowMatrix.h"

///////////////////////////////////
// Epetra_BasicRowMatrix support //
///////////////////////////////////
%teuchos_rcp(Epetra_BasicRowMatrix)
%warnfilter(473)        Epetra_BasicRowMatrix;
%feature("director")    Epetra_BasicRowMatrix;
%rename(BasicRowMatrix) Epetra_BasicRowMatrix;
%ignore Epetra_BasicRowMatrix::ExtractMyEntryView(int,const double*&,int&,int&) const;
// Typemap for double * & Value
%typemap(directorin) double *&Value
%{
  npy_intp dims$argnum[ ] = { (npy_intp) 1 };
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
%teuchos_rcp(Epetra_CrsMatrix)
%teuchos_rcp_epetra_argout(Epetra_CrsMatrix)
%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, int numEntriesPerRow, 
    bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with implicit column map and constant number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    numEntriesPerRow  - constant number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&, int,
				   bool);
%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, PySequence numEntriesPerRow, 
    bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with implicit column map and variable number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    numEntriesPerRow  - variable number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&,
				   const int*, int, bool);
%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, Map colMap, int numEntriesPerRow, 
    bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with specified column map and constant number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    colMap            - describes distribution of columns across processors
    numEntriesPerRow  - constant number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&,
				   const Epetra_Map&, int, bool);
%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, Map rowMap, Map colMap, PySequence
    numEntriesPerRow, bool staticProfile=False) -> CrsMatrix

  CrsMatrix constructor with specified column map and variable number
  of entries per row.  Arguments:

    CV                - Epetra.Copy or Epetra.View
    rowMap            - describes distribution of rows across processors
    colMap            - describes distribution of columns across processors
    numEntriesPerRow  - variable number of entries per row
    staticProfile     - static profile flag
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_Map&,
				   const Epetra_Map&, const int*, int, bool);
%feature("autodoc",
"__init__(self, Epetra_DataAccess CV, CrsGraph graph) -> CrsMatrix

  CrsMatrix constructor with CrsGraph.  Arguments:

    CV     - Epetra.Copy or Epetra.View
    graph  - CrsGraph describing structure of matrix
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess, const Epetra_CrsGraph&);
%feature("autodoc",
"__init__(self, CrsMatrix matrix) -> CrsMatrix

  CrsMatrix copy constructor.  Argument:

    matrix  - source CrsMatrix
"
)
Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix&);
%feature("autodoc",
"ExtractGlobalRowCopy(self, int globalRow) -> (numpy.ndarray,numpy.ndarray)

  Returns a two-tuple of numpy arrays of the same size; the first is
  an array of integers that represent the nonzero columns on the
  matrix; the second is an array of doubles that represent the values
  of the matrix entries.  The input argument is a global row index."
)
Epetra_CrsMatrix::ExtractGlobalRowCopy;
%feature("autodoc",
"ExtractMyRowCopy(self, int myRow) -> (numpy.ndarray,numpy.ndarray)

  Returns a two-tuple of numpy arrays of the same size; the first is
  an array of integers that represent the nonzero columns on the
  matrix; the second is an array of doubles that represent the values
  of the matrix entries.  The input argument is a local row index."
)
Epetra_CrsMatrix::ExtractMyRowCopy;
%feature("autodoc",
"InsertGlobalValues(self, int globalRow, PySequence values, PySequence
    indices) -> int

  Arguments:

    globalRow - global row index
    values    - a sequence of doubles that represent the values to insert
    indices   - a sequence of integers that represent the indices to insert
")
Epetra_CrsMatrix::InsertGlobalValues(int,double*,int,int*,int);
%feature("autodoc",
"InsertGlobalValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to insert
    cols    - a sequence of integers that represent the column indices to
              insert
    values  - a sequence of doubles that represent the values to insert
")
Epetra_CrsMatrix::InsertGlobalValues(PyObject*,PyObject*,PyObject*);
%feature("autodoc",
"InsertMyValues(self, int myRow, PySequence values, PySequence indices) -> int

  Arguments:

    myRow     - local row index
    values    - a sequence of doubles that represent the values to insert
    indices   - a sequence of integers that represent the indices to insert
")
Epetra_CrsMatrix::InsertMyValues(int,double*,int,int*,int);
%feature("autodoc",
"InsertMyValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to insert
    cols    - a sequence of integers that represent the column indices to
              insert
    values  - a sequence of doubles that represent the values to insert
")
Epetra_CrsMatrix::InsertMyValues(PyObject*,PyObject*,PyObject*);
%feature("autodoc",
"ReplaceGlobalValues(self, int globalRow, PySequence values, PySequence
    indices) -> int

  Arguments:

    globalRow - global row index
    values    - a sequence of doubles that represent the values to replace
    indices   - a sequence of integers that represent the indices to replace
")
Epetra_CrsMatrix::ReplaceGlobalValues(int,double*,int,int*,int);
%feature("autodoc",
"ReplaceGlobalValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to replace
    cols    - a sequence of integers that represent the column indices to
              replace
    values  - a sequence of doubles that represent the values to replace
")
Epetra_CrsMatrix::ReplaceGlobalValues(PyObject*,PyObject*,PyObject*);
%feature("autodoc",
"ReplaceMyValues(self, int myRow, PySequence values, PySequence indices) -> int

  Arguments:

    myRow     - local row index
    values    - a sequence of doubles that represent the values to replace
    indices   - a sequence of integers that represent the indices to replace
")
Epetra_CrsMatrix::ReplaceMyValues(int,double*,int,int*,int);
%feature("autodoc",
"ReplaceMyValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to replace
    cols    - a sequence of integers that represent the column indices to
              replace
    values  - a sequence of doubles that represent the values to replace
")
Epetra_CrsMatrix::ReplaceMyValues(PyObject*,PyObject*,PyObject*);
%feature("autodoc",
"SumIntoGlobalValues(self, int globalRow, PySequence values, PySequence
    indices) -> int

  Arguments:

    globalRow - global row index
    values    - a sequence of doubles that represent the values to sum into
    indices   - a sequence of integers that represent the indices to sum into
")
Epetra_CrsMatrix::SumIntoGlobalValues(int,double*,int,int*,int);
%feature("autodoc",
"SumIntoGlobalValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to sum into
    cols    - a sequence of integers that represent the column indices to
              sum into
    values  - a sequence of doubles that represent the values to sum into
")
Epetra_CrsMatrix::SumIntoGlobalValues(PyObject*,PyObject*,PyObject*);
%feature("autodoc",
"SumIntoMyValues(self, int myRow, PySequence values, PySequence indices) -> int

  Arguments:

    myRow     - local row index
    values    - a sequence of doubles that represent the values to sum into
    indices   - a sequence of integers that represent the indices to sum into
")
Epetra_CrsMatrix::SumIntoMyValues(int,double*,int,int*,int);
%feature("autodoc",
"SumIntoMyValues(self, PySequence rows, PySequence cols, PySequence
    values) -> int

  Arguments:

    rows    - a sequence of integers that represent the row indices to sum into
    cols    - a sequence of integers that represent the column indices to
              sum into
    values  - a sequence of doubles that represent the values to sum into
")
Epetra_CrsMatrix::SumIntoValues(PyObject*,PyObject*,PyObject*);
%feature("autodoc",
"__setitem__(self, PyTuple index, double val)

The __setitem__() method is called when square-bracket indexing is
used to set a value of the matrix.  For example, the last line of::

    comm = Epetra.SerialComm()
    m = Epetra.CrsMatrix(9,0,comm)
    m[0,0] = 3.14

calls::

    m.__setitem__((0,0), 3.14)

Thus, argument 'index' is a tuple filled with whatever indices you
give the square-bracket operator when setting.  For __setitem__(),
this raises an IndexError unless 'index' is a two-tuple of integers.
Argument 'val' must be convertible to a double.  Under the covers,
__setitem__() calls ReplaceGlobalValues() or InsertGlobalValues() as
necessary, so the indices are expected to be global IDs.  Note that if
you use __setitem__() to insert a new matrix element, you will need to
call FillComplete() again, whether or not you have called it before.
")
Epetra_CrsMatrix::__setitem__;
%feature("autodoc",
"
__getitem__(self, PyTuple index) -> double
__getitem__(self, int row) -> numpy.ndarray

The __getitem__() method is called when square-bracket indexing is
used to get a value from the matrix.  For example, the last two lines
of::

    comm = Epetra.SerialComm()
    m = Epetra.CrsMatrix(9,0,comm)
    m.InsertGlobalValues(0, [0.0, 1.0, 2.0], [0,1,2])
    diag = m[0,0]
    row  = m[0]

call::

    m.__getitem__((0,0))
    m.__getitem__(0)

The __getitem__() method behaves according to the following table:

                    FillComplete()    #    
    Index               called      procs  Return value
    --------------  --------------  -----  ---------------------------
    single integer       true        any   numpy array of doubles
    single integer       false        1    numpy array of doubles
    single integer       false       >1    raise IndexError
    two integers         either      any   double

You should provide global IDs as the integer indices if FillComplete()
has been called.  If not, you should provide local IDs.  If you
reference a matrix element that is off-processor, __getitem__() will
raise an IndexError.

Under the covers, __getitem__() will call ExtractGlobalRowView() if
FillComplete() has been called, or ExtractMyRowView() if it has not.
If either of these return a non-zero return code, this is converted to
a python RuntimeError.  The resulting data is copied to the output
array.
")
Epetra_CrsMatrix::__getitem__;
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
    npy_intp   dimensions[ ] = { 0 };
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
    npy_intp   dimensions[ ] = { 0 };
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
    npy_intp   dimensions[ ] = { 0 };
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
%teuchos_rcp(Epetra_FECrsMatrix)
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

  int InsertGlobalValues(const int row,
                         const int size, 
                         const Epetra_SerialDenseVector & values,
                         const Epetra_IntSerialDenseVector & entries)
  {
    return self->InsertGlobalValues(1, &row, size, (int*)entries.Values(),
                                    values.Values());
  }

  int InsertGlobalValues(const Epetra_IntSerialDenseVector & rows,
                         const Epetra_IntSerialDenseVector & cols,
                         const Epetra_SerialDenseMatrix & values,
                         int format=Epetra_FECrsMatrix::COLUMN_MAJOR)
  {
    return self->InsertGlobalValues(rows, cols, values, format);
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
%teuchos_rcp(Epetra_VbrMatrix)
%teuchos_rcp_epetra_argout(Epetra_VbrMatrix)
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, int
    numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with implicit column map and constant number
  of block entries per row.
")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
			           const Epetra_BlockMap&,
				   int);
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, PySequence
    numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with implicit column map and variable number
  of block entries per row.

")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
				   const Epetra_BlockMap&,
				   int*, int);
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    int numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with specified column map and constant number
  of block entries per row.

")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
			           const Epetra_BlockMap&,
			           const Epetra_BlockMap&,
				   int);
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, BlockMap rowMap, BlockMap colMap,
    PySequence numBlockEntriesPerRow) -> VbrMatrix

  VbrMatrix constructor with specified column map and variable number
  of block entries per row.

")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
				   const Epetra_BlockMap&,
				   const Epetra_BlockMap&,
				   int*, int);
%feature("autodoc",
"
__init__(self, Epetra_DataAccess CV, CrsGraph graph) -> VbrMatrix

  CrsGraph constructor.
")
Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess,
				   const Epetra_CrsGraph&);
%feature("autodoc",
"
__init__(self, VbrMatrix matrix) -> VbrMatrix

  Copy constructor.
")
Epetra_VbrMatrix::Epetra_VbrMatrix(const Epetra_VbrMatrix&);
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
%teuchos_rcp(Epetra_FEVbrMatrix)
%rename(FEVbrMatrix) Epetra_FEVbrMatrix;
%include "Epetra_FEVbrMatrix.h"

//////////////////////////////
// Epetra_JadMatrix support //
//////////////////////////////
%teuchos_rcp(Epetra_JadMatrix)
%ignore Epetra_JadMatrix::ExtractMyEntryView(int,double*&,int&,int&);
%rename(JadMatrix) Epetra_JadMatrix;
%include "Epetra_JadMatrix.h"

//////////////////////////////////
// Epetra_LinearProblem support //
//////////////////////////////////
%rename(LinearProblem) Epetra_LinearProblem;
%include "Epetra_LinearProblem.h"

/////////////////////////
// Epetra_Util support //
/////////////////////////
%feature("docstring")
Epetra_Util
"Epetra Util Wrapper Class.

The Epetra.Util class is a collection of useful functions that cut
across a broad set of other classes.  A random number generator is
provided, along with methods to set and retrieve the random-number
seed.

The random number generator is a multiplicative linear congruential
generator, with multiplier 16807 and modulus 2^31 - 1. It is based on
the algorithm described in 'Random Number Generators: Good Ones Are
Hard To Find', S. K. Park and K. W. Miller, Communications of the ACM,
vol. 31, no. 10, pp. 1192-1201.

The C++ Sort() method is not supported in python.

A static function is provided for creating a new Epetra.Map object
with 1-to-1 ownership of entries from an existing map which may have
entries that appear on multiple processors.

Epetra.Util is a serial interface only.  This is appropriate since the
standard utilities are only specified for serial execution (or shared
memory parallel)."
%ignore Epetra_Util::Sort;
%rename(Util) Epetra_Util;
%include "Epetra_Util.h"

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

%define EPETRA_DOCSTRING
"The Epetra module allows access to The Trilinos package Epetra.  Note
that the 'Epetra_' prefix has been stripped from all Epetra objects,
but that if imported with 'from PyTrilinos import Epetra', these
objects exist in the 'Epetra' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.


*) Brief Description

The Epetra module offers a suite of serial and distributed linear algebra
objects, like vectors, multivectors, sparse matrices. Serial dense vectors and
matrices are also supported. Using Epetra, several BLAS and LAPACK routines
can be easily accessed through Python. The Epetra communicator is the basic
intra-processor class of Trilinos.

The most important classes of the Epetra module is:
- PyComm: basic serial or parallel communicator
- SerialDenseVector: serial vector, to be used with a serial dense matrix
- SerialDenseMatrix: serial matrix, allocated as dense
- SerialDenseProblem: to solve serial dense problems
- Map: layout of distributed objects across processors
- Vector: distributed vector
- MultiVector: series of Vector's with the same Map
- CrsMatrix: distributed matrix with row access
- LinearProblem: to solve distributed linear systems
- Time: utility class to track CPU times
- Import: to easily redistribute distributed objects


*) Example of usage: Creating a Distributed Matrix

The following example builds a distributed tridiagonal matrix. The
matrix has size 10 and the elements are distributed linearly among the
processors.  It can be run is serial or parallel, depending on how
Trilinos was configured.

from PyTrilinos import Epetra
Epetra.Init()
Comm  = Epetra.PyComm()
NumGlobalRows = 10
Map   = Epetra.Map(NumGlobalRows, 0, Comm)
A     = Epetra.CrsMatrix(Epetra.Copy, Map, 0);
NumMyRows = Map.NumMyElements()
# Loop over all local rows to create a tridiagonal matrix
for ii in xrange(0, NumMyRows):
  # `i' is the global ID of local ID `ii'
  i = Map.GID(ii)
  if i != NumGlobalRows - 1:
    Indices = [i, i + 1]
    Values = [2.0, -1.0]
  else:
    Indices = [i]
    Values = [2.0];
  A.InsertGlobalValues(i, Values, Indices);
# transform the matrix into local representation -- no entries
# can be added after this point. However, existing entries can be
# modified.
ierr = A.FillComplete();
Epetra.Finalize()

*) Notes

An Epetra.Vector object is at the same time a Numeric vector, and it
can therefore be used everywhere Numeric vectors are accepted.
"
%enddef

%module(package="PyTrilinos", directors="1", docstring=EPETRA_DOCSTRING) Epetra

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "Epetra_ConfigDefs.h"

// Epetra includes
#include "Epetra_Version.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_LAPACK.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"
#include "Epetra_DataAccess.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseOperator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"

// Local includes
#include "FileStream.h"
#include "Epetra_NumPyVector.h"
#include "NumPyArray.h"
#include "NumPyWrapper.h"

#ifdef HAVE_MPI
#include "mpi.h"
PyObject* Init_Argv(PyObject *args) 
{  
  /*
  int flag = 0;
  // check that MPI_Init has not already been called 
  MPI_Initialized(&flag);
  if (flag)
    return Py_BuildValue("");
   */

  int i, error, myid, size;
  int argc = 0;  
  char **argv;   

  /* Reconstruct C-commandline */     
  /*                           */ 
  argc = PyList_Size(args); //Number of commandline arguments
  argv = (char**) malloc((argc+1)*sizeof(char*)); 
  
  for (i=0; i<argc; i++)  
    argv[i] = PyString_AsString( PyList_GetItem(args, i) );
    
  argv[i] = NULL; //Lam 7.0 requires last arg to be NULL  
  
  error = MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
  MPI_Comm_size(MPI_COMM_WORLD, &size);    

  if (error != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
    PyErr_SetString(PyExc_RuntimeError, "error");   
    return NULL;
  }  

  return Py_BuildValue("");  
} 
  
PyObject* Finalize() {  
  int error, myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  
  
  error = MPI_Finalize();
  if (error != 0) {
    PyErr_SetString(PyExc_RuntimeError, "error");    //raise ValueError, errmsg
    return NULL;
  }  
    
  return Py_BuildValue("");
} 

MPI_Comm CommWorld()
{
  return(MPI_COMM_WORLD);
}
#else
PyObject* Init_Argv(PyObject *args) 
{
  return Py_BuildValue("");
}

PyObject* Finalize()
{
  return Py_BuildValue("");
}
#endif

#define SUMALL    0
#define MINALL    1
#define MAXALL    2
#define GATHERALL 3
#define SCANSUM   4
#define BROADCAST 5

PyObject* Epetra_Comm_InterfaceInt(const Epetra_Comm& Comm,
                                   int action, PyObject* input,
                                   int root)
{
  PyObject* result;

  int len = PyList_Size(input);
  vector<int> orig(len);
  vector<int> dest(len);

  for (int i = 0 ; i < len ; ++i)
  {
    PyObject* value;
    value = PyList_GetItem(input, i);

    if (PyInt_Check(value) == 0)
    {
      cerr << "Indices must be integers" << endl;
      Py_INCREF(Py_None);
      return Py_None;
    }
    orig[i] = PyLong_AsLong(value);
  }

  // communicate here
  switch (action) {
  case SUMALL:
    Comm.SumAll(&orig[0], &dest[0], len);
    break;
  case MINALL:
    Comm.MinAll(&orig[0], &dest[0], len);
    break;
  case MAXALL:
    Comm.MaxAll(&orig[0], &dest[0], len);
    break;
  case GATHERALL:
    Comm.GatherAll(&orig[0], &dest[0], len);
    break;
  case SCANSUM:
    Comm.ScanSum(&orig[0], &dest[0], len);
    break;
  case BROADCAST:
    Comm.Broadcast(&orig[0], len, root);
    for (int i = 0 ; i < len ; ++i)
      dest[i] = orig[i];
    break;
  }

  // now insert in the result
  result = PyList_New(len);
  for (int i = 0 ; i < len ; ++i)
  {
    PyObject* value;
    value = PyLong_FromLong(dest[i]);
    PyList_SetItem(result, i, value);
  }
  return(result);
}

PyObject* Epetra_Comm_InterfaceDouble(const Epetra_Comm& Comm,
                                      int action, PyObject* input,
                                      int root)
{
  PyObject* result;

  int len = PyList_Size(input);
  vector<double> orig(len);
  vector<double> dest(len);

  for (int i = 0 ; i < len ; ++i)
  {
    PyObject* value;
    value = PyList_GetItem(input, i);

    if (PyFloat_Check(value) == 0)
    {
      cerr << "Indices must be doubles" << endl;
      Py_INCREF(Py_None);
      return Py_None;
    }
    orig[i] = PyFloat_AsDouble(value);
  }

  // communicate here
  switch (action) {
  case SUMALL:
    Comm.SumAll(&orig[0], &dest[0], len);
    break;
  case MINALL:
    Comm.MinAll(&orig[0], &dest[0], len);
    break;
  case MAXALL:
    Comm.MaxAll(&orig[0], &dest[0], len);
    break;
  case GATHERALL:
    Comm.GatherAll(&orig[0], &dest[0], len);
    break;
  case SCANSUM:
    Comm.ScanSum(&orig[0], &dest[0], len);
    break;
  case BROADCAST:
    Comm.Broadcast(&orig[0], len, root);
    for (int i = 0 ; i < len ; ++i)
      dest[i] = orig[i];
    break;
  }

  // now insert in the result
  result = PyList_New(len);
  for (int i = 0 ; i < len ; ++i)
  {
    PyObject* value;
    value = PyFloat_FromDouble(dest[i]);
    PyList_SetItem(result, i, value);
  }
  return(result);
}

PyObject * Epetra_RowMatrix_GetEntries(const Epetra_RowMatrix& Matrix,
                                       int GlobalRow) 
{
  if (!Matrix.Filled())
  {
    cout << "Matrix not FillComplete()'d" << endl;
    Py_INCREF(Py_None);
    return Py_None;
  }
  int ierr;
  PyObject* res;
  PyObject* PyIndices;
  PyObject* PyValues;
  int NumEntries, Length;
  Length = Matrix.MaxNumEntries();
  vector<int>    Indices(Length);
  vector<double> Values(Length);
  int MyRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  ierr = Matrix.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0],
                                 &Indices[0]);
  if (ierr < 0)
  {
    Py_INCREF(Py_None);
    return Py_None;
  }

  PyIndices = PyList_New(NumEntries);
  PyValues  = PyList_New(NumEntries);

  // return global indices
  for (int i = 0 ; i < NumEntries ; ++i)
  {
    int GID = Matrix.RowMatrixColMap().GID(Indices[i]);
    PyList_SetItem(PyIndices, i, PyInt_FromLong(GID));
    PyList_SetItem(PyValues,  i, PyFloat_FromDouble(Values[i]));
  }
  res = PyTuple_New(2);
  PyTuple_SetItem(res, 0, PyIndices);
  PyTuple_SetItem(res, 1, PyValues);
  return(res);
}

PyObject* Epetra_RowMatrix_GetEntry(Epetra_RowMatrix& Matrix, 
                                    int GlobalRow, int GlobalCol)
{
  int NumEntries, Length, ierr;
  double val = 0.0;
  Length = Matrix.MaxNumEntries();
  vector<int>    Indices(Length);
  vector<double> Values(Length);
  int MyRow = Matrix.RowMatrixRowMap().LID(GlobalRow);
  int MyCol = Matrix.RowMatrixColMap().LID(GlobalCol);

  ierr = Matrix.ExtractMyRowCopy(MyRow, Length, NumEntries, &Values[0], 
                                 &Indices[0]);
  if (ierr) NumEntries = 0;
  for (int i = 0 ; i < Length ; ++i)
  {
    if (Indices[i] == MyCol)
    {
      val = Values[i];
      break;
    }
  }
  return(PyFloat_FromDouble(val));
}

extern "C" {
  // on some MAC OS X with LAM/MPI _environ() is not found,
  // need to specify -Wl,-i_environ:_fake_environ as LDFLAGS
  void fake_environ()
  {
    exit(EXIT_FAILURE);
  }
}

#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"
%}

%feature("director") PyOperator;
%feature("director") PyRowMatrix;

// Ignore directives
%ignore operator<<(ostream &, const Epetra_Object &);// From python, use __str__
%ignore Epetra_Object::Print(ostream &) const;
%ignore Epetra_SerialComm::operator=(const Epetra_SerialComm &);
#ifdef HAVE_MPI
%ignore Epetra_MpiComm::operator=(const Epetra_MpiComm &);
#endif
%ignore Epetra_CompObject::operator=(const Epetra_CompObject &);
%ignore Epetra_CompObject::UpdateFlops(int) const;   // Use long int version
%ignore Epetra_CompObject::UpdateFlops(float) const; // Use double version
%ignore Epetra_BlockMap::operator=(const Epetra_BlockMap &);
%ignore Epetra_Map::operator=(const Epetra_Map &);
%ignore Epetra_LocalMap::operator=(const Epetra_LocalMap &);
%ignore Epetra_MultiVector::operator=(const Epetra_MultiVector &);
%ignore Epetra_MultiVector::operator[](int);         // See %extend MultiVector
%ignore Epetra_MultiVector::operator[](int) const;   //       __getitem__()
%ignore Epetra_MultiVector::operator()(int) const;
%ignore Epetra_IntVector::operator=(const Epetra_IntVector &);
%ignore Epetra_IntVector::operator[](int);           // See %extend IntVector
%ignore Epetra_IntVector::operator[](int) const;     //       __getitem__()
%ignore Epetra_CrsGraph::operator=(const Epetra_CrsGraph &);
%ignore Epetra_CrsGraph::operator[](int);            // See %extend CrsGraph
%ignore Epetra_CrsGraph::operator[](int) const;      //       __getitem__()
%ignore Epetra_MapColoring::operator[](int);         // See %extend Mapcoloring
%ignore Epetra_MapColoring::operator[](int) const;   //       __getitem__()
%ignore Epetra_MapColoring::operator()(int) const;
%ignore Epetra_CrsMatrix::operator=(const Epetra_CrsMatrix &);
%ignore Epetra_CrsMatrix::operator[](int);           // See %extend CrsMatrix
%ignore Epetra_CrsMatrix::operator[](int) const;     //       __getitem__()
%ignore Epetra_FECrsMatrix::operator=(const Epetra_FECrsMatrix &);
%ignore Epetra_FECrsMatrix::operator[](int);           // See %extend FECrsMatrix
%ignore Epetra_FECrsMatrix::operator[](int) const;     //       __getitem__()
// Epetra_VbrMatrix member function
// Solve(bool,bool,bool,Epetra_Vector const&,Epetra_Vector&) const does not
// appear to be implemented.  Apparently overridden by
// Solve(bool,bool,bool,Epetra_MultiVector const&,Epetra_MultiVector&) const
%ignore Epetra_VbrMatrix::operator=(const Epetra_VbrMatrix &);
%ignore Epetra_VbrMatrix::Solve(bool, bool, bool,
				Epetra_Vector const&, Epetra_Vector&) const;
%ignore Epetra_SerialDenseMatrix::operator=(const Epetra_SerialDenseMatrix &);
%ignore Epetra_SerialDenseMatrix::operator[](int);
%ignore Epetra_SerialDenseMatrix::operator[](int) const;
%ignore Epetra_SerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_SerialDenseMatrix::A() const;
%ignore Epetra_SerialDenseVector::operator=(const Epetra_SerialDenseVector &);
%ignore Epetra_SerialDenseVector::operator[](int);
%ignore Epetra_SerialDenseVector::operator[](int) const;
%ignore Epetra_SerialDenseVector::operator()(int);
%ignore Epetra_SerialDenseVector::operator()(int) const;
%ignore Epetra_IntSerialDenseMatrix::operator=(const Epetra_IntSerialDenseMatrix &);
%ignore Epetra_IntSerialDenseMatrix::operator[](int);
%ignore Epetra_IntSerialDenseMatrix::operator[](int) const;
%ignore Epetra_IntSerialDenseMatrix::operator()(int,int) const;
%ignore Epetra_IntSerialDenseMatrix::A() const;
%ignore Epetra_IntSerialDenseVector::operator=(const Epetra_IntSerialDenseVector &);
%ignore Epetra_IntSerialDenseVector::operator[](int);
%ignore Epetra_IntSerialDenseVector::operator[](int) const;
%ignore Epetra_IntSerialDenseVector::operator()(int);
%ignore Epetra_IntSerialDenseVector::operator()(int) const;
%ignore NumPyArrayBase::print(std::ostream &) const;       // faciltated by __str__
%ignore NumPyArray::print(std::ostream &) const;           // faciltated by __str__
%ignore NumPyArrayContiguous::print(std::ostream &) const; // faciltated by __str__
%ignore NumPyArrayBase::getDataArray() const;
%ignore NumPyArrayBase::getArrayObject() const;
%ignore Epetra_LinearProblem::SetOperator()(int) const;
%ignore Epetra_Time::operator=(const Epetra_Time &);

// Rename directives
%rename(Version             ) Epetra_Version;
%rename(Object              ) Epetra_Object;
%rename(Comm                ) Epetra_Comm;
%rename(SerialComm          ) Epetra_SerialComm;
#ifdef HAVE_MPI
%rename(MpiComm             ) Epetra_MpiComm;
#endif
%rename(BlockMap            ) Epetra_BlockMap;
%rename(Map                 ) Epetra_Map;
%rename(LocalMap            ) Epetra_LocalMap;
%rename(SrcDistObject       ) Epetra_SrcDistObject;
%rename(DistObject          ) Epetra_DistObject;
%rename(CompObject          ) Epetra_CompObject;
%rename(BLAS                ) Epetra_BLAS;
%rename(LAPACK              ) Epetra_LAPACK;
%rename(MultiVector         ) Epetra_MultiVector;
%rename(IntVector           ) Epetra_IntVector;
%rename(CrsGraph            ) Epetra_CrsGraph;
%rename(MapColoring         ) Epetra_MapColoring;
%rename(Operator            ) Epetra_Operator;
%rename(RowMatrix           ) Epetra_RowMatrix;
%rename(CrsMatrix           ) Epetra_CrsMatrix;
%rename(FECrsMatrix         ) Epetra_FECrsMatrix;
%rename(VbrMatrix           ) Epetra_VbrMatrix;
%rename(SerialDenseSolver   ) Epetra_SerialDenseSolver;
%rename(SerialDenseOperator ) Epetra_SerialDenseOperator;
%rename(SerialDenseMatrix   ) Epetra_SerialDenseMatrix;
%rename(SerialDenseVector   ) Epetra_SerialDenseVector;
%rename(IntSerialDenseMatrix) Epetra_IntSerialDenseMatrix;
%rename(IntSerialDenseVector) Epetra_IntSerialDenseVector;
%rename(SerialSymDenseMatrix) Epetra_SerialSymDenseMatrix;
%rename(LinearProblem       ) Epetra_LinearProblem;
%rename(Import              ) Epetra_Import;
%rename(Export              ) Epetra_Export;
%rename(Time                ) Epetra_Time;
%rename(NumPyVector         ) Epetra_NumPyVector;

// Auto-documentation feature
%feature("autodoc", "1");

// Typemap directives
%typemap(in) (int * Indices)
{
  PyArrayObject * array = NumPyWrapper::contiguous_typed_array($input, PyArray_INT,
							       0, 0);
  assert (array && "Function should have already checked this.");
  $1 = (int *) array->data;
}

%typemap(argout) (int & NumIndices, int *& Indices)
{
  // Decrement the current result object
  Py_XDECREF($result);
  // Check for bad result
  if (result == -1)
    PyErr_SetString(PyExc_RuntimeError, "Invalid row index");
  if (result == -2)
    PyErr_SetString(PyExc_RuntimeError, "Graph not completed");
  // Copy the indices into a python tuple
  $result = PyTuple_New(*$1);
  for (int i=0; i<*$1; ++i)
    PyTuple_SetItem($result, i, PyInt_FromLong((long) *$2[i])); 
}

// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

// Epetra interface includes
using namespace std;
%include "Epetra_Version.h"
%include "Epetra_Object.h"
%include "Epetra_Comm.h"
%include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
%include "Epetra_MpiComm.h"
#endif
%include "Epetra_BlockMap.h"
%include "Epetra_Map.h"
%include "Epetra_LocalMap.h"
%include "Epetra_SrcDistObject.h"
%include "Epetra_DistObject.h"
%include "Epetra_CompObject.h"
%include "Epetra_BLAS.h"
%include "Epetra_LAPACK.h"
%include "Epetra_MultiVector.h"
%include "Epetra_Vector.h"
%include "Epetra_IntVector.h"
%include "Epetra_CrsGraph.h"
%include "Epetra_MapColoring.h"
%include "Epetra_Operator.h"
%include "Epetra_RowMatrix.h"
%include "Epetra_CrsMatrix.h"
%include "Epetra_FECrsMatrix.h"
%include "Epetra_VbrMatrix.h"
%include "Epetra_DataAccess.h"
%include "Epetra_SerialDenseSolver.h"
%include "Epetra_SerialDenseOperator.h"
%include "Epetra_SerialDenseMatrix.h"
%include "Epetra_SerialDenseVector.h"
%include "Epetra_SerialSymDenseMatrix.h"
%include "Epetra_IntSerialDenseMatrix.h"
%include "Epetra_IntSerialDenseVector.h"
%include "Epetra_LinearProblem.h"
%include "Epetra_CombineMode.h"
%include "Epetra_Import.h"
%include "Epetra_Export.h"
%include "Epetra_Time.h"

// Local interface includes
%include "NumPyArray.h"
%include "Epetra_NumPyVector.h"

// Extensions
%extend Epetra_Object {

  // Define the __str__() method, used by the python str() operator on any
  // object given to the python print command.
  string __str__() {
    stringstream os;
    self->Print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }

  // The Epetra_Object::Print(ostream) method is ignored and replaced by a
  // Print() method here that takes a python file as its argument.  If no
  // argument is given, then output is to standard out.
  PyObject * Print(PyObject*pf=NULL) const {
    if (pf == NULL) {
      self->Print(std::cout);
    } else {
      if (!PyFile_Check(pf)) return NULL;

      std::FILE*   f = PyFile_AsFile(pf);
      FileStream   buffer(f);
      std::ostream os(&buffer);
      self->Print(os);
    }
    Py_INCREF(Py_None);
    return Py_None;
  }
}

%extend Epetra_MultiVector {
  double * & __getitem__(int i) {
    return self->operator[](i);
  }

  PyObject * Norm1() {
    int n = self->NumVectors();
    double result[n];
    int numVectors[1] = {n};
    int status        = self->Norm1(result);
    PyObject * output = Py_BuildValue("(iO)", status, 
				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
							      (char *)result));
    return output;
  }

  PyObject * __setitem__(PyObject * args, double val) {
    int i, j;
    if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    (*self)[i][j] = val;
    Py_INCREF(Py_None);
    return Py_None;
  }

  void Set(const int vector, const int element, const double value) {
    (*self)[vector][element] = value;
  }

  PyObject * __getitem__(PyObject * args) {
    int i, j;
    if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    double val = (*self)[i][j];
    return PyFloat_FromDouble(val);
  }

  double Get(const int vector, const int element) {
    return((*self)[vector][element]);
  }

  PyObject * Norm2() {
    int n = self->NumVectors();
    double result[n];
    int numVectors[1] = {n};
    int status        = self->Norm2(result);
    PyObject * output = Py_BuildValue("(iO)", status,
				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
							      (char *)result));
    return output;
  }

  PyObject * NormInf() {
    int n = self->NumVectors();
    double result[n];
    int numVectors[1] = {n};
    int status        = self->NormInf(result);
    PyObject * output = Py_BuildValue("(iO)", status,
				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
							      (char *)result));
    return output;
  }

   PyObject * Dot(const Epetra_MultiVector& A){
    int n = self->NumVectors();
    double result[n];
    int numVectors[1] = {n};
    int status        = self->Dot(A, result);
    PyObject * output = Py_BuildValue("(iO)", status,
				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
							      (char *)result));
    return output;
   }
}

%extend Epetra_Vector{

  PyObject * Norm1() {
    double result[1];
    int status        = self->Norm1(result);
    PyObject * output = Py_BuildValue("(id)", status, result[0]);
    return output;
  }

 PyObject * Norm2() {
    double result[1];
    int status        = self->Norm2(result);
    PyObject * output = Py_BuildValue("(id)", status, result[0]);
    return output;
  }

  PyObject * NormInf() {
    double result[1];
    int status        = self->NormInf(result);
    PyObject * output = Py_BuildValue("(id)", status, result[0]);    
    return output;
  }
   PyObject * Dot(const Epetra_MultiVector& A){
    double result[1];
    int status        = self->Dot(A, result);
    PyObject * output = Py_BuildValue("(id)", status, result[0]);
    return output;
   }
}

%extend Epetra_IntVector {
  int & __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_CrsGraph {
  int * __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_MapColoring {
  int & __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_CrsMatrix {
  void __setitem__(PyObject* args, double val) 
  {
    int Row, Col;
    if (!PyArg_ParseTuple(args, "ii", &Row, &Col)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return;
    }

    if (self->ReplaceGlobalValues(Row, 1, &val, &Col))
      self->InsertGlobalValues(Row, 1, &val, &Col);
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
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  int InsertGlobalValues(const int Row, const int Size, 
                         const Epetra_SerialDenseVector& Values,
                         const Epetra_IntSerialDenseVector& Entries)
  {
    return self->InsertGlobalValues(Row, Size, Values.Values(), (int*)Entries.Values());
  }

  int InsertGlobalValue(int i, int j, double val) {
    double val2 = val;
    int j2 = j;
    return self->InsertGlobalValues(i, 1, &val2, &j2);
  }


  int InsertGlobalValues(const int row, PyObject* Values, PyObject* Indices)
  {
    if (row < 0)
      return(-1);

    if (PyList_Check(Values) == 0 || PyList_Check(Indices) == 0) 
    {
      cerr << "Input object is not a list" << endl;
      return(-1);
    }

    int len = PyList_Size(Values);
    if (len != PyList_Size(Indices))
    {
      cerr << "Length of input lists differ" << endl;
      return(-1);
    }

    for (int i = 0 ; i < len ; ++i)
    {
      PyObject* Value,* Index;
      Value = PyList_GetItem(Values, i);
      Index = PyList_GetItem(Indices, i);

      if (PyInt_Check(Index) == 0)
      {
        cerr << "Indices must be integers" << endl;
        return(-1);
      }

      if (PyFloat_Check(Value) == 0)
      {
        cerr << "Values must be doubles" << endl;
        return(-1);
      }

      int cIndex = PyLong_AsLong(Index);
      double cValue = PyFloat_AsDouble(Value);
      if (self->InsertGlobalValues(row, 1, &cValue, &cIndex) < 0)
        return(-1);
    }
    return(0);
  }
}

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
      Py_INCREF(Py_None);
      return Py_None;
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

  int InsertGlobalValues(const int row, PyObject* Values, PyObject* Indices)
  {
    if (row < 0)
      return(-1);

    if (PyList_Check(Values) == 0 || PyList_Check(Indices) == 0) 
    {
      cerr << "Input object is not a list" << endl;
      return(-1);
    }

    int len = PyList_Size(Values);
    if (len != PyList_Size(Indices))
    {
      cerr << "Length of input lists differ" << endl;
      return(-1);
    }

    for (int i = 0 ; i < len ; ++i)
    {
      PyObject* Value,* Index;
      Value = PyList_GetItem(Values, i);
      Index = PyList_GetItem(Indices, i);

      if (PyInt_Check(Index) == 0)
      {
        cerr << "Indices must be integers" << endl;
        return(-1);
      }

      if (PyFloat_Check(Value) == 0)
      {
        cerr << "Values must be doubles" << endl;
        return(-1);
      }

      int cIndex = PyLong_AsLong(Index);
      double cValue = PyFloat_AsDouble(Value);
      if (self->InsertGlobalValues(1, &row, 1, &cIndex, &cValue) < 0)
        return(-1);
    }
    return(0);
  }
}

%extend Epetra_SerialDenseMatrix {

  double * __getitem__(int i) {
    return self->operator[](i);
  }

  PyObject * __getitem__(PyObject * args) {
    int i, j;
    if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    double * column = self->operator[](j);
    return PyFloat_FromDouble(column[i]);
  }

  PyObject * __setitem__(PyObject * args, double val) {
    int i, j;
    if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
      PyErr_SetString(PyExc_IndexError, "Invalid index");
      return NULL;
    }
    double * column = self->operator[](j);
    column[i] = val;
    Py_INCREF(Py_None);
    return Py_None;
  }

}

%extend Epetra_SerialDenseVector {

  double __call__(int i) {
    return self->operator()(i);
  }

  double __getitem__(int i) {
    return self->operator[](i);
  }

  void __setitem__(int i, const double val) {
    double * column = self->Values();
    column[i] = val;
  }
}

%extend Epetra_IntSerialDenseVector {

  int __call__(int i) {
    return self->operator()(i);
  }

  int __getitem__(int i) {
    return self->operator[](i);
  }

  void __setitem__(int i, const int val) {
    int * column = self->Values();
    column[i] = val;
  }
}

%extend Epetra_Map 
{
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

%extend NumPyArrayBase {
  using namespace std;
  string __str__() {
    stringstream os;
    self->print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

// must be as defined above
#define SUMALL    0
#define MINALL    1
#define MAXALL    2
#define GATHERALL 3
#define SCANSUM   4
#define BROADCAST 5

#ifndef HAVE_MPI
%extend Epetra_SerialComm {
#else
%extend Epetra_MpiComm {
#endif
  // works for INTs and DOUBLEs, like Epetra_Comm
  PyObject* GlobalOp(int action, PyObject* input, int root = 0)
  {
    char what;

    if (PyList_Check(input) == 0)
    {
      cerr << "Input object is not a list" << endl;
      Py_INCREF(Py_None);
      return Py_None;
    }

    // get the first object to understand the type
    PyObject* value = PyList_GetItem(input, 0);

    if (PyInt_Check(value))
      what = 'I';
    else if (PyFloat_Check(value))
      what = 'D';
    else
    {
      cerr << "Only integers and doubles are supported" << endl;
      Py_INCREF(Py_None);
      return Py_None;
    }

    PyObject* result;
    if (what == 'I')
      result = Epetra_Comm_InterfaceInt(*self, action, input, root);
    else
      result = Epetra_Comm_InterfaceDouble(*self, action, input, root);
    return(result);
  }
}

// MPI stuff
PyObject* Init_Argv(PyObject *args);
PyObject* Finalize();
#ifdef HAVE_MPI
MPI_Comm CommWorld();
#endif

// Python code.  Here we set the __version__ string, call MPI_Init()
// and arrange for MPI_Finalize to be called (if appropriate), and
// declare classes that inherit both from Epetra objects and
// UserArrays, to give users additional functionality.
%pythoncode %{

__version__ = Version().split()[2]

# Call MPI_Init if appropriate
import sys
Init_Argv(sys.argv)
del sys

# Arrange for MPI_Finalize to be called at exit, if appropriate
import atexit
atexit.register(Finalize)

from UserArray import *

class Vector(UserArray,NumPyVector):
    def __init__(self, *args):
        NumPyVector.__init__(self, *args)
        UserArray.__init__(self,self.getArray(),'d',0,1)
    def __str__(self):
        return str(self.array)
    def __setattr__(self, name, value):
        "Protect the 'array' attribute"
        if name == "array":
            if name in self.__dict__:
                raise AttributeError, "Cannot change Epetra.Vector array attribute"
        UserArray.__setattr__(self, name, value)

%}

#ifndef HAVE_MPI
%pythoncode %{
def PyComm():
  return SerialComm();
%}
#else
%pythoncode %{
def PyComm():
  return MpiComm(CommWorld());
%}
#endif

%include "Epetra_PyOperator.h"
%include "Epetra_PyRowMatrix.h"

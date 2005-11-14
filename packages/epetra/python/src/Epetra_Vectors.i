%{
// Epetra includes
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

// Local includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
%}

// Ignore directives
%ignore Epetra_MultiVector::operator=(const Epetra_MultiVector &);
%ignore Epetra_MultiVector::operator[](int);         // See %extend MultiVector
%ignore Epetra_MultiVector::operator[](int) const;   //       __getitem__()
%ignore Epetra_MultiVector::operator()(int) const;

// Rename directives
%rename(NumPyMultiVector) Epetra_NumPyMultiVector;
%rename(NumPyVector     ) Epetra_NumPyVector;

// Import directives for Epetra
%include "Epetra_MultiVector.h"
%include "Epetra_Vector.h"

// Local interface includes
%include "Epetra_NumPyMultiVector.h"
%include "Epetra_NumPyVector.h"

// Extend directives
// %extend Epetra_MultiVector {
//   double * & __getitem__(int i) {
//     return self->operator[](i);
//   }

//   PyObject * Norm1() {
//     int n = self->NumVectors();
//     double result[n];
//     int numVectors[1] = {n};
//     int status        = self->Norm1(result);
//     PyObject * output = Py_BuildValue("(iO)", status, 
// 				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
// 							      (char *)result));
//     return output;
//   }

//   PyObject * __setitem__(PyObject * args, double val) {
//     int i, j;
//     if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
//       PyErr_SetString(PyExc_IndexError, "Invalid index");
//       return NULL;
//     }
//     (*self)[i][j] = val;
//     Py_INCREF(Py_None);
//     return Py_None;
//   }

//   void Set(const int vector, const int element, const double value) {
//     (*self)[vector][element] = value;
//   }

//   PyObject * __getitem__(PyObject * args) {
//     int i, j;
//     if (!PyArg_ParseTuple(args, "ii", &i, &j)) {
//       PyErr_SetString(PyExc_IndexError, "Invalid index");
//       return NULL;
//     }
//     double val = (*self)[i][j];
//     return PyFloat_FromDouble(val);
//   }

//   double Get(const int vector, const int element) {
//     return((*self)[vector][element]);
//   }

//   PyObject * Norm2() {
//     int n = self->NumVectors();
//     double result[n];
//     int numVectors[1] = {n};
//     int status        = self->Norm2(result);
//     PyObject * output = Py_BuildValue("(iO)", status,
// 				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
// 							      (char *)result));
//     return output;
//   }

//   PyObject * NormInf() {
//     int n = self->NumVectors();
//     double result[n];
//     int numVectors[1] = {n};
//     int status        = self->NormInf(result);
//     PyObject * output = Py_BuildValue("(iO)", status,
// 				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
// 							      (char *)result));
//     return output;
//   }

//    PyObject * Dot(const Epetra_MultiVector& A){
//     int n = self->NumVectors();
//     double result[n];
//     int numVectors[1] = {n};
//     int status        = self->Dot(A, result);
//     PyObject * output = Py_BuildValue("(iO)", status,
// 				      PyArray_FromDimsAndData(1,numVectors, PyArray_DOUBLE,
// 							      (char *)result));
//     return output;
//    }
// }

// %extend Epetra_Vector{

//   PyObject * Norm1() {
//     double result[1];
//     int status        = self->Norm1(result);
//     PyObject * output = Py_BuildValue("(id)", status, result[0]);
//     return output;
//   }

//  PyObject * Norm2() {
//     double result[1];
//     int status        = self->Norm2(result);
//     PyObject * output = Py_BuildValue("(id)", status, result[0]);
//     return output;
//   }

//   PyObject * NormInf() {
//     double result[1];
//     int status        = self->NormInf(result);
//     PyObject * output = Py_BuildValue("(id)", status, result[0]);    
//     return output;
//   }
//    PyObject * Dot(const Epetra_MultiVector& A){
//     double result[1];
//     int status        = self->Dot(A, result);
//     PyObject * output = Py_BuildValue("(id)", status, result[0]);
//     return output;
//    }
// }

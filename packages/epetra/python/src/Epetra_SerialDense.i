// -*- C++ -*-

%{
// Epetra includes
#include "Epetra_SerialDenseOperator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
%}

// Ignore directives
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

// Rename directives
%rename(SerialDenseOperator ) Epetra_SerialDenseOperator;
%rename(SerialDenseMatrix   ) Epetra_SerialDenseMatrix;
%rename(SerialSymDenseMatrix) Epetra_SerialSymDenseMatrix;
%rename(SerialDenseVector   ) Epetra_SerialDenseVector;
%rename(SerialDenseSolver   ) Epetra_SerialDenseSolver;
%rename(IntSerialDenseMatrix) Epetra_IntSerialDenseMatrix;
%rename(IntSerialDenseVector) Epetra_IntSerialDenseVector;

// Include directives
%include "Epetra_SerialDenseOperator.h"
%include "Epetra_SerialDenseMatrix.h"
%include "Epetra_SerialSymDenseMatrix.h"
%include "Epetra_SerialDenseVector.h"
%include "Epetra_SerialDenseSolver.h"
%include "Epetra_IntSerialDenseMatrix.h"
%include "Epetra_IntSerialDenseVector.h"

// Extend directives
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

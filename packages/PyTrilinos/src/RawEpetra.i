// -*- c++ -*-

%module RawEpetra

%{
// System includes
#include <sstream>
#include <vector>

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
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
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseOperator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialSymDenseMatrix.h"
// Local includes
#include "Callback.h"
#include "CallbackVectorLoadTest.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_VectorHelper.h"
#include "NumPyArray.h"
#include "NumPyWrapper.h"
#include "PyObjectHolder.h"
%}

// Ignore directives
%ignore operator<<(ostream &, const Epetra_Object &);// From python, use print
%ignore Epetra_Object::Print(ostream &) const;       // faciltated by __str__
%ignore Epetra_SerialComm::operator=(const Epetra_SerialComm &);
%ignore Epetra_CompObject::operator=(const Epetra_CompObject &);
%ignore Epetra_CompObject::UpdateFlops(int) const;   // Use long int version
%ignore Epetra_CompObject::UpdateFlops(float) const; // Use double version
%ignore Epetra_BlockMap::operator=(const Epetra_BlockMap &);
%ignore Epetra_Map::operator=(const Epetra_Map &);
%ignore Epetra_LocalMap::operator=(const Epetra_LocalMap &);
%ignore Epetra_MultiVector::operator=(const Epetra_MultiVector &);
%ignore Epetra_MultiVector::operator[](int);         // See %extend MultiVector
%ignore Epetra_MultiVector::operator[](int) const;   //       __getitem__()
%ignore Epetra_IntVector::operator=(const Epetra_IntVector &);
%ignore Epetra_IntVector::operator[](int);           // See %extend IntVector
%ignore Epetra_IntVector::operator[](int) const;     //       __getitem__()
%ignore Epetra_CrsGraph::operator=(const Epetra_CrsGraph &);
%ignore Epetra_CrsGraph::operator[](int);            // See %extend CrsGraph
%ignore Epetra_CrsGraph::operator[](int) const;      //       __getitem__()
%ignore Epetra_MapColoring::operator[](int);         // See %extend Mapcoloring
%ignore Epetra_MapColoring::operator[](int) const;   //       __getitem__()
%ignore Epetra_CrsMatrix::operator=(const Epetra_CrsMatrix &);
%ignore Epetra_CrsMatrix::operator[](int);           // See %extend CrsMatrix
%ignore Epetra_CrsMatrix::operator[](int) const;     //       __getitem__()
// Epetra_VbrMatrix member function
// Solve(bool,bool,bool,Epetra_Vector const&,Epetra_Vector&) const doesn't
// appear to be implemented.  Apparently overridden by
// Solve(bool,bool,bool,Epetra_MultiVector const&,Epetra_MultiVector&) const
%ignore Epetra_VbrMatrix::operator=(const Epetra_VbrMatrix &);
%ignore Epetra_VbrMatrix::Solve(bool, bool, bool,
				Epetra_Vector const&, Epetra_Vector&) const;
%ignore Epetra_SerialDenseMatrix::operator=(const Epetra_SerialDenseMatrix &);
//%ignore Epetra_SerialDenseMatrix::operator+=(const Epetra_SerialDenseMatrix &);
%ignore Epetra_SerialDenseMatrix::operator[](int);
%ignore Epetra_SerialDenseMatrix::operator[](int) const;
%ignore Epetra_SerialDenseVector::operator=(const Epetra_SerialDenseVector &);
%ignore Epetra_SerialDenseVector::operator[](int);
%ignore Epetra_SerialDenseVector::operator[](int) const;
%ignore Epetra_SerialDenseVector::operator()(int);
%ignore Epetra_SerialDenseVector::operator()(int) const;
%ignore PyObjectHolder::operator PyObject * ();
%ignore NumPyArrayBase::print(std::ostream &) const;       // faciltated by __str__
%ignore NumPyArray::print(std::ostream &) const;           // faciltated by __str__
%ignore NumPyArrayContiguous::print(std::ostream &) const; // faciltated by __str__

// Rename directives
%rename(Object              ) Epetra_Object;
%rename(Comm                ) Epetra_Comm;
%rename(SerialComm          ) Epetra_SerialComm;
%rename(BlockMap            ) Epetra_BlockMap;
%rename(Map                 ) Epetra_Map;
%rename(LocalMap            ) Epetra_LocalMap;
%rename(SrcDistObject       ) Epetra_SrcDistObject;
%rename(DistObject          ) Epetra_DistObject;
%rename(CompObject          ) Epetra_CompObject;
%rename(BLAS                ) Epetra_BLAS;
%rename(LAPACK              ) Epetra_LAPACK;
%rename(MultiVector         ) Epetra_MultiVector;
//%rename(Vector              ) Epetra_Vector;
%rename(IntVector           ) Epetra_IntVector;
%rename(CrsGraph            ) Epetra_CrsGraph;
%rename(MapColoring         ) Epetra_MapColoring;
%rename(Operator            ) Epetra_Operator;
%rename(RowMatrix           ) Epetra_RowMatrix;
%rename(CrsMatrix           ) Epetra_CrsMatrix;
%rename(VbrMatrix           ) Epetra_VbrMatrix;
%rename(SerialDenseSolver   ) Epetra_SerialDenseSolver;
%rename(SerialDenseOperator ) Epetra_SerialDenseOperator;
%rename(SerialDenseMatrix   ) Epetra_SerialDenseMatrix;
%rename(SerialDenseVector   ) Epetra_SerialDenseVector;
%rename(SerialSymDenseMatrix) Epetra_SerialSymDenseMatrix;
%rename(NumPyVector         ) Epetra_NumPyVector;

// Typemap directives
%typemap(in) (int * Indices)
{
  PyArrayObject * array = NumPyWrapper::contiguous_typed_array($input, PyArray_INT,
							       0, 0);
  assert (array && "Function should have already checked this.");
  $1 = (int *) array->data;
}

%typemap(in,numinputs=0) (int & NumIndices, int *& Indices)
{
  int   temp1;
  int * temp2;
  $1 = &temp1;
  $2 = &temp2;
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
%include "Epetra_Object.h"
%include "Epetra_Comm.h"
%include "Epetra_SerialComm.h"
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
%include "Epetra_VbrMatrix.h"
%include "Epetra_DataAccess.h"
%include "Epetra_SerialDenseSolver.h"
%include "Epetra_SerialDenseOperator.h"
%include "Epetra_SerialDenseMatrix.h"
%include "Epetra_SerialDenseVector.h"
%include "Epetra_SerialSymDenseMatrix.h"

// Local interface includes
%include "PyObjectHolder.h"
%include "NumPyArray.h"
%include "Epetra_NumPyVector.h"

// Extensions
%extend Epetra_Object {
  using namespace std;
  string __str__() {
    stringstream os;
    self->Print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

%extend Epetra_MultiVector {
  double * & __getitem__(int i) {
    return self->operator[](i);
  }
}

%extend Epetra_Vector {
  void load(PyObject * p_pyObject) {
    Epetra_VectorHelper::loadViaCopy(self, p_pyObject);
  }

  void unload(PyObject * p_pyObject) {
    Epetra_VectorHelper::unloadViaCopy(self, p_pyObject);
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
  double * & __getitem__(int i) {
    return self->operator[](i);
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

//   Epetra_SerialDenseMatrix & __iadd__(const Epetra_SerialDenseMatrix & other) {
//     return self->operator+=(other);
//   }
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

%extend NumPyArrayBase {
  using namespace std;
  string __str__() {
    stringstream os;
    self->print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}

%extend Epetra_Vector {

  // HACK: Constructor that takes a wrapped PyObject pointer to get around
  //       a SWIG bug that has been fixed in SWIG 1.3.19.
  // NOTE: This must come before the constructor that takes a PyObject
  Epetra_Vector(Epetra_BlockMap & map, PyObjectHolder & holder)  {
    return Epetra_VectorHelper::new_Epetra_Vector(map, holder);
  }

  Epetra_Vector(Epetra_BlockMap & map, PyObject * p_pyObject)  {
    return Epetra_VectorHelper::new_Epetra_Vector(map, p_pyObject);
  }

  // HACK: See note for Epetra_Vector constructor.
  // NOTE: This must come before the load method that takes a PyObject
  void loadViaCopy (PyObjectHolder & holder) {
    Epetra_VectorHelper::loadViaCopy(self, holder);
  }

  void loadViaCopy (PyObject * p_pyObject) {
    Epetra_VectorHelper::loadViaCopy(self, p_pyObject);
  }

  // HACK: See note for Epetra_Vector constructor.
  // NOTE: This must come before the unload method that takes a PyObject
  void unloadViaCopy (PyObjectHolder & holder) {
    Epetra_VectorHelper::unloadViaCopy(self, holder);
  }

  void unloadViaCopy (PyObject * p_pyObject) {
    Epetra_VectorHelper::unloadViaCopy(self, p_pyObject);
  }
}

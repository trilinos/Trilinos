// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef EPETRA_NUMPYINTSERIALDENSEMATRIX_HPP
#define EPETRA_NUMPYINTSERIALDENSEMATRIX_HPP

#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

#include "PyTrilinos_PythonException.hpp"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Epetra_IntSerialDenseMatrix.h"

namespace PyTrilinos
{

class Epetra_NumPyIntSerialDenseMatrix : public Epetra_IntSerialDenseMatrix
{
public:

  // Constructors
  Epetra_NumPyIntSerialDenseMatrix();
  Epetra_NumPyIntSerialDenseMatrix(int numRows, int numCols);
  Epetra_NumPyIntSerialDenseMatrix(PyObject * pyObject);
  Epetra_NumPyIntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix & src);

  // Destructor
  ~Epetra_NumPyIntSerialDenseMatrix();

  // Overridden Epetra_IntSerialDenseMatrix methods.  These are
  // overriden for one of two reasons: (1) to provide a more
  // python-like signature, or (2) to maintain synchronization between
  // the Epetra_SerialDenseMatrix and the numpy array.
  int        operator() (int rowIndex, int colIndex);
  int        Shape(  int numRows, int numCols);
  int        Reshape(int numRows, int numCols);
  PyObject * A();

  // Static cleanup function, to be called when python exceptions are
  // encountered
  static void cleanup();

private:

  // This private static pointer is for use with constructors only.
  // Outside of a constructor, it should be NULL.  It should only be
  // set by the static helper functions below, and when a constructor
  // sets this pointer, it should reset it to NULL when done.
  static PyArrayObject * tmp_array;

  // Static helper functions.  These are intended to be called from
  // the constructors, specifically to compute arguments in the
  // Epetra_SerialDenseMatrix constructors called in the constructor
  // initialization lists.  They all assume that if tmp_array is
  // already set, it has been set by the same PyObject.
  static int * getArray(  PyObject *);
  static int   getNumCols(PyObject *);
  static int   getNumRows(PyObject *);

  // Private method.  This method is typically called after an
  // Epetra_IntSerialDenseMatrix constructor has been called to
  // synchronize the internal PyArrayObject.
  void setArray(bool copy=false);

  // Private data
  PyArrayObject * array;

};

}  // Namespace PyTrilinos

#endif

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

#ifndef EPETRA_NUMPYSERIALDENSEVECTOR_H
#define EPETRA_NUMPYSERIALDENSEVECTOR_H

#define NO_IMPORT_ARRAY
#include "numpy_include.h"

#include "PythonException.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Epetra_SerialDenseVector.h"

namespace PyTrilinos
{

class Epetra_NumPySerialDenseVector : public Epetra_SerialDenseVector
{
public:

  // Constructors
  Epetra_NumPySerialDenseVector();
  Epetra_NumPySerialDenseVector(int length);
  Epetra_NumPySerialDenseVector(PyObject * pyObject);
  Epetra_NumPySerialDenseVector(const Epetra_SerialDenseVector & src);

  // Destructor
  ~Epetra_NumPySerialDenseVector();

  // Overridden Epetra_SerialDenseVector methods.  These are overriden
  // for one of two reasons: (1) to provide a more python-like
  // signature, or (2) to maintain synchronization between the
  // Epetra_SerialDenseVector and the numpy array.
  int        Size(int length);
  int        Resize(int length);
  PyObject * Values() const;

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
  // Epetra_SerialDenseVector constructors called in the constructor
  // initialization lists.  They all assume that if tmp_array is
  // already set, it has been set by the same PyObject.
  static double * getArray(     PyObject *);
  static int      getVectorSize(PyObject *);

  // Private method.  This method is typically called after an
  // Epetra_SerialDenseVector constructor has been called to
  // synchronize the internal PyArrayObject.
  void setArray();

  // Private data
  PyArrayObject * array;

};

}  // Namespace PyTrilinos

#endif

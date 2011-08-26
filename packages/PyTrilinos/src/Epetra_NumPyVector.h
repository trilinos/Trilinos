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

#ifndef EPETRA_NUMPYVECTOR_H
#define EPETRA_NUMPYVECTOR_H

#define NO_IMPORT_ARRAY
#include "numpy_include.h"

#include "PyTrilinos_PythonException.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

namespace PyTrilinos
{

class Epetra_NumPyVector : public Epetra_Vector
{
public:

  // Constructors
  Epetra_NumPyVector(const Epetra_BlockMap & blockMap, bool zeroOut=true);
  Epetra_NumPyVector(const Epetra_Vector & source);
  Epetra_NumPyVector(const Epetra_BlockMap & blockMap, PyObject * pyObject);
  Epetra_NumPyVector(Epetra_DataAccess CV, const Epetra_Vector & source);
  Epetra_NumPyVector(Epetra_DataAccess CV, const Epetra_MultiVector & source, int index);
  Epetra_NumPyVector(PyObject * pyObject);

  // Destructor
  ~Epetra_NumPyVector();

  // Overridden Epetra_Vector methods with more python-like signatures
  PyObject * ExtractCopy() const;
  PyObject * ExtractView() const;
  double     Dot(const Epetra_Vector & A) const;
  double     Norm1() const;
  double     Norm2() const;
  double     NormInf() const;
  double     NormWeighted(const Epetra_Vector & weights) const;
  double     MinValue() const;
  double     MaxValue() const;
  double     MeanValue() const;
  int        ReplaceGlobalValues(PyObject * values, PyObject * indices);
  int        ReplaceGlobalValues(int blockOffset, PyObject * values, PyObject * indices);
  int        ReplaceMyValues(PyObject * values, PyObject * indices);
  int        ReplaceMyValues(int blockOffset, PyObject * values, PyObject * indices);
  int        SumIntoGlobalValues(PyObject * values, PyObject * indices);
  int        SumIntoGlobalValues(int blockOffset, PyObject * values, PyObject * indices);
  int        SumIntoMyValues(PyObject * values, PyObject * indices);
  int        SumIntoMyValues(int blockOffset, PyObject * values, PyObject * indices);

  // Static cleanup function, to be called when python exceptions are
  // encountered
  static void cleanup();

private:

  // Private method thus not callable
  Epetra_NumPyVector();

  // This static private constant is the default communicator for any
  // Epetra_NumPyVector constructed without specifying an
  // Epetra_BlockMap
  static const Epetra_SerialComm defaultComm;

  // These private static pointers are for use with constructors only.
  // Outside of a constructor, they should be NULL.  They should only
  // be set by the static helper functions below, and when a
  // constructor sets any of these pointers, it should reset it to
  // NULL when done.
  static PyArrayObject * tmp_array;
  static Epetra_Map    * tmp_map;

  // Static helper functions.  These are intended to be called from
  // the constructors, specifically to compute arguments in the
  // Epetra_Vector constructors called in the constructor
  // initialization lists.  They all assume that if tmp_array or
  // tmp_map is already set, they have been set by the same PyObject.
  static double     * getArray(     PyObject *);
  static Epetra_Map & getMap(       PyObject *);
  static int          getVectorSize(PyObject *);
  static double     * getArray(const Epetra_BlockMap &, PyObject *);

  // Private data
  Epetra_BlockMap * map;
  PyArrayObject   * array;

};

}  // Namespace PyTrilinos

#endif

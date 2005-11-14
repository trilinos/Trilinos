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

#ifndef EPETRA_NUMPYMULTIVECTOR_H
#define EPETRA_NUMPYMULTIVECTOR_H

#define NO_IMPORT_ARRAY
#include "numeric_include.h"

#include "Epetra_SerialComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

class Epetra_NumPyMultiVector : public Epetra_MultiVector {

public:

  Epetra_NumPyMultiVector(const Epetra_BlockMap & blockMap, int numVectors, bool zeroOut=true);
  Epetra_NumPyMultiVector(const Epetra_NumPyMultiVector & source);
  Epetra_NumPyMultiVector(const Epetra_BlockMap & blockMap, PyObject * pyObject);
  Epetra_NumPyMultiVector(Epetra_DataAccess CV, const Epetra_NumPyMultiVector & source,
			  PyObject * range);
  Epetra_NumPyMultiVector(PyObject * pyObject);

  ~Epetra_NumPyMultiVector();

  PyObject * getArray() const;
  PyObject * Norm1() const;
  PyObject * Norm2() const;
  PyObject * NormInf() const;
  PyObject * Dot(const Epetra_MultiVector & A) const;

private:

  // Private methods thus not callable
  Epetra_NumPyMultiVector();

  // Static helper functions
  static Epetra_Map & getEpetraMap(PyObject *);
  static int        * getRange(    PyObject *);
  static double     * getArrayFromObject(PyObject *);
  static double     * getArrayFromMapAndObject(const Epetra_BlockMap &, PyObject *);

  // Static private data
  static const Epetra_SerialComm defaultComm;
  static       PyArrayObject *   tmp_array;
  static       Epetra_Map    *   tmp_map;
  static       PyArrayObject *   tmp_range;
  static       int               tmp_range_len;

  // Private data
  Epetra_BlockMap * map;
  PyArrayObject   * array;

};

#endif

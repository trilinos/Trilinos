// @HEADER
// ***********************************************************************
//
//                PyTrilinos: Rapid Prototyping Package
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

#ifndef NUMPY_INCLUDE_H
#define NUMPY_INCLUDE_H

// This include file takes care of three of the four things necessary
// when including the numpy header file arrayobject.h.  First, the
// Python.h header file is included.  Second, the
// PY_ARRAY_UNIQUE_SYMBOL is defined.  Third, the
// numpy/arrayobject.h header file is included.

// The user is responsible for defining the macro NO_IMPORT_ARRAY in
// those source files that do not call the numpy routine
// import_array().

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL PyTrilinos
#include <numpy/arrayobject.h>

#endif // NUMPY_INCLUDE_H

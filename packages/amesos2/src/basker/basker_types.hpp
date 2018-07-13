// @HEADER
// ***********************************************************************
//
//                   Basker: A Direct Linear Solver package
//                    Copyright 2011 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Mike A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#ifndef BASKER_TYPES_HPP
#define BASKER_TYPES_HPP

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"

#define BASKERASSERT(a) mxAssert(a,"")
#define BASKERREALLOC(ptr,size) mxRealloc(ptr, size)
#define BASKERCALLOC(num, size) mxCalloc(num, size)
#define BASKERFREE(ptr)

#else
#include <cassert>
#include <cstdlib>

#define BASKERASSERT(a) assert(a)
#define BASKERREALLOC(ptr, size) realloc(ptr, size)
#define BASKERCALLOC(num, size)  calloc(num,size)
#define BASKERFREE(ptr)          free(ptr)

#endif



template <class Int, class Entry>
struct basker_matrix
{
  Int nrow, ncol, nnz;
  Int *col_ptr, *row_idx;
  Entry *val;

};

#endif

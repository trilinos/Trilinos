// @HEADER
// *****************************************************************************
//                   Basker: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the Basker contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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

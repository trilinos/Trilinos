#ifndef BASKER_TYPES_HPP
#define BASKER_TYPES_HPP

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"

#define ASSERT(a) mxAssert(a,"")
#define REALLOC(ptr,size) mxRealloc(ptr, size)
#define CALLOC(num, size) mxCalloc(num, size)
#define FREE(ptr)         

#else
#include <assert.h>
#include <stdlib.h>

#define ASSERT(a) assert(a)
#define REALLOC(ptr, size) realloc(ptr, size)
#define CALLOC(num, size)  calloc(num,size)
#define FREE(ptr)          free(ptr) 

#endif



template <class Int, class Entry>
struct basker_matrix
{
  Int nrow, ncol, nnz;
  Int *col_ptr, *row_idx;
  Entry *val;
  
};

#endif

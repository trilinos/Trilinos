#ifndef BASKER_HPP
#define BASKER_HPP

#include "basker_types.hpp"
#include "basker_decl.hpp"
#include "basker_def.hpp"

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
#define CALLOC(num, size) calloc(num,size)
#define FREE(ptr)         free(ptr) 

#endif


#endif


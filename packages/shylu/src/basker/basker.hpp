#ifndef BASKER__HPP
#define BASKER_HPP

#include "basker_decl.hpp"
#include "basker_def.hpp"

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"

#define ASSERT(a) mxAssert(a,"")
#define REALLOC(ptr,size) mxRealloc(ptr,size)

#else
#include <assert.h>
#define ASSERT(a) assert(a)
#define REALLOC(ptr,size) realloc(ptr,size)


#endif






#endif


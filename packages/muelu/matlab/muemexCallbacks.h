//Muemex callbacks
//Brian Kelley

#ifndef MUEMEX_CALLBACKS_H
#define MUEMEX_CALLBACKS_H

#include "stdexcept"
#include "mex.h"

#if !defined(HAVE_MUELU_MATLAB)
#error "Muemex callbacks require MATLAB."
#else

enum DATA_TYPES
{
	REAL_MATRIX,
	COMPLEX_MATRIX,
	REAL_SCALAR,
	COMPLEX_SCALAR,
	STRING,
	REAL_MULTIVECTOR,
	COMPLEX_MULTIVECTOR
};

typedef void (*MuemexCallback) ();

#endif //HAVE_MUELU_MATLAB
#endif //MUEMEX_CALLBACKS_H guard

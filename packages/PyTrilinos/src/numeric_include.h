#ifndef numeric_include_h
#define numeric_include_h

// When including arrayobject.h to use Numeric objects,
// care must be taken w.r.t the PY_ARRAY_UNIQUE_SYMBOL
// symbol.  This include file addresses this and should be
// included by PyNOX code instead of including
// Numeric/arrayobject.h directly.

#include <Python.h>

//#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_NumPyWrapper
#include <Numeric/arrayobject.h>

#endif // numeric_include

#ifndef NUMERIC_INCLUDE_H
#define NUMERIC_INCLUDE_H

// This include file takes care of three of the four things necessary
// when including the Numeric header file arrayobject.h.  First, the
// Python.h header file is included.  Second, the
// PY_ARRAY_UNIQUE_SYMBOL is defined.  Third, the
// Numeric/arrayobject.h header file is included.

// The user is responsible for defining the macro NO_IMPORT_ARRAY in
// those source files that do not call the Numeric routine
// import_array().

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL PyTrilinos
#include <Numeric/arrayobject.h>

#endif // numeric_include

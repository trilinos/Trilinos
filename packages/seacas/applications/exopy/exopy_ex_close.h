#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"

static PyObject * exopy_ex_close(PyObject *self, PyObject *args) {
  /* ex_close(exoid) */

  int exoid, error;

  if ( !PyArg_ParseTuple(args, "i", &exoid) ) {
    return NULL;
  }

  error = ex_close(exoid);

  if ( error < 0 ) {
    PyErr_SetString(PyExc_RuntimeError, "error in exopy_ex_close()");
    return NULL;
  }

  return Py_BuildValue("i", error);
}

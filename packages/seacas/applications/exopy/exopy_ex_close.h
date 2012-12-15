#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"
#include "exopy_ref.h"

static PyObject * exopy_ex_close(PyObject *self, PyObject *args) {
  /* ex_close(exoid) */

  int exoid, error;

  if ( !PyArg_ParseTuple(args, "i", &exoid) ) {
    return NULL;
  }

  error = ex_close(exoid);

  // Do this so error is raised pointing to this function
  if ( PyErr_Occurred() ) { return NULL; }

  return Py_BuildValue("i", error);
}

#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"

static PyObject * exopy_ex_open_test(PyObject *self, PyObject *args) {
  /* (exoid, comp_ws, io_ws, version) = ex_open(path, mode) */

  int exoid, mode, comp_ws=0, io_ws=0; 
  float version=0.0;
  const char *path;

  if ( !PyArg_ParseTuple(args, "si", &path, &mode ) ) {
    return NULL;
  }

  exoid = ex_open(path, mode, &comp_ws, &io_ws, &version);

  return Py_BuildValue("(iiif)", exoid, comp_ws, io_ws, version);
}

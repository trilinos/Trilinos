#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"

static PyObject * exopy_ex_open(PyObject *self, PyObject *args) {
  /* (exoid,comp_ws,io_ws,version) = ex_open(path,mode,comp_ws_ref,io_ws_ref,version_ref) */

  int exoid, mode, comp_ws_in, io_ws_in;
  const char *path;

  if ( !PyArg_ParseTuple(args, "siii", &path, &mode, &comp_ws_in, &io_ws_in) ) {
    return NULL;
  }

  int comp_ws = comp_ws_in;
  int io_ws = io_ws_in;
  float version;

  exoid = ex_open(path, mode, &comp_ws, &io_ws, &version);

  if ( exoid < 0 ) {
    PyErr_SetString(PyExc_RuntimeError, "error in exopy_ex_open()");
    return NULL;
  }

  return Py_BuildValue("(iiif)", exoid, comp_ws, io_ws, version);
}

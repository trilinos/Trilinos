#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"
#include "exopy_ref.h"

static PyObject * exopy_ex_open(PyObject *self, PyObject *args) {
  /* ex_open(path,mode,comp_ws_ref,io_ws_ref,version_ref) */

  int exoid, mode;
  const char *path;

  ref *comp_ws_ref, *io_ws_ref, *version_ref;

  if ( !PyArg_ParseTuple(args, "siOOO", &path, &mode,
                         &comp_ws_ref, &io_ws_ref, &version_ref) ) {
    return NULL;
  }

  int comp_ws = Ref_AsInt(comp_ws_ref);
  int io_ws = Ref_AsInt(io_ws_ref);
  float version = Ref_AsFloat(version_ref);

  exoid = ex_open(path, mode, &comp_ws, &io_ws, &version);

  comp_ws_ref->value = PyInt_FromLong(comp_ws);
  io_ws_ref->value   = PyInt_FromLong(io_ws);
  version_ref->value = PyFloat_FromDouble(version);

  // Do this so error is raised pointing to this function
  if ( PyErr_Occurred() ) { return NULL; }

  return Py_BuildValue("i", exoid);
}

#include "Python.h" // must be included first
#include "exodusII.h"
#include "netcdf.h"
#include "exopy_ref.h"

#define is_ref(r) ((r)->ob_type == &refType)

static PyObject * exopy_ex_open(PyObject *self, PyObject *args) {
  /* ex_open(path,mode,comp_ws_ref,io_ws_ref,version_ref) */

  int exoId, mode, comp_ws, io_ws;
  const char *path;
  float version;
  ref *comp_ws_ref, *io_ws_ref, *version_ref;

  if ( !PyArg_ParseTuple(args, "siOOO", &path, &mode,
                         &comp_ws_ref, &io_ws_ref, &version_ref) ) {
    return NULL;
  }

  if ( !is_ref(comp_ws_ref) ) {
    PyErr_SetString(PyExc_TypeError,
      "argument 3 should be an exopy.ref object");
    return NULL;
  }

  if ( !is_ref(io_ws_ref) ) {
    PyErr_SetString(PyExc_TypeError,
      "argument 4 should be an exopy.ref object");
    return NULL;
  }

  if ( !is_ref(version_ref) ) {
    PyErr_SetString(PyExc_TypeError,
      "argument 5 should be an exopy.ref object");
    return NULL;
  }

// MGV: the following checks restrict the ref.value to be specific
//      PyObject types; this is really only neccessary for references
//      whose passed-in values are used by the Exodus API method, but
//      we are always doing the check here for the sake of code quality
  if ( PyInt_Check(comp_ws_ref->value) ) {
    comp_ws = PyInt_AsLong(comp_ws_ref->value);
  } else if ( comp_ws_ref->value == Py_None ) {
    comp_ws = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,
      "ref.value of argument 3 should be an int object");
    return NULL;
  }

  if ( PyInt_Check(io_ws_ref->value) ) {
    io_ws = PyInt_AsLong(io_ws_ref->value);
  } else if ( io_ws_ref->value == Py_None ) {
    io_ws = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,
      "ref.value of argument 4 should be an int object");
    return NULL;
  }

  if ( PyFloat_Check(version_ref->value) ) {
    version = PyFloat_AsDouble(version_ref->value);
  } else if ( version_ref->value == Py_None ) {
    version = 0.0;
  } else {
    PyErr_SetString(PyExc_TypeError,
      "ref.value of argument 5 should be a float object");
    return NULL;
  }

  exoId = ex_open(path, mode, &comp_ws, &io_ws, &version);

  comp_ws_ref->value = PyInt_FromLong(comp_ws);
  io_ws_ref->value   = PyInt_FromLong(io_ws);
  version_ref->value = PyFloat_FromDouble(version);

  return Py_BuildValue("i", exoId);
}

#include "Python.h"
#include "exodusII.h"
#include "netcdf.h"

static PyObject * py_ex_open(PyObject *self, PyObject *args) {
  /* ex_open(path,mode,comp_ws,io_ws,version) */

  const char *path;
  int mode;
  int comp_ws;
  int io_ws;
  float version;

  int exoId;

  if ( !PyArg_ParseTuple(args, "siiif", &path, &mode,
                         &comp_ws, &io_ws, &version) ) {
    return NULL;
  }

  exoId = ex_open(path, mode, &comp_ws, &io_ws, &version);

  //args[5] = PyFloat(version);
  //args[4] = Py_BuildValue("f",version);

  return Py_BuildValue("i", exoId);
}

static PyMethodDef ExoPyMethods[] = {
  {"ex_open", py_ex_open, METH_VARARGS, "Open an ExodusII file."},
  {     NULL,       NULL,            0,                     NULL}
};

PyMODINIT_FUNC initexopy(void) {
  (void) Py_InitModule("exopy", ExoPyMethods);
}


#include "Python.h" // must be included first
#include "exopy_ex_open.h"
#include "exopy_ex_close.h"
#include "exopy_ex_get_init.h"

static PyMethodDef exopy_methods[] = {
  {"ex_open",  exopy_ex_open,  METH_VARARGS, "open an ExodusII file"},
  {"ex_close", exopy_ex_close, METH_VARARGS, "close an ExodusII file"},
  {"ex_get_init", exopy_ex_get_init, METH_VARARGS, "read parameters in an ExodusII file"},
  {     NULL,  NULL,           0,            NULL}
};

PyMODINIT_FUNC initexopy(void) {
  PyObject *m;
  
  if ( PyType_Ready(&refType) < 0 ) return;

  m = Py_InitModule3("exopy", exopy_methods, "A Python module for the ExodusII API");

  if ( m == NULL ) return;

  Py_INCREF(&refType);
  PyModule_AddObject(m, "ref", (PyObject*)&refType);
}


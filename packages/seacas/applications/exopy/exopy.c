#include "Python.h" // must be included first
#include "structmember.h"
#include "exodusII.h"
#include "netcdf.h"

typedef struct {
  PyObject_HEAD
  PyObject *value; /* referenced value can be any PyObject */
} ref;

#define is_ref(r) ((r)->ob_type == &refType)

static void ref_dealloc(ref* self) {
  Py_XDECREF(self->value);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject * ref_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  ref *self;
  self = (ref *)type->tp_alloc(type, 0);

  if (self != NULL) {
    // initialize member variables
    
    // ref.value should return None until assigned (?? so use Py_None ??)
    self->value = Py_None;
    if (self->value == NULL) {
      Py_DECREF(self);
      return NULL;
    }

  }

  return (PyObject *)self;
}

static int ref_init(ref *self, PyObject *args, PyObject *kwds) {
  PyObject *value=NULL, *tmp;
  
  // allow ref to be initialized by any of the following calls:
  //   ref()                   // ref.value initializes as <None>
  //   ref(<object>)           // ref.value initializes as <object>
  //   ref("value" = <object>) // (ditto)
  static char *kwlist[] = {"value", NULL};
  if ( !PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &value) ) {
    return -1;
  }

  if (value) {
    tmp = self->value;
    Py_INCREF(value);
    self->value = value;
    Py_XDECREF(tmp);
  }

  return 0;
}

static PyMemberDef ref_members[] = {
  {"value", T_OBJECT_EX, offsetof(ref, value), 0, "referenced value"},
  {NULL}
};

static PyMethodDef ref_methods[] = {
  {NULL}
};

static PyTypeObject refType = {
    PyObject_HEAD_INIT(NULL)
    0,                                        /*ob_size*/
    "exopy.ref",                              /*tp_name*/
    sizeof(ref),                              /*tp_basicsize*/
    0,                                        /*tp_itemsize*/
    (destructor)ref_dealloc,                  /*tp_dealloc*/
    0,                                        /*tp_print*/
    0,                                        /*tp_getattr*/
    0,                                        /*tp_setattr*/
    0,                                        /*tp_compare*/
    0,                                        /*tp_repr*/
    0,                                        /*tp_as_number*/
    0,                                        /*tp_as_sequence*/
    0,                                        /*tp_as_mapping*/
    0,                                        /*tp_hash*/
    0,                                        /*tp_call*/
    0,                                        /*tp_str*/
    0,                                        /*tp_getattro*/
    0,                                        /*tp_setattro*/
    0,                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "pass by reference actor",                /*tp_doc*/
    0,		                              /*tp_traverse*/
    0,		                              /*tp_clear*/
    0,		                              /*tp_richcompare*/
    0,		                              /*tp_weaklistoffset*/
    0,		                              /*tp_iter*/
    0,		                              /*tp_iternext*/
    ref_methods,                              /*tp_methods*/
    ref_members,                              /*tp_members*/
    0,                                        /*tp_getset*/
    0,                                        /*tp_base*/
    0,                                        /*tp_dict*/
    0,                                        /*tp_descr_get*/
    0,                                        /*tp_descr_set*/
    0,                                        /*tp_dictoffset*/
    (initproc)ref_init,                       /*tp_init*/
    0,                                        /*tp_alloc*/
    ref_new,                                  /*tp_new*/
};

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
      "ref.value of argument 5 should be an int object");
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

static PyMethodDef exopy_methods[] = {
  {"ex_open", exopy_ex_open, METH_VARARGS, "open an ExodusII file"},
  {     NULL,          NULL,            0,                     NULL}
};

PyMODINIT_FUNC initexopy(void) {
  PyObject *m;
  
  if ( PyType_Ready(&refType) < 0 ) return;

  m = Py_InitModule3("exopy", exopy_methods,
                     "A Python module for the ExodusII API");

  if ( m == NULL ) return;

  Py_INCREF(&refType);
  PyModule_AddObject(m, "ref", (PyObject*)&refType);
}


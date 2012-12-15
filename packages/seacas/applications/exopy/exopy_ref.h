#ifndef EXOPY_REF_H
#define EXOPY_REF_H

#include "Python.h" // must be included first
#include "structmember.h"

#define is_ref(r) ((r)->ob_type == &refType)

typedef struct {
  PyObject_HEAD
  PyObject *value; /* referenced value can be any PyObject */
} ref;

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

int Ref_AsInt( ref * r ) { 
  if ( !is_ref(r) ) {
    PyErr_SetString(PyExc_TypeError, "argument should be an exopy.ref object");
    return -1;
  }
  if ( PyInt_Check(r->value) ) {
    return PyInt_AsLong(r->value);
  } else if ( r->value == Py_None ) {
    return 0;
  } else {
    PyErr_SetString(PyExc_TypeError, "ref.value of argument should be an int object");
    return -1;
  }
}

float Ref_AsFloat( ref * r ) { 
  if ( !is_ref(r) ) {
    PyErr_SetString(PyExc_TypeError, "argument should be an exopy.ref object");
    return -1.0;
  }
  if ( PyFloat_Check(r->value) ) {
    return PyFloat_AsDouble(r->value);
  } else if ( r->value == Py_None ) {
    return 0.0;
  } else {
    PyErr_SetString(PyExc_TypeError, "ref.value of argument should be a float object");
    return -1.0;
  }
}

double Ref_AsDouble( ref * r ) { 
  if ( !is_ref(r) ) {
    PyErr_SetString(PyExc_TypeError, "argument should be an exopy.ref object");
    return -1.0;
  }
  if ( PyFloat_Check(r->value) ) {
    return PyFloat_AsDouble(r->value);
  } else if ( r->value == Py_None ) {
    return 0.0;
  } else {
    PyErr_SetString(PyExc_TypeError, "ref.value of argument should be a float object");
    return -1.0;
  }
}

#endif

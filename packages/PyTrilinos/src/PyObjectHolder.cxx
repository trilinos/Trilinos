#include <iostream>
#include "PyObjectHolder.h"

PyObjectHolder::PyObjectHolder(PyObject * p_pyObject) :
  mp_pyObject(p_pyObject)
{
}

PyObjectHolder::~PyObjectHolder() 
{
}

PyObjectHolder::operator PyObject * ()
{
  return mp_pyObject;
}

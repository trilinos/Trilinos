#ifndef PYOBJECTHOLDER_H
#define PYOBJECTHOLDER_H

#include "Python.h"

class PyObjectHolder
{
public:
  PyObjectHolder(PyObject * p_pyObject);
  ~PyObjectHolder();

  operator PyObject * ();

private:
  // Private and not implemented
  PyObjectHolder();
  PyObjectHolder(const PyObjectHolder &);
  PyObjectHolder & operator=(const PyObjectHolder &);

private:
  PyObject * mp_pyObject;
};

#endif // PYOBJECTHOLDER_H

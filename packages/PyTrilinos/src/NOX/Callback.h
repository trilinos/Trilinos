#ifndef Callback_h
#define Callback_h

#include "Python.h"

class Callback
{
public:
  Callback();
  ~Callback();

  PyObject       * setFunction(PyObject * args);
  PyObject       * getFunction();
  const PyObject * getFunction() const;

private:
  // Private and not implemented so never can be called
  Callback(const Callback & a_ref);
  const Callback & operator = (const Callback & a_ref); 

  PyObject * mp_callback;
};

#endif //Callback_h

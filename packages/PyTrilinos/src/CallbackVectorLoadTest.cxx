#include "CallbackVectorLoadTest.h"

#ifndef PYEPETRA_VECTORHELPER_H
#  include "Epetra_VectorHelper.h"
#endif

#include "Epetra_Vector.h"

#include <iostream>

CallbackVectorLoadTest::CallbackVectorLoadTest(Epetra_Vector * p_epetraVector1,
					       Epetra_Vector * p_epetraVector2,
					       PyObject * p_callbackFunction)
: mp_epetraVector1(p_epetraVector1),
  mp_epetraVector2(p_epetraVector2)
{
  assert (0 != p_epetraVector1    && "NULL pointer passed in argument list");
  assert (0 != p_epetraVector2    && "NULL pointer passed in argument list");
  assert (0 != p_callbackFunction && "NULL pointer passed in argument list");

  // Note, we assume callback function takes two arguments that are PyObjects of some-sort
  m_callback.setFunction(p_callbackFunction);
}

CallbackVectorLoadTest::~CallbackVectorLoadTest()
{
  // Reset pointers.
  // NOTE: We do not own the objects pointed to.
  mp_epetraVector1 = 0;
  mp_epetraVector2 = 0;
}

void CallbackVectorLoadTest::execute(int numTimes)
{
  // Copy vector 2 to vector 1 and call the callback function

  assert (numTimes >= 0         && "Function expects to be asked to perform callback at least once.");
  assert (0 != mp_epetraVector1 && "Function expects non-NULL Epetra_Vector pointer."               );
  assert (0 != mp_epetraVector2 && "Function expects non-NULL Epetra_Vector pointer."               );

  PyObject * p_arglist = 0;
  PyObject * p_result  = 0;

  
  register int counter = 0;
  for (counter = 0; counter < numTimes; ++counter)
  {
    std::cout << "Copying vector 2 to vector 1" << std::endl;
    *mp_epetraVector1 = *mp_epetraVector2;

    
    // Here we assume callback function takes two arguments that are PyObjects of some-sort
    p_arglist = 0;
    p_result  = 0;
    p_arglist = Py_BuildValue("(OO)",mp_epetraVector1, mp_epetraVector2); // Build argument list
    assert (0 != p_arglist && "Argument list for callback function not built correctly");

    std::cout << "Calling callback function..." << std::endl;
    p_result  = PyEval_CallObject(m_callback.getFunction(), p_arglist);
    Py_DECREF(p_arglist);  // All done with argument list
    
    // Handle bad function call
    // From Extending and Embedding the Python Interpreter
    // http://www.python.org/doc/current/ext/callingPython.html
    // ... it is important to check that the return
    // value isn't NULL. If it is, the Python function terminated by raising an
    // exception. If the C code that called PyEval_CallObject() is called from
    // Python, it should now return an error indication to its Python caller, so
    // the interpreter can print a stack trace, or the calling Python code can
    // handle the exception. If this is not possible or desirable, the exception
    // should be cleared by calling PyErr_Clear()

    PyObject * p_pyErr_Occurred = PyErr_Occurred();
    if (p_pyErr_Occurred ) {
      std::cout << "There was an exeption: " << PyString_AsString(PyObject_Str(p_pyErr_Occurred)) << std::endl;
      assert (0 && "An exception occurred.  Sorry, no stack trace yet...");
    }
    
    if (0 == p_result)
    {
      PyErr_SetString(PyExc_ValueError,
                      "Bad result returned when attempting to execute callback function.");
      assert (0 && "Bad result returned when attempting to execute callback function.  Sorry, no stack trace yet...");
    }
    
    Py_DECREF(p_result); // All done with returned result object

    std::cout << "Back from callback." << std::endl;
  }
}

void CallbackVectorLoadTest::reset1 (Epetra_Vector * p_epetraVector1)
{
  assert (0 != p_epetraVector1 && "NULL pointer passed in argument list");
  mp_epetraVector1 = p_epetraVector1;
}

void CallbackVectorLoadTest::reset2 (Epetra_Vector * p_epetraVector2)
{
  assert (0 != p_epetraVector2 && "NULL pointer passed in argument list");
  mp_epetraVector2 = p_epetraVector2;
}

Epetra_Vector * CallbackVectorLoadTest::getEpetraVector1()
{
  assert (0 != mp_epetraVector1 && "Found NULL pointer.  Object in bad state.");
  return mp_epetraVector1;
}

Epetra_Vector * CallbackVectorLoadTest::getEpetraVector2()
{
  assert (0 != mp_epetraVector2 && "Found NULL pointer.  Object in bad state.");
  return mp_epetraVector2;
}

void CallbackVectorLoadTest::load1ViaCopy (PyObject * p_pyObject)
{
  assert (0 != mp_epetraVector1 && "Found NULL pointer.  Object in bad state.");
  Epetra_VectorHelper::loadViaCopy (mp_epetraVector1, p_pyObject);
}

void CallbackVectorLoadTest::load2ViaCopy (PyObject * p_pyObject)
{
  assert (0 != mp_epetraVector2 && "Found NULL pointer.  Object in bad state.");
  Epetra_VectorHelper::loadViaCopy (mp_epetraVector1, p_pyObject);
}

void CallbackVectorLoadTest::unload1ViaCopy (PyObject * p_pyObject)
{
  assert (0 != mp_epetraVector1 && "Found NULL pointer.  Object in bad state.");
  Epetra_VectorHelper::unloadViaCopy (mp_epetraVector1, p_pyObject);
}

void CallbackVectorLoadTest::unload2ViaCopy (PyObject * p_pyObject)
{
  assert (0 != mp_epetraVector2 && "Found NULL pointer.  Object in bad state.");
  Epetra_VectorHelper::unloadViaCopy (mp_epetraVector2, p_pyObject);
}



#ifndef CallbackVectorLoadTest_h
#define CallbackVectorLoadTest_h

#include "Python.h"
#include "Callback.h"

class Epetra_Vector;

class CallbackVectorLoadTest
{
public:
  CallbackVectorLoadTest(Epetra_Vector * p_epetraVector1,
			 Epetra_Vector * p_epetraVector2,
			 PyObject      * p_callbackFunction);
  ~CallbackVectorLoadTest();

  void execute        (int numTimes);

  void reset1         (Epetra_Vector * p_epetraVector1);
  void reset2         (Epetra_Vector * p_epetraVector2);

  void load1ViaCopy   (PyObject * p_pyObject);
  void load2ViaCopy   (PyObject * p_pyObject);

  void unload1ViaCopy (PyObject * p_pyObject);
  void unload2ViaCopy (PyObject * p_pyObject);

  Epetra_Vector * getEpetraVector1();
  Epetra_Vector * getEpetraVector2();
  

private:
  // Private and not implemented so never can be called
  CallbackVectorLoadTest();
  CallbackVectorLoadTest(const CallbackVectorLoadTest & a_ref);
  const CallbackVectorLoadTest & operator = (const CallbackVectorLoadTest & a_ref);

  Epetra_Vector * mp_epetraVector1;
  Epetra_Vector * mp_epetraVector2;

  Callback m_callback;
};

#endif //CallbackVectorLoadTest_h

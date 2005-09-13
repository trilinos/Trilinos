#ifndef PYINTERFACE_H
#define PYINTERFACE_H

#include "Python.h"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "Callback.h"

// Forward declarations
class Epetra_Vector;
class Epetra_BlockMap;

class PyInterface : 
  public NOX::Epetra::Interface::Required,
  public NOX::Epetra::Interface::Jacobian,
  public NOX::Epetra::Interface::Preconditioner
{
public:
  PyInterface(PyObject * = Py_None);
  ~PyInterface();

  bool computeF(             const Epetra_Vector & x,
			     Epetra_Vector       & RHS,
			     FillType ft);

  bool computeJacobian(      const Epetra_Vector & x,
			     Epetra_Operator     & Jac);

  bool computePreconditioner(const Epetra_Vector & x,
			     Epetra_Operator     & M,
			     NOX::Parameter::List* precParams = 0);

  PyObject * setComputeF(             PyObject *);

  PyObject * setComputeJacobian(      PyObject *);

  PyObject * setComputePreconditioner(PyObject *);

  void unloadX  (PyObject *);
  void loadRHS  (PyObject *);
  void unloadRHS(PyObject *);

private:
  // Private and not implemented
  //PyInterface();
  //PyInterface(const PyInterface &);
  //PyInterface & operator=(const PyInterface &);

private:
  PyObject            * mp_problem;
  const Epetra_Vector * mp_x;
  Epetra_Vector       * mp_rhs;
  Epetra_Operator     * mp_jac;
  Epetra_Operator     * mp_precOp;
  Callback              m_computeF;
  Callback              m_computeJacobian;
  Callback              m_computePreconditioner;
};

#endif // PYINTERFACE_H

#ifndef PYINTERFACE_H
#define PYINTERFACE_H

#include "Python.h"
#include "NOX_Epetra_Interface.H"
#include "Callback.h"

// Forward declarations
class Epetra_Vector;
class Epetra_BlockMap;

class PyInterface:
  public NOX::Epetra::Interface
{
public:
  PyInterface(PyObject * = Py_None);
  ~PyInterface();

  bool computeF(             const Epetra_Vector & x,
	                     Epetra_Vector       & RHS,
	                     FillType = NOX::Epetra::Interface::F);

  bool computeJacobian(      const Epetra_Vector & x,
		             Epetra_Operator     & Jac);

  bool computePrecMatrix(    const Epetra_Vector & x,
			     Epetra_RowMatrix    & M);

  bool computePreconditioner(const Epetra_Vector & x,
			     Epetra_Operator     & M);

  PyObject * setComputeF(             PyObject *);

  PyObject * setComputeJacobian(      PyObject *);

  PyObject * setComputePrecMatrix(    PyObject *);

  PyObject * setComputePreconditioner(PyObject *);

  void unloadX  (PyObject *);
  void loadRHS  (PyObject *);
  void unloadRHS(PyObject *);

private:
  // Private and not implemented
  PyInterface();
  PyInterface(const PyInterface &);
  PyInterface & operator=(const PyInterface &);

private:
  PyObject            * mp_problem;
  const Epetra_Vector * mp_x;
  Epetra_Vector       * mp_rhs;
  Epetra_Operator     * mp_jac;
  Epetra_RowMatrix    * mp_precMx;
  Epetra_Operator     * mp_precOp;
  Callback              m_computeF;
  Callback              m_computeJacobian;
  Callback              m_computePrecMatrix;
  Callback              m_computePreconditioner;
};

#endif // PYINTERFACE_H

#include "PyInterface.h"
#include "Epetra_VectorHelper.h"

// Constructor
PyInterface::PyInterface(PyObject * problem) :
  mp_problem(problem)
{
  // Increment the reference count
  Py_XINCREF(mp_problem);

  // If the python object has a "computeF" attribute, set it
  // to be the PyInterface computeF callback function
  if (PyObject_HasAttrString (mp_problem,
			      "computeF")) {
    setComputeF(PyObject_GetAttrString(mp_problem,
				       "computeF"));
  }

  // If the python object has a "computeJacobian" attribute, set it
  // to be the PyInterface computeJacobian callback function
  if (PyObject_HasAttrString (mp_problem,
			      "computeJacobian")) {
    setComputeJacobian(PyObject_GetAttrString(mp_problem,
					      "computeJacobian"));
  }

  // If the python object has a "computePrecMatrix" attribute, set it
  // to be the PyInterface computePrecMatrix callback function
  if (PyObject_HasAttrString (mp_problem,
			      "computePrecMatrix")) {
    setComputePrecMatrix(PyObject_GetAttrString(mp_problem,
						"computePrecMatrix"));
  }

  // If the python object has a "computePeconditioner" attribute, set it
  // to be the PyInterface computePreconditioner callback function
  if (PyObject_HasAttrString (mp_problem,
			      "computePreconditioner")) {
    setComputePreconditioner(PyObject_GetAttrString(mp_problem,
						    "computePreconditioner"));
  }
}

// Destructor
PyInterface::~PyInterface() {
  // Decrement the reference count
  Py_XDECREF(mp_problem);
}

// Compute F
bool PyInterface::computeF(const Epetra_Vector & x,
				  Epetra_Vector       & RHS,
				  FillType              flag)
{
  mp_x   = &x;
  mp_rhs = &RHS;

  PyObject * arglist;
  PyObject * result;

  arglist = Py_BuildValue("()");
  result = PyEval_CallObject(m_computeF.getFunction(), arglist);
  Py_DECREF(arglist);  // All done with argument list

  if (0 == result) {
    PyErr_Clear();
    return false;
  }
  Py_DECREF(result); // All done with returned result object
  return (bool) PyObject_IsTrue(result);
}

// Compute Jacobian operator
bool PyInterface::computeJacobian(const Epetra_Vector & x,
					 Epetra_Operator     & Jac)
{
  mp_x   = &x;
  mp_jac = &Jac;

  PyObject * arglist;
  PyObject * result;

  arglist = Py_BuildValue("()");
  result  = PyEval_CallObject(m_computeJacobian.getFunction(), arglist);
  Py_DECREF(arglist);  // All done with argument list

  if (0 == result) {
    PyErr_Clear();
    return false;
  }
  Py_DECREF(result); // All done with returned result object
  return (bool) PyObject_IsTrue(result);
}

// Compute the preconditioning matrix
bool PyInterface::computePrecMatrix(const Epetra_Vector & x,
					   Epetra_RowMatrix    & M)
{
  mp_x      = &x;
  mp_precMx = &M;

  PyObject * arglist;
  PyObject * result;

  arglist = Py_BuildValue("()");
  result  = PyEval_CallObject(m_computePrecMatrix.getFunction(), arglist);
  Py_DECREF(arglist);  // All done with argument list

  if (0 == result) {
    PyErr_Clear();
    return false;
  }
  Py_DECREF(result); // All done with returned result object
  return (bool) PyObject_IsTrue(result);
}

// Compute the preconditioning operator
bool PyInterface::computePreconditioner(const Epetra_Vector & x,
					       Epetra_Operator     & M)
{
  mp_x      = &x;
  mp_precOp = &M;

  PyObject * arglist;
  PyObject * result;

  arglist = Py_BuildValue("()");
  result  = PyEval_CallObject(m_computePreconditioner.getFunction(), arglist);
  Py_DECREF(arglist);  // All done with argument list

  if (0 == result) {
    PyErr_Clear();
    return false;
  }
  Py_DECREF(result); // All done with returned result object
  return (bool) PyObject_IsTrue(result);
}

// Set the compute F callback function
PyObject * PyInterface::setComputeF(PyObject * p_pyObject)
{
  return m_computeF.setFunction(p_pyObject);
}

// Set the compute Jacobian callback function
PyObject * PyInterface::setComputeJacobian(PyObject * p_pyObject)
{
  return m_computeJacobian.setFunction(p_pyObject);
}

// Set the compute preconditioning matrix callback function
PyObject * PyInterface::setComputePrecMatrix(PyObject * p_pyObject)
{
  return m_computePrecMatrix.setFunction(p_pyObject);
}

// Set the compute preconditioning operator callback function
PyObject * PyInterface::setComputePreconditioner(PyObject * p_pyObject)
{
  return m_computePreconditioner.setFunction(p_pyObject);
}

// Unload from the Epetra_Vector x to a NumPy array
void PyInterface::unloadX  (PyObject * p_pyObject)
{
  Epetra_VectorHelper::unloadViaCopy(mp_x, p_pyObject);
}

// Load from a NumPy array to the Epetra_Vector RHS
void PyInterface::loadRHS  (PyObject * p_pyObject)
{
  Epetra_VectorHelper::loadViaCopy(mp_rhs, p_pyObject);
}

// Unload from the Epetra_Vector RHS to a NumPy array
void PyInterface::unloadRHS(PyObject * p_pyObject)
{
  Epetra_VectorHelper::unloadViaCopy(mp_rhs, p_pyObject);
}

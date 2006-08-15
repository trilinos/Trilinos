// @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
#                 Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#  CVS Information
#  $Source$
#  $Author$
#  $Date$
#  $Revision$
# ************************************************************************
// @HEADER

#include "PyInterface.h"
#include "Epetra_Vector.h"
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
			   Epetra_Vector & RHS,
			   NOX::Epetra::Interface::Required::FillType flag)
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
bool PyInterface::computeJacobian(const Epetra_Vector   & x,
					Epetra_Operator & Jac)
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

// Compute the preconditioning operator
bool PyInterface::computePreconditioner(const Epetra_Vector   & x,
					Epetra_Operator & M,
					NOX::Parameter::List* precParams)
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

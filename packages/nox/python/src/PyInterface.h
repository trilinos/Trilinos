// @HEADER
// ***********************************************************************
//
//               PyTrilinos.NOX: Python Interface to NOX
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

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

// @HEADER
// ***********************************************************************
//
//                 PyTrilinos: Rapid Prototyping Package
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

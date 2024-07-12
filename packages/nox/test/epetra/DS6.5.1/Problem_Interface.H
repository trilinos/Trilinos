// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to the NLS_PetraGroup to provide for
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "NOX_Epetra_Interface_Required.H" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class

// ---------- Forward Declarations ----------
class DennisSchnabel;

class  Problem_Interface : public NOX::Epetra::Interface::Required,
                           public NOX::Epetra::Interface::Jacobian
{
public:
  Problem_Interface(DennisSchnabel& Problem);
  ~Problem_Interface();

  //! Compute and return F.  Returns true if computation was successful.
  bool computeF(const Epetra_Vector& x, Epetra_Vector& FVec,
                NOX::Epetra::Interface::Required::FillType fillType);

  //! Compute an explicit Jacobian.  Returns true if computation was successful.
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is a type of
   *  passdown class design by the application code.
   */
  DennisSchnabel& problem;
};

#endif


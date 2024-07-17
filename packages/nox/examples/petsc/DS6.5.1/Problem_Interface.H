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
#include "petscksp.h"
#include "NOX_Petsc_Interface.H"

// ---------- Forward Declarations ----------
class DennisSchnabel;

class  Problem_Interface : public NOX::Petsc::Interface
{
public:
  Problem_Interface(DennisSchnabel& Problem);
  ~Problem_Interface();

  //! Compute and return RHS.  Returns true if computation was successful.
  bool computeF(const Vec& x, Vec& RHS);

  //! Compute an explicit Jacobian.  Returns true if computation was successful.
  bool computeJacobian(const Vec& x, Mat& Jac);

  //! Compute the matrix M that will be used as the preconditioner.  Returns true if computation was successful.
  bool computePreconditioner(Mat& M);

  //! Return the action of the preconditioner on a vector.  Returns true if computation was successful.
  bool preconditionVector(Vec& y);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of
   *  passdown class design by the application code.
   */
  DennisSchnabel& problem;
};

#endif


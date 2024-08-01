// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to the NOX_Epetra_Group to provide for
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Epetra_Interface_Required.H" // base class
#include "NOX_Epetra_Interface_Jacobian.H" // base class
#include "LOCA_Epetra_Interface_TimeDependent.H" // base class

// ---------- Forward Declarations ----------
class Brusselator;

class Problem_Interface : public virtual LOCA::Epetra::Interface::Required,
                          public NOX::Epetra::Interface::Jacobian,
                          public LOCA::Epetra::Interface::TimeDependent
{
public:
  Problem_Interface(Brusselator& Problem);
  ~Problem_Interface();

  //! Compute and return F (LOCA::Epetra::Interface::Required)
  bool computeF(const Epetra_Vector& x, Epetra_Vector& FVec,
                NOX::Epetra::Interface::Required::FillType fillType);

  //! Compute an explicit Jacobian (NOX::Epetra::Interface::Jacobian)
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Set a parameter in the user's code. (LOCA::Epetra::Interface::Required)
  void setParameters(const LOCA::ParameterVector& params);

  /*!
   * \brief Call user routine for computing the shifted matrix
   * \f$\alpha J + \beta M\f$ where \f$J\f$ is the Jacobian matrix
   * and \f$M\f$ is the mass matrix.
   */
  virtual bool
  computeShiftedMatrix(double alpha, double beta, const Epetra_Vector& x,
               Epetra_Operator& A);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of
   *  passdown class design by the application code.
   */
  Brusselator& problem;
};

#endif


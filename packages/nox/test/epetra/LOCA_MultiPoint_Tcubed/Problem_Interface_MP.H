// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//-----------------------------------------------------------------------------
#ifndef Problem_Interface_MP_H
#define Problem_Interface_MP_H

// Interface to the NLS_PetraGroup to provide for
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>
#include "Epetra_Vector.h"
#include "LOCA_Epetra.H"

// ---------- Forward Declarations ----------
class FiniteElementProblem;
namespace LOCA {
  class ParameterVector;
}


class  Problem_Interface_MP :
  public LOCA::Epetra::Interface::Required,
  public NOX::Epetra::Interface::Jacobian
{
public:
  Problem_Interface_MP(FiniteElementProblem& Problem);
  ~Problem_Interface_MP();

  //! Compute and return F
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag);

  //! Compute an explicit Jacobian
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Set a parameter in the user's code.
  void setParameters(const LOCA::ParameterVector& params);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of
   *  passdown class design by the application code.
   */
  FiniteElementProblem& problem;

  //! Number of residual fills
  int numFillsF;


  /*
    These two functions are added for the multipoint problem.
    All above functions are the same as for regular problem.
   */

  //! Provides data to application for output files.
  virtual void dataForPrintSolution(const int conStep,
                                    const int timeStep,
                                    const int totalTimeSteps);

  //! Set multipoint parameter in the user's application.
  virtual void setMultiPointParameter(const int stepNum);

};

#endif


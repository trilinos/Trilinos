//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
#include "NOX_Epetra_Interface_Preconditioner.H" // base class
#include "LOCA_Epetra_Interface_TimeDependent.H" // base class

// ---------- Forward Declarations ----------
class Brusselator;

class Problem_Interface : public NOX::Epetra::Interface::Jacobian,
                          public NOX::Epetra::Interface::Preconditioner,
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

  //! Call user's own print routine for vector-parameter pair (LOCA::Epetra::Interface::Required)
  virtual void printSolution(const Epetra_Vector& x_, double conParam);

  //! Routine to set data for file indexing in subsequent printSolution call
  virtual void dataForPrintSolution(const int conStep, const int timeStep,
                            const int numTimeSteps);

  //! Calculation of preconditioner (NOX::Epetra::Interface::Preconditioner)
  bool computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M);

  //! Computes a user supplied preconditioner based on input vector x.  Returns true if computation was successful.  (NOX::Epetra::Interface::Preconditioner)
  bool computePreconditioner(const Epetra_Vector& x,
                 Epetra_Operator& Prec,
                 Teuchos::ParameterList* p);

  //! Compute a Time Dependent Shifted Matrix (LOCA::Epetra::Interface::TimeDependent)
  bool computeShiftedMatrix(double alpha, double beta, const Epetra_Vector& x,
                            Epetra_Operator& shMat);

  //! Set the time derivative of the solution for use in the computeF routine
  //  (LOCA::Epetra::Interface::TimeDependent), and absolut time as well
  void setXdot(const Epetra_Vector& xDot, const double time);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of
   *  passdown class design by the application code.
   */
  Brusselator& problem;

private:
  // temporarily keep this routine around while incrementally refactoring
  bool computeMassMatrix(const Epetra_Vector& x);

  int conStep;  // Output continuation step index
  int timeStep;  // Output time step index
};

#endif


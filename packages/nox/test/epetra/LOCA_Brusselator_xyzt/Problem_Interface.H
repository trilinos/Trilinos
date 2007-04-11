//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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


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
//  $Source: /space/CVS/Trilinos/packages/nox/examples/epetra/LOCA_Tcubed/Problem_Interface.H,v $
//  $Author: etphipp $
//  $Date: 2008/10/24 21:04:40 $
//  $Revision: 1.5 $
// ************************************************************************
//@HEADER
                                                                                
//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H

// Interface to the NLS_PetraGroup to provide for 
// residual and matrix fill routines.

// ---------- Standard Includes ----------
#include <iostream>

#include "Epetra_Vector.h"
#include "LOCA_Epetra.H"

#include "LOCA_Parameter_Library.H"
#include "LOCA_Epetra_Interface_FreeEnergy.H"

// ---------- Forward Declarations ----------
class FiniteElementProblem;
namespace LOCA {
  class ParameterVector;
}


class  Problem_Interface : 
  public LOCA::Epetra::Interface::TimeDependent,
  public LOCA::Epetra::Interface::FreeEnergy,
  public NOX::Epetra::Interface::Jacobian
  
{
public:
  Problem_Interface(FiniteElementProblem& Problem);
  ~Problem_Interface();
  
  //! Compute and return F
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType flag);

  //! Compute an explicit Jacobian
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

 /*! 
  * \brief Call user routine for computing the shifted matrix
  * \f$\alpha J + \beta M\f$ where \f$J\f$ is the Jacobian matrix
  * and \f$M\f$ is the mass matrix.
  */
  bool computeShiftedMatrix(double alpha, double beta, 
			    const Epetra_Vector& x,
			     Epetra_Operator& A); 

  //! Set a parameter in the user's code.
  void setParameters(const LOCA::ParameterVector& params);

  //! Print solution to output file
  virtual void printSolution(const Epetra_Vector& x_, double conParam);

  //! Interface for Phase Transition alg
  double computeFreeEnergy(const Epetra_Vector& x_);

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   */ 
  FiniteElementProblem& problem;

  LOCA::Parameter::Library paramLib;

  
};

#endif


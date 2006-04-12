//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
// Finite Element Problem Class
/* Provides function (F) and Jacobian evaluations for the following equation 
 * via a 1D linear finite element discretization with Epetra objects.
 *
 * d2u 
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
 */

#ifndef _NOX_TEST_CONVDIFF_PDE_EQUATION_H
#define _NOX_TEST_CONVDIFF_PDE_EQUATION_H

#include "GenericEpetraProblem.H" // Base class

// Single PDE Equation Problem Class for HMX Cook-off problem

class ConvDiff_PDE : public GenericEpetraProblem
{ 

public:

  // Constructor
  ConvDiff_PDE(
          Epetra_Comm& Comm, 
          double peclet,
          double radiation,
          double kappa,
	  double xmin,
	  double xmax,
	  double Tleft,
	  double Tright,
          int NumGlobalUnknowns = 0,
          string name = "");
  
  // Destructor
  ~ConvDiff_PDE();

  // Set initial condition for solution vector
  void initializeSolution(double val = 0.0);

  // Compute and output analytic slution here
  void initialize();

  // Evaluates the function (F) and/or the Jacobian using the solution 
  // values in solnVector.
  bool evaluate( NOX::Epetra::Interface::Required::FillType fill,
                const Epetra_Vector *solnVector, 
		Epetra_Vector *rhsVector);


protected:

  // inserts the global column indices into the Graph
  void generateGraph();

protected:

  double xmin   ;
  double xmax   ;
  double Tleft  ;
  double Tright ;
  double dx     ;

  double peclet         ;
  double radiation      ;
  double kappa          ;

  Epetra_Import *ColumnToOverlapImporter        ;
  Epetra_Vector *oldSolution                    ;
  Epetra_Vector *rhs                            ;
};
#endif

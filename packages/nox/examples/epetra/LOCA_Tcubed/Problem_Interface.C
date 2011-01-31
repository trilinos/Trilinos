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
//  $Source: /space/CVS/Trilinos/packages/nox/examples/epetra/LOCA_Tcubed/Problem_Interface.C,v $
//  $Author: etphipp $
//  $Date: 2008/10/24 21:04:40 $
//  $Revision: 1.6 $
// ************************************************************************
//@HEADER
                                                                                
// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "Problem_Interface.H"
#include "LOCA_Parameter_Vector.H"

// ----------   User Defined Includes   ----------
#include "FiniteElementProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteElementProblem& Problem) :
  problem(Problem),
  paramLib()
{ 
  paramLib.addParameterEntry("Nonlinear Factor", problem, 
			     &FiniteElementProblem::factor);
  paramLib.addParameterEntry("Left BC", problem, 
			     &FiniteElementProblem::leftBC);
  paramLib.addParameterEntry("Right BC", problem, 
			     &FiniteElementProblem::rightBC);
}

Problem_Interface::~Problem_Interface()
{ 
}

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType flag)
{
  return problem.evaluate(F_ONLY, &x, &F, NULL);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
					Epetra_Operator& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,&problem.getJacobian(),
			  1.0, 0.0);
}

bool Problem_Interface::computeShiftedMatrix(double alpha, double beta, 
					     const Epetra_Vector& x,
					     Epetra_Operator& A)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,&problem.getJacobian(),
			  alpha, beta);
}

void Problem_Interface::setParameters(const LOCA::ParameterVector& params)
{
  for (int i = 0; i < params.length(); i++ ) {
    //problem.set(params.getLabel(i), params.getValue(i)); 
    paramLib.setValue(params.getLabel(i), params.getValue(i));
  }
}

void Problem_Interface::printSolution(const Epetra_Vector& x, double conParam)
{
  double ave; x.MeanValue(&ave);
  cout << "   FreeEnergy Mock Function (param, soln norm, FreeEnergy)= " 
       << conParam << "   " <<  ave << "   " << computeFreeEnergy(x) << endl;

  problem.printSolution(x, conParam);
}

double Problem_Interface::computeFreeEnergy(const Epetra_Vector& x)
{
    double ave; x.MeanValue(&ave);
    double fe =  abs(ave - 1.2);
    return abs(ave - 1.2);
}
//-----------------------------------------------------------------------------


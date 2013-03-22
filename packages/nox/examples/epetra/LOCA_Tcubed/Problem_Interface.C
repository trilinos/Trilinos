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
  std::cout << "   FreeEnergy Mock Function (param, soln norm, FreeEnergy)= " 
       << conParam << "   " <<  ave << "   " << computeFreeEnergy(x) << std::endl;

  problem.printSolution(x, conParam);
}

double Problem_Interface::computeFreeEnergy(const Epetra_Vector& x)
{
    double ave; x.MeanValue(&ave);
    double fe =  abs(ave - 1.2);
    return abs(ave - 1.2);
}
//-----------------------------------------------------------------------------


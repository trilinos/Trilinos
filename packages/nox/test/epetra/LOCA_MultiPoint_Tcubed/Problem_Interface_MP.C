// $Id$
// $Source$

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

// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "Problem_Interface_MP.H"
#include "LOCA_Parameter_Vector.H"

// ----------   User Defined Includes   ----------
#include "FiniteElementProblem.H"

//-----------------------------------------------------------------------------
Problem_Interface_MP::Problem_Interface_MP(FiniteElementProblem& Problem) :
  problem(Problem),
  numFillsF(0)
{ 
}

Problem_Interface_MP::~Problem_Interface_MP()
{ }

bool Problem_Interface_MP::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag)
{
  ++numFillsF;
  return problem.evaluate(F_ONLY, &x, &F, NULL);
}

bool Problem_Interface_MP::computeJacobian(const Epetra_Vector& x,
					Epetra_Operator& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,&problem.getJacobian());
}

void Problem_Interface_MP::setParameters(const LOCA::ParameterVector& params)
{
  for (int i = 0; i < params.length(); i++ ) {
    problem.setParameter(params.getLabel(i), params.getValue(i)); 
  }
}

void Problem_Interface_MP::dataForPrintSolution(const int conStep,
                          const int timeStep,
                          const int totalTimeSteps)
{
  // If this example had a printSolution method, the conStep and
  // timeStep indicies could be used to generate unique file names.
}

// The parameter that is varied in the multipoint run is set here
void Problem_Interface_MP::setMultiPointParameter(const int stepNum)
{
  problem.setParameter("Nonlinear Factor", 1.0 + 0.2*stepNum); 
}

//-----------------------------------------------------------------------------


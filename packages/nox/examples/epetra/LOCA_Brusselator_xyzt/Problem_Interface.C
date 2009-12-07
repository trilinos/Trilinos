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
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "Brusselator.H"
#include "Epetra_CrsMatrix.h"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(Brusselator& Problem) :
  problem(Problem),
  conStep(0),
  timeStep(1)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, 
                       NOX::Epetra::Interface::Required::FillType fillType) 
{
  return problem.evaluate(fillType, &x, &FVec, NULL, 0.0, 0.0);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
					Epetra_Operator& Jac)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0, 0, 1.0, 0.0);
}

void Problem_Interface::setParameters(const LOCA::ParameterVector& params)
{
  problem.setParameters( params.getValue("alpha"), params.getValue("beta") );
}
 
void Problem_Interface::printSolution(const Epetra_Vector& x, double conParam)
{
   if (timeStep==1)
      cout << "Writing solution at continuation step " << conStep
           << "  for parameter = " << conParam << endl;
   char file_name[25];
   FILE *ifp;
   Epetra_Vector& xMesh = problem.getMesh();
   int NumMyNodes = xMesh.Map().NumMyElements();
   (void) sprintf(file_name, "output.p%02d_t%03d_s%03d", xMesh.Comm().MyPID(),
		  timeStep, conStep);
   ifp = fopen(file_name, "w");
   for (int i=0; i<NumMyNodes; i++)
     fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i, xMesh[i],
                       x[2*i], x[2*i+1]);
   fclose(ifp);
   timeStep++; //for time integration runs
}

void Problem_Interface::dataForPrintSolution(const int conStep_,
                       	const int timeStep_, const int numTimeSteps_)
{
  conStep = conStep_; // Zero based
  timeStep = timeStep_ + 1; // One based
}

bool Problem_Interface::computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M)
{
  cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << endl;
  throw 1;
}

bool Problem_Interface::computePreconditioner(const Epetra_Vector& x, 
					      Epetra_Operator& Prec,
					      Teuchos::ParameterList* p)
{
  cout << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jacobian only for this test problem!" << endl;
  throw 1;
}

bool Problem_Interface::computeMassMatrix(const Epetra_Vector& x)
{
  return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0, 0, 0.0, 1.0);
}

bool Problem_Interface::computeShiftedMatrix(double alpha, double beta,
                                    const Epetra_Vector& x, Epetra_Operator& shMat)
{
    return problem.evaluate(NOX::Epetra::Interface::Required::Jac, &x, 0, 0, alpha, beta);
}

void Problem_Interface::setXdot(const Epetra_Vector& xDot, const double time_)
{
   problem.setxDot(xDot);

  // Can set time for nonautonomous systems
}
//-----------------------------------------------------------------------------


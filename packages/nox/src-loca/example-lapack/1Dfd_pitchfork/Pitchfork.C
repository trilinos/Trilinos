// $Id$ 
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include <iomanip>
#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "LOCA_Bifurcation_PitchforkBordGroup.H"
#include "PitchforkProblemInterface.H"

int main()
{
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(14);

  int n = 100;
  double alpha = 1.0;
  double beta = 0.0;
  double lambda = -2.25;
  double pi = 4.0*atan(1.0);
  double h = 2.0/(n-1);

  // Set up the problem interface
  PitchforkProblemInterface pf(n, alpha, beta, lambda);
  LOCA::ParameterVector p;
  p.addParameter("alpha",alpha);
  p.addParameter("beta",beta);
  p.addParameter("lambda",lambda);
  
  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  LOCA::LAPACK::Group grp(pf);
  grp.setParams(p);

  // Create initial guess for the null vector of jacobian
  NOX::LAPACK::Vector nullVec(n);  // length n
  nullVec.init(1.0);             // initial value 1.0

  // Create asymmetry vector
  NOX::LAPACK::Vector asymVec(n);  // length n
  for (int i=0; i<n; i++)
    asymVec(i) = sin( pi/2.0 * (-1.0 + h*i) );

  // Create a turning point group that uses the lapack group
  LOCA::DerivUtils du;
  LOCA::Bifurcation::PitchforkBordGroup pfgrp(grp, asymVec, asymVec, 2, du);

  // Set up the status tests
  NOX::StatusTest::NormF statusTestA(1.0e-8);
  NOX::StatusTest::MaxIters statusTestB(10);
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, 
					  statusTestA, statusTestB);

  // Create the list of solver parameters
  NOX::Parameter::List solverParameters;

  // Set the level of output (this is the default)
  solverParameters.setParameter("Output Information", 
		     NOX::Utils::Warning + 
		     NOX::Utils::OuterIteration + 
		     NOX::Utils::InnerIteration + 
		     NOX::Utils::Parameters);

  // Set the solver (this is the default)
  solverParameters.setParameter("Nonlinear Solver", "Line Search Based");

  // Create the line search parameters sublist
  NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");

  // Set the line search method
  //lineSearchParameters.setParameter("Method","More'-Thuente");
  lineSearchParameters.setParameter("Method","Full Step");

  // Create the newton and  linear solver parameters sublist
  NOX::Parameter::List& directionParameters = solverParameters.sublist("Direction");
  NOX::Parameter::List& newtonParameters = directionParameters.sublist("Newton");
  NOX::Parameter::List& linearSolverParameters = newtonParameters.sublist("Linear Solver");
  linearSolverParameters.setParameter("Bifurcation Solve", "Default");
// The following non-Default choices are still unproven research ideas
//  linearSolverParameters.setParameter("Bifurcation Solve", "Nic");
//  linearSolverParameters.setParameter("Bifurcation Solve", "NicDay");
//  linearSolverParameters.setParameter("Bifurcation Solve", "Iterative Refinement");

  // Create the solver
  NOX::Solver::Manager solver(pfgrp, statusTestsCombo, solverParameters);

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver.solve();

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getParameterList().print(cout);

  // Get the answer
  pfgrp = solver.getSolutionGroup();
  

  // Print the answer
  cout << "\n" << "-- Final Solution From Solver --" << "\n";
  pfgrp.printSolution(0.0);
}

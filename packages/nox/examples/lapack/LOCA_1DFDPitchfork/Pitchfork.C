// $Id$ 
// $Source$

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include <iomanip>
#include "LOCA.H"
#include "LOCA_LAPACK.H"
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
  Teuchos::RefCountPtr<LOCA::Bifurcation::PitchforkBord::ExtendedGroup> pfgrp = Teuchos::rcp(new LOCA::Bifurcation::PitchforkBord::ExtendedGroup(grp, asymVec, asymVec, asymVec,2));

  // Set up the status tests
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> statusTestA = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> statusTestB = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(10));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> statusTestsCombo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					    statusTestA, statusTestB));

  // Create the list of solver parameters
  Teuchos::RefCountPtr<NOX::Parameter::List> solverParametersPtr = 
    Teuchos::rcp(new NOX::Parameter::List);
  NOX::Parameter::List& solverParameters = *solverParametersPtr.get();

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
  NOX::Solver::Manager solver(pfgrp, statusTestsCombo, solverParametersPtr);

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver.solve();

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getParameterList().print(cout);

  // Get the answer
  const LOCA::Bifurcation::PitchforkBord::ExtendedGroup& soln_pfgrp = dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedGroup&>(solver.getSolutionGroup());
  
  // Print the answer
  cout << "\n" << "-- Final Solution From Solver --" << "\n";
  soln_pfgrp.printSolution(0.0);
}

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
#include "ChanProblemInterface.H"

int main()
{
  // cout.setf(ios::scientific, ios::floatfield);
//   cout.precision(14);

  int n = 100;
  double alpha = 3.0;
  double beta = 0.5;
  double scale = 1.0;

  // Set up the problem interface
  ChanProblemInterface tp(n, alpha, beta, scale);
  LOCA::ParameterVector p;
  p.addParameter("alpha",alpha);
  p.addParameter("beta",beta);
  p.addParameter("scale",scale);

  // Create a group which uses that problem interface. The group will
  // be initialized to contain the default initial guess for the
  // specified problem.
  LOCA::LAPACK::Group grp(tp);
  grp.setParams(p);

  // Create initial guess for the null vector of jacobian
  NOX::LAPACK::Vector nullVec(n);  // length 1
  nullVec.init(1.0);             // initial value 1.0

  // Set up the status tests
  //NOX::StatusTest::NormUpdate statusTestA(tpgrp, 1.0e-7);
  NOX::StatusTest::NormF statusTestA(1.0e-11);
  NOX::StatusTest::MaxIters statusTestB(20);
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

  // Create the list of solver parameters
  NOX::Parameter::List solverParameters;

  // Create the NOX printing parameter list
  NOX::Parameter::List& nlPrintParams = solverParameters.sublist("Printing");
  nlPrintParams.setParameter("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning);

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

  // Create a turning point group that uses the lapack group
  LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup tpgrp(grp, nullVec, 
							  nullVec, 0);

  // Create the solver
  NOX::Solver::Manager solver(tpgrp, statusTestsCombo, solverParameters);

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver.solve();

  // Print the answer
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver.getParameterList().print(cout);

  // Get the answer
  tpgrp = solver.getSolutionGroup();

  // Print the answer
  cout << "\n" << "-- Final Solution From Solver --" << "\n";
  tpgrp.printSolution(0.0);
}

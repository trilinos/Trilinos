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
#include "BrusselatorProblemInterface.H"

int main()
{
  double pi = 4.0*atan(1.0);
  int n = 100;
  double alpha = 0.33;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  double beta = 2.0*pi*pi*D1 + 1.0 + alpha*alpha;
  int maxNewtonIters = 20;

  try {

    // Create output file to save solutions
    ofstream outFile("BrusselatorHopf.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << endl;

    // Set up the problem interface
    BrusselatorProblemInterface brus(n, alpha, beta, D1, D2, outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("D1",D1);
    p.addParameter("D2",D2);

    // Create a group which uses that problem interface.  For the Hopf
    // right now we have to set the mass matrix flag as true.
    LOCA::LAPACK::Group grp(brus, true);
    grp.setParams(p);

    // Create initial guess for the real and imaginary eigenvectors
    NOX::LAPACK::Vector y(2*n), z(2*n);
    double h = 1.0 / double(n-1);
    double lambda_real = (beta - 1.0 - alpha*alpha)/2.0;
    double lambda_imag = sqrt(alpha*alpha - lambda_real*lambda_real);
    double v1_real = -alpha*alpha;
    double v1_imag = 0.0;
    double v2_real = beta - 1.0 - lambda_real;
    double v2_imag = -lambda_imag;
    double x;
    for (int i=0; i<n; i++) {
      x = sin(pi*h*i);
      y(i) = v1_real*x;
      z(i) = v1_imag*x;

      y(i+n) = v2_real*x;
      z(i+n) = v2_imag*x;
    }

    // Initial guess for frequency (valid for |alpha| > (pi^2)*|D1|)
    double w = lambda_imag;

    // Create length scaling vector (phi)
    NOX::LAPACK::Vector phi(2*n);
    phi.init(1.0);

    grp.setParam("alpha",0.35);

    // Create a Hopf point group that uses the lapack group
    LOCA::Bifurcation::HopfBord::ExtendedGroup hopfgrp(grp, y, z, phi, w, 1);

    // Set up the status tests
    NOX::StatusTest::NormF statusTestA(1.0e-8, NOX::StatusTest::NormF::Scaled);
    NOX::StatusTest::MaxIters statusTestB(maxNewtonIters);
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
    
    // Create the solver
    NOX::Solver::Manager solver(hopfgrp, statusTestsCombo, solverParameters);

    // Solve the nonlinear system
    NOX::StatusTest::StatusType status = solver.solve();

    // Print the answer
    cout << "\n" << "-- Parameter List From Solver --" << "\n";
    solver.getParameterList().print(cout);
    
    // Get the answer
    hopfgrp = solver.getSolutionGroup();
    grp = hopfgrp.getUnderlyingGroup();

    // Print the answer
    cout << "\n" << "-- Final Solution From Solver --" << "\n";
    hopfgrp.printSolution(grp.getParam("alpha"));

    // Close output file
    outFile.close();

    cout << "beta = " << grp.getParam("beta") << endl;

    grp.computeEigenvalues(solverParameters);
  }

  catch (char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }
}

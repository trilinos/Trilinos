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

#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "BrusselatorProblemInterface.H"

int main()
{
  double pi = 4.0*atan(1.0);
  int n = 100;
  double alpha = 0.25;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  double beta = 2.0*pi*pi*D1 + 1.0 + alpha*alpha;
  int maxNewtonIters = 20;

  try {
    // Create output file to save solutions
    ofstream outFile("BrusselatorHopfContinuation.dat");
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

    // Create a Hopf point group that uses the lapack group
    LOCA::Bifurcation::HopfBord::ExtendedGroup hopfgrp(grp, y, z, phi, w, 1);

    // Create parameter list
    NOX::Parameter::List paramList;

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList.sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.setParameter("Continuation Method", "Natural");
    stepperList.setParameter("Continuation Method", "Arc Length");
    stepperList.setParameter("Continuation Parameter", "alpha");
    stepperList.setParameter("Initial Value", alpha);
    stepperList.setParameter("Max Value", 1.0);
    stepperList.setParameter("Min Value", 0.24);
    stepperList.setParameter("Max Steps", 2);
    stepperList.setParameter("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.setParameter("Goal g", 0.5);
    stepperList.setParameter("Max g", 0.7);
    stepperList.setParameter("Initial Scale Factor", 1.0);
    stepperList.setParameter("Min Scale Factor", 1.0e-8);
    stepperList.setParameter("Min Tangent Factor", -1.0);
    stepperList.setParameter("Tangent Factor Exponent",1.0);
    stepperList.setParameter("Compute Eigenvalues",true);

    // Create predictor sublist
    NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.setParameter("Method", "Constant");
    predictorList.setParameter("Method", "Secant");
    //predictorList.setParameter("Method", "Random");
    //predictorList.setParameter("Epsilon", 1.0e-3);

    NOX::Parameter::List& firstStepPredictor 
      = predictorList.sublist("First Step Predictor");
    firstStepPredictor.setParameter("Method", "Random");
    firstStepPredictor.setParameter("Epsilon", 1.0e-3);

    // Create step size sublist
    NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.setParameter("Method", "Adaptive");
    stepSizeList.setParameter("Initial Step Size", 0.02);
    stepSizeList.setParameter("Min Step Size", 1.0e-3);
    stepSizeList.setParameter("Max Step Size", 0.02);
    stepSizeList.setParameter("Aggressiveness", 0.5);

    // Set the LOCA Utilities
    NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
    locaUtilsList.setParameter("Output Information", 
			       LOCA::Utils::Warning +
			       LOCA::Utils::StepperIteration +
			       LOCA::Utils::StepperDetails +
			       LOCA::Utils::Solver +
			       LOCA::Utils::SolverDetails);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information", 
			       //NOX::Utils::OuterIteration + 
			       //NOX::Utils::OuterIterationStatusTest + 
			       //NOX::Utils::InnerIteration +
			       //NOX::Utils::Parameters +
			       //NOX::Utils::Details + 
			       NOX::Utils::Warning);

    // Create the "Line Search" sublist for the "Line Search Based" solver
    NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
    searchParams.setParameter("Method", "Full Step");

    // Set up the status tests
    NOX::StatusTest::NormF statusTestA(1.0e-8, NOX::StatusTest::NormF::Scaled);
    NOX::StatusTest::MaxIters statusTestB(maxNewtonIters);
    NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

    // Create the stepper  
    LOCA::Stepper stepper(hopfgrp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished)
      cout << "Stepper failed to converge!" << endl;

    // Output the parameter list
    if (NOX::Utils::doPrint(NOX::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      stepper.getParameterList().print(cout);
      cout << endl;
    }
  }

  catch (char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}

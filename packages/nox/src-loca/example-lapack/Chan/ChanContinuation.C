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
#include "LOCA_LAPACK_DataOutput.H"
#include "ChanProblemInterface.H"

int main()
{
  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int maxNewtonIters = 20;

  alpha = alpha / scale;

  try {

    // Set up the problem interface
    ChanProblemInterface chan(n, alpha, beta, scale);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);

    // Create scaling vector
    NOX::LAPACK::Vector s(n);
    s.init(1.0);
    //s.init(1.0/ sqrt((double) n));
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    LOCA::LAPACK::Group grp(chan, s);
    grp.setParams(p);

    // Create DataOutput object to save intermediate continuation points
    fstream file("chan_eqcont.dat", ios::out); 
    file.setf(ios::scientific, ios::floatfield);
    file.precision(15);
    LOCA::LAPACK::DataOutput dataOut(file);

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
    stepperList.setParameter("Max Value", 5.0/scale);
    stepperList.setParameter("Min Value", 0.0/scale);
    stepperList.setParameter("Max Steps", 100);
    stepperList.setParameter("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.setParameter("Goal g", 0.5);
    stepperList.setParameter("Max g", 0.7);
    stepperList.setParameter("Initial Scale Factor", 1.0);
    stepperList.setParameter("Min Scale Factor", 1.0e-8);
    stepperList.setParameter("Min Tangent Factor", -1.0);
    stepperList.setParameter("Tangent Factor Exponent",1.0);

    // Create predictor sublist
    NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.setParameter("Method", "Constant");
    //predictorList.setParameter("Method", "Tangent");
    predictorList.setParameter("Method", "Secant");

    // Create step size sublist
    NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.setParameter("Method", "Constant");
    stepSizeList.setParameter("Method", "Adaptive");
    stepSizeList.setParameter("Initial Step Size", 0.1/scale);
    stepSizeList.setParameter("Min Step Size", 1.0e-3/scale);
    stepSizeList.setParameter("Max Step Size", 10.0/scale);
    //stepSizeList.setParameter("Max Step Size", 1.0);
    stepSizeList.setParameter("Aggressiveness", 0.5);

    // Set the LOCA Utilities
    NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
    locaUtilsList.setParameter("Output Information", 
			       LOCA::Utils::Warning +
			       LOCA::Utils::StepperIteration +
			       LOCA::Utils::StepperDetails +
			       LOCA::Utils::Solver +
			       LOCA::Utils::Parameters +
			       LOCA::Utils::SolverDetails);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information", 
			 //  NOX::Utils::Details +
// 			  NOX::Utils::OuterIteration + 
// 			  NOX::Utils::InnerIteration + 
			  NOX::Utils::Warning);

    NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
    NOX::Parameter::List& lsParams = dirParams.sublist("Linear Solver");

    // Set up the status tests
    NOX::StatusTest::NormF statusTestA(grp, 1.0e-8);
    NOX::StatusTest::MaxIters statusTestB(maxNewtonIters);
    NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);

    // Create the stepper  
    LOCA::Stepper stepper(grp, combo, paramList, dataOut);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished)
      cout << "Stepper failed to converge!" << endl;

    // Get the final solution from the solver
    const LOCA::LAPACK::Group& finalGroup = dynamic_cast<const LOCA::LAPACK::Group&>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution = dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup.getX());

    // Output the parameter list
    if (LOCA::Utils::doPrint(LOCA::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      stepper.getParameterList().print(cout);
      cout << endl;
    }

    file.close();

  }

  catch (string& s) {
    cout << s << endl;
  }
  catch (exception& e) {
    cout << e.what() << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}

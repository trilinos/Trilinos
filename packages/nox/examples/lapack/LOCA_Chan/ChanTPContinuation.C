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

#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "ChanProblemInterface.H"

int main()
{

  try {
    int n = 100;
    double alpha = 4.0;
    double beta = 0.0;
    double scale = 1.0;
    int maxNewtonIters = 10;
    //int maxNewtonIters = 6;
    bool use_initial_guess = false;

    // Create output file to save solutions
    ofstream outFile("ChanTPContinuation.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << endl;

    // Create initial guess for the null vector of jacobian
    Teuchos::RefCountPtr<NOX::Abstract::Vector> nullVec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
   
    double conparam;
    double bifparam;
    NOX::LAPACK::Vector initial_x(n);
    NOX::LAPACK::Vector initial_v(n);
    if (use_initial_guess) {
      // read in initial guess
      ifstream inFile("ChanTPContinuation.dat.at_tp0");
      
      // get size
      int n2;
      inFile >> n2;

      if (n2 != n) {
	cerr << "Error!  Size of discretization in initial guess file" << endl
	     << "does not match specifed size!" << endl;
	throw "NOX::Error";
      }
      
      // get conparam
      inFile >> conparam;
      
      // get x
      for (int i=0; i<n; i++)
	inFile >> initial_x(i);

      // get bifparam
      inFile >> bifparam;

      // get v
      for (int i=0; i<n; i++)
	inFile >> initial_v(i);

      alpha = bifparam;
      beta = conparam;
      *nullVec = initial_v;
    }
    else {
       nullVec->init(1.0);               // initial value 1.0
    }

    // Create initial values for a and b
    Teuchos::RefCountPtr<NOX::Abstract::Vector> a_vec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    a_vec->init(1.0);

    Teuchos::RefCountPtr<NOX::Abstract::Vector> b_vec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));
    b_vec->init(1.0);

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Method", "Arc Length");
    stepperList.set("Bordered Solver Method", "Nested");
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value", 0.0);
    stepperList.set("Max Steps", 20);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.set("Enable Arc Length Scaling", true);
    stepperList.set("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.set("Max Arc Length Parameter Contribution", 0.7);
    stepperList.set("Initial Scale Factor", 1.0);
    stepperList.set("Min Scale Factor", 1.0e-8);
    stepperList.set("Enable Tangent Factor Step Size Scaling",false);
    stepperList.set("Min Tangent Factor", -1.0);
    stepperList.set("Tangent Factor Exponent",1.0);

    Teuchos::ParameterList& nestedList = 
      stepperList.sublist("Nested Bordered Solver");
    nestedList.set("Bordered Solver Method", "LAPACK Direct Solve");

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Turning Point");
    //bifurcationList.set("Formulation", "Moore-Spence");
    bifurcationList.set("Formulation", "Minimally Augmented");
    //bifurcationList.set("Solver Method", "Salinger Bordering");
    bifurcationList.set("Solver Method", "Phipps Bordering");
    bifurcationList.set("Bifurcation Parameter", "alpha");
    bifurcationList.set("Length Normalization Vector", nullVec);
    bifurcationList.set("Initial Null Vector", nullVec);
    bifurcationList.set("Initial A Vector", a_vec);
    bifurcationList.set("Initial B Vector", b_vec);
    bifurcationList.set("Bordered Solver Method", 
				 "LAPACK Direct Solve");

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.set("Method", "Constant");
    predictorList.set("Method", "Secant");
    //predictorList.set("Method", "Random");
    //predictorList.set("Epsilon", 1.0e-3);

    Teuchos::ParameterList& firstStepPredictor 
      = predictorList.sublist("First Step Predictor");
    firstStepPredictor.set("Method", "Random");
    //firstStepPredictor.set("Method", "Constant");
    firstStepPredictor.set("Epsilon", 1.0e-3);

    Teuchos::ParameterList& lastStepPredictor 
      = predictorList.sublist("Last Step Predictor");
    lastStepPredictor.set("Method", "Random");
    //firstStepPredictor.set("Method", "Constant");
    lastStepPredictor.set("Epsilon", 1.0e-3);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 1.0);
    stepSizeList.set("Aggressiveness", 0.5);
    stepSizeList.set("Failed Step Reduction Factor", 0.5);
    stepSizeList.set("Successful Step Increase Factor", 1.26); // for constant

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    nlParams.set("Nonlinear Solver", "Line Search Based");

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information", 
			       NOX::Utils::OuterIteration + 
			       NOX::Utils::OuterIterationStatusTest + 
			       NOX::Utils::InnerIteration +
			       NOX::Utils::Parameters +
			       NOX::Utils::Details + 
			       NOX::Utils::Warning + 
			       NOX::Utils::StepperIteration +
			       NOX::Utils::StepperDetails);
    //nlPrintParams.set("Output Precision", 12);

    // Create the "Line Search" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Create LAPACK Factory
    Teuchos::RefCountPtr<LOCA::LAPACK::Factory> lapackFactory = 
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    ChanProblemInterface chan(globalData, n, alpha, beta, scale, outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grp = 
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, chan));
    
    grp->setParams(p);
    if (use_initial_guess) {
      grp->setX(initial_x);
    }

    // Set up the status tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> statusTestA = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-10, 
					      NOX::StatusTest::NormF::Scaled));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> statusTestB = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      statusTestA, statusTestB));

    // Create the stepper  
    LOCA::NewStepper stepper(globalData, grp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    outFile.close();

    destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    cout << e.what() << endl;
  }
  catch (const char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}

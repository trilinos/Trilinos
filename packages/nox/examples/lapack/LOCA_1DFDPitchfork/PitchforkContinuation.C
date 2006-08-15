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

#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "PitchforkProblemInterface.H"

int main()
{

  try {
    int n = 100;
    double alpha = 1.0;
    double beta = 0.0;
    double lambda = -2.25;
    double pi = 4.0*atan(1.0);
    double h = 2.0/(n-1);
    int maxNewtonIters = 20;

    // Create asymmetry vector
    Teuchos::RefCountPtr<NOX::LAPACK::Vector> asymLapackVec = 
      Teuchos::rcp(new NOX::LAPACK::Vector(n));  // length n
    for (int i=0; i<n; i++)
      (*asymLapackVec)(i) = sin( pi/2.0 * (-1.0 + h*i) );
    Teuchos::RefCountPtr<NOX::Abstract::Vector> asymVec = asymLapackVec;

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");         // Default
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value", -0.01);
    stepperList.set("Max Steps", 20);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Pitchfork");
    bifurcationList.set("Formulation", "Moore-Spence");           // Default
    //bifurcationList.set("Solver Method", "Salinger Bordering"); // Default
    bifurcationList.set("Solver Method", "Phipps Bordering");
    bifurcationList.set("Bordered Solver Method", 
			"LAPACK Direct Solve");  // For Phipps Bordering
    bifurcationList.set("Bifurcation Parameter", "lambda");       // Must set
    bifurcationList.set("Length Normalization Vector", asymVec);  // Must set
    bifurcationList.set("Initial Null Vector", asymVec);          // Must set
    bifurcationList.set("Antisymmetric Vector", asymVec);         // Must set

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");                       // Default

    // Should use for Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& firstStepPredictor 
    //  = predictorList.sublist("First Step Predictor");
    //firstStepPredictor.set("Method", "Random");
    //firstStepPredictor.set("Epsilon", 1.0e-3);

    // Should use for Salinger Bordering & Secant predictor
    //Teuchos::ParameterList& lastStepPredictor 
    //  = predictorList.sublist("Last Step Predictor");
    //lastStepPredictor.set("Method", "Random");
    //lastStepPredictor.set("Epsilon", 1.0e-3);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");                      // Default
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 1.0);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information", 
		      NOX::Utils::OuterIteration + 
		      NOX::Utils::OuterIterationStatusTest + 
		      NOX::Utils::InnerIteration +
		      NOX::Utils::Details + 
		      NOX::Utils::Warning + 
		      NOX::Utils::StepperIteration +
		      NOX::Utils::StepperDetails +
		      NOX::Utils::StepperParameters);

    // Create LAPACK Factory
    Teuchos::RefCountPtr<LOCA::LAPACK::Factory> lapackFactory = 
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    PitchforkProblemInterface pf(n, alpha, beta, lambda);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("lambda",lambda);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grp = 
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, pf));
    
    grp->setParams(p);

    // Set up the status tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> statusTestA = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8, 
					      NOX::StatusTest::NormF::Scaled));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> statusTestB = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
					      statusTestA, statusTestB));

    // Create the stepper  
    LOCA::Stepper stepper(globalData, grp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

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

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
#include "ChanProblemInterface.H"
#include "LOCA_MultiStepper.H"

int main()
{
  int n = 100;
  double alpha = 0.01;
  double beta = 0.01;
  double scale = 1.0;
  int maxNewtonIters = 10;

  alpha = alpha / scale;

  try {

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");         // Default
    stepperList.set("Number of Continuation Parameters", 2);      // Must set
    stepperList.set("Epsilon", 0.1);                              // Default
    stepperList.set("Max Charts", 10000);                         // Must set
    stepperList.set("Verbosity", 1);                              // Default
    stepperList.set("Page Charts", 1);                            // Default
    stepperList.set("Dump Polyhedra", true);                      // Default
    stepperList.set("Dump Centers", false);                       // Default
    stepperList.set("Filename", "MFresults");                     // Default
    stepperList.set("Enable Arc Length Scaling", false);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);  // Should set
    stepperList.set("Aggressiveness", 0.01);
    stepperList.set("Max Solution Component", 30.0);

    // Create sublist for each continuation parameter
    Teuchos::ParameterList& paramList1 = 
      stepperList.sublist("Continuation Parameter 1");
    paramList1.set("Parameter Name", "alpha");
    paramList1.set("Initial Value", alpha);
    paramList1.set("Max Value", 5.0/scale);
    paramList1.set("Min Value", 0.0/scale);
    paramList1.set("Initial Step Size", 0.1/scale);
    paramList1.set("Max Step Size", 0.2/scale);
    paramList1.set("Min Step Size", 1.0e-3/scale);
    
    Teuchos::ParameterList& paramList2 = 
      stepperList.sublist("Continuation Parameter 2");
    paramList2.set("Parameter Name", "beta");
    paramList2.set("Initial Value", beta);
    paramList2.set("Max Value", 2.0);
    paramList2.set("Min Value", 0.0);
    paramList2.set("Initial Step Size", 0.1);
    paramList2.set("Max Step Size", 0.2);
    paramList2.set("Min Step Size", 1.0e-3);


    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");                         // Default

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    //predictorList.set("Method", "Constant");
    predictorList.set("Method", "Tangent");                      // Best choice
    //predictorList.set("Method", "Secant");

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Information", 
		      NOX::Utils::Details +
		      NOX::Utils::OuterIteration + 
		      NOX::Utils::InnerIteration + 
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
    ChanProblemInterface chan(globalData, n, alpha, beta, scale);
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

    // Set up the status tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> normF = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxIters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RefCountPtr<NOX::StatusTest::Generic> comboOR = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      normF, 
					      maxIters));

    // Create the stepper  
    LOCA::MultiStepper stepper(globalData, grp, comboOR, paramList);

    // Perform continuation run
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

    LOCA::destroyGlobalData(globalData);

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

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
#include "BrusselatorProblemInterface.H"

int main()
{
  int n = 100;
  double alpha = 0.25;
  double beta = 0.0;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  int maxNewtonIters = 10;

  try {

    // Create output file to save solutions
    ofstream outFile("BrusselatorContinuation.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << endl;

    // Create parameter list
    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");    // Default
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 2.0);
    stepperList.set("Min Value", 0.0);
    stepperList.set("Max Steps", 100);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);
    
    // Use LAPACK solver for arclength bordered solves for efficiency
    stepperList.set("Bordered Solver Method", "LAPACK Direct Solve");

    // Create Eigensolver list using LAPACK eigensolver
    stepperList.set("Compute Eigenvalues",true);
    Teuchos::ParameterList& eigenList = stepperList.sublist("Eigensolver");
    eigenList.set("Method", "DGGEV");     // LAPACK eigensolver
    eigenList.set("NEV", 10);             // Show 10 eigenvalues
    eigenList.set("Sorting Order", "LR"); // Sort by largest real component

    // Create predictor sublist.  We use the constant predictor because the
    // solution is linear in the continuation parameter.  Using Secant or
    // Tangent is not very interesting, but are perfectly valid
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Constant");
    //predictorList.set("Method", "Tangent");
    //predictorList.set("Method", "Secant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");    // Default
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 10.0);
    stepSizeList.set("Aggressiveness", 0.1);

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

    // Create LAPACK Factory (necessary for LAPACK eigensolver)
    Teuchos::RefCountPtr<LOCA::LAPACK::Factory> lapackFactory = 
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, lapackFactory);

    // Set up the problem interface
    BrusselatorProblemInterface brus(globalData, n, alpha, beta, D1, D2, 
				     outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("D1",D1);
    p.addParameter("D2",D2);
  
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grp = 
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, brus));

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
    LOCA::Stepper stepper(globalData, grp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    // Check for convergence
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

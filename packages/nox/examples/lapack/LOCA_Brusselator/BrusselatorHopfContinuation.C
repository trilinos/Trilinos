// $Id$ 
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
  double pi = 4.0*atan(1.0);
  int n = 100;
  double alpha = 0.25;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  double beta = 2.0*pi*pi*D1 + 1.0 + alpha*alpha;
  int maxNewtonIters = 10;

  try {
    // Create output file to save solutions
    std::ofstream outFile("BrusselatorHopfContinuation.dat");
    outFile.setf(std::ios::scientific, std::ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << std::endl;

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

    Teuchos::RCP<NOX::Abstract::Vector> y_vec =
      Teuchos::rcp(&y,false);
    Teuchos::RCP<NOX::Abstract::Vector> z_vec =
      Teuchos::rcp(&z,false);
    Teuchos::RCP<NOX::Abstract::Vector> phi_vec =
      Teuchos::rcp(&phi,false);

    // Create initial values for a and b for minimally augmented method
    Teuchos::RCP<NOX::Abstract::Vector> a_vec_real = 
      Teuchos::rcp(new NOX::LAPACK::Vector(2*n));
    Teuchos::RCP<NOX::Abstract::Vector> a_vec_imag = 
      Teuchos::rcp(new NOX::LAPACK::Vector(2*n));
    Teuchos::RCP<NOX::Abstract::Vector> b_vec_real = 
      Teuchos::rcp(new NOX::LAPACK::Vector(2*n));
    Teuchos::RCP<NOX::Abstract::Vector> b_vec_imag = 
      Teuchos::rcp(new NOX::LAPACK::Vector(2*n));
    *a_vec_real = *y_vec;
    *a_vec_imag = *z_vec;
    *b_vec_real = *y_vec;
    *b_vec_imag = *z_vec;

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Method", "Arc Length");   // Default
    //stepperList.set("Continuation Method", "Natural");   // Default
    stepperList.set("Continuation Parameter", "alpha");
    stepperList.set("Initial Value", alpha);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value", 0.24);
    stepperList.set("Max Steps", 100);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Hopf");
    bifurcationList.set("Bifurcation Parameter", "beta");        // Must set
    bifurcationList.set("Initial Frequency", w);                 // Must set

//     // For Moore-Spence Formulation
//     bifurcationList.set("Formulation", "Moore-Spence");          // Default
//     bifurcationList.set("Solver Method", "Salinger Bordering");  // Default    
//     bifurcationList.set("Length Normalization Vector", phi_vec); // Must set
//     bifurcationList.set("Initial Real Eigenvector", y_vec);      // Must set
//     bifurcationList.set("Initial Imaginary Eigenvector", z_vec); // Must set

    // For minimally augmented formulation
    bifurcationList.set("Formulation", "Minimally Augmented");
    bifurcationList.set("Initial Real A Vector", a_vec_real);       // Must set
    bifurcationList.set("Initial Imaginary A Vector", a_vec_imag);  // Must set
    bifurcationList.set("Initial Real B Vector", b_vec_real);       // Must set
    bifurcationList.set("Initial Imaginary B Vector", b_vec_imag);  // Must set
    bifurcationList.set("Update Null Vectors Every Continuation Step", true);

    // For minimally augmented method, should set these for good performance
    // Direct solve of bordered equations
    bifurcationList.set("Bordered Solver Method",  "LAPACK Direct Solve");
    // Combine arc-length and turning point bordered rows & columns
    stepperList.set("Bordered Solver Method", "Nested");
    Teuchos::ParameterList& nestedList = 
      stepperList.sublist("Nested Bordered Solver");
    // Direct solve of combined bordered system
    nestedList.set("Bordered Solver Method", "LAPACK Direct Solve");

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");     // Default

//     // Should use w/Secant predictor & Moore-Spence formulation
//     Teuchos::ParameterList& firstStepPredictor 
//       = predictorList.sublist("First Step Predictor");
//     firstStepPredictor.set("Method", "Random");
//     firstStepPredictor.set("Epsilon", 1.0e-3);

//     // Should use w/Secant predictor & Moore-Spence fomulation
//     Teuchos::ParameterList& lastStepPredictor 
//       = predictorList.sublist("Last Step Predictor");
//     lastStepPredictor.set("Method", "Random");
//     lastStepPredictor.set("Epsilon", 1.0e-3);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");      // Default
    stepSizeList.set("Initial Step Size", 0.02);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 0.1);
    stepSizeList.set("Aggressiveness", 0.5);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("Output Precision", 3);
    nlPrintParams.set("Output Information", 
		      NOX::Utils::OuterIteration + 
		      NOX::Utils::OuterIterationStatusTest + 
		      NOX::Utils::InnerIteration +
		      NOX::Utils::Details + 
		      NOX::Utils::Warning + 
		      NOX::Utils::StepperIteration +
		      NOX::Utils::StepperDetails +
		      NOX::Utils::StepperParameters);

    // Create the "Line Search" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Create LAPACK Factory
    Teuchos::RCP<LOCA::LAPACK::Factory> lapackFactory = 
      Teuchos::rcp(new LOCA::LAPACK::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
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
    Teuchos::RCP<LOCA::LAPACK::Group> grp = 
      Teuchos::rcp(new LOCA::LAPACK::Group(globalData, brus));

    grp->setParams(p);

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-10, 
					      NOX::StatusTest::NormF::Scaled));
    Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					      statusTestA, statusTestB));

    // Create the stepper  
    LOCA::Stepper stepper(globalData, grp, combo, paramList);

    // Solve the nonlinear system
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status == LOCA::Abstract::Iterator::Finished) 
      std::cout << "All examples passed" << std::endl;
    else {
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

    LOCA::destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
  }

  return 0;
}

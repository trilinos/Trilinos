// @HEADER
// ************************************************************************
// @HEADER

/*! \file  test_02.cpp
    \brief Test Rosenbrock.
*/

#define USE_HESSVEC 1

#include "ROL_TestObjectives.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Test body.

  try {

    Teuchos::ParameterList parlist;
    // Step Information
    parlist.set("Secant Type",                          ROL::SECANT_LBFGS);
    // Secant Information
    parlist.set("Use Secant Preconditioning",           false); 
    parlist.set("Maximum Secant Storage",               10);
    parlist.set("Barzilai-Borwein Type",                1);
    // Inexactness Information
    parlist.set("Use Inexact Objective Function",       false);
    parlist.set("Use Inexact Gradient",                 false);
    parlist.set("Use Inexact Hessian-Times-A-Vector",   false);
    // Trust-Region Parameters
    parlist.set("Initial Trust-Region Radius",          10.0);
    parlist.set("Minimum Trust-Region Radius",          1.e-8);
    parlist.set("Maximum Trust-Region Radius",          5000.0);
    parlist.set("Step Acceptance Parameter",            0.05);
    parlist.set("Radius Shrinking Threshold",           0.05);
    parlist.set("Radius Growing Threshold",             0.9);
    parlist.set("Radius Shrinking Rate (Negative rho)", 0.0625);
    parlist.set("Radius Shrinking Rate (Positive rho)", 0.25);
    parlist.set("Radius Growing Rate",                  2.5);
    parlist.set("Trust-Region Safeguard",               1.0);
    // CG Parameters
    parlist.set("Absolute CG Tolerance",                1.e-4);
    parlist.set("Relative CG Tolerance",                1.e-2);
    

    // Define Status Test
    ROL::StatusTest<RealT> status(1.e-10,1.e-12,1000);    

    // Loop Through Test Objectives
    for ( ROL::ETestObjectives objFunc = ROL::TESTOBJECTIVES_ROSENBROCK; objFunc < ROL::TESTOBJECTIVES_LAST; objFunc++ ) {
      *outStream << "\n\n" << ROL::ETestObjectivesToString(objFunc) << "\n\n";

      // Initial Guess Vector 
      Teuchos::RCP<std::vector<RealT> > x0_rcp = Teuchos::rcp( new std::vector<RealT> );
      ROL::StdVector<RealT> x0(x0_rcp);

      // Exact Solution Vector
      Teuchos::RCP<std::vector<RealT> > z_rcp = Teuchos::rcp( new std::vector<RealT> );
      ROL::StdVector<RealT> z(z_rcp);

      // Get Objective Function
      Teuchos::RCP<ROL::Objective<RealT> > obj = Teuchos::null;
      ROL::getTestObjectives<RealT>(obj,x0,z,objFunc);

      // Get Dimension of Problem
      int dim = 
        Teuchos::rcp_const_cast<std::vector<RealT> >((Teuchos::dyn_cast<ROL::StdVector<RealT> >(x0)).getVector())->size();
      parlist.set("Maximum Number of CG Iterations", 2*dim);

      // Iteration Vector
      Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
      ROL::StdVector<RealT> x(x_rcp);
      x.set(x0);

      // Error Vector
      Teuchos::RCP<std::vector<RealT> > e_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
      ROL::StdVector<RealT> e(e_rcp);
      e.zero();

      for ( ROL::ETrustRegion tr = ROL::TRUSTREGION_CAUCHYPOINT; tr < ROL::TRUSTREGION_LAST; tr++ ) {
        *outStream << "\n\n" << ROL::ETrustRegionToString(tr) << "\n\n";
        parlist.set("Trust-Region Subproblem Solver Type", tr);
        if ( tr == ROL::TRUSTREGION_DOGLEG || tr == ROL::TRUSTREGION_DOUBLEDOGLEG ) {
          parlist.set("Use Secant Hessian-Times-A-Vector", false);
        } 
        else {
          parlist.set("Use Secant Hessian-Times-A-Vector", false);
        }

        // Define Step
        ROL::TrustRegionStep<RealT> step(parlist);

        // Define Algorithm
        ROL::DefaultAlgorithm<RealT> algo(step,status,false);

        // Run Algorithm
        x.set(x0);
        std::vector<std::string> output = algo.run(x, *obj);
        for ( unsigned i = 0; i < output.size(); i++ ) {
          std::cout << output[i];
        }

        // Compute Error 
        e.set(x);
        e.axpy(-1.0,z);
        *outStream << "\nNorm of Error: " << e.norm() << "\n";
      }
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}


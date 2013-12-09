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
    // Enumerations
    parlist.set("Linesearch Type",                        ROL::LINESEARCH_CUBICINTERP);
    parlist.set("Linesearch Curvature Condition",         ROL::CURVATURECONDITION_STRONGWOLFE);
    parlist.set("Secant Type",                            ROL::SECANT_LBFGS);
    // Inexactness Information
    parlist.set("Use Inexact Objective Function",         false);
    parlist.set("Use Inexact Gradient",                   false);
    parlist.set("Use Inexact Hessian-Times-A-Vector",     true);
    // Secant Information
    parlist.set("Maximum Secant Storage",                 10);
    parlist.set("Barzilai-Borwein Type",                  1);
    // Linesearch Parameters
    parlist.set("Maximum Number of Function Evaluations", 20);
    parlist.set("Sufficient Decrease Parameter",          1.e-4);
    parlist.set("Curvature Conditions Parameter",         0.9);
    parlist.set("Bracketing Tolerance",                   1.e-8);
    parlist.set("Backtracking Rate",                      0.5);
    parlist.set("Initial Linesearch Parameter",           1.0);
    parlist.set("User Defined Linesearch Parameter",      false);
    // Krylov Parameters
    parlist.set("Absolute Krylov Tolerance",              1.e-4);
    parlist.set("Relative Krylov Tolerance",              1.e-2);                    

    // Define Status Test
    ROL::StatusTest<RealT> status(1.e-6,1.e-12,1000);    

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
      parlist.set("Maximum Number of Krylov Iterations", 2*dim);

      // Iteration Vector
      Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
      ROL::StdVector<RealT> x(x_rcp);
      x.set(x0);

      // Error Vector
      Teuchos::RCP<std::vector<RealT> > e_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
      ROL::StdVector<RealT> e(e_rcp);
      e.zero();

      for ( ROL::EDescent desc = ROL::DESCENT_STEEPEST; desc < ROL::DESCENT_LAST; desc++ ) {
        parlist.set("Descent Type", desc);
        if ( desc == ROL::DESCENT_NEWTON && 
             ((objFunc == ROL::TESTOBJECTIVES_LEASTSQUARES) || objFunc == ROL::TESTOBJECTIVES_POISSONCONTROL) ) {
          parlist.set("Descent Type", ROL::DESCENT_NEWTONKRYLOV);
        }
        *outStream << "\n\n" << ROL::EDescentToString(desc) << "\n\n";

        // Define Step
        ROL::LineSearchStep<RealT> step(parlist);
      
        // Define Algorithm
        ROL::DefaultAlgorithm<RealT> algo(step,status);

        // Run Algorithm
        x.set(x0);
        algo.run(x, *obj);

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


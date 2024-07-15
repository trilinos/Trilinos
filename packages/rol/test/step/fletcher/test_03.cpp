// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_03.cpp
*/

#define USE_HESSVEC 1

#include "ROL_GetTestProblems.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"


#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Test body.

  try {

    std::string filename = "input_ex03.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",true);
#if USE_HESSVEC
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
#endif

    // Krylov parameters.
    parlist->sublist("Step").set("Type","Fletcher");

    for ( ROL::ETestOptProblem prob = ROL::TESTOPTPROBLEM_ROSENBROCK; prob < ROL::TESTOPTPROBLEM_LAST; prob++ ) { 
      if ((prob != ROL::TESTOPTPROBLEM_HS55) && (prob != ROL::TESTOPTPROBLEM_HS14) && (prob != ROL::TESTOPTPROBLEM_HS63) && (prob != ROL::TESTOPTPROBLEM_HS41)) {
        // Get Objective Function
        ROL::Ptr<ROL::Vector<RealT> > x0;
        std::vector<ROL::Ptr<ROL::Vector<RealT> > > z;
        ROL::Ptr<ROL::OptimizationProblem<RealT> > problem;
        ROL::GetTestProblem<RealT>(problem,x0,z,prob);

        if ((problem->getProblemType() == ROL::TYPE_E || problem->getProblemType() == ROL::TYPE_EB)
            && (prob != ROL::TESTOPTPROBLEM_CANTILEVERBEAM && prob != ROL::TESTOPTPROBLEM_QUARTIC && prob != ROL::TESTOPTPROBLEM_CANTILEVER && prob != ROL::TESTOPTPROBLEM_CYLINDERHEAD)) {
          *outStream << std::endl << std::endl << ROL:: ETestOptProblemToString(prob)  << std::endl << std::endl;

          // Get Dimension of Problem
          // int dim = x0->dimension(); 
          // parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);

          // Error Vector
          ROL::Ptr<ROL::Vector<RealT> > e = x0->clone();
          e->zero();
          
          // Define Solver
          ROL::OptimizationSolver<RealT> solver(*problem,*parlist);

          // Run Solver
          solver.solve(*outStream);

          // Compute Error
          RealT err(0);
          for (int i = 0; i < static_cast<int>(z.size()); ++i) {
            e->set(*x0);
//std::cout << "\n\n  e dim =" << e->dimension() << "  z[i] dim =" << z[i]->dimension() << "\n\n";
            e->axpy(-1.0,*z[i]);
            if (i == 0) {
              err = e->norm();
            }
            else {
              err = std::min(err,e->norm());
            }
          }
          *outStream << std::endl << "Norm of Error: " << err << std::endl;

          // Update error flag
          ROL::Ptr<const ROL::AlgorithmState<RealT> > state = solver.getAlgorithmState();
          errorFlag += ((err < std::max(1.e-6*z[0]->norm(),1.e-8) || (state->gnorm < 1.e-6)) ? 0 : 1);
        }
      }
   }
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}

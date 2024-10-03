// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_04.cpp
    \brief Test bound constrained trust-region steps.
*/

#define USE_HESSVEC 1

#include "ROL_GetTestProblems.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_TypeB_LinMoreAlgorithm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <unordered_map>
//#include <fenv.h>

typedef double RealT;

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

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

    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("General").set("Output Level",iprint);
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",true);
#if USE_HESSVEC
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
#endif

    std::unordered_map<ROL::ETestOptProblem, int, std::hash<int>> iters;

    int totProb = 0;
    for ( ROL::ETestOptProblem prob = ROL::TESTOPTPROBLEM_ROSENBROCK; prob < ROL::TESTOPTPROBLEM_LAST; prob++ ) { 
      // Get Objective Function
      ROL::Ptr<ROL::Vector<RealT>> x0;
      std::vector<ROL::Ptr<ROL::Vector<RealT>>> z;
      ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
      ROL::GetTestProblem<RealT>(problem,x0,z,prob);
      if (problem->getProblemType() == ROL::TYPE_B) {
        totProb++;
        if ( prob == ROL::TESTOPTPROBLEM_HS2 || prob == ROL::TESTOPTPROBLEM_BVP ) {
          parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",-1.e1);
          parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size",1.e-4);
          parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
        }
        else if ( prob == ROL::TESTOPTPROBLEM_HS25 ) {
          parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",1.e3);
          parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size",1.e4);
          parlist->sublist("Status Test").set("Gradient Tolerance",1.e-8);
        }
        else {
          parlist->sublist("Step").sublist("Trust Region").set("Initial Radius",-1.e1);
          parlist->sublist("Step").sublist("Trust Region").set("Safeguard Size",1.e4);
          parlist->sublist("Status Test").set("Gradient Tolerance",1.e-6);
        }
        *outStream << std::endl << std::endl << ROL:: ETestOptProblemToString(prob)  << std::endl << std::endl;

        // Get Dimension of Problem
        int dim = x0->dimension();
        parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);

        // Error Vector
        ROL::Ptr<ROL::Vector<RealT>> e = problem->getSolutionVector()->clone();
        e->zero();

        // Define Solver
        ROL::TypeB::LinMoreAlgorithm<RealT> algo(*parlist);
        algo.run(*problem->getSolutionVector(),
                 *problem->getObjective(),
                 *problem->getBoundConstraint());

        // Compute Error
        RealT err(ROL::ROL_INF<RealT>());
        for (int i = 0; i < static_cast<int>(z.size()); ++i) {
          e->set(*problem->getSolutionVector());
          e->axpy(static_cast<RealT>(-1),*z[i]);
          err = std::min(err,e->norm());
        }
        *outStream << std::endl << "Norm of Error: " << err << std::endl;
        iters.insert({prob,algo.getState()->iter});
      }
    }
    *outStream << std::endl << std::string(80,'=') << std::endl;
    *outStream << "Performance Summary" << std::endl;
    *outStream << "  "
               << std::setw(45) << std::left << "Problem"
               << std::setw(20) << std::left << "Iteration Count"
               << std::endl; 
    *outStream << std::string(80,'-') << std::endl;
    int totIter = 0;
    for ( ROL::ETestOptProblem prob = ROL::TESTOPTPROBLEM_ROSENBROCK; prob < ROL::TESTOPTPROBLEM_LAST; prob++ ) { 
      if (iters.count(prob)>0) {
        *outStream << "  "
                   << std::setw(45) << std::left << ROL::ETestOptProblemToString(prob)
                   << std::setw(20) << std::left << iters[prob]
                   << std::endl;
        totIter += iters[prob];
      }
    }
    *outStream << std::string(80,'-') << std::endl;
    *outStream << "  Total Iterations:    " << totIter << std::endl;
    *outStream << "  Average Iterations:  " << static_cast<RealT>(totIter)/static_cast<RealT>(totProb) << std::endl;
    *outStream << std::string(80,'=') << std::endl << std::endl;
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


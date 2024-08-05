// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test line search.
*/

#define USE_HESSVEC 1

#include "ROL_GetTestProblems.hpp"
#include "ROL_TypeU_BundleAlgorithm.hpp"
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

    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("Step").sublist("Bundle").set("Epsilon Solution Tolerance",1e-4);
    parlist->sublist("General").set("Output Level", iprint);

    ROL::Ptr<ROL::TypeU::Algorithm<RealT>> algo;
    ROL::Ptr<ROL::Vector<RealT>> e, x, x0;

    for ( ROL::ETestOptProblem objFunc = ROL::TESTOPTPROBLEM_ROSENBROCK;
          objFunc < ROL::TESTOPTPROBLEM_LAST; objFunc++ ) {
      std::vector<ROL::Ptr<ROL::Vector<RealT>>> z;
      ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
      ROL::GetTestProblem<RealT>(problem,x0,z,objFunc);
      e = x0->clone();
      x = x0->clone();
      if (problem->getProblemType() == ROL::TYPE_U) {
        if (objFunc == ROL::TESTOPTPROBLEM_MINIMAX1
            || objFunc == ROL::TESTOPTPROBLEM_MINIMAX2
            || objFunc == ROL::TESTOPTPROBLEM_MINIMAX3) {
          *outStream << std::endl << std::endl
                     << ROL::ETestOptProblemToString(objFunc)
                     << std::endl << std::endl;
          // Define Line Search Solver
          algo = ROL::makePtr<ROL::TypeU::BundleAlgorithm<RealT>>(*parlist);

          // Run Solver
          x->set(*x0);
          algo->run(*x,
                    *problem->getObjective(),
                    *outStream);

          // Compute Error
          e->zero();
          RealT err(0);
          for (int i = 0; i < static_cast<int>(z.size()); ++i) {
            e->set(*x);
            e->axpy(-1.0,*z[i]);
            if (i == 0) {
              err = e->norm();
            }
            else {
              err = std::min(err,e->norm());
            }
          }
          *outStream << std::endl << "Norm of Error: " << err << std::endl;
          //errorFlag += (int)(e.norm() < std::sqrt(ROL::ROL_EPSILON<RealT>()));
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

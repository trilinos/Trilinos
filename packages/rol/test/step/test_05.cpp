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

    std::string filename = "input.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",true);
#if USE_HESSVEC
    parlist->sublist("General").set("Inexact Hessian-Times-A-Vector",false);
#endif

    // Krylov parameters.
    parlist->sublist("General").sublist("Krylov").set("Type", "Conjugate Residuals");
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance", 1.e-8);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance", 1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 50);
    parlist->sublist("Step").set("Type","Primal Dual Active Set");

    for ( ROL::ETestOptProblem prob = ROL::TESTOPTPROBLEM_ROSENBROCK; prob < ROL::TESTOPTPROBLEM_LAST; prob++ ) { 
      // Get Objective Function
      ROL::Ptr<ROL::Vector<RealT> > x0;
      std::vector<ROL::Ptr<ROL::Vector<RealT> > > z;
      ROL::Ptr<ROL::OptimizationProblem<RealT> > problem;
      ROL::GetTestProblem<RealT>(problem,x0,z,prob);

      if (problem->getProblemType() == ROL::TYPE_B) {
        if ( prob != ROL::TESTOPTPROBLEM_HS5 ) {
          // PDAS parameters.
          if (prob == ROL::TESTOPTPROBLEM_HS1 ||
              prob == ROL::TESTOPTPROBLEM_HS2 ||
              prob == ROL::TESTOPTPROBLEM_HS3 ||
              prob == ROL::TESTOPTPROBLEM_HS4 ||
              prob ==  ROL::TESTOPTPROBLEM_HS45) {
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-10);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-8);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit",1);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",1.e8);
          }
          else if (prob == ROL::TESTOPTPROBLEM_HS5) {
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-10);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-8);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit",10);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",1.e-2);
          }
          else if (prob == ROL::TESTOPTPROBLEM_HS25) {
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-10);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-8);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit",10);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",1.e10);
          }
          else if (prob == ROL::TESTOPTPROBLEM_HS38) {
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-10);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-8);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit",1);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",1.e-3);
          }
          else if (prob == ROL::TESTOPTPROBLEM_BVP) {
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-10);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-8);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit",1);
            parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",1.e0);
          }
          *outStream << std::endl << std::endl << ROL:: ETestOptProblemToString(prob)  << std::endl << std::endl;
  
          // Get Dimension of Problem
          int dim = x0->dimension(); 
          parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);
  
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

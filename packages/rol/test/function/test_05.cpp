// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_05.cpp
    \brief Shows how to use the nonlinear least squares interface
           to find a feasible point for the equality constrained NLP
           from Nocedal/Wright, 2nd edition, page 574, example 18.2.
*/

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_NonlinearLeastSquaresObjective.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
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

  // *** Example body.

  try {

    ROL::Ptr<ROL::Objective<RealT>> obj;
    ROL::Ptr<ROL::Constraint<RealT>> constr;
    ROL::Ptr<ROL::Vector<RealT>> x;
    ROL::Ptr<ROL::Vector<RealT>> sol;

    // Retrieve objective, constraint, iteration vector, solution vector.
    ROL::ZOO::getSimpleEqConstrained<RealT> SEC;
    obj    = SEC.getObjective();
    constr = SEC.getEqualityConstraint();
    x      = SEC.getInitialGuess();
    sol    = SEC.getSolution();

    // Inititalize vectors
    int dim = 5;
    int nc = 3;
    RealT left = -1e0, right = 1e0;
    ROL::Ptr<std::vector<RealT>> xtest_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT>> g_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT>> d_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT>> v_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT>> vc_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
    ROL::Ptr<std::vector<RealT>> vl_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
    ROL::StdVector<RealT> xtest(xtest_ptr);
    ROL::StdVector<RealT> g(g_ptr);
    ROL::StdVector<RealT> d(d_ptr);
    ROL::StdVector<RealT> v(v_ptr);
    ROL::StdVector<RealT> vc(vc_ptr);
    ROL::StdVector<RealT> vl(vl_ptr);
    // set xtest, d, v
    for (int i=0; i<dim; i++) {
      (*xtest_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*d_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*v_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }
    // set vc, vl
    for (int i=0; i<nc; i++) {
      (*vc_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*vl_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }

    xtest.set(*x);

    // Initialize nonlinear least squares objectives
    ROL::NonlinearLeastSquaresObjective<RealT> nlls(constr,*x,vc,false);
    ROL::NonlinearLeastSquaresObjective<RealT> gnnlls(constr,*x,vc,true);

    // Check derivatives
    constr->checkApplyJacobian(xtest, v, vc, true, *outStream);                 *outStream << "\n";
    constr->checkApplyAdjointJacobian(xtest, vl, vc, xtest, true, *outStream);  *outStream << "\n";
    constr->checkApplyAdjointHessian(xtest, vl, d, xtest, true, *outStream);    *outStream << "\n";
    nlls.checkGradient(xtest, d, true, *outStream);                             *outStream << "\n"; 
    nlls.checkHessVec(xtest, v, true, *outStream);                              *outStream << "\n";
    nlls.checkHessSym(xtest, d, v, true, *outStream);                           *outStream << "\n";
    
    // Define algorithm.
    ROL::ParameterList parlist;
    parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-10);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-10);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-18);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Ptr<ROL::StatusTest<RealT>> status = ROL::makePtr<ROL::StatusTest<RealT>>(parlist);
    ROL::Ptr<ROL::Step<RealT>> step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Run Algorithm
    *outStream << "\nSOLVE USING FULL HESSIAN\n";
    x->set(xtest);
    algo.run(*x, nlls, true, *outStream);
    algo.reset();
    *outStream << "\nSOLVE USING GAUSS-NEWTON HESSIAN\n";
    x->set(xtest);
    algo.run(*x, gnnlls, true, *outStream);
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}


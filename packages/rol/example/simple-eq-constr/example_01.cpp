// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the equality constrained NLP
           from Nocedal/Wright, 2nd edition, page 574, example 18.2.
*/

#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"
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

    int dim = 5;
    int nc = 3;
    ROL::Ptr<ROL::Objective<RealT> > obj;
    ROL::Ptr<ROL::Constraint<RealT> > constr;
    ROL::Ptr<std::vector<RealT>> x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT>> sol_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::StdVector<RealT> x(x_ptr);      // Iteration vector.
    ROL::StdVector<RealT> sol(sol_ptr);  // Reference solution vector.

    // Retrieve objective, constraint, iteration vector, solution vector.
    ROL::ZOO::getSimpleEqConstrained <RealT> SEC;
    obj = SEC.getObjective();
    constr = SEC.getEqualityConstraint();
    x.set(*SEC.getInitialGuess());
    sol.set(*SEC.getSolution());

    // Run derivative checks, etc.
    RealT left = -1e0, right = 1e0;
    ROL::Ptr<std::vector<RealT> > xtest_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > g_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > d_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > v_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > vc_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
    ROL::Ptr<std::vector<RealT> > vl_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
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
    obj->checkGradient(xtest, d, true, *outStream);                             *outStream << "\n"; 
    obj->checkHessVec(xtest, v, true, *outStream);                              *outStream << "\n";
    obj->checkHessSym(xtest, d, v, true, *outStream);                           *outStream << "\n";
    constr->checkApplyJacobian(xtest, v, vc, true, *outStream);                 *outStream << "\n";
    constr->checkApplyAdjointJacobian(xtest, vl, vc, xtest, true, *outStream);  *outStream << "\n";
    constr->checkApplyAdjointHessian(xtest, vl, d, xtest, true, *outStream);    *outStream << "\n";

    ROL::Ptr<std::vector<RealT> > v1_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > v2_ptr = ROL::makePtr<std::vector<RealT>>(nc, 0.0);
    ROL::StdVector<RealT> v1(v1_ptr);
    ROL::StdVector<RealT> v2(v2_ptr);
    RealT augtol = 1e-8;
    constr->solveAugmentedSystem(v1, v2, d, vc, xtest, augtol);
    
    // Define algorithm.
    ROL::ParameterList parlist;
    std::string stepname = "Composite Step";
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Nominal Relative Tolerance",1.e-4);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist.sublist("Step").sublist(stepname).set("Output Level",0);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-18);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(parlist);
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::CompositeStep<RealT>>(parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    // Run Algorithm
    vl.zero();
    //(*x_ptr)[0] = 3.0; (*x_ptr)[1] = 2.0; (*x_ptr)[2] = 2.0; (*x_ptr)[3] = 1.0; (*x_ptr)[4] = 1.0;
    //(*x_ptr)[0] = -5.0; (*x_ptr)[1] = -5.0; (*x_ptr)[2] = -5.0; (*x_ptr)[3] = -6.0; (*x_ptr)[4] = -6.0;

    algo.run(x, g, vl, vc, *obj, *constr, true, *outStream);

    // Compute Error
    *outStream << "\nReference solution x_r =\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[0] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[1] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[2] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[3] << "\n";
    *outStream << std::scientific << "  " << (*sol_ptr)[4] << "\n";
    *outStream << "\nOptimal solution x =\n";
    *outStream << std::scientific << "  " << (*x_ptr)[0] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[1] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[2] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[3] << "\n";
    *outStream << std::scientific << "  " << (*x_ptr)[4] << "\n";
    x.axpy(-1.0, sol);
    RealT abserr = x.norm();
    RealT relerr = abserr/sol.norm();
    *outStream << std::scientific << "\n   Absolute Error: " << abserr;
    *outStream << std::scientific << "\n   Relative Error: " << relerr << "\n";
    if ( relerr > sqrt(ROL::ROL_EPSILON<RealT>()) ) {
      errorFlag += 1;
    }
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


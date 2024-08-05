// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_23.cpp
    \brief Validate projected gradient algorithm.
*/

#include "ROL_HS14.hpp"
#include "ROL_HS32.hpp"
#include "ROL_HS63.hpp"
#include "ROL_TypeG_AugmentedLagrangianAlgorithm.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a
  // (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {
    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());

    ROL::ParameterList list;
    list.sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Iteration Limit",20);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling",false);
    list.sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    //list.sublist("Step").sublist("Augmented Lagrangian").set("Print Intermediate Optimization History",true);
    list.sublist("Status Test").set("Gradient Tolerance",1e-8);
    list.sublist("Status Test").set("Constraint Tolerance",1e-8);
    list.sublist("Status Test").set("Step Tolerance",1e-12);
    list.sublist("Status Test").set("Iteration Limit", 250);
    list.sublist("General").set("Output Level",iprint);

    ROL::Ptr<ROL::Vector<RealT>>     sol, emul, imul;
    ROL::Ptr<ROL::Objective<RealT>>  obj;
    ROL::Ptr<ROL::Constraint<RealT>> econ, icon;
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd, ibnd;
    ROL::Ptr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>> algo;
    std::vector<RealT> data;
    RealT e1, e2, err;

    *outStream << std::endl << "Hock and Schittkowski Problem #14" << std::endl << std::endl;
    ROL::ZOO::getHS14<RealT> HS14;
    obj  = HS14.getObjective();
    sol  = HS14.getInitialGuess();
    econ = HS14.getEqualityConstraint();
    emul = HS14.getEqualityMultiplier();
    icon = HS14.getInequalityConstraint();
    imul = HS14.getInequalityMultiplier();
    ibnd = HS14.getSlackBoundConstraint();

    algo = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algo->run(*sol,*obj,*icon,*imul,*ibnd,*econ,*emul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0.5 *(std::sqrt(7.0)-1.0)));
    e2 = (data[1]-static_cast<RealT>(0.25*(std::sqrt(7.0)+1.0)));
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    RealT e3;
    *outStream << std::endl << "Hock and Schittkowski Problem #32" << std::endl << std::endl;
    ROL::ZOO::getHS32<RealT> HS32;
    obj  = HS32.getObjective();
    sol  = HS32.getInitialGuess();
    econ = HS32.getEqualityConstraint();
    emul = HS32.getEqualityMultiplier();
    icon = HS32.getInequalityConstraint();
    imul = HS32.getInequalityMultiplier();
    ibnd = HS32.getSlackBoundConstraint();

    algo = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algo->run(*sol,*obj,*icon,*imul,*ibnd,*econ,*emul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(-0.5));
    e2 = (data[1]-static_cast<RealT>(-0.5));
    e3 = (data[2]-static_cast<RealT>(2));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #32" << std::endl << std::endl;
    ROL::ZOO::getHS32<RealT> HS32a;
    obj  = HS32a.getObjective();
    bnd  = HS32a.getBoundConstraint();
    sol  = HS32a.getInitialGuess();
    econ = HS32a.getEqualityConstraint();
    emul = HS32a.getEqualityMultiplier();
    icon = HS32a.getInequalityConstraint();
    imul = HS32a.getInequalityMultiplier();
    ibnd = HS32a.getSlackBoundConstraint();

    algo = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algo->run(*sol,*obj,*bnd,*icon,*imul,*ibnd,*econ,*emul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0));
    e2 = (data[1]-static_cast<RealT>(0));
    e3 = (data[2]-static_cast<RealT>(1));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #63" << std::endl << std::endl;
    ROL::ZOO::getHS63<RealT> HS63;
    obj  = HS63.getObjective();
    sol  = HS63.getInitialGuess();
    bnd  = HS63.getBoundConstraint();
    econ = ROL::makePtr<ROL::ZOO::Constraint_HS63b<RealT>>();
    emul = ROL::makePtr<ROL::StdVector<RealT>>(1);
    icon = ROL::makePtr<ROL::ZOO::Constraint_HS63a<RealT>>();
    imul = ROL::makePtr<ROL::StdVector<RealT>>(1);

    algo = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algo->run(*sol,*obj,*bnd,*econ,*emul,*icon,*imul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(3.512118414));
    e2 = (data[1]-static_cast<RealT>(0.2169881741));
    e3 = (data[2]-static_cast<RealT>(3.552174034));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > 1e3*tol ? 1 : 0);
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


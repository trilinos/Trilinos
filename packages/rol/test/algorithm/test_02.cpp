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

#include "ROL_TypeU_TrustRegionAlgorithm.hpp"
#include "ROL_TypeB_LinMoreAlgorithm.hpp"
#include "ROL_TypeE_AugmentedLagrangianAlgorithm.hpp"
#include "ROL_TypeG_AugmentedLagrangianAlgorithm.hpp"
#include "ROL_Rosenbrock.hpp"
#include "ROL_HS3.hpp"
#include "ROL_HS9.hpp"
#include "ROL_HS14.hpp"
#include "ROL_HS21.hpp"
#include "ROL_HS32.hpp"
#include "ROL_HS41.hpp"
#include "ROL_HS42.hpp"
#include "ROL_HS63.hpp"

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
    list.sublist("Step").set("Type","Trust Region");
    list.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
    list.sublist("Status Test").set("Gradient Tolerance",1e-8);
    list.sublist("Status Test").set("Constraint Tolerance",1e-8);
    list.sublist("Status Test").set("Step Tolerance",1e-12);
    list.sublist("Status Test").set("Iteration Limit", 250);
    list.sublist("General").set("Output Level",iprint);

    ROL::Ptr<ROL::Vector<RealT>>     x, sol, emul, imul, el, il;
    ROL::Ptr<ROL::Objective<RealT>>  obj;
    ROL::Ptr<ROL::Constraint<RealT>> econ, icon;
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd, ibnd;
    ROL::Ptr<ROL::Problem<RealT>> problem;
    ROL::Ptr<ROL::TypeU::Algorithm<RealT>> algoU;
    ROL::Ptr<ROL::TypeB::Algorithm<RealT>> algoB;
    ROL::Ptr<ROL::TypeE::Algorithm<RealT>> algoE;
    ROL::Ptr<ROL::TypeG::Algorithm<RealT>> algoG;
    std::vector<RealT> data;
    RealT e1, e2, e3, e4, err;

    *outStream << std::endl << std::endl << "Type U: Rosenbrock" << std::endl << std::endl;
    ROL::ZOO::getRosenbrock<RealT> Ros;
    obj  = Ros.getObjective();
    sol  = Ros.getInitialGuess();
    x    = sol->clone();

    x->set(*sol);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    ROL::Ptr<ROL::Vector<RealT>> x0 = x->clone(); x->randomize(static_cast<RealT>(-2),static_cast<RealT>(2));
    problem->check(true,*outStream,x0,static_cast<RealT>(0.1));
    problem->finalize(false,true,*outStream);
    algoU = ROL::makePtr<ROL::TypeU::TrustRegionAlgorithm<RealT>>(list);
    algoU->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    e1 = (data[0]-static_cast<RealT>(1));
    e2 = (data[1]-static_cast<RealT>(1));
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << std::endl << "Type U: Hock and Schittkowski Problem #9" << std::endl << std::endl;
    ROL::ZOO::getHS9<RealT> HS9;
    obj  = HS9.getObjective();
    sol  = HS9.getInitialGuess();
    econ = HS9.getEqualityConstraint(); 
    emul = HS9.getEqualityMultiplier();
    x    = sol->clone();
    el   = emul->clone();

    x->set(*sol); el->set(*emul);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addLinearConstraint("linear_econ1",econ,el);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoU = ROL::makePtr<ROL::TypeU::TrustRegionAlgorithm<RealT>>(list);
    algoU->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]+static_cast<RealT>(3))/static_cast<RealT>(12);
    e2 = (data[1]+static_cast<RealT>(4))/static_cast<RealT>(16);
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    err = std::max(std::abs(e1-std::round(e1)),std::abs(e2-std::round(e2)));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag = (err > tol ? 1 : 0);

    *outStream << std::endl << "Type B: Hock and Schittkowski Problem #3" << std::endl << std::endl;
    ROL::ZOO::getHS3<RealT> HS3;
    obj  = HS3.getObjective();
    sol  = HS3.getInitialGuess();
    bnd  = HS3.getBoundConstraint();
    x    = sol->clone();

    x->set(*sol);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addBoundConstraint(bnd);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoB = ROL::makePtr<ROL::TypeB::LinMoreAlgorithm<RealT>>(list);
    algoB->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0));
    e2 = (data[1]-static_cast<RealT>(0));
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Type B: Hock and Schittkowski Problem #21" << std::endl << std::endl;
    ROL::ZOO::getHS21<RealT> HS21;
    obj  = HS21.getObjective();
    sol  = HS21.getInitialGuess();
    bnd  = HS21.getBoundConstraint();
    icon = HS21.getInequalityConstraint();
    imul = HS21.getInequalityMultiplier();
    ibnd = HS21.getSlackBoundConstraint();
    x    = sol->clone();
    il   = imul->clone();

    x->set(*sol); il->set(*imul);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addBoundConstraint(bnd);
    problem->addLinearConstraint("linear_icon1",icon,il,ibnd);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoB = ROL::makePtr<ROL::TypeB::LinMoreAlgorithm<RealT>>(list);
    algoB->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    e1 = (data[0]-static_cast<RealT>(2));
    e2 = (data[1]-static_cast<RealT>(0));
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Type B: Hock and Schittkowski Problem #41" << std::endl << std::endl;
    ROL::ZOO::getHS41<RealT> HS41;
    obj  = HS41.getObjective();
    sol  = HS41.getInitialGuess();
    bnd  = HS41.getBoundConstraint();
    econ = HS41.getEqualityConstraint();
    emul = HS41.getEqualityMultiplier();
    x    = sol->clone();
    el   = emul->clone();

    x->set(*sol); el->set(*emul);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addBoundConstraint(bnd);
    problem->addLinearConstraint("linear_econ1",econ,el);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoB = ROL::makePtr<ROL::TypeB::LinMoreAlgorithm<RealT>>(list);
    algoB->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3] << std::endl;
    e1 = (data[0]-static_cast<RealT>(2.0/3.0));
    e2 = (data[1]-static_cast<RealT>(1.0/3.0));
    e3 = (data[2]-static_cast<RealT>(1.0/3.0));
    e4 = (data[3]-static_cast<RealT>(2.0));
    err = std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << std::endl << "Type E: Hock and Schittkowski Problem #42" << std::endl << std::endl;
    ROL::ZOO::getHS42<RealT> HS42;
    obj  = HS42.getObjective();
    sol  = HS42.getInitialGuess();
    econ = ROL::makePtr<ROL::ZOO::Constraint_HS42b<RealT>>();
    emul = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
    icon = ROL::makePtr<ROL::ZOO::Constraint_HS42a<RealT>>();
    imul = ROL::makePtr<ROL::StdVector<RealT>>(1,0.0);
    x    = sol->clone();
    el   = emul->clone();
    il   = imul->clone();

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addConstraint("econ1",econ,el);
    problem->addLinearConstraint("linear_econ1",icon,il);
    problem->check(true,*outStream);
    problem->finalize(true,true,*outStream);
    algoE = ROL::makePtr<ROL::TypeE::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoE->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(2));
    e2 = (data[1]-static_cast<RealT>(2));
    e3 = (data[2]-static_cast<RealT>(0.6*std::sqrt(2)));
    e4 = (data[3]-static_cast<RealT>(0.8*std::sqrt(2)));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3] << std::endl;
    err = std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag = (err > tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoE = ROL::makePtr<ROL::TypeE::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoE->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(2));
    e2 = (data[1]-static_cast<RealT>(2));
    e3 = (data[2]-static_cast<RealT>(0.6*std::sqrt(2)));
    e4 = (data[3]-static_cast<RealT>(0.8*std::sqrt(2)));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3] << std::endl;
    err = std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag = (err > tol ? 1 : 0);

    *outStream << std::endl << std::endl << "Type G: Hock and Schittkowski Problem #14" << std::endl << std::endl;
    ROL::ZOO::getHS14<RealT> HS14;
    obj  = HS14.getObjective();
    sol  = HS14.getInitialGuess();
    econ = HS14.getEqualityConstraint();
    emul = HS14.getEqualityMultiplier();
    icon = HS14.getInequalityConstraint();
    imul = HS14.getInequalityMultiplier();
    ibnd = HS14.getSlackBoundConstraint();
    x    = sol->clone();
    el   = emul->clone();
    il   = imul->clone();

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addConstraint("econ1",econ,el);
    problem->addConstraint("icon1",icon,il,ibnd);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0.5 *(std::sqrt(7.0)-1.0)));
    e2 = (data[1]-static_cast<RealT>(0.25*(std::sqrt(7.0)+1.0)));
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->removeConstraint("econ1");
    problem->addLinearConstraint("linear_econ1",econ,el);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0.5 *(std::sqrt(7.0)-1.0)));
    e2 = (data[1]-static_cast<RealT>(0.25*(std::sqrt(7.0)+1.0)));
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->check(true,*outStream);
    problem->finalize(true,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0.5 *(std::sqrt(7.0)-1.0)));
    e2 = (data[1]-static_cast<RealT>(0.25*(std::sqrt(7.0)+1.0)));
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << std::endl << "Type G: Hock and Schittkowski Problem #32" << std::endl << std::endl;
    ROL::ZOO::getHS32<RealT> HS32;
    obj  = HS32.getObjective();
    sol  = HS32.getInitialGuess();
    bnd  = HS32.getBoundConstraint();
    econ = HS32.getEqualityConstraint();
    emul = HS32.getEqualityMultiplier();
    icon = HS32.getInequalityConstraint();
    imul = HS32.getInequalityMultiplier();
    ibnd = HS32.getSlackBoundConstraint();
    x    = sol->clone();
    el   = emul->clone();
    il   = imul->clone();

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addBoundConstraint(bnd);
    problem->addConstraint("econ1",econ,el);
    problem->addConstraint("icon1",icon,il,ibnd);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0));
    e2 = (data[1]-static_cast<RealT>(0));
    e3 = (data[2]-static_cast<RealT>(1));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->removeConstraint("econ1");
    problem->addLinearConstraint("linear_econ1",econ,el);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0));
    e2 = (data[1]-static_cast<RealT>(0));
    e3 = (data[2]-static_cast<RealT>(1));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->check(true,*outStream);
    problem->finalize(true,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(0));
    e2 = (data[1]-static_cast<RealT>(0));
    e3 = (data[2]-static_cast<RealT>(1));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->removeBoundConstraint();
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(-0.5));
    e2 = (data[1]-static_cast<RealT>(-0.5));
    e3 = (data[2]-static_cast<RealT>(2));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->check(true,*outStream);
    problem->finalize(true,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(-0.5));
    e2 = (data[1]-static_cast<RealT>(-0.5));
    e3 = (data[2]-static_cast<RealT>(2));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << std::endl << "Type G: Hock and Schittkowski Problem #63" << std::endl << std::endl;
    ROL::ZOO::getHS63<RealT> HS63;
    obj  = HS63.getObjective();
    sol  = HS63.getInitialGuess();
    bnd  = HS63.getBoundConstraint();
    econ = ROL::makePtr<ROL::ZOO::Constraint_HS63b<RealT>>();
    emul = ROL::makePtr<ROL::StdVector<RealT>>(1);
    icon = ROL::makePtr<ROL::ZOO::Constraint_HS63a<RealT>>();
    imul = ROL::makePtr<ROL::StdVector<RealT>>(1);
    x    = sol->clone();
    el   = emul->clone();
    il   = imul->clone();

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addBoundConstraint(bnd);
    problem->addConstraint("econ1",econ,el);
    problem->addConstraint("econ2",icon,il);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(3.512118414));
    e2 = (data[1]-static_cast<RealT>(0.2169881741));
    e3 = (data[2]-static_cast<RealT>(3.552174034));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > 1e3*tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->removeConstraint("econ2");
    problem->addLinearConstraint("linear_econ1",icon,il);
    problem->check(true,*outStream);
    problem->finalize(false,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    e1 = (data[0]-static_cast<RealT>(3.512118414));
    e2 = (data[1]-static_cast<RealT>(0.2169881741));
    e3 = (data[2]-static_cast<RealT>(3.552174034));
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > 1e3*tol ? 1 : 0);

    x->set(*sol); el->set(*emul); il->set(*imul);
    problem->edit();
    problem->check(true,*outStream);
    problem->finalize(true,true,*outStream);
    algoG = ROL::makePtr<ROL::TypeG::AugmentedLagrangianAlgorithm<RealT>>(list);
    algoG->run(*problem,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
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


// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_14.cpp
    \brief Validate ExplicitLinearConstraint infrastructure.

*/

#include "ROL_TypeU_AlgorithmFactory.hpp"
#include "ROL_HS9.hpp"
#include "ROL_HS28.hpp"
#include "ROL_HS48.hpp"
#include "ROL_HS49.hpp"
#include "ROL_HS50.hpp"
#include "ROL_HS51.hpp"
#include "ROL_HS52.hpp"

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
    list.sublist("Status Test").set("Gradient Tolerance",1e-12);
    list.sublist("Status Test").set("Constraint Tolerance",1e-12);
    list.sublist("Status Test").set("Step Tolerance",1e-14);
    list.sublist("Step").set("Type","Trust Region");
    list.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
    ROL::Ptr<ROL::Vector<RealT>>      sol, mul, x;
    ROL::Ptr<ROL::Objective<RealT>>   obj;
    ROL::Ptr<ROL::Constraint<RealT>>  con;
    ROL::Ptr<ROL::TypeU::Algorithm<RealT>> algo;
    std::vector<RealT> data;
    RealT e1(0), e2(0), e3(0), e4(0), e5(0), err(0);

    *outStream << "Hock and Schittkowski Problem #9" << std::endl << std::endl;
    ROL::ZOO::getHS9<RealT> HS9;
    obj = HS9.getObjective();
    sol = HS9.getInitialGuess();
    con = HS9.getEqualityConstraint(); 
    mul = HS9.getEqualityMultiplier();
    x   = sol->clone(); x->set(*sol);

    algo = ROL::TypeU::AlgorithmFactory<RealT>(list);
    algo->run(*x,*obj,*con,*mul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]+static_cast<RealT>(3))/static_cast<RealT>(12);
    e2 = (data[1]+static_cast<RealT>(4))/static_cast<RealT>(16);
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    err = std::max(std::abs(e1-std::round(e1)),std::abs(e2-std::round(e2)));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag = (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #28" << std::endl << std::endl;
    ROL::ZOO::getHS28<RealT> HS28;
    obj = HS28.getObjective();
    sol = HS28.getInitialGuess();
    con = HS28.getEqualityConstraint(); 
    mul = HS28.getEqualityMultiplier();
    x   = sol->clone(); x->set(*sol);

    algo = ROL::TypeU::AlgorithmFactory<RealT>(list);
    algo->run(*x,*obj,*con,*mul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(0.5));
    e2 = (data[1]+static_cast<RealT>(0.5));
    e3 = (data[2]-static_cast<RealT>(0.5));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1] << "  x3 = " << data[2] << std::endl;
    err = std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #48" << std::endl << std::endl;
    ROL::ZOO::getHS48<RealT> HS48;
    obj = HS48.getObjective();
    sol = HS48.getInitialGuess();
    con = HS48.getEqualityConstraint(); 
    mul = HS48.getEqualityMultiplier();
    x   = sol->clone(); x->set(*sol);

    algo = ROL::TypeU::AlgorithmFactory<RealT>(list);
    algo->run(*x,*obj,*con,*mul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(1));
    e2 = (data[1]-static_cast<RealT>(1));
    e3 = (data[2]-static_cast<RealT>(1));
    e4 = (data[3]-static_cast<RealT>(1));
    e5 = (data[4]-static_cast<RealT>(1));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << std::endl;
    err = std::max(std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4)),std::abs(e5));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #49" << std::endl << std::endl;
    ROL::ZOO::getHS49<RealT> HS49;
    obj = HS49.getObjective();
    sol = HS49.getInitialGuess();
    con = HS49.getEqualityConstraint(); 
    mul = HS49.getEqualityMultiplier();
    x   = sol->clone(); x->set(*sol);

    algo = ROL::TypeU::AlgorithmFactory<RealT>(list);
    algo->run(*x,*obj,*con,*mul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(1));
    e2 = (data[1]-static_cast<RealT>(1));
    e3 = (data[2]-static_cast<RealT>(1));
    e4 = (data[3]-static_cast<RealT>(1));
    e5 = (data[4]-static_cast<RealT>(1));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << std::endl;
    err = std::max(std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4)),std::abs(e5));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > 1e4*tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #50" << std::endl << std::endl;
    ROL::ZOO::getHS50<RealT> HS50;
    obj = HS50.getObjective();
    sol = HS50.getInitialGuess();
    con = HS50.getEqualityConstraint(); 
    mul = HS50.getEqualityMultiplier();
    x   = sol->clone(); x->set(*sol);

    algo = ROL::TypeU::AlgorithmFactory<RealT>(list);
    algo->run(*x,*obj,*con,*mul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(1));
    e2 = (data[1]-static_cast<RealT>(1));
    e3 = (data[2]-static_cast<RealT>(1));
    e4 = (data[3]-static_cast<RealT>(1));
    e5 = (data[4]-static_cast<RealT>(1));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << std::endl;
    err = std::max(std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4)),std::abs(e5));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #51" << std::endl << std::endl;
    ROL::ZOO::getHS51<RealT> HS51;
    obj = HS51.getObjective();
    sol = HS51.getInitialGuess();
    con = HS51.getEqualityConstraint(); 
    mul = HS51.getEqualityMultiplier();
    x   = sol->clone(); x->set(*sol);

    algo = ROL::TypeU::AlgorithmFactory<RealT>(list);
    algo->run(*x,*obj,*con,*mul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(1));
    e2 = (data[1]-static_cast<RealT>(1));
    e3 = (data[2]-static_cast<RealT>(1));
    e4 = (data[3]-static_cast<RealT>(1));
    e5 = (data[4]-static_cast<RealT>(1));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << std::endl;
    err = std::max(std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4)),std::abs(e5));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #52" << std::endl << std::endl;
    ROL::ZOO::getHS52<RealT> HS52;
    obj = HS52.getObjective();
    sol = HS52.getInitialGuess();
    con = HS52.getEqualityConstraint(); 
    mul = HS52.getEqualityMultiplier();
    x   = sol->clone(); x->set(*sol);

    algo = ROL::TypeU::AlgorithmFactory<RealT>(list);
    algo->run(*x,*obj,*con,*mul,*outStream);

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    e1 = (data[0]-static_cast<RealT>(-33.0/349.0));
    e2 = (data[1]-static_cast<RealT>(11.0/349.0));
    e3 = (data[2]-static_cast<RealT>(180.0/349.0));
    e4 = (data[3]-static_cast<RealT>(-158.0/349.0));
    e5 = (data[4]-static_cast<RealT>(11.0/349.0));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << std::endl;
    err = std::max(std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4)),std::abs(e5));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);
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


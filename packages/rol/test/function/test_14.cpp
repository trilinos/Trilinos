// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_15.cpp
    \brief Validate polyhedral projection infrastructure.

*/

#include "ROL_HS41.hpp"
#include "ROL_HS53.hpp"
#include "ROL_HS55.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_PolyhedralProjectionFactory.hpp"

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
    RealT cnorm(0);
    ROL::Ptr<ROL::Vector<RealT>>     sol, mul, x, lam, l, u, c;
    ROL::Ptr<ROL::Objective<RealT>>  obj;
    ROL::Ptr<ROL::Constraint<RealT>> con;
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd;
    ROL::Ptr<ROL::PolyhedralProjection<RealT>> proj;
    ROL::ParameterList list;
    list.sublist("General").set("Output Level",2);
    std::vector<RealT> data;

    *outStream << std::endl << "Hock and Schittkowski Problem #41" << std::endl << std::endl;
    ROL::ZOO::getHS41<RealT> HS41;
    obj = HS41.getObjective();
    sol = HS41.getInitialGuess();
    con = HS41.getEqualityConstraint();
    mul = HS41.getEqualityMultiplier();
    bnd = HS41.getBoundConstraint();

    lam = mul->clone(); lam->set(*mul);
    x   = sol->clone(); x->set(*sol);
    l   = sol->clone(); l->zero();
    u   = sol->clone(); u->setScalar(static_cast<RealT>(1));
    c   = mul->dual().clone();

    list.sublist("General").sublist("Polyhedral Projection").set("Type","Dai-Fletcher");
    proj = ROL::PolyhedralProjectionFactory<RealT>(*sol,sol->dual(),bnd,con,*lam,*c,list);
    proj->project(*x,*outStream);

    con->value(*c,*x,tol);
    cnorm = c->norm();

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Initial:    x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(lam)->getVector();
    *outStream << "  Multiplier: l1 = " << data[0] << std::endl;

    *outStream << std::endl;
    *outStream << "  is equality feasible = " << (cnorm<=tol)        << std::endl
               << "  are bounds feasible  = " << bnd->isFeasible(*x) << std::endl;

    errorFlag += !bnd->isFeasible(*x);
    errorFlag += (cnorm > tol);

    *outStream << std::endl << "Hock and Schittkowski Problem #41" << std::endl << std::endl;
    ROL::ZOO::getHS41<RealT> HS41a;
    obj = HS41a.getObjective();
    sol = HS41a.getInitialGuess();
    con = HS41a.getEqualityConstraint();
    mul = HS41a.getEqualityMultiplier();
    bnd = HS41a.getBoundConstraint();

    lam = mul->clone(); lam->set(*mul);
    x   = sol->clone(); x->set(*sol);
    l   = sol->clone(); l->zero();
    u   = sol->clone(); u->setScalar(static_cast<RealT>(1));
    c   = mul->dual().clone();

    list.sublist("General").sublist("Polyhedral Projection").set("Type","Ridders");
    proj = ROL::PolyhedralProjectionFactory<RealT>(*sol,sol->dual(),bnd,con,*lam,*c,list);
    proj->project(*x,*outStream);

    con->value(*c,*x,tol);
    cnorm = c->norm();

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Initial:    x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(lam)->getVector();
    *outStream << "  Multiplier: l1 = " << data[0] << std::endl;

    *outStream << std::endl;
    *outStream << "  is equality feasible = " << (cnorm<=tol)        << std::endl
               << "  are bounds feasible  = " << bnd->isFeasible(*x) << std::endl;

    errorFlag += !bnd->isFeasible(*x);
    errorFlag += (cnorm > tol);

    *outStream << std::endl << "Hock and Schittkowski Problem #53" << std::endl << std::endl;
    ROL::ZOO::getHS53<RealT> HS53;
    obj = HS53.getObjective();
    sol = HS53.getInitialGuess();
    con = HS53.getEqualityConstraint();
    mul = HS53.getEqualityMultiplier();
    bnd = HS53.getBoundConstraint();

    lam = mul->clone(); lam->set(*mul);
    x   = sol->clone(); x->set(*sol);
    l   = sol->clone(); l->zero();
    u   = sol->clone(); u->setScalar(static_cast<RealT>(1));
    c   = mul->dual().clone();

    list.sublist("General").sublist("Polyhedral Projection").set("Type","Dykstra");
    proj = ROL::PolyhedralProjectionFactory<RealT>(*sol,sol->dual(),bnd,con,*lam,*c,list);
    proj->project(*x,*outStream);

    con->value(*c,*x,tol);
    cnorm = c->norm();

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Initial:    x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(lam)->getVector();
    *outStream << "  Multiplier: l1 = " << data[0] << "  l2 = " << data[1]
               << "  l3 = " << data[2] << std::endl;

    *outStream << std::endl;
    *outStream << "  is equality feasible = " << (cnorm<=tol)        << std::endl
               << "  are bounds feasible  = " << bnd->isFeasible(*x) << std::endl;

    errorFlag += !bnd->isFeasible(*x);
    errorFlag += (cnorm > tol);

    *outStream << std::endl << "Hock and Schittkowski Problem #55" << std::endl << std::endl;
    ROL::ZOO::getHS55<RealT> HS55;
    obj = HS55.getObjective();
    sol = HS55.getInitialGuess();
    con = HS55.getEqualityConstraint();
    mul = HS55.getEqualityMultiplier();
    bnd = HS55.getBoundConstraint();

    //ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
    //ROL::Ptr<ROL::Vector<RealT>> xt;
    //std::vector<ROL::Ptr<ROL::Vector<RealT>>> xv;
    //HS55.get(problem,xt,xv);
    //problem->check(*outStream);

    lam = mul->clone(); lam->set(*mul);
    x   = sol->clone(); x->set(*sol);
    l   = sol->clone(); l->zero();
    u   = sol->clone(); u->setScalar(static_cast<RealT>(1));
    c   = mul->dual().clone();

    list.sublist("General").sublist("Polyhedral Projection").set("Type","Semismooth Newton");
    proj = ROL::PolyhedralProjectionFactory<RealT>(*sol,sol->dual(),bnd,con,*lam,*c,list);
    proj->project(*x,*outStream);

    con->value(*c,*x,tol);
    cnorm = c->norm();

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(sol)->getVector();
    *outStream << "  Initial:    x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << "  x6 = " << data[5] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(x)->getVector();
    *outStream << "  Result:     x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3]
               << "  x5 = " << data[4] << "  x6 = " << data[5] << std::endl;
    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(lam)->getVector();
    *outStream << "  Multiplier: l1 = " << data[0] << "  l2 = " << data[1]
               << "  l3 = " << data[2] << "  l4 = " << data[3]
               << "  l5 = " << data[4] << "  l6 = " << data[5] << std::endl;

    *outStream << std::endl;
    *outStream << "  is equality feasible = " << (cnorm<=tol)        << std::endl
               << "  are bounds feasible  = " << bnd->isFeasible(*x) << std::endl;
    *outStream << std::endl;

    errorFlag += !bnd->isFeasible(*x);
    errorFlag += (cnorm > tol);
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


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

/*! \file  test_14.cpp
    \brief Validate ExplicitLinearConstraint infrastructure.

*/

#include "ROL_OptimizationSolver.hpp"
#include "ROL_ExplicitLinearConstraint.hpp"
#include "ROL_HS9.hpp"
#include "ROL_HS14.hpp"
#include "ROL_HS28.hpp"
#include "ROL_HS42.hpp"
#include "ROL_HS48.hpp"
#include "ROL_HS49.hpp"
#include "ROL_HS50.hpp"
#include "ROL_HS51.hpp"
#include "ROL_HS52.hpp"

#include "ROL_BinaryConstraint.hpp"
#include "ROL_DiagonalOperator.hpp"
#include "ROL_QuadraticObjective.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"

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
    RealT tol = 1e4*std::sqrt(ROL::ROL_EPSILON<RealT>());
    ROL::ParameterList list;
    list.sublist("Status Test").set("Gradient Tolerance",1e-12);
    list.sublist("Status Test").set("Constraint Tolerance",1e-12);
    list.sublist("Status Test").set("Step Tolerance",1e-14);
    list.sublist("Step").set("Type","Trust Region");
    list.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
    list.sublist("Step").sublist("Augmented Lagrangian").set("Print Intermediate Optimization History",true);
    ROL::Ptr<ROL::Vector<RealT>>     sol, mul, s, ds, Ps, imul, c;
    ROL::Ptr<ROL::Objective<RealT>>  obj, tobj;
    ROL::Ptr<ROL::Constraint<RealT>> con, icon, tcon;
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd;
    ROL::Ptr<ROL::ExplicitLinearConstraint<RealT>> elc;
    ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
    ROL::Ptr<ROL::OptimizationSolver<RealT>>  solver;
    std::vector<RealT> data;
    RealT e1(0), e2(0), e3(0), e4(0), e5(0), err(0);

    *outStream << "Hock and Schittkowski Problem #9" << std::endl << std::endl;
    ROL::ZOO::getHS9<RealT> HS9;
    obj = HS9.getObjective();
    sol = HS9.getInitialGuess();
    con = HS9.getEqualityConstraint(); 
    mul = HS9.getEqualityMultiplier();

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,sol,mul);
    tobj = elc->getTransformedObjective();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
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

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,sol,mul);
    tobj = elc->getTransformedObjective();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
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

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,sol,mul);
    tobj = elc->getTransformedObjective();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
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

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,sol,mul);
    tobj = elc->getTransformedObjective();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
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

    *outStream << std::endl << "Hock and Schittkowski Problem #50" << std::endl << std::endl;
    ROL::ZOO::getHS50<RealT> HS50;
    obj = HS50.getObjective();
    sol = HS50.getInitialGuess();
    con = HS50.getEqualityConstraint(); 
    mul = HS50.getEqualityMultiplier();

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,sol,mul);
    tobj = elc->getTransformedObjective();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
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

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,sol,mul);
    tobj = elc->getTransformedObjective();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
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

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,sol,mul);
    tobj = elc->getTransformedObjective();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
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

    *outStream << std::endl << "Hock and Schittkowski Problem #14" << std::endl << std::endl;
    list.sublist("Step").set("Type","Moreau-Yosida Penalty");
    ROL::ZOO::getHS14<RealT> HS14;
    obj  = HS14.getObjective();
    sol  = HS14.getInitialGuess();
    con  = HS14.getEqualityConstraint(); 
    mul  = HS14.getEqualityMultiplier();
    icon = HS14.getInequalityConstraint();
    imul = HS14.getInequalityMultiplier();
    ibnd = HS14.getSlackBoundConstraint();

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,icon,sol,mul);
    tobj = elc->getTransformedObjective();
    tcon = elc->getTransformedConstraint();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();
    c    = imul->clone(); c->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);
    tcon->checkApplyJacobian(*s,*ds,*c,true,*outStream);
    tcon->checkApplyAdjointJacobian(*s,*c,*c,*ds,true,*outStream);
    tcon->checkApplyAdjointHessian(*s,*c,*ds,*s,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s,tcon,imul,ibnd);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
    e1 = (data[0]-static_cast<RealT>(0.5 *(std::sqrt(7)-1.0)));
    e2 = (data[1]-static_cast<RealT>(0.25*(std::sqrt(7)+1.0)));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1] << std::endl;
    err = std::max(std::abs(e1),std::abs(e2));
    *outStream << "  Max-Error = " << err << std::endl;
    errorFlag += (err > tol ? 1 : 0);

    *outStream << std::endl << "Hock and Schittkowski Problem #42" << std::endl << std::endl;
    list.sublist("Step").set("Type","Composite Step");
    ROL::ZOO::getHS42<RealT> HS42;
    obj  = HS42.getObjective();
    sol  = HS42.getInitialGuess();
    con  = ROL::staticPtrCast<ROL::Constraint_Partitioned<RealT>>(HS42.getEqualityConstraint())->get(0);
    mul  = ROL::staticPtrCast<ROL::PartitionedVector<RealT>>(HS42.getEqualityMultiplier())->get(0);
    icon = ROL::staticPtrCast<ROL::Constraint_Partitioned<RealT>>(HS42.getEqualityConstraint())->get(1);
    imul = ROL::staticPtrCast<ROL::PartitionedVector<RealT>>(HS42.getEqualityMultiplier())->get(1);

    elc  = ROL::makePtr<ROL::ExplicitLinearConstraint<RealT>>(con,obj,icon,sol,mul);
    tobj = elc->getTransformedObjective();
    tcon = elc->getTransformedConstraint();
    s    = sol->clone();  s->randomize();
    ds   = sol->clone(); ds->randomize();
    c    = imul->clone(); c->randomize();

    obj->checkGradient(*s,*ds,true,*outStream);
    obj->checkHessVec(*s,*ds,true,*outStream);
    tobj->checkGradient(*s,*ds,true,*outStream);
    tobj->checkHessVec(*s,*ds,true,*outStream);
    tcon->checkApplyJacobian(*s,*ds,*c,true,*outStream);
    tcon->checkApplyAdjointJacobian(*s,*c,*c,*ds,true,*outStream);
    tcon->checkApplyAdjointHessian(*s,*c,*ds,*s,true,*outStream);

    s->zero();
    problem = ROL::makePtr<ROL::OptimizationProblem<RealT>>(tobj,s,tcon,imul);
    solver  = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*problem,list);
    solver->solve(*outStream);

    Ps = s->clone();
    elc->project(Ps,s);
    Ps->plus(*elc->getFeasibleVector());

    data = *ROL::staticPtrCast<ROL::StdVector<RealT>>(Ps)->getVector();
    e1 = (data[0]-static_cast<RealT>(2.0));
    e2 = (data[1]-static_cast<RealT>(2.0));
    e3 = (data[2]-static_cast<RealT>(0.6*std::sqrt(2.0)));
    e4 = (data[3]-static_cast<RealT>(0.8*std::sqrt(2.0)));
    *outStream << "  x1 = " << data[0] << "  x2 = " << data[1]
               << "  x3 = " << data[2] << "  x4 = " << data[3] << std::endl;
    err = std::max(std::max(std::max(std::abs(e1),std::abs(e2)),std::abs(e3)),std::abs(e4));
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


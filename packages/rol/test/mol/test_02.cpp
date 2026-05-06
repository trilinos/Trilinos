// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_02.cpp
    \brief Test ParetoSampler interface.
*/

#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_MultiObjectiveFactory.hpp"
#include "ROL_ParetoSampler.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "ROL_TeuchosBatchManager.hpp"

#include <iostream>

template<typename Real>
class Objective0_Test01 : public ROL::StdObjective<Real> {
public:
  Objective0_Test01() {}
  Real value(const std::vector<Real>& x, Real& tol) {
    return x[0]*x[0]+x[1]*x[1];
  }
  void gradient(std::vector<Real>& g, const std::vector<Real>& x, Real& tol) {
    g[0] = static_cast<Real>(2)*x[0];
    g[1] = static_cast<Real>(2)*x[1];
  }
  void hessVec(std::vector<Real>& hv, const std::vector<Real>& v, const std::vector<Real>& x, Real& tol) {
    hv[0] = static_cast<Real>(2)*v[0];
    hv[1] = static_cast<Real>(2)*v[1];
  }
};

template<typename Real>
class Objective1_Test01 : public ROL::StdObjective<Real> {
public:
  Objective1_Test01() {}
  Real value(const std::vector<Real>& x, Real& tol) {
    const Real diff = x[0]-static_cast<Real>(1);
    return diff*diff+x[1]*x[1];
  }
  void gradient(std::vector<Real>& g, const std::vector<Real>& x, Real& tol) {
    const Real diff = x[0]-static_cast<Real>(1);
    g[0] = static_cast<Real>(2)*diff;
    g[1] = static_cast<Real>(2)*x[1];
  }
  void hessVec(std::vector<Real>& hv, const std::vector<Real>& v, const std::vector<Real>& x, Real& tol) {
    hv[0] = static_cast<Real>(2)*v[0];
    hv[1] = static_cast<Real>(2)*v[1];
  }
};

template<typename Real>
class Constraint_Test01 : public ROL::StdConstraint<Real> {
public:
  Constraint_Test01() {}
  void value(std::vector<Real>& c, const std::vector<Real>& x, Real& tol) {
    c[0] = static_cast<Real>(2)*(x[0]-static_cast<Real>(0.1))*(x[0]-static_cast<Real>(0.9))/static_cast<Real>(0.18);
    c[1] = static_cast<Real>(-20)*(x[0]-static_cast<Real>(0.4))*(x[0]-static_cast<Real>(0.6))/static_cast<Real>(4.8);
  }
  void applyJacobian(std::vector<Real>& jv, const std::vector<Real>& v, const std::vector<Real>& x, Real& tol) {
    jv[0] = static_cast<Real>(2)*((x[0]-static_cast<Real>(0.1))+(x[0]-static_cast<Real>(0.9)))/static_cast<Real>(0.18) * v[0];
    jv[1] = static_cast<Real>(-20)*((x[0]-static_cast<Real>(0.4))+(x[0]-static_cast<Real>(0.6)))/static_cast<Real>(4.8) * v[0];
  }
  void applyAdjointJacobian(std::vector<Real>& ajv, const std::vector<Real>& v, const std::vector<Real>& x, Real& tol) {
    ajv[0] = static_cast<Real>(2)*((x[0]-static_cast<Real>(0.1))+(x[0]-static_cast<Real>(0.9)))/static_cast<Real>(0.18) * v[0]
           + static_cast<Real>(-20)*((x[0]-static_cast<Real>(0.4))+(x[0]-static_cast<Real>(0.6)))/static_cast<Real>(4.8) * v[1];
    ajv[1] = static_cast<Real>(0);
  }
  void applyAdjointHessian(std::vector<Real>& ahuv, const std::vector<Real>& u, const std::vector<Real>& v, const std::vector<Real>& x, Real& tol) {
    ahuv[0] = static_cast<Real>(4)/static_cast<Real>(0.18) * u[0]*v[0]
            + static_cast<Real>(-40)/static_cast<Real>(4.8) * u[1]*v[0];
    ahuv[1] = static_cast<Real>(0);
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  ROL::GlobalMPISession mpiSession(&argc, &argv);
  auto comm = ROL::toPtr(Teuchos::DefaultComm<int>::getComm());

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && comm->getRank()==0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  RealT errtol = static_cast<RealT>(1e2)*std::sqrt(ROL::ROL_EPSILON<RealT>());

  // *** Test body.

  try {

    // Create multi-objective problem
    const int xdim = 2, cdim = 2; 
    std::vector<RealT> lb(xdim,-2), ub(xdim,2), iub(cdim,0);
    auto x    = ROL::makePtr<ROL::StdVector<RealT>>(xdim);
    auto obj0 = ROL::makePtr<Objective0_Test01<RealT>>();
    auto obj1 = ROL::makePtr<Objective1_Test01<RealT>>();
    auto bnd  = ROL::makePtr<ROL::StdBoundConstraint<RealT>>(lb,ub);
    auto icon = ROL::makePtr<Constraint_Test01<RealT>>();
    auto imul = ROL::makePtr<ROL::StdVector<RealT>>(cdim);
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd = ROL::makePtr<ROL::StdBoundConstraint<RealT>>(iub,false);

    auto prob0 = ROL::makePtr<ROL::Problem<RealT>>(obj0,x);
    prob0->addBoundConstraint(bnd);
    prob0->addConstraint("Quad",icon,imul,ibnd);
    prob0->finalize(false,true,*outStream);
    prob0->check(true,*outStream);

    auto prob1 = ROL::makePtr<ROL::Problem<RealT>>(obj1,x);
    prob1->addBoundConstraint(bnd);
    prob1->addConstraint("Quad",icon,imul,ibnd);
    prob1->finalize(false,true,*outStream);
    prob1->check(true,*outStream);

    auto mof = ROL::makePtr<ROL::MultiObjectiveFactory<RealT>>(x);
    mof->addObjective("Obj0",obj0);
    mof->addObjective("Obj1",obj1);
    mof->addBoundConstraint(bnd);
    mof->addConstraint("Quad",icon,imul,ibnd);

    ROL::ParameterList parlist;
    parlist.sublist("General").set("Output Level", 1);
    parlist.sublist("Step").sublist("Augmented Lagrangian").set("Initial Penalty Parameter",1.0);
    parlist.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Initial Penalty Parameter",false);
    parlist.sublist("Step").sublist("Augmented Lagrangian").set("Use Default Problem Scaling",false);
    parlist.sublist("Mult-Objective").set("Scalarization Type","Convex Combination");
    parlist.sublist("Mult-Objective").sublist("Pareto Sampler").set("Number of Points",10);
    parlist.sublist("Mult-Objective").sublist("Pareto Sampler").set("Warm Start",false);

    auto bman = ROL::makePtr<ROL::TeuchosBatchManager<RealT,int>>(comm);
    auto ps = ROL::makePtr<ROL::ParetoSampler<RealT>>(bman);
    ps->run(mof,parlist,*outStream);
    ps->print("output_01.txt");

    // Pareto front is given by f1 in [0.01,0.16] or [0.36,0.81]
    // and f2 = (sqrt(f1)-1)^2.  The optimal solutions are
    // given by x1 in [0.1,0.4] or [0.6,0.9] and x2 = 0.
    auto sol = ps->getParetoData();
    for (const auto& pd : sol) {
      bool feas1 = (pd.values[0] <= static_cast<RealT>(0.16)+errtol && pd.values[0] >= static_cast<RealT>(0.01)-errtol);
      bool feas2 = (pd.values[0] <= static_cast<RealT>(0.81)+errtol && pd.values[0] >= static_cast<RealT>(0.36)-errtol);
      bool paret = (std::abs(pd.values[1]-std::pow(std::sqrt(pd.values[0])-static_cast<RealT>(1),2))<=errtol);
      if ((!feas1&&!feas2)||!paret) errorFlag++;
    }
  }

  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";

  return 0;

}


// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_BINARYDESIGNALGORITHM_HPP
#define ROL_OED_BINARYDESIGNALGORITHM_HPP

#include <iostream>

#include "ROL_Ptr.hpp"
#include "ROL_Problem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_UnaryFunctions.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_SlacklessObjective.hpp"
#include "ROL_OED_Factory.hpp"
#include "ROL_OED_DesignVector.hpp"

namespace ROL::OED {

// General Greedy Algorithm Implementation
template<typename Real>
class BinaryDesignAlgorithm {
private:
  Ptr<Factory<Real>> factory_;
  bool budgetSet_, equality_, eflag_;

  Ptr<Vector<Real>> getDesign(Ptr<Vector<Real>>& x) const {
    if (eflag_) {
      return staticPtrCast<DesignVector<Real>>(staticPtrCast<PartitionedVector<Real>>(x)->get(0))->getVector();
    }
    else {
      if (!equality_)
        return staticPtrCast<PartitionedVector<Real>>(x)->get(0);
      else
        return x;
    }
  }

  bool runBase(const Ptr<Problem<Real>>& problem0,
               ParameterList& list,
               std::ostream& os=std::cout) {
    if (!budgetSet_) {
      throw Exception::NotImplemented(">>> ROL::OED::BinaryDesignAlgorithm::run: Budget constraint not set!");
    }
    std::ios_base::fmtflags osFlags(os.flags());
    const Real zero(0), one(1); //two(2), three(3);
    Real tol(1e-2*std::sqrt(ROL_EPSILON<Real>()));
    // Finalize unpenalized optimization problem
    auto problem = problem0;
    problem->setProjectionAlgorithm(list);
    problem->finalize(false,true,os);
    // Parse problem
    auto obj  = problem0->getObjective();
    auto x    = problem0->getPrimalOptimizationVector();
    auto proj = problem0->getPolyhedralProjection();
    auto bnd  = proj->getBoundConstraint();
    auto lcon = proj->getLinearConstraint();
    auto lmul = proj->getMultiplier();
    auto lres = proj->getResidual();
    auto econ = problem0->getConstraint();
    auto emul = problem0->getMultiplierVector();
    auto eres = problem0->getResidualVector();
    eflag_    = (econ!=nullPtr) && (emul!=nullPtr) && (eres!=nullPtr);
    // Initialize penalty function
    auto dwpen0 = makePtr<DoubleWellPenalty<Real>>(1u);
    dwpen0->forBinaryDesignAlgorithm(eflag_);
    Ptr<Objective<Real>> dwpen = dwpen0;
    if (!equality_||eflag_) dwpen = makePtr<SlacklessObjective<Real>>(dwpen0);
    std::vector<Ptr<Objective<Real>>> ovec = {obj,dwpen};
    std::vector<Real> wvec = {one,zero};
    // Setup optimization solver and solve unpenalized problem
    os << std::endl << "  Solving unpenalized problem" << std::endl << std::endl;
    std::clock_t timer = std::clock();
    auto solver = makePtr<Solver<Real>>(problem,list);
    solver->solve(os);
    os << "  OED optimization time:      "
       << static_cast<Real>(std::clock()-timer)/static_cast<Real>(CLOCKS_PER_SEC)
       << " seconds" << std::endl;
    factory_->profile(os);
    factory_->reset();
    // Compute objective function information
    obj->update(*x,UpdateType::Trial);
    auto c  = obj->value(*x,tol);
    auto xb = x->clone();
    xb->set(*x);
    getDesign(xb)->applyUnary(Elementwise::Round<Real>());
    obj->update(*xb,UpdateType::Trial);
    auto cb = obj->value(*xb,tol);
    getDesign(xb)->axpy(-one,*getDesign(x));
    auto db = getDesign(xb)->norm();
    os << std::string(80,'=') << std::endl;
    os << "  Fractional optimality criterion value: " << c  << std::endl;
    os << "  Encumbent optimality criterion value:  " << cb << std::endl;
    os << "  Distance to binary design:             " << db << std::endl;
    os << std::string(80,'=') << std::endl;
    // Grab initial penalty parameters
    Real pen     = list.sublist("OED").sublist("Binary Design Algorithm").get("Initial Penalty Parameter",1e0);
    Real penscal = list.sublist("OED").sublist("Binary Design Algorithm").get("Penalty Update Scaling",1e1);
    int contype  = list.sublist("OED").sublist("Binary Design Algorithm").get("Continuation Type",1);
    int maxit    = list.sublist("OED").sublist("Binary Design Algorithm").get("Iteration Limit",10);
    Real etol    = list.sublist("OED").sublist("Binary Design Algorithm").get("Tolerance",tol);
    bool rndws   = list.sublist("OED").sublist("Binary Design Algorithm").get("Rounded Warm Start",false);
    // Prepare path-following loop
    Real tau(one/penscal), A(0), B(0), C(0), v(0); //, vtmp(0);
    dwpen->update(*x,UpdateType::Trial);
    auto p  = dwpen->value(*x,tol);
    auto v0 = c;
    pen = contype==1 ? std::max(pen,(cb-c)/p) : pen;
    // Path-following loop
    for (int i=0; i<maxit; ++i) {
      os << std::endl << "  Solving penalized problem with penalty parameter " << pen << std::endl << std::endl;
      // Generate optimization problem
      wvec[1] = pen;
      auto pobj = makePtr<LinearCombinationObjective<Real>>(wvec,ovec);
      if (rndws) x->applyUnary(Elementwise::Round<Real>());
      problem = makePtr<Problem<Real>>(pobj,x);
      problem->addBoundConstraint(bnd);
      problem->addLinearConstraint("Budget",lcon,lmul,lres);
      if (eflag_) problem->addConstraint("Robust",econ,emul,eres);
      problem->setProjectionAlgorithm(list);
      problem->finalize(false,true,os);
      // Setup ROL solver
      timer = std::clock();
      solver = makePtr<Solver<Real>>(problem,list);
      solver->solve(os);
      os << "  OED optimization time:      "
         << static_cast<Real>(std::clock()-timer)/static_cast<Real>(CLOCKS_PER_SEC)
         << " seconds" << std::endl;
      factory_->profile(os);
      factory_->reset();
      // Compute objective function information
      obj->update(*x,UpdateType::Trial);
      c = obj->value(*x,tol);
      xb->set(*x);
      getDesign(xb)->applyUnary(Elementwise::Round<Real>());
      obj->update(*xb,UpdateType::Trial);
      cb = obj->value(*xb,tol);
      getDesign(xb)->axpy(-one,*getDesign(x));
      db = getDesign(xb)->norm();
      os << std::string(80,'=') << std::endl;
      os << "  Fractional optimality criterion value: " << c  << std::endl;
      os << "  Rounded optimality criterion value:    " << cb << std::endl;
      os << "  Distance to binary design:             " << db << std::endl;
      os << std::string(80,'=') << std::endl;
      // Check distance to binary solution 
      if (db < etol) break;
      // Perform path following
      if (contype==1) {
        // Exact path following
        dwpen->update(*x,UpdateType::Trial);
        p = dwpen->value(*x,tol);
        v = c + pen*p;
        C = pen*pen*p/(c-v0);
        B = (C*(C+pen)*(v-v0))/pen;
        A = v0 + B/C;
        pen = B/(tau * std::abs(A-v))-C;
        tau /= penscal;
      }
      //else if (contype==2) {
      //  // Inexact path following
      //  vtmp = v;
      //  crit->update(*x,UpdateType::Trial);
      //  c = crit->value(*x,tol);
      //  v = c + pen*p;
      //  C = pen*pen*p/(c-v0);
      //  B = (C*(C+pen)*(v-v0))/pen;
      //  A = v0 + B/C;
      //  pen = std::max(pen * std::max(tau,two*p),std::pow(two*p,-(three/two)));
      //  t = v + p*(pen'-pen)-m(pen') = tau|vtmp-v|
      //}
      else {
        pen *= penscal;
      }
    }
    os.flags(osFlags);
    return (db<etol);
  }

public:
  BinaryDesignAlgorithm(const Ptr<Objective<Real>>&       model,
                        const Ptr<SampleGenerator<Real>>& sampler,
                        const Ptr<Vector<Real>>&          theta,
                        const Ptr<MomentOperator<Real>>&  cov,
                        ParameterList&                    list)
    : budgetSet_(false), equality_(false), eflag_(false) {
    list.sublist("OED").set("Use Design Upper Bound",true);
    list.sublist("OED").set("Use Double-Well Penalty",false);
    factory_ = makePtr<Factory<Real>>(model,sampler,theta,cov,list);
  }
  BinaryDesignAlgorithm(const Ptr<Constraint<Real>>&      model,
                        const Ptr<SampleGenerator<Real>>& sampler,
                        const Ptr<Vector<Real>>&          theta,
                        const Ptr<Vector<Real>>&          obs,
                        const Ptr<MomentOperator<Real>>&  cov,
                        ParameterList&                    list)
    : budgetSet_(false), equality_(false), eflag_(false) {
    list.sublist("OED").set("Use Design Upper Bound",true);
    list.sublist("OED").set("Use Double-Well Penalty",false);
    factory_ = makePtr<Factory<Real>>(model,sampler,theta,obs,cov,list);
  }

  Ptr<Vector<Real>> createDesignVector() const {
    return factory_->createDesignVector();
  }

  Ptr<Vector<Real>> createCostVector() const {
    return factory_->createDesignVector()->dual().clone();
  }

  void setBudgetConstraint(const Ptr<Vector<Real>> &cost, Real budget, bool equality=false) {
    if (!budgetSet_) {
      factory_->setBudgetConstraint(cost,budget,equality);
      equality_  = equality;
      budgetSet_ = true;
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OED::BinaryDesignAlgorithm::setBudgetConstraint: Constraint already set!");
    }
  }

  void printDesign(const std::string &name, const std::string &ext = ".txt") const {
    factory_->printDesign(name,ext);
  }

  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           const std::vector<Real>& weights,
           const Ptr<Vector<Real>>& c,
           ParameterList& list,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,weights,c,list);
    return runBase(problem,list,os);
  }

  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           const Ptr<Vector<Real>>& c,
           ParameterList& list,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,c,list);
    return runBase(problem,list,os);
  }

  bool run(const Ptr<Vector<Real>>& c,
           ParameterList& list,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(c,list);
    return runBase(problem,list,os);
  }

  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           const std::vector<Real>& weights,
           ParameterList& list,
           const Ptr<SampleGenerator<Real>>& sampler=nullPtr,
           const Ptr<Objective<Real>>& predFun=nullPtr,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,weights,list,sampler,predFun);
    return runBase(problem,list,os);
  }
  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           const std::vector<Real>& weights,
           ParameterList& list,
           const Ptr<SampleGenerator<Real>>& sampler,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,weights,list,sampler);
    return runBase(problem,list,os);
  }
  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           const std::vector<Real>& weights,
           ParameterList& list,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,weights,list);
    return runBase(problem,list,os);
  }

  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           ParameterList& list,
           const Ptr<SampleGenerator<Real>>& sampler=nullPtr,
           const Ptr<Objective<Real>>& predFun=nullPtr,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,list,sampler,predFun);
    return runBase(problem,list,os);
  }
  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           ParameterList& list,
           const Ptr<SampleGenerator<Real>>& sampler,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,list,sampler);
    return runBase(problem,list,os);
  }
  bool run(const std::vector<Ptr<Vector<Real>>>& thetas,
           ParameterList& list,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(thetas,list);
    return runBase(problem,list,os);
  }

  bool run(ParameterList& list,
           const Ptr<SampleGenerator<Real>>& sampler,
           const Ptr<Objective<Real>>& predFun,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(list,sampler,predFun);
    return runBase(problem,list,os);
  }
  bool run(ParameterList& list,
           const Ptr<SampleGenerator<Real>>& sampler,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(list,sampler);
    return runBase(problem,list,os);
  }
  bool run(ParameterList& list,
           std::ostream& os=std::cout) {
    auto problem = factory_->get(list);
    return runBase(problem,list,os);
  }
};

} // End ROL::OED Namespace

#endif

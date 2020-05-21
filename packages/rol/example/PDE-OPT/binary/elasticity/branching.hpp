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

#ifndef ROL_PDEOPT_MULTIMAT_BRANCHING_PEBBL_H
#define ROL_PDEOPT_MULTIMAT_BRANCHING_PEBBL_H

#include "ROL_PEBBL_Interface.hpp"
#include "../../TOOLS/pdevector.hpp"

template<class Real>
class MultiMatBranchSub;

template<class Real>
class MultiMatBranching : public ROL::PEBBL::Branching<Real> {
private:
  ROL::Ptr<ROL::Vector<Real>> z0_;

  using ROL::PEBBL::Branching<Real>::verbosity_;
  using ROL::PEBBL::Branching<Real>::outStream_;
  using ROL::PEBBL::Branching<Real>::parlist_;

public:
  MultiMatBranching(const ROL::Ptr<ElasticityFactory<Real>>        &factory,
                    const ROL::Ptr<ROL::ParameterList>             &parlist,
                    const ROL::Ptr<ROL::PEBBL::BranchHelper<Real>> &bHelper,
                    int                                             verbosity = 0,
                    const ROL::Ptr<std::ostream>                   &outStream = ROL::nullPtr)
    : ROL::PEBBL::Branching<Real>(factory,parlist,bHelper,verbosity,outStream) {
    z0_ = factory->buildSolutionVector();
  }

  pebbl::branchSub* blankSub() {
    return new MultiMatBranchSub<Real>(parlist_,ROL::makePtrFromRef<MultiMatBranching<Real>>(*this),verbosity_,outStream_);
  }

//  pebbl::solution* iniitalGuess() {
//
//  }
}; // MultiMatBranching

template <class Real>
class MultiMatBranchSub : public ROL::PEBBL::BranchSub<Real> {
private:
  int method_;
  Real ctol_;
  int T_;
  std::string methodName_;
  ROL::Ptr<ROL::Vector<Real>> c_;
  ROL::Ptr<ROL::Constraint<Real>> con_;

  using ROL::PEBBL::BranchSub<Real>::anyChild;
  using ROL::PEBBL::BranchSub<Real>::index_;
  using ROL::PEBBL::BranchSub<Real>::branching_;
  using ROL::PEBBL::BranchSub<Real>::problem0_;
  using ROL::PEBBL::BranchSub<Real>::solution_;
  using ROL::PEBBL::BranchSub<Real>::rndSolution_;
  using ROL::PEBBL::BranchSub<Real>::verbosity_;
  using ROL::PEBBL::BranchSub<Real>::outStream_;

  void round(ROL::Vector<Real> &rx, const ROL::Vector<Real> &x, Real t) const {
    rx.set(x);
    Teuchos::ArrayView<Real> data = getData(rx);
    Real val(0);
    for (auto it = data.begin(); it != data.end(); ++it) {
      val = *it;
      *it = (val < t ? std::floor(val) : std::ceil(val));
    }
  }

  ROL::Ptr<ROL::Vector<Real>> get(ROL::Vector<Real> &x, int ind) const {
    return dynamic_cast<ROL::PartitionedVector<Real>&>(x).get(ind);
  }

  ROL::Ptr<const ROL::Vector<Real>> get(const ROL::Vector<Real> &x, int ind) const {
    return dynamic_cast<const ROL::PartitionedVector<Real>&>(x).get(ind);
  }

  Teuchos::ArrayView<Real> getData(ROL::Vector<Real> &x) const {
    try {
      return (ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(get(x,0))->getVector()->getDataNonConst(0))();
    }
    catch (std::exception &e) {
      return (ROL::dynamicPtrCast<PDE_OptVector<Real>>(get(x,0))->getField()->getVector()->getDataNonConst(0))();
    }
  }

  void zeroSlack(ROL::Vector<Real> &x) const {
    size_t nv = dynamic_cast<ROL::PartitionedVector<Real>&>(x).numVectors();
    for (size_t i = 1; i < nv; ++i) {
      get(x,i)->zero();
    }
  }

  void setSlack(ROL::Vector<Real> &x, const ROL::Vector<Real> &c) const {
    size_t nv = dynamic_cast<ROL::PartitionedVector<Real>&>(x).numVectors();
    for (size_t i = 1; i < nv; ++i) {
      get(x,i)->set(*get(c,i-1));
    }
    problem0_->getBoundConstraint()->project(x);
  }

public:
  MultiMatBranchSub(const ROL::Ptr<ROL::ParameterList> &parlist,
                    const ROL::Ptr<ROL::PEBBL::Branching<Real>> &branching,
                    int verbosity = 0,
                    const ROL::Ptr<std::ostream> &outStream = ROL::nullPtr)
    : ROL::PEBBL::BranchSub<Real>(branching, verbosity, outStream) {
    method_ = parlist->sublist("Problem").get("Incumbent Heuristic",0);
    ctol_   = parlist->sublist("Status Test").get("Constraint Tolerance",1e-8);
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(parlist->sublist("Problem"), "Young's Modulus");
    T_ = ym.size();
    methodName_ = "Default";
    if (problem0_->getConstraint()==ROL::nullPtr) {
      c_   = problem0_->getPolyhedralProjection()->getResidual()->clone();
      con_ = problem0_->getPolyhedralProjection()->getLinearConstraint();
    }
    else {
      c_   = problem0_->getResidualVector()->clone();
      con_ = problem0_->getConstraint();
    }
  }

  MultiMatBranchSub(const MultiMatBranchSub &rpbs)
    : ROL::PEBBL::BranchSub<Real>(rpbs),
      method_(rpbs.method_), ctol_(rpbs.ctol_), T_(rpbs.T_), methodName_(rpbs.methodName_),
      c_(rpbs.c_->clone()), con_(rpbs.con_) {}

  void incumbentHeuristic() {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    const Real zero(0), one(1);
    Teuchos::ArrayView<Real> data = getData(*rndSolution_);
    const int nc = data.size()/T_;
    Real r(0), sum(0), cnorm(0);
    int cnt(0);
    bool oneSet(false);
    while (true) {
      rndSolution_->set(*solution_);
      for (int i = 0; i < nc; ++i) { 
        r = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
        sum = zero;
        oneSet = false;
        for (int t = 0; t < T_; ++t) {
          sum += data[t+i*T_];
          if (r <= sum && !oneSet) {
            data[t+i*T_] = one;
            oneSet = true;
          }
          else {
            data[t+i*T_] = zero;
          }
        }
      }
      zeroSlack(*rndSolution_);
      con_->value(*c_,*rndSolution_,tol);
      setSlack(*rndSolution_,*c_);
      con_->value(*c_,*rndSolution_,tol);
      cnorm = c_->norm();
      cnt++;
      if (cnorm < ctol_) break;
      if (verbosity_ > 1) {
        *outStream_ << "  cnt = " << cnt << "  infeasibility = " << cnorm << std::endl;
      }
    }
    problem0_->getObjective()->update(*rndSolution_);
    Real val = problem0_->getObjective()->value(*rndSolution_,tol);
    branching_->foundSolution(new ROL::PEBBL::IntegerSolution<Real>(*rndSolution_,val));
    if (verbosity_ > 0) {
      *outStream_ << "MultiMatBranchSub::incumbentHeuristic: " << methodName_ << std::endl;
      *outStream_ << "  Incumbent Value:       " <<   val << std::endl;
      *outStream_ << "  Incumbent Feasibility: " << cnorm << std::endl;
      *outStream_ << "  Number of Samples:     " <<   cnt << std::endl;
    }
  }

  pebbl::branchSub* makeChild(int whichChild = anyChild) override {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild, std::logic_error,
      ">>> MultiMatBranchSub::makeChild: whichChild is equal to anyChild!");
    MultiMatBranchSub<Real>* child = new MultiMatBranchSub<Real>(*this);
    child->updateFixed(index_,
      (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

}; // class MultiMatBranchSub

#endif

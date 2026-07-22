// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STEFANBOLTZMANN_BRANCHING_H
#define ROL_STEFANBOLTZMANN_BRANCHING_H

#include "ROL_PEBBL_Interface.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "hilbert.hpp"

template<class Real>
class TpetraStefanBoltzmannBranchSub;

template<class Real>
class StdStefanBoltzmannBranchSub;

template<class Real>
class TpetraStefanBoltzmannBranching : public ROL::PEBBL::Branching<Real> {
private:
  const int method_;
  ROL::Ptr<ROL::Vector<Real>> z0_;

  using ROL::PEBBL::Branching<Real>::verbosity_;
  using ROL::PEBBL::Branching<Real>::outStream_;
  using ROL::PEBBL::Branching<Real>::parlist_;

public:
  TpetraStefanBoltzmannBranching(const ROL::Ptr<BinaryStefanBoltzmannFactory<Real>> &factory,
                                 const ROL::Ptr<ROL::ParameterList>                 &parlist,
                                 const ROL::Ptr<ROL::PEBBL::BranchHelper<Real>>     &bHelper,
                                 int                                                 verbosity = 0,
                                 const ROL::Ptr<std::ostream>                       &outStream = ROL::nullPtr,
                                 int                                                 method = 0)
    : ROL::PEBBL::Branching<Real>(factory,parlist,bHelper,verbosity,outStream),
      method_(method) {
    z0_ = factory->buildSolutionVector();
  }

  pebbl::branchSub* blankSub() {
    return new TpetraStefanBoltzmannBranchSub<Real>(ROL::makePtrFromRef<TpetraStefanBoltzmannBranching<Real>>(*this),verbosity_,outStream_,method_);
  }

//  pebbl::solution* iniitalGuess() {
//
//  }
}; // StefanBoltzmannBranching

template <class Real>
class TpetraStefanBoltzmannBranchSub : public ROL::PEBBL::BranchSub<Real> {
private:
  const int method_;
  std::string methodName_;

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

  Teuchos::ArrayView<Real> getData(ROL::Vector<Real> &x) const {
    try {
      return (dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector()->getDataNonConst(0))();
    }
    catch (std::exception &e) {
      return (dynamic_cast<PDE_OptVector<Real>&>(x).getField()->getVector()->getDataNonConst(0))();
    }
  }

public:
  TpetraStefanBoltzmannBranchSub(const ROL::Ptr<ROL::PEBBL::Branching<Real>> &branching,
                                 int verbosity = 0,
                                 const ROL::Ptr<std::ostream> &outStream = ROL::nullPtr,
                                 int method = 0)
    : ROL::PEBBL::BranchSub<Real>(branching, verbosity, outStream), method_(method) {
    switch (method_) {
      case 1:  methodName_ = "Mass Preserving Rounding"; break;
      case 2:  methodName_ = "Objective Gap Rounding";   break;
      case 3:  methodName_ = "Sum Up Rounding";          break;
      default: methodName_ = "Naive Rounding";
    }
  }

  TpetraStefanBoltzmannBranchSub(const TpetraStefanBoltzmannBranchSub &rpbs)
    : ROL::PEBBL::BranchSub<Real>(rpbs),
      method_(rpbs.method_), methodName_(rpbs.methodName_) {}

  void incumbentHeuristic() {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    Real t(0), val(0);
    if (method_==1) { // mass preserving rounding
      const Real zero(0), one(1);
      rndSolution_->set(*solution_);
      Teuchos::ArrayView<Real> data = getData(*rndSolution_);
      std::map<Real,size_t> sdata;
      Real sum(0);
      for (size_t i = 0; i < static_cast<size_t>(data.size()); ++i) {
        sdata.insert(std::pair<Real,size_t>(data[i],i));
        sum += data[i];
      }
      int rsum = static_cast<int>(std::round(sum)), psum(0);
      for (auto it = sdata.rbegin(); it != sdata.rend(); ++it) {
        data[it->second] = (psum < rsum ? one : zero);
        t = (psum < rsum ? it->first : t);
        psum++;
      }
    }
    else if (method_==2) { // objective gap rounding
      const Real invphi(0.5*(std::sqrt(5.0)-1.0));
      const Real invphi2(0.5*(3.0-std::sqrt(5.0)));
      const Real itol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      Real a(0), b(1), c(0), d(1), h(b-a), fc(0), fd(0);
      // Evaluate at c
      c = a + invphi2 * h;
      round(*rndSolution_,*solution_,c);
      problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
      fc = problem0_->getObjective()->value(*rndSolution_,tol);
      // Evaluate at d
      d = a + invphi * h;
      round(*rndSolution_,*solution_,d);
      problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
      fd = problem0_->getObjective()->value(*rndSolution_,tol);
      while (std::abs(c-d) > itol) {
        h *= invphi;
        if (fc < fd) {
          b  = d;
          d  = c;
          fd = fc;
          c  = a + invphi2 * h;
          round(*rndSolution_,*solution_,c);
          problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
          fc = problem0_->getObjective()->value(*rndSolution_,tol);
        }
        else {
          a  = c;
          c  = d;
          fc = fd;
          d  = a + invphi * h;
          round(*rndSolution_,*solution_,d);
          problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
          fd = problem0_->getObjective()->value(*rndSolution_,tol);
        }
      }
      t = static_cast<Real>(0.5)*(a+b);
      round(*rndSolution_,*solution_,t);
    }
    else if (method_==3) { // Sum up rounding
      const Real zero(0), one(1);
      rndSolution_->set(*solution_);
      Teuchos::ArrayView<Real> data = getData(*rndSolution_);
      const int n = data.size();
      const int m = (std::log2(n)-1)/2;
      Real g1(0), g2(0);
      if (n == std::pow(2,2*m+1)) { // Use Hilbert curves if n = 2^{2m+1}
        const int nx = std::pow(2,m+1);
        const int ny = std::pow(2,m);
        int x(0), y(0), ind(0);
        // Domain [0,1]x[0,1]
        for (int j = 0; j < 2; ++j) {
          for (int i = 0; i < ny*ny; ++i) {
            hilbert::d2xy(m,i,x,y);
            ind = x + y*nx + j*ny;
            g1 += data[ind];
            g2 += one - data[ind];
            if (g1 >= g2) {
              data[ind] = one;
              g1 -= one;
            }
            else {
              data[ind] = zero;
              g2 -= one;
            }
          }
        }
      }
      else { // Otherwise use native ordering of cells
        for (int i = 0; i < n; ++i) {
          g1 += data[i];
          g2 += one - data[i];
          if (g1 >= g2) {
            data[i] = one;
            g1 -= one;
          }
          else {
            data[i] = zero;
            g2 -= one;
          }
        }
      }
    }
    else { // naive rounding
      t = static_cast<Real>(0.5);
      round(*rndSolution_,*solution_,t);
    }
    problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
    val = problem0_->getObjective()->value(*rndSolution_,tol);
    branching_->foundSolution(new ROL::PEBBL::IntegerSolution<Real>(*rndSolution_,val));
    if (verbosity_ > 0) {
      *outStream_ << "StefanBoltzmannBranchSub::incumbentHeuristic: " << methodName_ << std::endl;
      *outStream_ << "  Incumbent Value:    " << val  << std::endl;
      *outStream_ << "  Rounding Threshold: " << t    << std::endl;
    }
  }

  pebbl::branchSub* makeChild(int whichChild = anyChild) override {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild,std::logic_error,
      ">>> StefanBoltzmannBranchSub::makeChild: whichChild is equal to anyChild!");
    TpetraStefanBoltzmannBranchSub<Real>* child
      = new TpetraStefanBoltzmannBranchSub<Real>(*this);
    child->updateFixed(index_,
      (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

}; // class StefanBoltzmannBranchSub

template<class Real>
class StdStefanBoltzmannBranching : public ROL::PEBBL::Branching<Real> {
private:
  const int method_;
  ROL::Ptr<ROL::Vector<Real>> z0_;

  using ROL::PEBBL::Branching<Real>::verbosity_;
  using ROL::PEBBL::Branching<Real>::outStream_;
  using ROL::PEBBL::Branching<Real>::parlist_;

public:
  StdStefanBoltzmannBranching(const ROL::Ptr<BinaryStefanBoltzmannFactory<Real>> &factory,
                              const ROL::Ptr<ROL::ParameterList>                 &parlist,
                              const ROL::Ptr<ROL::PEBBL::BranchHelper<Real>>     &bHelper,
                              int                                                 verbosity = 0,
                              const ROL::Ptr<std::ostream>                       &outStream = ROL::nullPtr,
                              int                                                 method = 0)
    : ROL::PEBBL::Branching<Real>(factory,parlist,bHelper,verbosity,outStream),
      method_(method) {
    z0_ = factory->buildSolutionVector();
  }

  pebbl::branchSub* blankSub() {
    return new StdStefanBoltzmannBranchSub<Real>(ROL::makePtrFromRef<StdStefanBoltzmannBranching<Real>>(*this),verbosity_,outStream_,method_);
  }

//  pebbl::solution* iniitalGuess() {
//
//  }
}; // StefanBoltzmannBranching

template <class Real>
class StdStefanBoltzmannBranchSub : public ROL::PEBBL::BranchSub<Real> {
private:
  const int method_;
  std::string methodName_;

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
    ROL::Ptr<std::vector<Real>> data = getData(rx);
    Real val(0);
    for (auto it = data->begin(); it != data->end(); ++it) {
      val = *it;
      *it = (val < t ? std::floor(val) : std::ceil(val));
    }
  }

  ROL::Ptr<std::vector<Real>> getData(ROL::Vector<Real> &x) const {
    try {
      return ROL::dynamicPtrCast<ROL::StdVector<Real>>(dynamic_cast<ROL::PartitionedVector<Real>&>(x).get(0))->getVector();
    }
    catch (std::exception &e) {
      return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
  }

public:
  StdStefanBoltzmannBranchSub(const ROL::Ptr<ROL::PEBBL::Branching<Real>> &branching,
                              int verbosity = 0,
                              const ROL::Ptr<std::ostream> &outStream = ROL::nullPtr,
                              int method = 0)
    : ROL::PEBBL::BranchSub<Real>(branching, verbosity, outStream), method_(method) {
    switch (method_) {
      case 1:  methodName_ = "Mass Preserving Rounding"; break;
      case 2:  methodName_ = "Objective Gap Rounding";   break;
      case 3:  methodName_ = "Sum Up Rounding";          break;
      default: methodName_ = "Naive Rounding";
    }
  }

  StdStefanBoltzmannBranchSub(const StdStefanBoltzmannBranchSub &rpbs)
    : ROL::PEBBL::BranchSub<Real>(rpbs),
      method_(rpbs.method_), methodName_(rpbs.methodName_) {}

  void incumbentHeuristic() {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    Real t(0), val(0);
    if (method_==1) { // mass preserving rounding
      const Real zero(0), one(1);
      rndSolution_->set(*solution_);
      ROL::Ptr<std::vector<Real>> data = getData(*rndSolution_);
      std::map<Real,size_t> sdata;
      Real sum(0);
      for (size_t i = 0; i < static_cast<size_t>(data->size()); ++i) {
        sdata.insert(std::pair<Real,size_t>((*data)[i],i));
        sum += (*data)[i];
      }
      int rsum = static_cast<int>(std::round(sum)), psum(0);
      for (auto it = sdata.rbegin(); it != sdata.rend(); ++it) {
        (*data)[it->second] = (psum < rsum ? one : zero);
        t = (psum < rsum ? it->first : t);
        psum++;
      }
    }
    else if (method_==2) { // objective gap rounding
      const Real invphi(0.5*(std::sqrt(5.0)-1.0));
      const Real invphi2(0.5*(3.0-std::sqrt(5.0)));
      const Real itol(std::sqrt(ROL::ROL_EPSILON<Real>()));
      Real a(0), b(1), c(0), d(1), h(b-a), fc(0), fd(0);
      // Evaluate at c
      c = a + invphi2 * h;
      round(*rndSolution_,*solution_,c);
      problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
      fc = problem0_->getObjective()->value(*rndSolution_,tol);
      // Evaluate at d
      d = a + invphi * h;
      round(*rndSolution_,*solution_,d);
      problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
      fd = problem0_->getObjective()->value(*rndSolution_,tol);
      while (std::abs(c-d) > itol) {
        h *= invphi;
        if (fc < fd) {
          b  = d;
          d  = c;
          fd = fc;
          c  = a + invphi2 * h;
          round(*rndSolution_,*solution_,c);
          problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
          fc = problem0_->getObjective()->value(*rndSolution_,tol);
        }
        else {
          a  = c;
          c  = d;
          fc = fd;
          d  = a + invphi * h;
          round(*rndSolution_,*solution_,d);
          problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
          fd = problem0_->getObjective()->value(*rndSolution_,tol);
        }
      }
      t = static_cast<Real>(0.5)*(a+b);
      round(*rndSolution_,*solution_,t);
    }
    else { // naive rounding
      t = static_cast<Real>(0.5);
      round(*rndSolution_,*solution_,t);
    }
    problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
    val = problem0_->getObjective()->value(*rndSolution_,tol);
    branching_->foundSolution(new ROL::PEBBL::IntegerSolution<Real>(*rndSolution_,val));
    if (verbosity_ > 0) {
      *outStream_ << "StefanBoltzmannBranchSub::incumbentHeuristic: " << methodName_ << std::endl;
      *outStream_ << "  Incumbent Value:    " << val  << std::endl;
      *outStream_ << "  Rounding Threshold: " << t    << std::endl;
    }
  }

  pebbl::branchSub* makeChild(int whichChild = anyChild) override {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild,std::logic_error,
      ">>> StefanBoltzmannBranchSub::makeChild: whichChild is equal to anyChild!");
    StdStefanBoltzmannBranchSub<Real>* child
      = new StdStefanBoltzmannBranchSub<Real>(*this);
    child->updateFixed(index_,
      (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

}; // class StefanBoltzmannBranchSub

#endif

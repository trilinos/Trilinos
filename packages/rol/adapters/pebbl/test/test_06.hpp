// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

// ROL_Types contains predefined constants and objects
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Solver.hpp"
#include "ROL_PEBBL_MixedVector.hpp"
#include "ROL_PEBBL_IntegerConstraint.hpp"
#include "ROL_PEBBL_IntegerProblemFactory.hpp"
#include "ROL_PEBBL_BranchAndBound.hpp"
#include "ROL_PEBBL_StdBranchHelper.hpp"

template<class Real>
class EqualityConstraint_FacilityLocation : public ROL::Constraint<Real> {
private:
  ROL::Ptr<ROL::Vector<Real>> getIntegerVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::PEBBL::MixedVector<Real>&>(x).getIntegerVariables();
  }

  ROL::Ptr<ROL::Vector<Real>> getContinuousVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::PEBBL::MixedVector<Real>&>(x).getContinuousVariables();
  }

  ROL::Ptr<const ROL::Vector<Real>> getIntegerVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::PEBBL::MixedVector<Real>&>(x).getIntegerVariables();
  }

  ROL::Ptr<const ROL::Vector<Real>> getContinuousVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::PEBBL::MixedVector<Real>&>(x).getContinuousVariables();
  }


  ROL::Ptr<std::vector<Real>> getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

public:
  EqualityConstraint_FacilityLocation() {}

  void value(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &x, 
             Real &tol) {
    ROL::Ptr<std::vector<Real>>       cp = getVector(c);
    ROL::Ptr<const std::vector<Real>> xi = getConstVector(*getIntegerVector(x));
    ROL::Ptr<const std::vector<Real>> xc = getConstVector(*getContinuousVector(x));
    const Real one(1);
    const int M = xi->size(), N = xc->size()/M;
    for (int j = 0; j < N; ++j) {
      (*cp)[j] = -one;
      for (int i = 0; i < M; ++i) {
        (*cp)[j] += (*xc)[i + j*M];
      }
    }
  }

  void applyJacobian(ROL::Vector<Real> &jv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &x, 
                     Real &tol) {
    ROL::Ptr<std::vector<Real>>      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>> vi = getConstVector(*getIntegerVector(v));
    ROL::Ptr<const std::vector<Real>> vc = getConstVector(*getContinuousVector(v));
    const Real zero(0);
    const int M = vi->size(), N = vc->size()/M;
    for (int j = 0; j < N; ++j) {
      (*jvp)[j] = zero;
      for (int i = 0; i < M; ++i) {
        (*jvp)[j] += (*vc)[i + j*M];
      }
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                      const ROL::Vector<Real> &v,
                      const ROL::Vector<Real> &x, 
                            Real &tol) {
    ROL::Ptr<std::vector<Real>>      jvi = getVector(*getIntegerVector(ajv));
    ROL::Ptr<std::vector<Real>>      jvc = getVector(*getContinuousVector(ajv));
    ROL::Ptr<const std::vector<Real>> vp = getConstVector(v);
    const int M = jvi->size(), N = jvc->size()/M;
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        (*jvc)[i + j*M] = (*vp)[j];
      }
    }
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv,
                     const ROL::Vector<Real> &w,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &x,
                           Real &tol) {
    ahwv.zero();
  }
};

template<class Real>
class InequalityConstraint_FacilityLocation : public ROL::Constraint<Real> {
private:
  ROL::Ptr<ROL::Vector<Real>> getIntegerVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::PEBBL::MixedVector<Real>&>(x).getIntegerVariables();
  }

  ROL::Ptr<ROL::Vector<Real>> getContinuousVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::PEBBL::MixedVector<Real>&>(x).getContinuousVariables();
  }

  ROL::Ptr<const ROL::Vector<Real>> getIntegerVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::PEBBL::MixedVector<Real>&>(x).getIntegerVariables();
  }

  ROL::Ptr<const ROL::Vector<Real>> getContinuousVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::PEBBL::MixedVector<Real>&>(x).getContinuousVariables();
  }


  ROL::Ptr<std::vector<Real>> getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

public:
  InequalityConstraint_FacilityLocation() {}

  void value(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &x, 
             Real &tol) {
    ROL::Ptr<std::vector<Real>>       cp = getVector(c);
    ROL::Ptr<const std::vector<Real>> xi = getConstVector(*getIntegerVector(x));
    ROL::Ptr<const std::vector<Real>> xc = getConstVector(*getContinuousVector(x));
    const int M = xi->size(), N = xc->size()/M;
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        (*cp)[i + j*M] = (*xc)[i + j*M] - (*xi)[i];
      }
    }
  }

  void applyJacobian(ROL::Vector<Real> &jv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &x, 
                     Real &tol) {
    ROL::Ptr<std::vector<Real>>      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>> vi = getConstVector(*getIntegerVector(v));
    ROL::Ptr<const std::vector<Real>> vc = getConstVector(*getContinuousVector(v));
    const int M = vi->size(), N = vc->size()/M;
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        (*jvp)[i + j*M] = (*vc)[i + j*M] - (*vi)[i];
      }
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                      const ROL::Vector<Real> &v,
                      const ROL::Vector<Real> &x, 
                            Real &tol) {
    ROL::Ptr<std::vector<Real>>      jvi = getVector(*getIntegerVector(ajv));
    ROL::Ptr<std::vector<Real>>      jvc = getVector(*getContinuousVector(ajv));
    ROL::Ptr<const std::vector<Real>> vp = getConstVector(v);
    const Real zero(0);
    const int M = jvi->size(), N = jvc->size()/M;
    for (int i = 0; i < M; ++i) {
        (*jvi)[i] = zero; 
      for (int j = 0; j < N; ++j) {
        (*jvc)[i + j*M] = (*vp)[i + j*M];
        (*jvi)[i] -= (*vp)[i + j*M];
      }
    }
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv,
                     const ROL::Vector<Real> &w,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &x,
                           Real &tol) {
    ahwv.zero();
  }
};

template<class Real>
class Objective_FacilityLocation : public ROL::Objective<Real> {
private:
  const std::vector<Real> c_, q_;
  const int M_, N_;

  ROL::Ptr<ROL::Vector<Real>> getIntegerVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::PEBBL::MixedVector<Real>&>(x).getIntegerVariables();
  }

  ROL::Ptr<ROL::Vector<Real>> getContinuousVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::PEBBL::MixedVector<Real>&>(x).getContinuousVariables();
  }

  ROL::Ptr<const ROL::Vector<Real>> getIntegerVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::PEBBL::MixedVector<Real>&>(x).getIntegerVariables();
  }

  ROL::Ptr<const ROL::Vector<Real>> getContinuousVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::PEBBL::MixedVector<Real>&>(x).getContinuousVariables();
  }

  ROL::Ptr<std::vector<Real>> getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

public:
  Objective_FacilityLocation(const std::vector<Real> &c, const std::vector<Real> &q)
    : c_(c), q_(q), M_(c_.size()), N_(q_.size()/M_) {}

  Real value(const ROL::Vector<Real> &x, Real &tol) {
    Real val(0);
    ROL::Ptr<const std::vector<Real>> xi = getConstVector(*getIntegerVector(x));
    ROL::Ptr<const std::vector<Real>> xc = getConstVector(*getContinuousVector(x));
    const Real two(2);
    for (int i = 0; i < M_; ++i) {
      val += c_[i] * (*xi)[i];
      for (int j = 0; j < N_; j++) {
        val += q_[i + j*M_]*std::pow((*xc)[i + j*M_],two);
      }
    }
    return val;
  }

  void gradient(ROL::Vector<Real> &g,
          const ROL::Vector<Real> &x,
                Real &tol) {
    g.zero();
    ROL::Ptr<std::vector<Real>>       gi = getVector(*getIntegerVector(g));
    ROL::Ptr<std::vector<Real>>       gc = getVector(*getContinuousVector(g));
    ROL::Ptr<const std::vector<Real>> xi = getConstVector(*getIntegerVector(x));
    ROL::Ptr<const std::vector<Real>> xc = getConstVector(*getContinuousVector(x));
    const Real two(2);
    for (int i = 0; i < M_; ++i) {
      (*gi)[i] = c_[i];
      for (int j = 0; j < N_; ++j) {
        (*gc)[i + j*M_] = two*q_[i + j*M_]*(*xc)[i + j*M_];
      }
    }
  }

  void hessVec(ROL::Vector<Real> &hv,
         const ROL::Vector<Real> &v,
         const ROL::Vector<Real> &x,
               Real &tol ) {
    hv.zero();
    ROL::Ptr<std::vector<Real>>       hc = getVector(*getContinuousVector(hv));
    ROL::Ptr<const std::vector<Real>> vc = getConstVector(*getContinuousVector(v));
    const Real two(2);
    for (int i = 0; i < M_; ++i) {
      for (int j = 0; j < N_; ++j) {
        (*hc)[i + j*M_] = two*q_[i + j*M_]*(*vc)[i + j*M_];
      }
    }
  }
};

template<class Real>
class FacilityLocationFactory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  ROL::Ptr<ROL::Vector<Real>> x_, emul_, imul_;
  ROL::Ptr<ROL::BoundConstraint<Real>> bnd_, ibnd_;
  ROL::Ptr<ROL::Objective<Real>> obj_;
  ROL::Ptr<ROL::Constraint<Real>> econ_, icon_;
  bool useLinearCon_;

public:
  FacilityLocationFactory(const std::vector<Real> &c, const std::vector<Real> &q,
                          ROL::ParameterList &pl) {
    const int M = c.size();
    const int N = q.size()/M;
    obj_  = ROL::makePtr<Objective_FacilityLocation<Real>>(c,q);

    ROL::Ptr<ROL::Vector<Real>> xc, xi, xl, xu, zl, zu, lo, up, iup;
    xc    = ROL::makePtr<ROL::StdVector<Real>>(M*N,0.0); 
    xi    = ROL::makePtr<ROL::StdVector<Real>>(M,0.0); 
    x_    = ROL::makePtr<ROL::PEBBL::MixedVector<Real>>(xc,xi);

    xl    = ROL::makePtr<ROL::StdVector<Real>>(M*N,0.0);
    xu    = ROL::makePtr<ROL::StdVector<Real>>(M*N,1.0);
    zl    = ROL::makePtr<ROL::StdVector<Real>>(M,0.0);
    zu    = ROL::makePtr<ROL::StdVector<Real>>(M,1.0);
    lo    = ROL::makePtr<ROL::PEBBL::MixedVector<Real>>(xl,zl);
    up    = ROL::makePtr<ROL::PEBBL::MixedVector<Real>>(xu,zu);
    bnd_  = ROL::makePtr<ROL::Bounds<Real>>(lo,up);

    econ_ = ROL::makePtr<EqualityConstraint_FacilityLocation<Real>>();
    emul_ = ROL::makePtr<ROL::StdVector<Real>>(N,0.0);

    icon_ = ROL::makePtr<InequalityConstraint_FacilityLocation<Real>>();
    imul_ = ROL::makePtr<ROL::StdVector<Real>>(M*N,0.0);
    iup   = ROL::makePtr<ROL::StdVector<Real>>(M*N,0.0);
    ibnd_ = ROL::makePtr<ROL::Bounds<Real>>(*iup,false);

    useLinearCon_ = pl.sublist("Problem").get("Maintain Linear Constraints",true);
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) {
    ROL::Ptr<ROL::Vector<Real>>    x = x_->clone();
    ROL::Ptr<ROL::Vector<Real>> emul = emul_->clone();
    ROL::Ptr<ROL::Vector<Real>> imul = imul_->clone();
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(obj_,x);
    problem->addBoundConstraint(bnd_);
    if (useLinearCon_) {
      problem->addLinearConstraint("Equality",econ_,emul);
      problem->addLinearConstraint("Inequality",icon_,imul,ibnd_);
    }
    else {
      problem->addConstraint("Equality",econ_,emul);
      problem->addConstraint("Inequality",icon_,imul,ibnd_);
    }
    return problem;
  }
};

template<class Real>
class FacilityLocationBranchSub;

template<class Real>
class FacilityLocationBranching : public ROL::PEBBL::Branching<Real> {
private:

  using ROL::PEBBL::Branching<Real>::verbosity_;
  using ROL::PEBBL::Branching<Real>::outStream_;
  using ROL::PEBBL::Branching<Real>::parlist_;

public:
  FacilityLocationBranching(const ROL::Ptr<FacilityLocationFactory<Real>>  &factory,
                            const ROL::Ptr<ROL::ParameterList>             &parlist,
                            const ROL::Ptr<ROL::PEBBL::BranchHelper<Real>> &bHelper,
                            int                                             verbosity = 0,
                            const ROL::Ptr<std::ostream>                   &outStream = ROL::nullPtr)
    : ROL::PEBBL::Branching<Real>(factory,parlist,bHelper,verbosity,outStream) {}

  pebbl::branchSub* blankSub() {
    return new FacilityLocationBranchSub<Real>(parlist_,ROL::makePtrFromRef<FacilityLocationBranching<Real>>(*this),verbosity_,outStream_);
  }
}; // FacilityLocationBranching

template <class Real>
class FacilityLocationBranchSub : public ROL::PEBBL::BranchSub<Real> {
private:
  Real ctol_;
  ROL::Ptr<ROL::Vector<Real>> x_, c_;
  ROL::Ptr<ROL::Constraint<Real>> con_;
  ROL::Ptr<ROL::PolyhedralProjection<Real>> proj_;

  using ROL::PEBBL::BranchSub<Real>::anyChild;
  using ROL::PEBBL::BranchSub<Real>::index_;
  using ROL::PEBBL::BranchSub<Real>::branching_;
  using ROL::PEBBL::BranchSub<Real>::problem0_;
  using ROL::PEBBL::BranchSub<Real>::solution_;
  using ROL::PEBBL::BranchSub<Real>::rndSolution_;
  using ROL::PEBBL::BranchSub<Real>::verbosity_;
  using ROL::PEBBL::BranchSub<Real>::outStream_;

  ROL::Ptr<std::vector<Real>> getData(ROL::Vector<Real> &x) const {
    return ROL::dynamicPtrCast<ROL::StdVector<Real>>(
             ROL::dynamicPtrCast<ROL::PEBBL::MixedVector<Real>>(get(x,0))->getIntegerVariables())->getVector();
  }

  ROL::Ptr<std::vector<Real>> getData(ROL::Vector<Real> &x, int comp) const {
    if (comp == 0)
      return ROL::dynamicPtrCast<ROL::StdVector<Real>>(
             ROL::dynamicPtrCast<ROL::PEBBL::MixedVector<Real>>(get(x,0))->getIntegerVariables())->getVector();
    else
      return ROL::dynamicPtrCast<ROL::StdVector<Real>>(
             ROL::dynamicPtrCast<ROL::PEBBL::MixedVector<Real>>(get(x,0))->getContinuousVariables())->getVector();
  }


  ROL::Ptr<const std::vector<Real>> getConstData(const ROL::Vector<Real> &x) const {
    return ROL::dynamicPtrCast<const ROL::StdVector<Real>>(
             ROL::dynamicPtrCast<const ROL::PEBBL::MixedVector<Real>>(get(x,0))->getContinuousVariables())->getVector();
  }

  ROL::Ptr<ROL::Vector<Real>> get(ROL::Vector<Real> &x, int ind) const {
    return dynamic_cast<ROL::PartitionedVector<Real>&>(x).get(ind);
  }

  ROL::Ptr<const ROL::Vector<Real>> get(const ROL::Vector<Real> &x, int ind) const {
    return dynamic_cast<const ROL::PartitionedVector<Real>&>(x).get(ind);
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
      get(x,i)->set(*get(c,i));
    }
    problem0_->getBoundConstraint()->project(x);
  }

  Real infeasibility(ROL::Vector<Real> &x) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    zeroSlack(x);
    con_->update(x,ROL::UpdateType::Temp);
    con_->value(*c_,x,tol);
    setSlack(x,*c_);
    x_->set(x);
    proj_->project(*x_);
    x_->axpy(static_cast<Real>(-1),x);
    Real infeas = x_->norm();
    return infeas;
  }

public:
  FacilityLocationBranchSub(const ROL::Ptr<ROL::ParameterList> &parlist,
                            const ROL::Ptr<ROL::PEBBL::Branching<Real>> &branching,
                            int verbosity = 0,
                            const ROL::Ptr<std::ostream> &outStream = ROL::nullPtr)
    : ROL::PEBBL::BranchSub<Real>(branching, verbosity, outStream) {
    ctol_ = parlist->sublist("Status Test").get("Constraint Tolerance",1e-8);
    x_    = solution_->clone();
    bool useLinearCon = parlist->sublist("Problem").get("Maintain Linear Constraints",true);
    if (useLinearCon) {
      proj_ = problem0_->getPolyhedralProjection();
      c_    = proj_->getResidual()->clone();
      con_  = proj_->getLinearConstraint();
    }
    else {
      c_    = problem0_->getResidualVector()->clone();
      con_  = problem0_->getConstraint();
      proj_ = ROL::PolyhedralProjectionFactory<Real>(
                       *problem0_->getPrimalOptimizationVector(),
                       *problem0_->getDualOptimizationVector(),
                       problem0_->getBoundConstraint(),
                       con_,*problem0_->getMultiplierVector(),*c_,*parlist);
    }
  }

  FacilityLocationBranchSub(const FacilityLocationBranchSub &rpbs)
    : ROL::PEBBL::BranchSub<Real>(rpbs), ctol_(rpbs.ctol_), x_(rpbs.x_->clone()),
      c_(rpbs.c_->clone()), con_(rpbs.con_), proj_(rpbs.proj_) {}

  void incumbentHeuristic() {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    const Real zero(0), one(1);
    // Max Rounding
    //rndSolution_->set(*solution_);
    //ROL::Ptr<std::vector<Real>>       idata = getData(*rndSolution_);
    //ROL::Ptr<const std::vector<Real>> cdata = getConstData(*rndSolution_);
    //const int M = idata->size();
    //const int N = cdata->size()/M;
    //Real maxj(0);
    //for (int i = 0; i < M; ++i) {
    //  maxj = (*cdata)[i];
    //  for (int j = 1; j < N; ++j) {
    //    maxj = std::max(maxj,(*cdata)[i + j*M]);
    //  }
    //  (*idata)[i] = (maxj > tol ? one : zero);
    //}
    // Real cnorm = infeasibility(*rndSolution_);
    // Randomized Rounding
    ROL::Ptr<std::vector<Real>> idata = getData(*rndSolution_,0);
    ROL::Ptr<std::vector<Real>> cdata = getData(*rndSolution_,1);
    const int M = idata->size();
    const int N = cdata->size()/M;
    Real z(0), r(0), cnorm(0), sum(0);
    int cnt(0);
    while (true) {
      rndSolution_->set(*solution_);
      std::vector<int> ind;
      for (int i = 0; i < M; ++i) {
        r = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
        z = (*idata)[i];
        if (r <= z) {
          (*idata)[i] = zero;
          for (int j = 0; j < N; ++j) (*cdata)[i + j*M] = zero;
        }
        else {
          (*idata)[i] = one;
          ind.push_back(i);
        }
      }
      for (int j = 0; j < N; ++j) {
        sum = zero;
        for (const auto i : ind) sum += (*cdata)[i + j*M];
        for (const auto i : ind) (*cdata)[i + j*M] /= sum;
      }
      cnorm = infeasibility(*rndSolution_);
      cnt++;
      if (cnorm < ctol_) break;
      if (verbosity_ > 1) {
        *outStream_ << "  cnt = " << cnt << "  infeasibility = " << cnorm << std::endl;
      }
    }

    problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
    Real val = problem0_->getObjective()->value(*rndSolution_,tol);
    branching_->foundSolution(new ROL::PEBBL::IntegerSolution<Real>(*rndSolution_,val));
    if (verbosity_ > 0) {
      *outStream_ << "FacilityLocationBranchSub::incumbentHeuristic" << std::endl;
      *outStream_ << "  Incumbent Value:         " << val  << std::endl;
      *outStream_ << "  Incumbent Infeasibility: " << cnorm << std::endl;
    }
  }

  pebbl::branchSub* makeChild(int whichChild = anyChild) override {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild,std::logic_error,
      ">>> FacilityLocationBranchSub::makeChild: whichChild is equal to anyChild!");
    FacilityLocationBranchSub<Real>* child
      = new FacilityLocationBranchSub<Real>(*this);
    child->updateFixed(index_,
      (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

}; // class FacilityLocationBranchSub

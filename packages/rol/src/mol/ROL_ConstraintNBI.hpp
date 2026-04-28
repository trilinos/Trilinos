
#ifndef ROL_CONSTRAINTNBI_HPP
#define ROL_CONSTRAINTNBI_HPP

#include "ROL_Constraint.hpp"
#include "ROL_Objective.hpp"
#include "ROL_StdVector.hpp"

namespace ROL {

template<typename Real>
class ConstraintNBI : public Constraint<Real> {
private:
  const std::vector<Ptr<Objective<Real>>> obj_;
  const unsigned nobj_;
  std::vector<Real> Pbeta_, normal_;

  // Need wdual_
  Ptr<Vector<Real>> wdual_;
  bool isinit_;

  void initialize(const Vector<Real>& w) {
    if (!isinit_) {
      wdual_ = w.dual().clone();
      isinit_ = true;
    }
  }

public:
  ConstraintNBI(const std::vector<Ptr<Objective<Real>>>& obj,
                const std::vector<Real>& beta,
                const std::vector<std::vector<Real>>& values)
    : obj_(obj), nobj_(obj.size()), isinit_(false) {
    Pbeta_.clear(); Pbeta_.resize(nobj_,static_cast<Real>(0));
    normal_.clear(); normal_.resize(nobj_,static_cast<Real>(0));
    Real Pij(0);
    for (unsigned i = 0u; i < nobj_; ++i) {
      for (unsigned j = 0u; j < nobj_; ++j) {
        Pij = (values[j][i] - values[i][i]);
        Pbeta_[i] += Pij * beta[j];
        normal_[i] += Pij;
      }
    }
  }

  void update(const ROL::Vector<Real> &x, ROL::UpdateType type, int iter = -1) {
    for (const auto& oi : obj_) oi->update(x,type,iter);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &x, Real &tol) {
    auto cs = static_cast<ROL::StdVector<Real>&>(c).getVector();
    Real F0 = (obj_[0]->value(x,tol)-Pbeta_[0])/normal_[0];
    for (unsigned i = 1u; i < nobj_; ++i)
      (*cs)[i-1] = (obj_[i]->value(x,tol)-Pbeta_[i])/normal_[i] - F0;
  }

  void applyJacobian(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    auto js = static_cast<ROL::StdVector<Real>&>(jv).getVector();
    initialize(x);
    obj_[0]->gradient(*wdual_,x,tol);
    Real F0 = wdual_->apply(v)/normal_[0];
    for (unsigned i = 1u; i < nobj_; ++i) {
      obj_[i]->gradient(*wdual_,x,tol);
      (*js)[i-1] = wdual_->apply(v)/normal_[i] - F0;
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    const auto vs = static_cast<const ROL::StdVector<Real>&>(v).getVector();
    initialize(x);
    ajv.zero();
    for (unsigned i = 1u; i < nobj_; ++i) {
      obj_[i]->gradient(*wdual_,x,tol);
      ajv.axpy((*vs)[i-1]/normal_[i],*wdual_);
    }
    obj_[0]->gradient(*wdual_,x,tol);
    for (unsigned i = 0u; i < nobj_-1; ++i) {
      ajv.axpy(-(*vs)[i]/normal_[0],*wdual_);
    }
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahuv, const ROL::Vector<Real> &u, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    const auto us = static_cast<const ROL::StdVector<Real>&>(u).getVector();
    initialize(x);
    ahuv.zero();
    for (unsigned i = 1u; i < nobj_; ++i) {
      obj_[i]->hessVec(*wdual_,v,x,tol);
      ahuv.axpy((*us)[i-1]/normal_[i],*wdual_);
    }
    obj_[0]->hessVec(*wdual_,v,x,tol);
    for (unsigned i = 0u; i < nobj_-1; ++i) {
      ahuv.axpy(-(*us)[i]/normal_[0],*wdual_);
    }
  }
};

}

#endif

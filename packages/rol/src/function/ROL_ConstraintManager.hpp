// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_MANAGER_H
#define ROL_CONSTRAINT_MANAGER_H

#include "ROL_Constraint.hpp"
#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_BoundConstraint_Partitioned.hpp"

/** @ingroup func_group
    \class ROL::ConstraintManager
    \brief Provides a wrapper for multiple constraints.

    ---
*/

namespace ROL {

template <class Real>
class ConstraintManager {
private:
  ROL::Ptr<Constraint<Real> >      con_;
  ROL::Ptr<Vector<Real> >          l_;
  ROL::Ptr<Vector<Real> >          x_;
  ROL::Ptr<BoundConstraint<Real> > bnd_;

  std::vector<ROL::Ptr<Constraint<Real> > >      cvec_;
  std::vector<ROL::Ptr<Vector<Real> > >          lvec_;
  std::vector<ROL::Ptr<Vector<Real> > >          svec_;
  std::vector<ROL::Ptr<BoundConstraint<Real> > > sbnd_;

  std::vector<bool> isInequality_;

  bool isNull_;
  bool hasInequality_;

  void initializeSlackVariable(const ROL::Ptr<Constraint<Real> >      &con,
                               const ROL::Ptr<BoundConstraint<Real> > &cbnd,
                               const ROL::Ptr<Vector<Real> >          &s,
                               const ROL::Ptr<Vector<Real> >          &x) const {
    // Set slack variable to s = proj(c(x))
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    con->value(*s,*x,tol);
    cbnd->project(*s);
  }

  void initialize(const std::vector<ROL::Ptr<Constraint<Real> > >      &cvec,
                  const std::vector<ROL::Ptr<Vector<Real> > >          &lvec,
                  const std::vector<ROL::Ptr<BoundConstraint<Real> > > &bvec,
                  const ROL::Ptr<Vector<Real> >                        &x,
                  const ROL::Ptr<BoundConstraint<Real> >               &bnd) {
    // Check size of multiplier vector and constraint vector
    int size = static_cast<int>(cvec.size());
    if ( size != static_cast<int>(lvec.size()) ) {
      throw Exception::NotImplemented(">>> ROL::ConstraintManager: Constraint and multiplier vectors are different sizes!");
    }
    if ( size != static_cast<int>(bvec.size()) ) {
      throw Exception::NotImplemented(">>> ROL::ConstraintManager: Constraint and BoundConstraint vectors are different sizes!");
    }
    // If bnd is null, then make a null BoundConstraint
    ROL::Ptr<BoundConstraint<Real> > bnd0;
    if ( bnd == ROL::nullPtr ) {
      bnd0 = ROL::makePtr<BoundConstraint<Real>>(*x);
      bnd0->deactivate();
    }
    else {
      bnd0 = bnd;
    }
    // Build slack variables
    svec_.clear(); svec_.push_back(x);
    sbnd_.clear(); sbnd_.push_back(bnd0);
    cvec_.clear(); lvec_.clear(); isInequality_.clear();
    int cnt = 1, cnt_con = 0;
    isNull_ = true;
    hasInequality_ = false;
    for (int i = 0; i < size; ++i) {
      ROL::Ptr<Constraint<Real> >      con  = cvec[i];
      ROL::Ptr<Vector<Real> >          l    = lvec[i];
      ROL::Ptr<BoundConstraint<Real> > cbnd = bvec[i];
      if (con != ROL::nullPtr) {
        if ( con->isActivated() ) {
          // Set default type to equality
          isInequality_.push_back(false);
          // Fill constraint and multiplier vectors
          cvec_.push_back(con);
          lvec_.push_back(l);
          if (cbnd != ROL::nullPtr) {
            if ( cbnd->isActivated() ) {
              // Set type to inequality
              isInequality_.back() = true;
              // Create slack variables
              svec_.push_back(l->dual().clone());
              initializeSlackVariable(con,cbnd,svec_[cnt],x);
              // Create slack bound
              sbnd_.push_back(cbnd);
              // Update inequality constraint counter
              cnt++;
              hasInequality_ = true;
            }
          }
          cnt_con++;
          isNull_ = false;
        }
      }
    }
    // Create partitioned constraint and multiplier vector
    if ( !isNull_ ) {
      if ( cnt_con > 1 || hasInequality_ ) {
        con_ = ROL::makePtr<Constraint_Partitioned<Real>>(cvec_,isInequality_);
        l_   = ROL::makePtr<PartitionedVector<Real>>(lvec_);
      }
      else {
        con_ = cvec_[0];
        l_   = lvec_[0];
      }
    }
    else {
      con_ = ROL::nullPtr;
      l_   = ROL::nullPtr;
    }
    // Create partitioned optimization vector and bound constraint
    if ( hasInequality_ ) {
      x_   = ROL::makePtr<PartitionedVector<Real>>(svec_);
      bnd_ = ROL::makePtr<BoundConstraint_Partitioned<Real>>(sbnd_,svec_);
    }
    else {
      x_   = x;
      bnd_ = bnd0;
    }
  }

public:
  virtual ~ConstraintManager(void) {}

  ConstraintManager(const std::vector<ROL::Ptr<Constraint<Real> > >      &cvec,
                    const std::vector<ROL::Ptr<Vector<Real> > >          &lvec,
                    const std::vector<ROL::Ptr<BoundConstraint<Real> > > &bvec,
                    const ROL::Ptr<Vector<Real> >                        &x,
                    const ROL::Ptr<BoundConstraint<Real> >               &bnd = ROL::nullPtr)
    : isNull_(true), hasInequality_(false) {
    initialize(cvec,lvec,bvec,x,bnd);
  }

  ConstraintManager(const std::vector<ROL::Ptr<Constraint<Real> > >      &cvec,
                    const std::vector<ROL::Ptr<Vector<Real> > >          &lvec,
                    const ROL::Ptr<Vector<Real> >                        &x,
                    const ROL::Ptr<BoundConstraint<Real> >               &bnd = ROL::nullPtr)
    : isNull_(true), hasInequality_(false) {
    std::vector<ROL::Ptr<BoundConstraint<Real> > > bvec(cvec.size(),ROL::nullPtr);
    initialize(cvec,lvec,bvec,x,bnd);
  }

  ConstraintManager(const ROL::Ptr<Constraint<Real> >                    &con,
                    const ROL::Ptr<Vector<Real> >                        &l,
                    const ROL::Ptr<BoundConstraint<Real> >               &cbnd,
                    const ROL::Ptr<Vector<Real> >                        &x,
                    const ROL::Ptr<BoundConstraint<Real> >               &bnd = ROL::nullPtr)
    : isNull_(true), hasInequality_(false) {
    std::vector<ROL::Ptr<Constraint<Real> > >      cvec(1,con);
    std::vector<ROL::Ptr<Vector<Real> > >          lvec(1,l);
    std::vector<ROL::Ptr<BoundConstraint<Real> > > bvec(1,cbnd);
    initialize(cvec,lvec,bvec,x,bnd);
  }

  ConstraintManager(const ROL::Ptr<Constraint<Real> >                    &con,
                    const ROL::Ptr<Vector<Real> >                        &l,
                    const ROL::Ptr<Vector<Real> >                        &x,
                    const ROL::Ptr<BoundConstraint<Real> >               &bnd = ROL::nullPtr)
    : isNull_(true), hasInequality_(false) {
    std::vector<ROL::Ptr<Constraint<Real> > >      cvec(1,con);
    std::vector<ROL::Ptr<Vector<Real> > >          lvec(1,l);
    std::vector<ROL::Ptr<BoundConstraint<Real> > > bvec(1,ROL::nullPtr);
    initialize(cvec,lvec,bvec,x,bnd);
  }

  const ROL::Ptr<Constraint<Real> > getConstraint(void) const {
    return con_;
  }

  const ROL::Ptr<Vector<Real> > getMultiplier(void) const {
    return l_;
  }

  const ROL::Ptr<Vector<Real> > getOptVector(void) const {
    return x_;
  }

  const ROL::Ptr<BoundConstraint<Real> > getBoundConstraint(void) const {
    return bnd_;
  }

  bool isNull(void) const {
    return isNull_;
  }

  bool hasInequality(void) const {
    return hasInequality_;
  }

  void resetSlackVariables(void) {
    if (hasInequality_) {
      int size = static_cast<int>(cvec_.size());
      int cnt = 1;
      for (int i = 0; i < size; ++i) {
        if (isInequality_[i]) {
          // Reset slack variables
          initializeSlackVariable(cvec_[i],sbnd_[cnt],svec_[cnt],svec_[0]);
          cnt++;
        }
      }
    }
  }

}; // class ConstraintManager

} // namespace ROL

#endif

// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_ASSEMBLER_DEF_H
#define ROL_CONSTRAINT_ASSEMBLER_DEF_H

namespace ROL {

template<typename Real>
void ConstraintAssembler<Real>::initializeSlackVariable(const Ptr<Constraint<Real>>      &con,
                                                         const Ptr<BoundConstraint<Real>> &cbnd,
                                                         const Ptr<Vector<Real>>          &s,
                                                         const Ptr<Vector<Real>>          &x) const {
  // Set slack variable to s = proj(c(x))
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  con->value(*s,*x,tol);
  cbnd->project(*s);
}

template<typename Real>
void ConstraintAssembler<Real>::initialize( const std::unordered_map<std::string,ConstraintData<Real>> &input_con,
                                            const Ptr<Vector<Real>>                                    &xprim,
                                            const Ptr<Vector<Real>>                                    &xdual,
                                            const Ptr<BoundConstraint<Real>>                           &bnd) {
  // If bnd is null, then make a null BoundConstraint
  Ptr<BoundConstraint<Real>> bnd0;
  if ( bnd == nullPtr ) {
    bnd0 = makePtr<BoundConstraint<Real>>(*xprim);
    bnd0->deactivate();
  }
  else {
    bnd0 = bnd;
  }
  // Build slack variables
  psvec_.clear(); psvec_.push_back(xprim);
  dsvec_.clear(); dsvec_.push_back(xdual);
  sbnd_.clear(); sbnd_.push_back(bnd0);
  cvec_.clear(); lvec_.clear(); rvec_.clear();
  isInequality_.clear();
  int cnt = 1, cnt_con = 0;
  isNull_ = true;
  hasInequality_ = false;
  for (auto it = input_con.begin(); it != input_con.end(); ++it) {
    Ptr<Constraint<Real>>      con  = it->second.constraint;
    Ptr<Vector<Real>>          l    = it->second.multiplier;
    Ptr<Vector<Real>>          r    = it->second.residual;
    Ptr<BoundConstraint<Real>> cbnd = it->second.bounds;
    if (con != nullPtr) {
      if ( con->isActivated() ) {
        // Set default type to equality
        isInequality_.push_back(false);
        // Fill constraint and multiplier vectors
        cvec_.push_back(con);
        lvec_.push_back(l);
        rvec_.push_back(r);
        if (cbnd != nullPtr) {
          if ( cbnd->isActivated() ) {
            // Set type to inequality
            isInequality_.back() = true;
            // Create slack variables
            psvec_.push_back(r->clone());
            dsvec_.push_back(l->clone());
            initializeSlackVariable(con,cbnd,psvec_[cnt],xprim);
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
      con_ = makePtr<Constraint_Partitioned<Real>>(cvec_,isInequality_);
      mul_ = makePtr<PartitionedVector<Real>>(lvec_);
      res_ = makePtr<PartitionedVector<Real>>(rvec_);
    }
    else {
      con_ = cvec_[0];
      mul_ = lvec_[0];
      res_ = rvec_[0];
    }
  }
  else {
    con_ = nullPtr;
    mul_ = nullPtr;
    res_ = nullPtr;
  }
  // Create partitioned optimization vector and bound constraint
  if ( hasInequality_ ) {
    xprim_ = makePtr<PartitionedVector<Real>>(psvec_);
    xdual_ = makePtr<PartitionedVector<Real>>(dsvec_);
    bnd_   = makePtr<BoundConstraint_Partitioned<Real>>(sbnd_,psvec_);
  }
  else {
    xprim_ = xprim;
    xdual_ = xdual;
    bnd_   = bnd0;
  }
}

template<typename Real>
void ConstraintAssembler<Real>::initialize( const std::unordered_map<std::string,ConstraintData<Real>> &input_con,
                                            const std::unordered_map<std::string,ConstraintData<Real>> &input_lcon,
                                            const Ptr<Vector<Real>>                                    &xprim,
                                            const Ptr<Vector<Real>>                                    &xdual,
                                            const Ptr<BoundConstraint<Real>>                           &bnd) {
  // If bnd is null, then make a null BoundConstraint
  Ptr<BoundConstraint<Real>> bnd0;
  if ( bnd == nullPtr ) {
    bnd0 = makePtr<BoundConstraint<Real>>(*xprim);
    bnd0->deactivate();
  }
  else {
    bnd0 = bnd;
  }
  // Build slack variables
  psvec_.clear(); psvec_.push_back(xprim);
  dsvec_.clear(); dsvec_.push_back(xdual);
  sbnd_.clear();  sbnd_.push_back(bnd0);
  cvec_.clear();  lvec_.clear();  rvec_.clear();
  lcvec_.clear(); llvec_.clear(); lrvec_.clear();
  isInequality_.clear(); isLinearInequality_.clear();
  int cnt = 1, cnt_con = 0, cnt_lcon = 0;
  isNull_ = true;
  hasInequality_ = false;
  for (auto it = input_con.begin(); it != input_con.end(); ++it) {
    Ptr<Constraint<Real>>      con  = it->second.constraint;
    Ptr<Vector<Real>>          l    = it->second.multiplier;
    Ptr<Vector<Real>>          r    = it->second.residual;
    Ptr<BoundConstraint<Real>> cbnd = it->second.bounds;
    if (con != nullPtr) {
      if ( con->isActivated() ) {
        // Set default type to equality
        isInequality_.push_back(false);
        // Fill constraint and multiplier vectors
        cvec_.push_back(con);
        lvec_.push_back(l);
        rvec_.push_back(r);
        if (cbnd != nullPtr) {
          if ( cbnd->isActivated() ) {
            // Set type to inequality
            isInequality_.back() = true;
            // Create slack variables
            psvec_.push_back(r->clone());
            dsvec_.push_back(l->clone());
            initializeSlackVariable(con,cbnd,psvec_[cnt],xprim);
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
  for (auto it = input_lcon.begin(); it != input_lcon.end(); ++it) {
    Ptr<Constraint<Real>>      con  = it->second.constraint;
    Ptr<Vector<Real>>          l    = it->second.multiplier;
    Ptr<Vector<Real>>          r    = it->second.residual;
    Ptr<BoundConstraint<Real>> cbnd = it->second.bounds;
    if (con != nullPtr) {
      if ( con->isActivated() ) {
        // Set default type to equality
        isLinearInequality_.push_back(false);
        // Fill constraint and multiplier vectors
        lcvec_.push_back(con);
        llvec_.push_back(l);
        lrvec_.push_back(r);
        if (cbnd != nullPtr) {
          if ( cbnd->isActivated() ) {
            // Set type to inequality
            isLinearInequality_.back() = true;
            // Create slack variables
            psvec_.push_back(r->clone());
            dsvec_.push_back(l->clone());
            initializeSlackVariable(con,cbnd,psvec_[cnt],xprim);
            // Create slack bound
            sbnd_.push_back(cbnd);
            // Update inequality constraint counter
            cnt++;
            hasInequality_ = true;
          }
        }
        cnt_lcon++;
        isNull_ = false;
      }
    }
  }
  // Create partitioned constraint and multiplier vector
  if ( !isNull_ ) {
    if ( cnt_con > 1 || hasInequality_ ) {
      con_  = makePtr<Constraint_Partitioned<Real>>(cvec_,isInequality_);
      mul_  = makePtr<PartitionedVector<Real>>(lvec_);
      res_  = makePtr<PartitionedVector<Real>>(rvec_);
    }
    else {
      con_  = cvec_[0];
      mul_  = lvec_[0];
      res_  = rvec_[0];
    }
    if ( cnt_lcon > 1 || hasInequality_ ) {
      linear_con_ = makePtr<Constraint_Partitioned<Real>>(lcvec_,isLinearInequality_,cnt_con);
      linear_mul_ = makePtr<PartitionedVector<Real>>(llvec_);
      linear_res_ = makePtr<PartitionedVector<Real>>(lrvec_);
    }
    else {
      linear_con_ = lcvec_[0];
      linear_mul_ = llvec_[0];
      linear_res_ = lrvec_[0];
    }
  }
  else {
    con_  = nullPtr;
    mul_  = nullPtr;
    res_  = nullPtr;
    linear_con_ = nullPtr;
    linear_mul_ = nullPtr;
    linear_res_ = nullPtr;
  }
  // Create partitioned optimization vector and bound constraint
  if ( hasInequality_ ) {
    xprim_ = makePtr<PartitionedVector<Real>>(psvec_);
    xdual_ = makePtr<PartitionedVector<Real>>(dsvec_);
    bnd_   = makePtr<BoundConstraint_Partitioned<Real>>(sbnd_,psvec_);
  }
  else {
    xprim_ = xprim;
    xdual_ = xdual;
    bnd_   = bnd0;
  }
}

template<typename Real>
ConstraintAssembler<Real>::ConstraintAssembler( const std::unordered_map<std::string,ConstraintData<Real>> &con,
                                                const Ptr<Vector<Real>>                                    &xprim,
                                                const Ptr<Vector<Real>>                                    &xdual,
                                                const Ptr<BoundConstraint<Real>>                           &bnd)
  : isNull_(true), hasInequality_(false) {
  initialize(con,xprim,xdual,bnd);
}

template<typename Real>
ConstraintAssembler<Real>::ConstraintAssembler( const std::unordered_map<std::string,ConstraintData<Real>> &con,
                                                const std::unordered_map<std::string,ConstraintData<Real>> &linear_con,
                                                const Ptr<Vector<Real>>                                    &xprim,
                                                const Ptr<Vector<Real>>                                    &xdual,
                                                const Ptr<BoundConstraint<Real>>                           &bnd)
  : isNull_(true), hasInequality_(false) {
  initialize(con,linear_con,xprim,xdual,bnd);
}

template<typename Real>
const Ptr<Constraint<Real>>& ConstraintAssembler<Real>::getConstraint() const {
  return con_;
}

template<typename Real>
const Ptr<Vector<Real>>& ConstraintAssembler<Real>::getMultiplier() const {
  return mul_;
}

template<typename Real>
const Ptr<Vector<Real>>& ConstraintAssembler<Real>::getResidual() const {
  return res_;
}

template<typename Real>
const Ptr<Constraint<Real>>& ConstraintAssembler<Real>::getLinearConstraint() const {
  return linear_con_;
}

template<typename Real>
const Ptr<Vector<Real>>& ConstraintAssembler<Real>::getLinearMultiplier() const {
  return linear_mul_;
}

template<typename Real>
const Ptr<Vector<Real>>& ConstraintAssembler<Real>::getLinearResidual() const {
  return linear_res_;
}

template<typename Real>
const Ptr<Vector<Real>>& ConstraintAssembler<Real>::getOptVector() const {
  return xprim_;
}

template<typename Real>
const Ptr<Vector<Real>>& ConstraintAssembler<Real>::getDualOptVector() const {
  return xdual_;
}

template<typename Real>
const Ptr<BoundConstraint<Real>>& ConstraintAssembler<Real>::getBoundConstraint() const {
  return bnd_;
}

template<typename Real>
bool ConstraintAssembler<Real>::isNull() const {
  return isNull_;
}

template<typename Real>
bool ConstraintAssembler<Real>::hasInequality() const {
  return hasInequality_;
}

template<typename Real>
void ConstraintAssembler<Real>::resetSlackVariables() {
  if (hasInequality_) {
    int size = static_cast<int>(cvec_.size());
    int cnt = 1;
    for (int i = 0; i < size; ++i) {
      if (isInequality_[i]) {
        // Reset slack variables
        initializeSlackVariable(cvec_[i],sbnd_[cnt],psvec_[cnt],psvec_[0]);
        cnt++;
      }
    }
    size = static_cast<int>(lcvec_.size());
    for (int i = 0; i < size; ++i) {
      if (isLinearInequality_[i]) {
        // Reset slack variables
        initializeSlackVariable(lcvec_[i],sbnd_[cnt],psvec_[cnt],psvec_[0]);
        cnt++;
      }
    }
  }
}

} // namespace ROL

#endif

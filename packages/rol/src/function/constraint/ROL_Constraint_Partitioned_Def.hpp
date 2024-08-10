// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_PARTITIONED_DEF_H
#define ROL_CONSTRAINT_PARTITIONED_DEF_H

namespace ROL {

template<typename Real>
Constraint_Partitioned<Real>::Constraint_Partitioned(const std::vector<Ptr<Constraint<Real>>> &cvec,
                                                     bool isInequality,
                                                     int offset)
 : cvec_(cvec), offset_(offset),
   scratch_(nullPtr), ncval_(0), initialized_(false) {
  isInequality_.clear(); isInequality_.resize(cvec.size(),isInequality);
}

template<typename Real>
Constraint_Partitioned<Real>::Constraint_Partitioned(const std::vector<Ptr<Constraint<Real>>> &cvec,
                                                     std::vector<bool> isInequality,
                                                     int offset)
 : cvec_(cvec), isInequality_(isInequality), offset_(offset),
   scratch_(nullPtr), ncval_(0), initialized_(false) {}

template<typename Real>
int Constraint_Partitioned<Real>::getNumberConstraintEvaluations(void) const {
  return ncval_;
}

template<typename Real>
Ptr<Constraint<Real>> Constraint_Partitioned<Real>::get(int ind) const {
  if (ind < 0 || ind > static_cast<int>(cvec_.size())) {
    throw Exception::NotImplemented(">>> Constraint_Partitioned::get : Index out of bounds!");
  }
  return cvec_[ind];
}

template<typename Real>
void Constraint_Partitioned<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  const int ncon = static_cast<int>(cvec_.size());
  for (int i = 0; i < ncon; ++i) {
    cvec_[i]->update(getOpt(x),type,iter);
  }
}

template<typename Real>
void Constraint_Partitioned<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  const int ncon = static_cast<int>(cvec_.size());
  for (int i = 0; i < ncon; ++i) {
    cvec_[i]->update(getOpt(x),flag,iter);
  }
}

template<typename Real>
void Constraint_Partitioned<Real>::value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
  PartitionedVector<Real> &cpv
    = dynamic_cast<PartitionedVector<Real>&>(c);

  const int ncon = static_cast<int>(cvec_.size());
  int cnt = offset_+1;
  for (int i = 0; i < ncon; ++i) {
    cvec_[i]->value(*cpv.get(i), getOpt(x), tol);
    if (isInequality_[i]) {
      cpv.get(i)->axpy(static_cast<Real>(-1),getSlack(x,cnt));
      cnt++;
    }
  }
  ++ncval_;
}

template<typename Real>
void Constraint_Partitioned<Real>::applyJacobian( Vector<Real> &jv,
                                            const Vector<Real> &v,
                                            const Vector<Real> &x,
                                                  Real &tol ) {
  PartitionedVector<Real> &jvpv
    = dynamic_cast<PartitionedVector<Real>&>(jv);

  const int ncon = static_cast<int>(cvec_.size());
  int cnt = offset_+1;
  for (int i = 0; i < ncon; ++i) {
    cvec_[i]->applyJacobian(*jvpv.get(i), getOpt(v), getOpt(x), tol);
    if (isInequality_[i]) {
      jvpv.get(i)->axpy(static_cast<Real>(-1),getSlack(v,cnt));
      cnt++;
    }
  }
}

template<typename Real>
void Constraint_Partitioned<Real>::applyAdjointJacobian( Vector<Real> &ajv,
                                                   const Vector<Real> &v,
                                                   const Vector<Real> &x,
                                                         Real &tol ) {
  if (!initialized_) {
    scratch_ = getOpt(ajv).clone();
    initialized_ = true;
  }

  const PartitionedVector<Real> &vpv
    = dynamic_cast<const PartitionedVector<Real>&>(v);

  const int ncon = static_cast<int>(cvec_.size());
  int cnt = offset_+1;
  getOpt(ajv).zero();
  for (int i = 0; i < ncon; ++i) {
    scratch_->zero();
    cvec_[i]->applyAdjointJacobian(*scratch_, *vpv.get(i), getOpt(x), tol);
    getOpt(ajv).plus(*scratch_);
    if (isInequality_[i]) {
      getSlack(ajv,cnt).set(*vpv.get(i));
      getSlack(ajv,cnt).scale(static_cast<Real>(-1));
      cnt++;
    }
  }
}

template<typename Real>
void Constraint_Partitioned<Real>::applyAdjointHessian( Vector<Real> &ahuv,
                                                  const Vector<Real> &u,
                                                  const Vector<Real> &v,
                                                  const Vector<Real> &x,
                                                        Real &tol ) {
  if (!initialized_) {
    scratch_ = getOpt(ahuv).clone();
    initialized_ = true;
  }

  const PartitionedVector<Real> &upv
    = dynamic_cast<const PartitionedVector<Real>&>(u);

  const int ncon = static_cast<int>(cvec_.size());
  int cnt = offset_+1;
  getOpt(ahuv).zero();
  for (int i = 0; i < ncon; ++i) {
    Ptr<const Vector<Real>> ui = upv.get(i);
    scratch_->zero();
    cvec_[i]->applyAdjointHessian(*scratch_, *ui, getOpt(v), getOpt(x), tol);
    getOpt(ahuv).plus(*scratch_);
    if (isInequality_[i]) {
      getSlack(ahuv,cnt).zero();
      cnt++;
    }
  }
}

template<typename Real>
void Constraint_Partitioned<Real>::applyPreconditioner(Vector<Real> &pv,
                                                 const Vector<Real> &v,
                                                 const Vector<Real> &x,
                                                 const Vector<Real> &g,
                                                       Real &tol) {
  PartitionedVector<Real> &pvpv
    = dynamic_cast<PartitionedVector<Real>&>(pv);
  const PartitionedVector<Real> &vpv
    = dynamic_cast<const PartitionedVector<Real>&>(v);
  
  const int ncon = static_cast<int>(cvec_.size());
  for (int i = 0; i < ncon; ++i) {
    cvec_[i]->applyPreconditioner(*pvpv.get(i), *vpv.get(i), getOpt(x), getOpt(g), tol);
  }
}

template<typename Real>
void Constraint_Partitioned<Real>::setParameter(const std::vector<Real> &param) {
  Constraint<Real>::setParameter(param);
  const int ncon = static_cast<int>(cvec_.size());
  for (int i = 0; i < ncon; ++i) {
    cvec_[i]->setParameter(param);
  }
}

template<typename Real>
Vector<Real>& Constraint_Partitioned<Real>::getOpt( Vector<Real> &xs ) const {
  try {
    return *dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
  }
  catch (std::exception &e) {
    return xs;
  }
}

template<typename Real>
const Vector<Real>& Constraint_Partitioned<Real>::getOpt( const Vector<Real> &xs ) const {
  try {
    return *dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
  }
  catch (std::exception &e) {
    return xs;
  }
}

template<typename Real>
Vector<Real>& Constraint_Partitioned<Real>::getSlack( Vector<Real> &xs, int ind ) const {
  return *dynamic_cast<PartitionedVector<Real>&>(xs).get(ind);
}

template<typename Real>
const Vector<Real>& Constraint_Partitioned<Real>::getSlack( const Vector<Real> &xs, int ind ) const {
  return *dynamic_cast<const PartitionedVector<Real>&>(xs).get(ind);
}

} // namespace ROL

#endif

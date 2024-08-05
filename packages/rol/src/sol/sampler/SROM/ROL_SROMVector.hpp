// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SROMVECTOR_H
#define ROL_SROMVECTOR_H

#include <algorithm>
#include <cstdlib>

#include "ROL_ProbabilityVector.hpp"
#include "ROL_AtomVector.hpp"

/** \class ROL::SROMVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real>
class SROMVector : public Vector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const ROL::Ptr<ProbabilityVector<Real> > pvec_;
  const ROL::Ptr<AtomVector<Real> > avec_;

  mutable ROL::Ptr<Vector<Real> > dual_pvec_;
  mutable ROL::Ptr<Vector<Real> > dual_avec_;
  mutable ROL::Ptr<SROMVector<Real> > dual_vec_;
  mutable bool isDualInitialized_;

public:

  SROMVector(const ROL::Ptr<ProbabilityVector<Real> >  &pvec,
             const ROL::Ptr<AtomVector<Real> >         &avec)
    : pvec_(pvec), avec_(avec), isDualInitialized_(false) {
    dual_pvec_ = (pvec_->dual()).clone();
    dual_avec_ = (avec_->dual()).clone();
  }

  void set( const Vector<Real> &x ) {
    const SROMVector &ex = dynamic_cast<const SROMVector&>(x);
    pvec_->set(*(ex.getProbabilityVector()));
    avec_->set(*(ex.getAtomVector()));
  }

  void plus( const Vector<Real> &x ) {
    const SROMVector &ex = dynamic_cast<const SROMVector&>(x);
    pvec_->plus(*(ex.getProbabilityVector()));
    avec_->plus(*(ex.getAtomVector()));
  }

  void scale( const Real alpha ) {
    pvec_->scale(alpha);
    avec_->scale(alpha);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const SROMVector &ex = dynamic_cast<const SROMVector&>(x);
    pvec_->axpy(alpha,*(ex.getProbabilityVector()));
    avec_->axpy(alpha,*(ex.getAtomVector()));
  }

  Real dot( const Vector<Real> &x ) const {
    const SROMVector & ex = dynamic_cast<const SROMVector&>(x);
    Real pval = pvec_->dot(*(ex.getProbabilityVector()));
    Real aval = avec_->dot(*(ex.getAtomVector()));
    return pval + aval;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  ROL::Ptr<Vector<Real> > clone(void) const {
    return ROL::makePtr<SROMVector>(
           ROL::staticPtrCast<ProbabilityVector<Real> >(pvec_->clone()),
           ROL::staticPtrCast<AtomVector<Real> >(avec_->clone()) );
  }

  const Vector<Real> & dual(void) const {
    if ( !isDualInitialized_ ) {
      dual_vec_ = ROL::makePtr<SROMVector>(
                  ROL::staticPtrCast<ProbabilityVector<Real> >(dual_pvec_),
                  ROL::staticPtrCast<AtomVector<Real> >(dual_avec_) );
      isDualInitialized_ = true;
    }
    dual_pvec_->set(pvec_->dual());
    dual_avec_->set(avec_->dual());
    return *dual_vec_;
  }

  int dimension(void) const {
    return avec_->dimension() + pvec_->dimension();
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    pvec_->applyUnary(f);
    avec_->applyUnary(f);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const SROMVector & ex = dynamic_cast<const SROMVector&>(x);
    pvec_->applyBinary(f,*(ex.getProbabilityVector()));
    avec_->applyBinary(f,*(ex.getAtomVector()));
  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    Real pval = pvec_->reduce(r);
    Real aval = avec_->reduce(r);
    r.reduce(pval,result);
    r.reduce(aval,result);
    return result;
  }

  const ROL::Ptr<const AtomVector<Real> > getAtomVector(void) const {
    return avec_;
  }
  
  const ROL::Ptr<const ProbabilityVector<Real> > getProbabilityVector(void) const {
    return pvec_;
  }

  ROL::Ptr<AtomVector<Real> > getAtomVector(void) {
    return avec_;
  }
  
  ROL::Ptr<ProbabilityVector<Real> > getProbabilityVector(void) {
    return pvec_;
  }

  void setAtomVector(const AtomVector<Real> &vec) {
    avec_->set(vec);
  }

  void setProbabilityVector(const ProbabilityVector<Real> &vec) {
    pvec_->set(vec);
  }
  
}; // class SROMVector

} // namespace ROL

#endif

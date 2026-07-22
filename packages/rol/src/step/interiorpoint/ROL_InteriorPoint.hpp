// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_INTERIORPOINT_H
#define ROL_INTERIORPOINT_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_ObjectiveFromBoundConstraint.hpp"

namespace ROL {

namespace InteriorPoint {

/** @ingroup func_group
 *  \class ROL::PenalizedObjective
 *  \brief Adds barrier term to generic objective
 */
template <class Real>
class PenalizedObjective : public ROL::Objective<Real> {
private:
  typedef typename PartitionedVector<Real>::size_type  size_type;

  ROL::Ptr<Objective<Real> >       obj_;
  ROL::Ptr<Objective<Real> >       barrier_;
  ROL::Ptr<Vector<Real> >          x_;
  ROL::Ptr<Vector<Real> >          g_;
  ROL::Ptr<Vector<Real> >          scratch_;

  Real mu_;
  Real fval_;
  Real gnorm_;
  int nfval_;
  int ngval_;

public:

  PenalizedObjective( const ROL::Ptr<Objective<Real> >       &obj,
                      const ROL::Ptr<BoundConstraint<Real> > &bnd,
                      const Vector<Real>                         &x,
                      ROL::ParameterList                     &parlist)
    : obj_(obj), x_(ROL::nullPtr), g_(ROL::nullPtr), scratch_(ROL::nullPtr),
      fval_(0), gnorm_(0), nfval_(0), ngval_(0) {
    ROL::ParameterList& IPlist = parlist.sublist("Step").sublist("Interior Point");
    barrier_ = ROL::makePtr<ObjectiveFromBoundConstraint<Real>>(*bnd,IPlist);
    x_       = x.clone();
    g_       = x.dual().clone();
    scratch_ = x.dual().clone();
    mu_      = IPlist.get("Initial Barrier Parameter",1.0);
  }

  void updatePenalty( Real mu ) {
    mu_ = mu;
  }

  int getNumberFunctionEvaluations(void) {
    return nfval_;
  }

  int getNumberGradientEvaluations(void) {
    return ngval_;
  }

  void reset(void) {
    nfval_ = 0; nfval_ = 0;
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    // Update original objective and bound penalty
    obj_->update(x,flag,iter);
    barrier_->update(x,flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    // Compute original objective value and bound penalty value
    fval_ = obj_->value(x,tol);
    Real val  = fval_;
    Real bval = barrier_->value(x,tol);
    val += mu_*bval;

    ++nfval_;
    return val;
  }

  Real getObjectiveValue(void) {
    return fval_;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Compute gradient of objective and bound penalty
    obj_->gradient(g,x,tol);
    barrier_->gradient(*scratch_,x,tol);
    scratch_->scale(mu_);
    g.plus(*scratch_);

    g_->set(g);
    gnorm_ = g.norm();
    ++ngval_;
  }

  void getObjectiveGradient( Vector<Real> &g ) {
    g.set(*g_);
  }

  Real getGradientNorm() {
    return gnorm_;
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                 const Vector<Real> &x, Real &tol ) {
    // Compute hessvec of objective and bound penalty
    obj_->hessVec(hv, v, x, tol);
    barrier_->hessVec(*scratch_,v,x,tol);
    scratch_->scale(mu_);
    hv.plus(*scratch_);
  }

}; // class InteriorPointObjective

} // namespace InteriorPoint
} // namespace ROL

#endif

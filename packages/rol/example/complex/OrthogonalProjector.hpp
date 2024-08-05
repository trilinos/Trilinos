// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ORTHOGONALPROJECTOR_HPP
#define ORTHOGONALPROJECTOR_HPP

#include <initializer_list>

#include "ROL_LinearOperator.hpp"
#include "ROL_ComplexStdVector.hpp"

/** \class HermitianMatrix
    \brief Implementation of a Hermitian Matrix
*/

template<typename Real>
class OrthogonalProjector: public ROL::LinearOperator<Real> {
public:

  OrthogonalProjector( const ROL::Vector<Real>& u ) : 
    u_(u.clone()), alpha_(1/u_->dot(*u_)) {
    u_->set(u);
  }

  void apply( ROL::Vector<Real>& Hv, 
              const ROL::Vector<Real>& v, 
              Real& tol ) const override {
    Hv.set(*u_);
    Hv.scale(alpha_*u_->dot(v));
  }

  void applyAdjoint( ROL::Vector<Real>& Hv, 
                     const ROL::Vector<Real>& v,
                     Real& tol ) const override {
    return apply(Hv,v,tol);
  }

private:

  ROL::Ptr<ROL::Vector<Real>> u_;
  Real alpha_;

};

#endif // ORTHOGONALPROJECTOR_HPP


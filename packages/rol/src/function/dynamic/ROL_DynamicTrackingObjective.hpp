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

#pragma once
#ifndef ROL_DYNAMICTRACKINGOBJECTIVE_HPP
#define ROL_DYNAMICTRACKINGOBJECTIVE_HPP

#include "ROL_DynamicObjective.hpp"
#include "ROL_VectorWorkspace.hpp"


/** @ingroup func_group
    \class ROL::DynamicTrackingObjective
    \brief Defines the time-local contribution to a quadratic tracking
           objective

           \f[ f_k(u,z) = \frac{1}{2} \int\limits_{t_{k-1}}^{t_k} \| u(t)-\tilde u(t)\|^2 
                                              + \alpha \|z(t)\^2\,\mathm{d}t \f]

           Currently approximates the integral with the trapezoidal rule.

*/


namespace ROL {

template<typename Real>
class DynamicTrackingObjective : public DynamicObjective<Real> {
public:

  using V         = Vector<Real>;
  using TS        = TimeStamp<Real>;
  using size_type = typename PartitionedVector<Real>::size_type;  

private:

  Ptr<PartitionedVector<Real>>  target_;

  size_type             Nt_;     // Number of time steps
  Real                  alpha_;  // Regularization parameter

  mutable VectorWorkspace<Real> workspace_;

public:

  DynamicTrackingObjective( const Ptr<PartitionedVector<Real>>& target, Real alpha=0.0 ) : 
    target_(target), Nt_(target_->numVectors()), alpha_(alpha) {}
  
  virtual ~DynamicTrackingObjective() {}  

  virtual Real value( const V& uo, const V& un, 
                      const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);

    size_type k = timeStamp.k;     

    auto udiff = workspace_.copy(un);
    auto utarg = target_->get(k);

    udiff->axpy(-1.0, *utarg);

    Real result = 0.5*dt*( 0.5*udiff->dot(*udiff) + alpha_*z.dot(z) );
    
    if( k>0 ) {
      utarg = target_->get(k-1);
      udiff->set(uo);
      udiff->axpy(-1.0,*utarg);
      result += 0.25*dt*(udiff->dot(*udiff));
    }
    return result;
  }

  //----------------------------------------------------------------------------
  // Gradient Terms
  virtual void gradient_uo( V& g, const V& uo, const V& un, 
                            const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    if( timeStamp.k>0 ) {
      g.set(uo);
      g.axpy(-1.0, *(target_->get(timeStamp.k-1)) );
      g.scale(0.5*dt);
    }
    else g.zero();
  }

  virtual void gradient_un( V& g, const V& uo, const V& un, 
                            const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    g.set(un);
    g.axpy(-1.0, *(target_->get(timeStamp.k)) );
    g.scale(0.5*dt);
  }

  virtual void gradient_z( V& g, const V& uo, const V& un, 
                           const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    g.set(z);
    g.scale(dt*alpha_);
  }

  //----------------------------------------------------------------------------
  // Hessian-Vector product terms
  virtual void hessVec_uo_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override {
    Real dt = timeStamp.t.at(1)-timeStamp.t.at(0);
    if( timeStamp.k>0 ) {
      hv.set(v);
      hv.axpy(-1.0, *(target_->get(timeStamp.k-1)) );
      hv.scale(0.5*dt);
    }
    else hv.zero();
  }

  virtual void hessVec_uo_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override { 
    hv.zero(); 
  }

  virtual void hessVec_uo_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const override { 
    hv.zero(); 
  }

  virtual void hessVec_un_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const override {
    hv.zero();
  }

  virtual void hessVec_un_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {
    hv.set(v);
    hv.scale(0.5*(timeStamp.t.at(1)-timeStamp.t.at(0)));
  }

  virtual void hessVec_un_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  virtual void hessVec_z_uo( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  virtual void hessVec_z_un( V& hv, const V& v, const V& uo, const V& un, 
                              const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  virtual void hessVec_z_z( V& hv, const V& v, const V& uo, const V& un, 
                             const V& z, const TS& timeStamp ) const {
    hv.set(v);
    hv.scale(alpha_*(timeStamp.t.at(1)-timeStamp.t.at(0)));
  }

}; // DynamicTrackingObjective


} // namespace ROL

#endif  // ROL_DYNAMICTRACKINGOBJECTIVE_HPP



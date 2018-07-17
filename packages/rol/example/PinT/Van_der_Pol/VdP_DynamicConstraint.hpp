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
#ifndef VDP_DYNAMICCONSTRAINT_HPP
#define VDP_DYNAMICCONSTRAINT_HPP

#include "ROL_DynamicConstraint.hpp"
#include "ROL_StdVector.hpp"

/** \class VdP::DynamicConstraint
    \brief Compute time-step for the Van der Pol optimal control problem

    The Van der Pol equation is

    \f[ \dot u_1 = u_2,\quad \dot u_2 = z(1-u_1^2)u_2-u_1 \f]
   
    The initial conditions are

    \f[ u_1(0) = 1,\quad u_2(0) = 0 \f]

*/

namespace VdP {

template<typename Real> 
class DynamicConstraint : public ROL::DynamicConstraint<Real> {

  using V  = ROL::Vector<Real>;
  using TS = ROL::TimeStamp<Real>;

private:

  // Jacobians
  Real Jo_[2][2];
  Real Jn_[2][2];
  Real Jz_[2];
  Real iJn_[2][2];
   
  // Hessians
  Real Hoo_[2][2][2];
  Real Hnn_[2][2][2];
  Real Hzo_[2][2];
  Real Hzn_[2][2];

  ROL::Ptr<std::vector<Real>> getVector( V& x ) const {
    return (static_cast<ROL::StdVector<Real>&>(x)).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getVector( const V& x ) const {
    return (static_cast<const ROL::StdVector<Real>&>(x)).getVector();
  }

public: 

  void update( const V& uo, const V& un, const V& z, const TS& ts ) override {

    auto uop  = getVector(uo);   
    auto unp  = getVector(un);   
    auto zp   = getVector(z);   

    Real dt = ts.t[1]-ts.t[0];

    Real uo0 = uop->at(0); 
    Real uo1 = uop->at(1); 
    Real un0 = unp->at(0); 
    Real un1 = unp->at(1); 

    Jo_[0][0] = -1.0;
    Jo_[0][1] = -0.5*dt;
    Jo_[1][0] = -0.5*dt + dt*zp->at(0)*uo0*uo1;
    Jo_[1][1] = -1.0 - 0.5*dt*zp->at(0)*(1-uo0*uo0);

    Jn_[0][0] =  1.0;
    Jn_[0][1] = -0.5*dt;
    Jn_[1][0] = -0.5*dt + dt*zp->at(0)*un0*un1;
    Jn_[1][1] =  1.0 - 0.5*dt*zp->at(0)*(1-un0*un0);
 
    Jz_[0] =  0.0;
    Jz_[1] = -0.5*dt*( (1.0-uo0*uo0)*uo1 + (1.0-un0*un0)*un1 );

    Real Jdet = Jn_[0][0]*Jn_[1][1]-Jn_[0][1]*Jn_[1][0];
    
    iJn_[0][0] =  Jn_[1][1]/Jdet;
    iJn_[0][1] = -Jn_[0][1]/Jdet;
    iJn_[1][0] = -Jn_[1][0]/Jdet;
    iJn_[1][1] =  Jn_[0][0]/Jdet;

    Hoo_[0][0][0] = 0.0;              Hoo_[0][0][1] = 0.0;
    Hoo_[0][1][0] = 0.0;              Hoo_[0][1][1] = 0.0;
    Hoo_[1][0][0] = dt*zp->at(0)*uo1; Hoo_[1][0][1] = dt*zp->at(0)*uo0; 
    Hoo_[1][1][0] = dt*zp->at(0)*uo0; Hoo_[1][1][1] = 0.0;

    Hnn_[0][0][0] = 0.0;              Hnn_[0][0][1] = 0.0;
    Hnn_[0][1][0] = 0.0;              Hnn_[0][1][1] = 0.0;
    Hnn_[1][0][0] = dt*zp->at(0)*un1; Hnn_[1][0][1] = dt*zp->at(0)*un0; 
    Hnn_[1][1][0] = dt*zp->at(0)*un0; Hnn_[1][1][1] = 0.0;

    Hzo_[0][0] = 0.0;
    Hzo_[0][1] = 0.0;
    Hzo_[1][0] = dt*uo0*uo1;
    Hzo_[1][1] = -0.5*dt*(1.0-uo0*uo0);

    Hzn_[0][0] = 0.0;
    Hzn_[0][1] = 0.0;
    Hzn_[1][0] = dt*un0*un1;
    Hzn_[1][1] = -0.5*dt*(1.0-un0*un0);

  } // update

  void value( V& c, const V& uo, const V& un, 
                    const V& z, const TS& ts ) const override {

    auto cp  = getVector(c);  
    auto uop = getVector(uo);   
    auto unp = getVector(un);
    auto zp  = getVector(z);

    Real dt = ts.t[1]-ts.t[0];

    Real uo0 = uop->at(0); 
    Real uo1 = uop->at(1); 
    Real un0 = unp->at(0); 
    Real un1 = unp->at(1); 

    cp->at(0) = un0 - uo0 - 0.5*dt*( un1 + uo1 );
    cp->at(1) = un1 - uo1 - 0.5*dt*( zp->at(0)*( (1-uo0*uo0)*uo1 + 
                                                 (1-un0*un0)*un1 ) + uo0 + un0 );
  }

  //----------------------------------------------------------------------------
  // Partial Jacobians
  void applyJacobian_uo( V& jv, const V& v,  const V& uo, 
                                const V& un, const V& z, 
                                const TS& ts ) const override {

    auto jvp = getVector(jv);   
    auto vp  = getVector(v);

    jvp->at(0) = Jo_[0][0]*vp->at(0) + Jo_[0][1]*vp->at(1);
    jvp->at(1) = Jo_[1][0]*vp->at(0) + Jo_[1][1]*vp->at(1);  
  }

  void applyJacobian_un( V& jv, const V& v,  const V& uo, 
                                const V& un, const V& z, 
                                const TS& ts ) const override {

    auto jvp = getVector(jv);   
    auto vp  = getVector(v);

    jvp->at(0) = Jn_[0][0]*vp->at(0) + Jn_[0][1]*vp->at(1);
    jvp->at(1) = Jn_[1][0]*vp->at(0) + Jn_[1][1]*vp->at(1);  
  }

  void applyJacobian_z( V& jv, const V& v,  const V& uo, 
                               const V& un, const V& z, 
                               const TS& ts ) const override {

    auto jvp = getVector(jv);   
    auto vp  = getVector(v);

    jvp->at(0) = 0.0;
    jvp->at(1) = Jz_[1]*vp->at(0);
  }

  //----------------------------------------------------------------------------
  // Adjoint partial Jacobians
  void applyAdjointJacobian_uo( V& ajv, const V& v, const V& uo, 
                                        const V& un, const V& z, 
                                        const TS& ts ) const override {

    auto ajvp = getVector(ajv);  
    auto vp   = getVector(v);

    ajvp->at(0) = Jo_[0][0]*vp->at(0) + Jo_[1][0]*vp->at(1);
    ajvp->at(1) = Jo_[0][1]*vp->at(0) + Jo_[1][1]*vp->at(1);  
  }

  void applyAdjointJacobian_un( V& ajv, const V& v,  const V& uo, 
                                        const V& un, const V& z, 
                                        const TS& ts ) const override {
    auto ajvp = getVector(ajv);   
    auto vp   = getVector(v);

    ajvp->at(0) = Jn_[0][0]*vp->at(0) + Jn_[1][0]*vp->at(1);
    ajvp->at(1) = Jn_[0][1]*vp->at(0) + Jn_[1][1]*vp->at(1);  
   }


  void applyAdjointJacobian_z( V& ajv, const V& v,  const V& uo, 
                                       const V& un, const V& z, 
                                       const TS& ts ) const override {

    auto ajvp = getVector(ajv);  auto vp  = getVector(v);
    
    ajvp->at(0) = Jz_[1]*vp->at(1);
  }

  //----------------------------------------------------------------------------
  // Inverses
  void applyInverseJacobian_un( V& ijv, const V& v,  const V& uo, 
                                        const V& un, const V& z, 
                                        const TS& ts ) const override {
    auto ijvp = getVector(ijv);   
    auto vp   = getVector(v);
    ijvp->at(0) = iJn_[0][0]*vp->at(0) + iJn_[0][1]*vp->at(1);      
    ijvp->at(1) = iJn_[1][0]*vp->at(0) + iJn_[1][1]*vp->at(1);    
  }
    
  void applyInverseAdjointJacobian_un( V& iajv, const V& v,  const V& uo, 
                                                const V& un, const V& z, 
                                                const TS& ts ) const override {
    auto iajvp = getVector(iajv);  
    auto vp    = getVector(v);
    iajvp->at(0) = iJn_[0][0]*vp->at(0) + iJn_[1][0]*vp->at(1);      
    iajvp->at(1) = iJn_[0][1]*vp->at(0) + iJn_[1][1]*vp->at(1);    

   }

  //----------------------------------------------------------------------------
  // Adjoint Hessian components

  void applyAdjointHessian_uo_uo( V& ahwv, const V& w,  const V& v,
                                           const V& uo, const V& un, 
                                           const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    Real v0 = vp->at(0);    Real v1 = vp->at(1);
    Real w0 = wp->at(0);    Real w1 = wp->at(1);

    ahwvp->at(0) = w0*Hoo_[0][0][0]*v0 + w1*Hoo_[1][0][0]*v0 +
                   w0*Hoo_[0][0][1]*v1 + w1*Hoo_[1][0][1]*v1;
    ahwvp->at(1) = w0*Hoo_[0][1][0]*v0 + w1*Hoo_[1][1][0]*v0 +
                   w0*Hoo_[0][1][1]*v1 + w1*Hoo_[1][1][1]*v1;
  }

  void applyAdjointHessian_uo_z( V& ahwv, const V& w,  const V& v,
                                          const V& uo, const V& un, 
                                          const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    Real v0 = vp->at(0);   
    Real v1 = vp->at(1);   
    Real w0 = wp->at(0); 
    Real w1 = wp->at(1); 

    ahwvp->at(0) = ( w0*Hzo_[0][0] + w1*Hzo_[1][0] )*v0 + ( w0*Hzo_[0][1] + w1*Hzo_[1][1] )*v1 ;
  }

  void applyAdjointHessian_un_un( V& ahwv, const V& w,  const V& v,
                                           const V& uo, const V& un, 
                                           const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    Real v0 = vp->at(0);    Real v1 = vp->at(1);
    Real w0 = wp->at(0);    Real w1 = wp->at(1);

    ahwvp->at(0) = w0*Hnn_[0][0][0]*v0 + w1*Hnn_[1][0][0]*v0 +
                   w0*Hnn_[0][0][1]*v1 + w1*Hnn_[1][0][1]*v1;
    ahwvp->at(1) = w0*Hnn_[0][1][0]*v0 + w1*Hnn_[1][1][0]*v0 +
                   w0*Hnn_[0][1][1]*v1 + w1*Hnn_[1][1][1]*v1;
  }

  void applyAdjointHessian_un_z( V& ahwv, const V& w,  const V& v,
                                          const V& uo, const V& un, 
                                          const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    Real v0 = vp->at(0);   
    Real v1 = vp->at(1);   
    Real w0 = wp->at(0); 
    Real w1 = wp->at(1); 

    ahwvp->at(0) = ( w0*Hzn_[0][0] + w1*Hzn_[1][0] )*v0 + ( w0*Hzn_[0][1] + w1*Hzn_[1][1] )*v1 ;
  }


  void applyAdjointHessian_z_uo( V& ahwv, const V& w,  const V& v,
                                          const V& uo, const V& un, 
                                          const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    Real v0 = vp->at(0);   
    Real w0 = wp->at(0); 
    Real w1 = wp->at(1); 

    ahwvp->at(0) = ( w0*Hzo_[0][0] + w1*Hzo_[1][0] )*v0 ;
    ahwvp->at(1) = ( w0*Hzo_[0][1] + w1*Hzo_[1][1] )*v0 ;
  }

  void applyAdjointHessian_z_un( V& ahwv, const V& w,  const V& v,
                                          const V& uo, const V& un, 
                                          const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    Real v0 = vp->at(0);   
    Real w0 = wp->at(0); 
    Real w1 = wp->at(1); 

    ahwvp->at(0) = ( w0*Hzn_[0][0] + w1*Hzn_[1][0] )*v0 ;
    ahwvp->at(1) = ( w0*Hzn_[0][1] + w1*Hzn_[1][1] )*v0 ;
  }
 
}; // VdP::DynamicConstraint

} // namespace VdP


#endif // VDP_DYNAMICCONSTRAINT_HPP

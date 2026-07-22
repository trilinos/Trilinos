// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  mutable Real uo0, uo1, un0, un1, z0, dt;

  // Jacobians
  mutable Real Jo[2][2];
  mutable Real Jn[2][2];
  mutable Real Jz_[2];
  mutable Real iJn[2][2];
   
  // Hessians
  mutable Real Hoo[2][2][2];
  mutable Real Hnn[2][2][2];
  mutable Real Hzo[2][2];
  mutable Real Hzn[2][2];

  ROL::Ptr<std::vector<Real>> getVector( V& x ) const {
    return (static_cast<ROL::StdVector<Real>&>(x)).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getVector( const V& x ) const {
    return (static_cast<const ROL::StdVector<Real>&>(x)).getVector();
  }

  void update_uz( const V& uo, const V& un, const V& z, const TS& ts ) const {
    auto uop  = getVector(uo);    
    auto unp  = getVector(un);   
    auto zp   = getVector(z);      
    dt  = ts.t[1]-ts.t[0];
    uo0 = uop->at(0); 
    un0 = unp->at(0); 
    uo1 = uop->at(1); 
    un1 = unp->at(1); 
    z0  = zp->at(0);
  }

  void update_Jo( const V& uo, const V& un, const V& z, const TS& ts ) const {
    update_uz(uo,un,z,ts);
    Jo[0][0] = -1.0;
    Jo[0][1] = -0.5*dt;
    Jo[1][0] = -0.5*dt + dt*z0*uo0*uo1;
    Jo[1][1] = -1.0 - 0.5*dt*z0*(1-uo0*uo0);
  }
  void update_Jn( const V& uo, const V& un, const V& z, const TS& ts ) const {
    update_uz(uo,un,z,ts);
    Jn[0][0] =  1.0;
    Jn[0][1] = -0.5*dt;
    Jn[1][0] = -0.5*dt + dt*z0*un0*un1;
    Jn[1][1] =  1.0 - 0.5*dt*z0*(1-un0*un0);

    Real Jdet = Jn[0][0]*Jn[1][1]-Jn[0][1]*Jn[1][0];
    
    iJn[0][0] =  Jn[1][1]/Jdet;
    iJn[0][1] = -Jn[0][1]/Jdet;
    iJn[1][0] = -Jn[1][0]/Jdet;
    iJn[1][1] =  Jn[0][0]/Jdet;
  }

  void update_Jz( const V& uo, const V& un, const V& z, const TS& ts ) const {
    update_uz(uo,un,z,ts);
    Jz_[0] =  0.0;
    Jz_[1] = -0.5*dt*( (1.0-uo0*uo0)*uo1 + (1.0-un0*un0)*un1 );
  }

  void update_Hoo( const V& uo, const V& un, const V& z, const TS& ts ) const {
    update_uz(uo,un,z,ts);
    Hoo[0][0][0] = 0.0;       Hoo[0][0][1] = 0.0;
    Hoo[0][1][0] = 0.0;       Hoo[0][1][1] = 0.0;
    Hoo[1][0][0] = dt*z0*uo1; Hoo[1][0][1] = dt*z0*uo0; 
    Hoo[1][1][0] = dt*z0*uo0; Hoo[1][1][1] = 0.0;
  }

  void update_Hnn( const V& uo, const V& un, const V& z, const TS& ts ) const {
    update_uz(uo,un,z,ts);
    Hnn[0][0][0] = 0.0;       Hnn[0][0][1] = 0.0;
    Hnn[0][1][0] = 0.0;       Hnn[0][1][1] = 0.0;
    Hnn[1][0][0] = dt*z0*un1; Hnn[1][0][1] = dt*z0*un0; 
    Hnn[1][1][0] = dt*z0*un0; Hnn[1][1][1] = 0.0;
  }

  void update_Hzo( const V& uo, const V& un, const V& z, const TS& ts ) const {
    update_uz(uo,un,z,ts);
    Hzo[0][0] = 0.0;
    Hzo[0][1] = 0.0;
    Hzo[1][0] = dt*uo0*uo1;
    Hzo[1][1] = -0.5*dt*(1.0-uo0*uo0);
  }

  void update_Hzn( const V& uo, const V& un, const V& z, const TS& ts ) const {
    update_uz(uo,un,z,ts);
    Hzn[0][0] = 0.0;
    Hzn[0][1] = 0.0;
    Hzn[1][0] = dt*un0*un1;
    Hzn[1][1] = -0.5*dt*(1.0-un0*un0);
  }

public: 

  void value( V& c, const V& uo, const V& un, 
                    const V& z, const TS& ts ) const override {

    auto cp  = getVector(c);  
    auto uop = getVector(uo);   
    auto unp = getVector(un);
    auto zp  = getVector(z);

    update_uz(uo,un,z,ts);

    cp->at(0) = un0 - uo0 - 0.5*dt*( un1 + uo1 );
    cp->at(1) = un1 - uo1 - 0.5*dt*( z0*( (1-uo0*uo0)*uo1 + 
                                          (1-un0*un0)*un1 ) + uo0 + un0 );
  }

  //----------------------------------------------------------------------------
  // Partial Jacobians
  void applyJacobian_uo( V& jv, const V& v,  const V& uo, 
                                const V& un, const V& z, 
                                const TS& ts ) const override {

    auto jvp = getVector(jv);   auto vp  = getVector(v);
    update_Jo(uo,un,z,ts);
    
    jvp->at(0) = Jo[0][0]*vp->at(0) + Jo[0][1]*vp->at(1);
    jvp->at(1) = Jo[1][0]*vp->at(0) + Jo[1][1]*vp->at(1);  
  }

  void applyJacobian_un( V& jv, const V& v,  const V& uo, 
                                const V& un, const V& z, 
                                const TS& ts ) const override {

    auto jvp = getVector(jv);   auto vp  = getVector(v);
    update_Jn(uo,un,z,ts);

    jvp->at(0) = Jn[0][0]*vp->at(0) + Jn[0][1]*vp->at(1);
    jvp->at(1) = Jn[1][0]*vp->at(0) + Jn[1][1]*vp->at(1);  
  }

  void applyJacobian_z( V& jv, const V& v,  const V& uo, 
                               const V& un, const V& z, 
                               const TS& ts ) const override {

    auto jvp = getVector(jv);  auto vp  = getVector(v);
    update_Jz(uo,un,z,ts);

    jvp->at(0) = 0.0;  jvp->at(1) = Jz_[1]*vp->at(0);
  }

  //----------------------------------------------------------------------------
  // Adjoint partial Jacobians
  void applyAdjointJacobian_uo( V& ajv, const V& v, const V& uo, 
                                        const V& un, const V& z, 
                                        const TS& ts ) const override {

    auto ajvp = getVector(ajv);  auto vp   = getVector(v);
    update_Jo(uo,un,z,ts);

    ajvp->at(0) = Jo[0][0]*vp->at(0) + Jo[1][0]*vp->at(1);
    ajvp->at(1) = Jo[0][1]*vp->at(0) + Jo[1][1]*vp->at(1);  
  }

  void applyAdjointJacobian_un( V& ajv, const V& v,  const V& uo, 
                                        const V& un, const V& z, 
                                        const TS& ts ) const override {
    auto ajvp = getVector(ajv);   auto vp   = getVector(v);
    update_Jn(uo,un,z,ts);

    ajvp->at(0) = Jn[0][0]*vp->at(0) + Jn[1][0]*vp->at(1);
    ajvp->at(1) = Jn[0][1]*vp->at(0) + Jn[1][1]*vp->at(1);  
   }


  void applyAdjointJacobian_z( V& ajv, const V& v,  const V& uo, 
                                       const V& un, const V& z, 
                                       const TS& ts ) const override {

    auto ajvp = getVector(ajv);  auto vp  = getVector(v);
    update_Jz(uo,un,z,ts);
    
    ajvp->at(0) = Jz_[1]*vp->at(1);
  }

  //----------------------------------------------------------------------------
  // Inverses
  void applyInverseJacobian_un( V& ijv, const V& v,  const V& uo, 
                                        const V& un, const V& z, 
                                        const TS& ts ) const override {
    auto ijvp = getVector(ijv);   auto vp   = getVector(v);
    update_Jn(uo,un,z,ts);

    ijvp->at(0) = iJn[0][0]*vp->at(0) + iJn[0][1]*vp->at(1);      
    ijvp->at(1) = iJn[1][0]*vp->at(0) + iJn[1][1]*vp->at(1);    
  }
    
  void applyInverseAdjointJacobian_un( V& iajv, const V& v,  const V& uo, 
                                                const V& un, const V& z, 
                                                const TS& ts ) const override {
    auto iajvp = getVector(iajv); auto vp = getVector(v);
    update_Jn(uo,un,z,ts);

    iajvp->at(0) = iJn[0][0]*vp->at(0) + iJn[1][0]*vp->at(1);      
    iajvp->at(1) = iJn[0][1]*vp->at(0) + iJn[1][1]*vp->at(1);    

   }

  //----------------------------------------------------------------------------
  // Adjoint Hessian components

  void applyAdjointHessian_uo_uo( V& ahwv, const V& w,  const V& v,
                                           const V& uo, const V& un, 
                                           const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    update_Hoo(uo,un,z,ts);

    Real v0 = vp->at(0);    Real v1 = vp->at(1);
    Real w0 = wp->at(0);    Real w1 = wp->at(1);

    ahwvp->at(0) = w0*Hoo[0][0][0]*v0 + w1*Hoo[1][0][0]*v0 +
                   w0*Hoo[0][0][1]*v1 + w1*Hoo[1][0][1]*v1;
    ahwvp->at(1) = w0*Hoo[0][1][0]*v0 + w1*Hoo[1][1][0]*v0 +
                   w0*Hoo[0][1][1]*v1 + w1*Hoo[1][1][1]*v1;
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

    update_Hzo(uo,un,z,ts);

    ahwvp->at(0) = ( w0*Hzo[0][0] + w1*Hzo[1][0] )*v0 + ( w0*Hzo[0][1] + w1*Hzo[1][1] )*v1 ;
  }

  void applyAdjointHessian_un_un( V& ahwv, const V& w,  const V& v,
                                           const V& uo, const V& un, 
                                           const V& z,  const TS& ts ) const override {
    auto ahwvp = getVector(ahwv);    
    auto wp    = getVector(w);
    auto vp    = getVector(v);     

    Real v0 = vp->at(0);    Real v1 = vp->at(1);
    Real w0 = wp->at(0);    Real w1 = wp->at(1);

    update_Hnn(uo,un,z,ts);

    ahwvp->at(0) = w0*Hnn[0][0][0]*v0 + w1*Hnn[1][0][0]*v0 +
                   w0*Hnn[0][0][1]*v1 + w1*Hnn[1][0][1]*v1;
    ahwvp->at(1) = w0*Hnn[0][1][0]*v0 + w1*Hnn[1][1][0]*v0 +
                   w0*Hnn[0][1][1]*v1 + w1*Hnn[1][1][1]*v1;

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

    update_Hzn(uo,un,z,ts);

    ahwvp->at(0) = ( w0*Hzn[0][0] + w1*Hzn[1][0] )*v0 + ( w0*Hzn[0][1] + w1*Hzn[1][1] )*v1 ;
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

    update_Hzo(uo,un,z,ts);

    ahwvp->at(0) = ( w0*Hzo[0][0] + w1*Hzo[1][0] )*v0 ;
    ahwvp->at(1) = ( w0*Hzo[0][1] + w1*Hzo[1][1] )*v0 ;
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

    update_Hzn(uo,un,z,ts);

    ahwvp->at(0) = ( w0*Hzn[0][0] + w1*Hzn[1][0] )*v0 ;
    ahwvp->at(1) = ( w0*Hzn[0][1] + w1*Hzn[1][1] )*v0 ;
  }
 
}; // VdP::DynamicConstraint

} // namespace VdP


#endif // VDP_DYNAMICCONSTRAINT_HPP

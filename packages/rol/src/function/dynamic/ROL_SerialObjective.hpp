// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_SERIALOBJECTIVE_HPP
#define ROL_SERIALOBJECTIVE_HPP

#include <type_traits>

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_DynamicObjective.hpp"
#include "ROL_SerialFunction.hpp"

/** @ingroup func_group
    \class ROL::SerialObjective
    \brief Evaluates ROL::DynamicObjective over a sequential set of time intervals


    \f[ f(u,z) = \sum\limits_{k=1}^n f_k(u_{k-1},u_k,z_k) \f]

    \f[ \frac{\partial f}{\partial u_j} = \frac{\partial f_j(u_{j-1},u_j,z_j}{\partial u_j} + 
                                          \frac{\partial f_{j+1}(u_j,u_{j+1},z_{j+1}}{\partial u_j} \f]

   
    ---
*/

namespace ROL {

template<typename Real>
class SerialObjective : public Objective_SimOpt<Real>,
                        public SerialFunction<Real> {
private:
  using PV = PartitionedVector<Real>;  
  using SerialFunction<Real>::ts;
  using SerialFunction<Real>::clone;

  Ptr<DynamicObjective<Real>>  obj_;        // Objective over a single time step

public:

  using size_type = typename std::vector<Real>::size_type;
  using SerialFunction<Real>::numTimeSteps;
  using SerialFunction<Real>::getZeroState;
  using SerialFunction<Real>::getInitialCondition;
  using SerialFunction<Real>::getSkipInitialCondition;

  SerialObjective( const Ptr<DynamicObjective<Real>>& obj,
                   const Vector<Real>& u_initial,
                   const TimeStampsPtr<Real> timeStampsPtr ) :
    SerialFunction<Real>::SerialFunction( u_initial, timeStampsPtr ),
    obj_(obj) {}

  using Objective_SimOpt<Real>::value;
  virtual Real value( const Vector<Real>& u, 
                      const Vector<Real>& z, 
                      Real& tol ) override {

    auto& up = partition(u);
    auto& zp = partition(z);
    Real result = 0;

    if( !getSkipInitialCondition() ) 
      result += obj_->value( getInitialCondition(), up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k ) 
      result += obj_->value( up[k-1], up[k], zp[k], ts(k) );  
   
    return result; 
  } // value

  virtual void gradient_1(       Vector<Real>& g, 
                           const Vector<Real>& u, 
                           const Vector<Real>& z, 
                           Real& tol ) override {

    auto& gp  = partition(g);
    auto& up  = partition(u);
    auto& zp  = partition(z);

    auto  tmp = clone(gp[0]);
    auto& x   = *tmp;

    // TODO: Implement skip initial condition

    obj_->gradient_un( gp[0], getInitialCondition(), up[0], zp[0], ts(0) );
    obj_->gradient_uo( x,     up[0],                 up[1], zp[1], ts(1) );
    gp[0].plus(x);

    for( size_type k=1; k<numTimeSteps()-1; ++k ) {
      obj_->gradient_un( gp[k], up[k-1], up[k],   zp[k],   ts(k)   );
      obj_->gradient_uo( x,     up[k],   up[k+1], zp[k+1], ts(k+1) );
      gp[k].plus(x);
    }
    
    size_t N = numTimeSteps()-1;

    obj_->gradient_un( gp[N], up[N-1], up[N], zp[N], ts(N) );

  } // gradient_1

  virtual void gradient_2(       Vector<Real>& g, 
                           const Vector<Real>& u, 
                           const Vector<Real>& z, 
                           Real& tol ) override {

    auto& gp = partition(g);
    auto& up = partition(u);
    auto& zp = partition(z);

    if( !getSkipInitialCondition() ) 
      obj_->gradient_z( gp[0], getInitialCondition(), up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k ) 
      obj_->gradient_z( gp[k], up[k-1], up[k], zp[k], ts(k) );    // df[k]/dz[k]
     
  } // gradient_2

  virtual void hessVec_11(       Vector<Real>& hv, 
                           const Vector<Real>& v,
                           const Vector<Real>& u, 
                           const Vector<Real>& z, 
                           Real& tol ) override {

    auto& hvp = partition(hv);   auto& vp  = partition(v);
    auto& up  = partition(u);    auto& zp  = partition(z);

    auto tmp = clone(hvp[0]);
    auto& x  = *tmp;

    // TODO: Implement skip initial condition

    obj_->hessVec_un_un( hvp[0], vp[0], getInitialCondition(), up[0], zp[0], ts(0) );
    obj_->hessVec_uo_uo( x,      vp[0], up[0],                 up[1], zp[1], ts(1) );
    hvp[0].plus(x);

    for( size_type k=1; k<numTimeSteps()-1; ++k ) {
      obj_->hessVec_un_un( hvp[k], vp[k], up[k-1], up[k],   zp[k],   ts(k)   );
      obj_->hessVec_uo_uo( x,      vp[k], up[k],   up[k+1], zp[k+1], ts(k+1) );
      hvp[k].plus(x);
   }
   
    size_t N = numTimeSteps()-1;

    obj_->hessVec_un_un( hvp[N], vp[N], up[N-1], up[N], zp[N], ts(N) );

  } // hessVec_11

  virtual void hessVec_12(       Vector<Real>& hv, 
                           const Vector<Real>& v,
                           const Vector<Real>& u,  
                           const Vector<Real>& z, 
                           Real& tol ) override {

    auto& hvp = partition(hv);   auto& vp  = partition(v);
    auto& up  = partition(u);    auto& zp  = partition(z);

    auto tmp = clone(hvp[0]);
    auto& x  = *tmp;

    // TODO: Implement skip initial condition

    obj_->hessVec_un_z( hvp[0], vp[0], getInitialCondition(), up[0], zp[0], ts(0) );
    obj_->hessVec_uo_z( x,      vp[0], up[0],                 up[1], zp[1], ts(1) );
    hvp[0].plus(x);

    for( size_type k=1; k<numTimeSteps()-1; ++k ) {
      obj_->hessVec_un_z( hvp[k], vp[k], up[k-1], up[k],   zp[k],   ts(k)   );
      obj_->hessVec_uo_z( x,      vp[k], up[k],   up[k+1], zp[k+1], ts(k+1) );
      hvp[k].plus(x);
   }
   
    size_t N = numTimeSteps()-1;

    obj_->hessVec_un_z( hvp[N], vp[N], up[N-1], up[N], zp[N], ts(N) );


  } // hessVec_22

  virtual void hessVec_21(       Vector<Real>& hv,
                           const Vector<Real>& v,
                           const Vector<Real>& u,
                           const Vector<Real>& z,
                           Real& tol ) override {

    auto& hvp = partition(hv);   auto& vp  = partition(v);
    auto& up  = partition(u);    auto& zp  = partition(z);

    auto tmp = clone(hvp[0]);
    auto& x  = *tmp;

    // TODO: Implement skip initial condition

    obj_->hessVec_z_un( hvp[0], vp[0], getInitialCondition(), up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k ) {
      obj_->hessVec_z_un( hvp[k], vp[k],   up[k-1], up[k], zp[k], ts(k) );
      obj_->hessVec_z_uo( x,      vp[k-1], up[k-1], up[k], zp[k], ts(k) );
      hvp[k].plus(x);
   }

  } // hessVec_21

  virtual void hessVec_22(       Vector<Real>& hv, 
                           const Vector<Real>& v,
                           const Vector<Real>& u,  
                           const Vector<Real>& z, 
                           Real& tol ) override {

    auto& hvp = partition(hv);   auto& vp  = partition(v);
    auto& up  = partition(u);    auto& zp  = partition(z);

    if( !getSkipInitialCondition() ) 
      obj_->hessVec_z_z( hvp[0], vp[0], getInitialCondition(), up[0], zp[0], ts(0) );

    for( size_type k=1; k<numTimeSteps(); ++k ) 
      obj_->hessVec_z_z( hvp[k], vp[k], up[k-1], up[k], zp[k], ts(k) );
     

  } // hessVec_22

}; // SerialObjective


// Helper function to create a new SerialObjective

template<typename DynObj, typename Real, typename P = Ptr<SerialObjective<Real>> >
inline typename std::enable_if<std::is_base_of<DynamicObjective<Real>,DynObj>::value,P>::type
make_SerialObjective( const Ptr<DynObj>& obj,
                      const Vector<Real>& u_initial,
                      const TimeStampsPtr<Real> timeStampsPtr ) {
  return makePtr<SerialObjective<Real>>(obj,u_initial,timeStampsPtr);
}

} // namespace ROL


#endif // ROL_SERIALOBJECTIVE_HPP

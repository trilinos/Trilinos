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
#ifndef TANKS_DYNAMICOBJECTIVE_HPP
#define TANKS_DYNAMICOBJECTIVE_HPP

#include "ROL_DynamicObjective.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_VectorWorkspace.hpp"

/** \class Tanks_DynamicObjective based on the new DynamicObjective interface
    \brief Compute contribution to the total objective from a given time step

    The objective is assumed to have the form:

    \f[ f(u,z) = \frac{\alpha}{2} | h(T)-\tilde h(T) |^2 +
                 \frac{\beta}{2 } \int\limits_0^T | h(t)-\tilde h(t) |^2 \,\mathrm{d}t
                 \frac{\gamma}{2} \int\limits_0^T z^2(t) \,\mathrm{d}t

    Integrals are approximated by the trapezoidal rule
*/

namespace Tanks {

template<typename Real>
class DynamicObjective : public ROL::DynamicObjective<Real> {

  using State     = StateVector<Real>;
  using Control   = ControlVector<Real>;
  using V         = ROL::Vector<Real>;
  using TS        = ROL::TimeStamp<Real>;
  using size_type = typename State::size_type; 
 
private:

  size_type rows_;  
  size_type cols_;
  size_type Nt_;     // Number of time steps
  Real      T_;      // Final time
  Real      htarg_;  // Target fluid level
  Real      tol_; 

  Real      alpha_;  // Penalty on final time error
  Real      beta_;   // Penalty on distributed time error
  Real      gamma_;  // Penalty on distributed time error

  ROL::VectorWorkspace<Real> workspace_;

  bool is_first_step( const TS& timeStamp ) const {
    return std::abs(timeStamp.t.at(0)) < tol_;
  }

  bool is_final_step( const TS& timeStamp ) const {
    return std::abs(timeStamp.t.at(1) - T_) < tol_;
  }

  void step_dependent_weights( Real& wo, Real& wn, Real& wz, const TS& timeStamp ) const {
    Real dt = timeStamp.t[1]-timeStamp.t[0];
    wz = dt*gamma_;
    if( std::abs(timeStamp.t.at(0)) < tol_ ) { // first step
      wo = 0.0;
      wn = 0.5*dt*beta_;
    } else if(  std::abs(timeStamp.t.at(1) - T_) < tol_ ) { // last step 
      wo = 0.5*dt*beta_;
      wn = 0.5*dt*beta_ + alpha_;
    } else { // Interior steps 
      wo = 0.5*dt*beta_;
      wn = 0.5*dt*beta_;
    }
  }

public:

  DynamicObjective( ROL::ParameterList& pl ) :
    rows_( static_cast<size_type>( pl.get( "Number of Rows",       3   ) ) ),
    cols_( static_cast<size_type>( pl.get( "Number of Columns",    3   ) ) ),
    Nt_  ( static_cast<size_type>( pl.get( "Number of Time Steps", 100 ) ) ),
    T_   ( pl.get("Total Time",20.0) ),
    tol_ ( T_*std::sqrt( ROL::ROL_EPSILON<Real>() ) ),
    target_height_( rows_*cols_ ),
    residual_n_( rows_*cols_ ) {
    residual_o_( rows_*cols_ ) {

    auto& penalty = pl.sublist( "Penalty Parameters" );
    alpha_ = penalty.get( "Final State",          0.0 );
    beta_  = penalty.get( "Distributed State",    1.0 );
    gamma_ = penalty.get( "Distributed Control" , 0.0 );
  }

  Real value( const V& uo, const V& un, 
              const V& z,  const TS& timeStamp ) const override {

    auto& uo_state = to_state(uo);
    auto& un_state = to_state(un);
    auto& z_ctrl   = to_control(z);

    Real result = 0.0;
    Real duo    = 0.0;    Real wo = 0.0;
    Real dun    = 0.0;    Real wn = 0.0;
    Real zij    = 0.0;    Real wz = 0.0;

    step_dependent_weights( wo, wn, wz, timeStamp );

    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) { 
        dun = un_state.h(i,j) - htarg_;
        duo = uo_state.h(i,j) - htarg_;
        zij = z_ctrl(i,j);
        result += wo*duo*duo + wn*dun*dun + wz*zij*zij;
      }
    }
    return 0.5*result;
  } 

  void gradient_uo( V& g, const V& uo, const V& un, 
                    const V& z, const TS& timeStamp ) const override {
    auto& uo_state = to_state(uo);
    auto& g_state  = to_state(g);
    Real wo = 0.0;
    Real wn = 0.0;
    Real wz = 0.0;
    step_dependent_weights( wo, wn, wz, timeStamp );
 
    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) { 
        g_state.h(i,j) = wo*(uo_state.h(i,j) - htarg_);     
      }
    }
  }

//  void gradient_un( V& g, const V& uo, const V& un, 
//                    const V& z, const TS& timeStamp ) const override {
//    auto& un_state = to_state(un);
//    auto& g_state  = to_state(g);
//    Real wo = 0.0;
//    Real wn = 0.0;
//    Real wz = 0.0;
//    step_dependent_weights( wo, wn, wz, timeStamp );
// 
//    for( size_type i=0; i<rows_; ++i ) {
//      for( size_type j=0; j<cols_; ++j ) { 
//        g_state.h(i,j) = wn*(un_state.h(i,j) - htarg_);     
//      }
//    }
//  }
//
//  void gradient_z( V& g, const V& uo, const V& un, 
//                   const V& z, const TS& timeStamp ) const override {
//    auto& z_ctrl   = to_control(z);
//    auto& g_state  = to_state(g);
//    Real wo = 0.0;
//    Real wn = 0.0;
//    Real wz = 0.0;
//    step_dependent_weights( wo, wn, wz, timeStamp );
// 
//    for( size_type i=0; i<rows_; ++i ) {
//      for( size_type j=0; j<cols_; ++j ) { 
//        g_state.h(i,j) = wn*(un_state.h(i,j) - htarg_);     
//      }
//    }
//
//  }
//  
//  void hessVec_uo_uo( V& hv, const V& vo, const V& uo, const V& un, 
//                      const V& z, const TS& timeStamp ) const {}
//
//  void hessVec_un_un( V& hv, const V& vo, const V& uo, const V& un, 
//                      const V& z, const TS& timeStamp ) const {}
//
//  void hessVec_z_z( V& hv, const V& vo, const V& uo, const V& un, 
//                     const V& z, const TS& timeStamp ) const {}
//  
// 
//


}; // Tanks::DynamicObjective


} // namespace Tanks





#endif  // TANKS_DYNAMICOBJECTIVE_HPP


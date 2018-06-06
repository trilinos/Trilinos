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
#ifndef TANKCONSTRAINTSERIAL_HPP
#define TANKCONSTRAINTSERIAL_HPP

#include "ROL_PartitionedVector.hpp"
#include "ROL_VectorClone.hpp"
#include "TankConstraint.hpp"

namespace details {

using namespace std;
using ROL::makePtr;
using ROL::VectorCloneMap;

template<typename Real> 
class TankConstraintSerial : ROL::Constraint_SimOpt<Real> {

  using V             = ROL::Vector<Real>;
  using PV            = ROL::PartitionedVector<Real>;
  using SV   = TankSV<Real>;
  using CV = TankCV<Real>;

  using size_type = typename vector<Real>::size_type;
  template<typename T> using Ptr = ROL::Ptr<T>;

private:

  Ptr<TankState<Real>>         state_;
  Ptr<TankConstraint<Real>>    con_;    // Constraint for a single time step

//  VectorCloneMap<Real>         clone_;  
  Ptr<SV> ui_;  // Initial condition
  Ptr<SV> u0_;  // Zero State
  Pre<CV> z0_;  // Zero Control

  int Nt_;      // Number of time steps
  
  PV& partition ( V& x )      { return static_cast<PV&>(x); }
  SV&   to_state( Vector& x ) { return static_cast<SV&>(x); }
  CV& to_control( Vector& x ) { return static_cast<CV&>(x); }

  const PV& partition ( const V& )        { return static_cast<const PV&>(x); }
  const SV& to_state  ( const Vector& x ) { return static_cast<const SV&>(x); }
  const CV& to_control( const Vector& x ) { return static_cast<const CV&>(x); }

public:

  TankConstraintSerial( ROL::ParameterList& pl, const Vector<Real>& u ) :
    ROL::Constraint_SimOpt<Real>(),
    state_( makePtr<TankState<Real>>( pl ) ),
    con_( state_, pl ), ui_( ROL::nullPtr ), u0_( ROL::nullPtr ), z0_( ROL::nullPtr ), 
    Nt_( static_cast<size_type>(pl.get( "Number of Time Steps", 100 ) ) ) { 

    Real h0 = pl.get( "Initial Fluid Level", 2.0 );
  
    size_type rows = static_cast<size_type>( pl.get("Number of Rows", 3) );
    size_type cols = static_cast<size_type>( pl.get("Number of Columns", 3) );

    ui_ = makePtr<SV>( rows, cols, "ui" );  ui_->zero();
    u0_ = makePtr<SV>( rows, cols, "u0" );  u0_->zero();
    z0_ = makePtr<CV>( rows, cols, "z0" );  z0_->zero();

    for( size_type i=0; i<rows; ++i ) 
      for( size_type j=0; j<cols; ++j ) 
        ui_->h(i,j) = h0;
  }

  // Update Sim Vector
  void update_1( const V& u, bool flag = true, int iter = -1 ) override {

  }

  // Update Opt Vector 
  void update_2( const V& z, bool flag = true, int iter = -1 ) override {
 
  }

  void value( V& c, const V& u, const V& z, Real& tol ) override {

    auto cp& = partition(c);
    auto up& = partition(u);
    auto zp& = partition(z);
  
    // First interval value uses initial condition
    con_->value( *(cp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)), tol );
    
    for( size_type k=1; k<Nt_; ++k ) 
      con_->value( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );

  } 

  void solve( V& c, V& u, const V& z, Real& tol ) override { 

    auto cp& = partition(c);
    auto up& = partition(u);
    auto zp& = partition(z);
 
    // First interval solve uses initial condition
    con_->solve( *(cp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)), tol );
    
    for( size_type k=1; k<Nt_; ++k ) 
      con_->solve( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), tol );

  }

  void applyJacobian_1( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol ) override {

    auto jvp& = partition(jv);   auto vp& = partition(v);
    auto up&  = partition(u);    auto zp& = partition(z);

    
  }

   void applyJacobian_2( V& jv, const V& v, const V& u,
                         const V &z, Real &tol ) override { 

   }

   void applyInverseJacobian_1( V& ijv, const V& v, const V& u,
                                const V& z, Real& tol) override {

   }


   void applyAdjointJacobian_1( V& ajv, const V& v, const V& u,
                                const V& z, const V& dualv, Real& tol) override {

   }

   void applyAdjointJacobian_2( V& ajv,  const V& v, const V& u,
                                const V& z, Real& tol ) override {

   }

   void applyAdjointJacobian_2( V& ajv, const V& v, const V& u,
                                const V& z, const V& dualv, Real& tol) override {
   }

   void applyInverseAdjointJacobian_1( V& iajv, const V& v, const V& u,
                                       const V& z, Real& tol) override {

   }

   void applyAdjointHessian_11( V& ahwv, const V& w, const V &v,
                                const V& u, const V& z, Real& tol) override {

    }

    void applyAdjointHessian_12(V &ahwv,
                                      const V &w,
                                      const V &v,
                                      const V &u,
                                      const V &z,
                                      Real &tol) override {

    }


    void applyAdjointHessian_21(V &ahwv, const V &w, const V &v,
                                const V &u, const V &z, Real &tol) override {
    }

    void applyAdjointHessian_22(V& ahwv, const V& w, const V& v,
                                const V& u, const V& z, Real& tol) override {
    }


}; // TankConstraintSerial


} // namespace details

using details::TankConstraintSerial;




#endif // TANKCONSTRAINTSERIAL_HPP


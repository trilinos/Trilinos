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
#ifndef TANKS_CONSTRAINTSERIAL_HPP
#define TANKS_CONSTRAINTSERIAL_HPP

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_VectorWorkspace.hpp"

namespace Tanks {

using namespace std;

template<typename Real> 
class ConstraintSerial : ROL::Constraint_SimOpt<Real> {

  using V  = ROL::Vector<Real>;
  using PV = ROL::PartitionedVector<Real>;
  using SV = StateVector<Real>;
  using CV = ControlVector<Real>;

  using size_type = typename vector<Real>::size_type;
  template<typename T> using Ptr = ROL::Ptr<T>;

private:

  Ptr<DynamicConstraint<Real>> con_;    // Constraint for a single time step

  ROL::VectorWorkspace<RealT>  workspace_;

  Ptr<SV> ui_;  // Initial condition
  Ptr<SV> u0_;  // Zero State
  Pre<CV> z0_;  // Zero Control

  int Nt_;      // Number of time steps

  ROL::TimeStep<Real> ts_; // placeholder

  PV& partition ( V& x )      { return static_cast<PV&>(x); }
  SV& to_state  ( Vector& x ) { return static_cast<SV&>(x); }
  CV& to_control( Vector& x ) { return static_cast<CV&>(x); }

  const PV& partition ( const V& )        { return static_cast<const PV&>(x); }
  const SV& to_state  ( const Vector& x ) { return static_cast<const SV&>(x); }
  const CV& to_control( const Vector& x ) { return static_cast<const CV&>(x); }

public:

  SerialConstraint( ROL::ParameterList& pl ) :
    con_( DynamicConstraint<Real>::create(pl) ), 
    ui_( SV::create(pl) ), 
    u0_( SV::create(pl) ), 
    z0_( CV::create(pl) ), 
    Nt_( static_cast<size_type>(pl.get( "Number of Time Steps", 100 ) ) ) { 

    Real h0 = pl.get( "Initial Fluid Level", 2.0 );
  
    size_type rows = static_cast<size_type>( pl.get("Number of Rows", 3) );
    size_type cols = static_cast<size_type>( pl.get("Number of Columns", 3) );

    for( size_type i=0; i<rows; ++i ) 
      for( size_type j=0; j<cols; ++j ) 
        ui_->h(i,j) = h0;
  }

  void value( V& c, const V& u, const V& z, Real& tol ) override {

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
  
    // First interval value uses initial condition
    con_->value( *(cp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)), ts_ );
    
    for( size_type k=1; k<Nt_; ++k ) 
      con_->value( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );

  } 

  void solve( V& c, V& u, const V& z, Real& tol ) override { 

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
 
    // First interval solve uses initial condition
    con_->solve( *(cp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)), ts_ );
    
    for( size_type k=1; k<Nt_; ++k ) 
      con_->solve( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );

  }

  void applyJacobian_1( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol ) override {

    auto& jvp = partition(jv);   auto& vp = partition(v);
    auto& up  = partition(u);    auto& zp = partition(z);

    auto tmp = workspace_.clone(jvp);

    con_->applyJacobian_uo( *(tmp->get(0)), *(vp.get(0)), *u0_, *(up.get(0)), *(zp.get(0)), ts_ );
    con_->applyJacobian_un( *(jvp.get(0)),  *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );
    
    for( size_type k=1; k<Nt_; ++k ) {
      con_->applyJacobian_uo( *(tmp->get(0)), *(vp.get(0)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
      con_->applyJacobian_un( *(jvp.get(0)), *(vp.get(0)),  *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
    }

    jvp.plus(*tmp);
    
  }

   void applyJacobian_2( V& jv, const V& v, const V& u,
                         const V &z, Real &tol ) override { 

    auto& jvp = partition(jv);   auto& vp = partition(v);
    auto& up  = partition(u);    auto& zp = partition(z);

    con_->applyJacobian_z( *(jvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );

    for( size_type k=1; k<Nt_; ++k ) 
      con_->applyJacobian_z( *(jvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
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


}; // Tanks::ConstraintSerial


} // namespace Tanks



#endif // TANKS_CONSTRAINTSERIAL_HPP


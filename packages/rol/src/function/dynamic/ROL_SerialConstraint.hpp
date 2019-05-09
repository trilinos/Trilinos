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
#ifndef ROL_SERIALCONSTRAINT_HPP
#define ROL_SERIALCONSTRAINT_HPP

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_DynamicConstraint.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_VectorWorkspace.hpp"

/** @ingroup func_group
    \class ROL::SerialConstraint
    \brief Evaluates ROL::DynamicConstraint over a sequential set of time intervals

    ---
*/

namespace ROL {

namespace details {

//using namespace std;

template<typename Real>
class SerialConstraint : public ROL::Constraint_SimOpt<Real> {

  using PV = PartitionedVector<Real>;  
  using size_type = typename std::vector<Real>::size_type;

private:
 
  Ptr<DynamicConstraint<Real>>       con_;        // Constraint over a single time step
  Ptr<Vector<Real>>                  u_initial_;  // Initial condition for state
  Ptr<Vector<Real>>                  u_zero_;     // Zero vector the size of u_initial_
  Ptr<std::vector<TimeStamp<Real>>>  timeStamp_; // Set of all time stamps
  VectorWorkspace<Real>              workspace_;  // For efficient cloning
  size_type                          Nt_;         // Number of time steps  
  bool                               skipInitialCond_; // skip the initial condition application

  PV& partition ( Vector<Real>& x ) const { return static_cast<PV&>(x); }
  const PV& partition ( const Vector<Real>& x ) const { return static_cast<const PV&>(x); }

public: 

  SerialConstraint( const Ptr<DynamicConstraint<Real>>& con, 
                    const Vector<Real>& u_initial, 
                    const Ptr<vector<TimeStamp<Real>>>& timeStamp ) : 
    con_(con), u_initial_(u_initial.clone()), 
    u_zero_(u_initial.clone()), 
    timeStamp_(timeStamp), Nt_(timeStamp_->size()),
    skipInitialCond_(false) {
    u_initial_->set(u_initial);
    u_zero_->zero();
  }

  size_type numTimeSteps(void) const { return Nt_; }

  const Vector<Real>& getInitialCondition() const { return *u_initial_; }
  void setInitialCondition( const Vector<Real>& u_initial ) { u_initial_->set(u_initial); }

  bool getSkipInitialCondition() const { return skipInitialCond_; }
  void setSkipInitialCondition(bool skip) { skipInitialCond_ = skip; }

  Ptr<vector<TimeStamp<Real>>> getTimeStamp() const { return timeStamp_; }
  void setTimeStamp(const Ptr<vector<TimeStamp<Real>>>& timeStamp ) 
  { timeStamp_ = timeStamp; Nt_ = timeStamp_->size(); }

  virtual void solve( Vector<Real>& c, Vector<Real>& u,
                      const Vector<Real>& z, Real& tol ) override {
    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);

    if(!skipInitialCond_)
      con_->solve( *(cp.get(0)), *u_initial_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );

    for( size_type k=1; k<Nt_; ++k ) 
       con_->solve( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
  } // solve

   
  virtual void value( Vector<Real>& c, const Vector<Real>& u, 
                      const Vector<Real>& z, Real& tol ) override {

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
  
    if(!skipInitialCond_)
      con_->value( *(cp.get(0)), *u_initial_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );
    
    for( size_type k=1; k<Nt_; ++k ) 
      con_->value( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
  } // value


  virtual void applyJacobian_1( Vector<Real>& jv, const Vector<Real>& v, 
                                const Vector<Real>& u, const Vector<Real>& z, 
                                Real& tol ) override {

    auto& jvp = partition(jv);   auto& vp = partition(v);
    auto& up  = partition(u);    auto& zp = partition(z);

    auto  tmp = workspace_.clone(jvp.get(0));
    auto& x   = *tmp; 

    if(!skipInitialCond_)
      con_->applyJacobian_un( *(jvp.get(0)),  *(vp.get(0)), *u_zero_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );
    
    for( size_type k=1; k<Nt_; ++k ) {
      con_->applyJacobian_uo( x, *(vp.get(k-1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
      con_->applyJacobian_un( *(jvp.get(k)),  *(vp.get(k)),   *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
      jvp.get(k)->plus(x);
    } // end for
  } // applyJacobian_1
  

  virtual void applyInverseJacobian_1( Vector<Real>& ijv, const Vector<Real>& v, const Vector<Real>& u,
                                       const Vector<Real>& z, Real& tol) override {

    auto& ijvp = partition(ijv);  auto& vp = partition(v);
    auto& up   = partition(u);    auto& zp = partition(z);

    auto  tmp  = workspace_.clone(ijvp.get(0));
    auto& x    = *tmp;

    if(!skipInitialCond_)
      con_->applyInverseJacobian_un( *(ijvp.get(0)), *(vp.get(0)), *u_zero_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );

    for( size_type k=1; k<Nt_; ++k ) {
      con_->applyJacobian_uo( x, *(ijvp.get(k-1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
      x.scale(-1.0);
      x.plus( *(vp.get(k)) );
      con_->applyInverseJacobian_un( *(ijvp.get(k)), x, *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
    } // end for
  } // applyInverseJacobian_1
  
 
  virtual void applyAdjointJacobian_1( Vector<Real>& ajv, const Vector<Real>& v, const Vector<Real>& u,
                                       const Vector<Real>& z, const Vector<Real>& dualv, Real& tol) override {

    auto& ajvp  = partition(ajv);   auto& vp = partition(v);
    auto& up    = partition(u);     auto& zp = partition(z);

    auto  tmp  = workspace_.clone(ajvp.get(0)); 
    auto& x    = *tmp; 

    if(!skipInitialCond_)
      con_->applyAdjointJacobian_un( *(ajvp.get(0)),  *(vp.get(0)), *u_zero_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );
   
    for( size_type k=1; k<Nt_; ++k ) {
      con_->applyAdjointJacobian_un( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
      con_->applyAdjointJacobian_uo( x, *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
      ajvp.get(k-1)->plus(x);
    } // end for
  } // applyAdjointJacobian_1

  void applyInverseAdjointJacobian_1( Vector<Real>& iajv, const Vector<Real>& v, const Vector<Real>& u,
                                      const Vector<Real>& z, Real& tol) override {

    auto& iajvp = partition(iajv);  auto& vp = partition(v);
    auto& up    = partition(u);     auto& zp = partition(z);

    auto  tmp  = workspace_.clone(iajvp.get(0));
    auto& x    = *tmp; 

    size_type k = Nt_-1;

    con_->applyInverseAdjointJacobian_un( *(iajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );

    for( size_type k=Nt_-2; k>0; --k ) {
      con_->applyAdjointJacobian_uo( x, *(iajvp.get(k+1)), *(up.get(k)), *(up.get(k+1)), *(zp.get(k+1)), timeStamp_->at(k+1) );
      x.scale(-1.0);
      x.plus( *(vp.get(k) ) );
      con_->applyInverseAdjointJacobian_un( *(iajvp.get(k)), x, *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );             
    } // end for

    con_->applyAdjointJacobian_uo( x, *(iajvp.get(1)), *(up.get(0)), *(up.get(1)), *(zp.get(1)), timeStamp_->at(1) );
    x.scale(-1.0);
    x.plus( *(vp.get(0) ) );
    if(!skipInitialCond_) {
      con_->applyInverseAdjointJacobian_un( *(iajvp.get(0)), x, *u_zero_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );           
    }
    else {
      // this weird condition places iajvp in the final vector slot
      iajvp.get(0)->set(x);
    }
     
  } // applyInverseAdjointJacobian_1


  virtual void applyJacobian_2( Vector<Real>& jv, const Vector<Real>& v, const Vector<Real>& u,
                                const Vector<Real>& z, Real &tol ) override { 
 
    auto& jvp = partition(jv);   auto& vp = partition(v);
    auto& up  = partition(u);    auto& zp = partition(z);
  
    if(!skipInitialCond_)
      con_->applyJacobian_z( *(jvp.get(0)), *(vp.get(0)), *u_zero_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );
  
    for( size_type k=1; k<Nt_; ++k ) 
      con_->applyJacobian_z( *(jvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
  } // applyJacobian_2

 
  virtual void applyAdjointJacobian_2( Vector<Real>& ajv, const Vector<Real>& v, const Vector<Real>& u,
                                       const Vector<Real>& z, Real& tol ) override {

    auto& ajvp  = partition(ajv);   auto& vp = partition(v);
    auto& up    = partition(u);     auto& zp = partition(z);

    if(!skipInitialCond_)
      con_->applyAdjointJacobian_z( *(ajvp.get(0)), *(vp.get(0)), *u_zero_, *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );

    for( size_type k=1; k<Nt_; ++k ) 
      con_->applyAdjointJacobian_z( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
  } // applyAdjointJacobian_2



  virtual void applyAdjointHessian_11( Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                       const Vector<Real> &u, const Vector<Real> &z, Real &tol) override {

    auto& ahwvp = partition(ahwv);    auto& wp = partition(w);
    auto& vp    = partition(v);       auto& up = partition(u);
    auto& zp    = partition(z);

    auto  tmp  = workspace_.clone(ahwvp.get(0)); 
    auto& x    = *tmp; 

   if( !skipInitialCond_ ) {
     con_->applyAdjointHessian_un_un( *(ahwvp.get(0)), *(wp.get(0)), *(vp.get(0)), *u_zero_,     *(up.get(0)), *(zp.get(0)), timeStamp_->at(0) );
     con_->applyAdjointHessian_un_uo(               x, *(wp.get(1)), *(vp.get(1)), *(up.get(0)), *(up.get(1)), *(zp.get(1)), timeStamp_->at(1) );
     ahwvp.get(0)->plus(x);
     con_->applyAdjointHessian_uo_uo( x, *(wp.get(1)), *(vp.get(0)), *(up.get(0)), *(up.get(1)), *(zp.get(1)), timeStamp_->at(1) );
     ahwvp.get(0)->plus(x);

   }

    for( size_type k=1; k<Nt_; ++k ) {
      con_->applyAdjointHessian_un_un( *(ahwvp.get(k)), *(wp.get(k)), *(vp.get(k)),   *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
      con_->applyAdjointHessian_uo_un(               x, *(wp.get(k)), *(vp.get(k-1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
      ahwvp.get(k)->plus(x);

      if( k < Nt_-1 ) {
        con_->applyAdjointHessian_un_uo( x, *(wp.get(k+1)), *(vp.get(k+1)), *(up.get(k)), *(up.get(k+1)), *(zp.get(k+1)), timeStamp_->at(k+1) );
        ahwvp.get(k)->plus(x);
        con_->applyAdjointHessian_uo_uo( x, *(wp.get(k+1)), *(vp.get(k)),   *(up.get(k)), *(up.get(k+1)), *(zp.get(k+1)), timeStamp_->at(k+1) );
        ahwvp.get(k)->plus(x);
      }
    }
  } // applyAdjointHessian_11

  
  virtual void applyAdjointHessian_12(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {}


  virtual void applyAdjointHessian_21(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {}

  virtual void applyAdjointHessian_22(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                      const Vector<Real> &u, const Vector<Real> &z, Real &tol) override {

    auto& ahwvp = partition(ahwv); auto& vp = partition(v);   auto& wp = partition(w);
    auto& up    = partition(u);    auto& zp = partition(z);
  
    con_->applyAdjointHessian_z_z( *(ahwvp.get(0)), *(wp.get(0)), *(vp.get(0)), *u_initial_, 
                                   *(up.get(0)),    *(zp.get(0)), timeStamp_->at(0) );
  
    for( size_type k=1; k<Nt_; ++k ) {
      con_->applyAdjointHessian_z_z( *(ahwvp.get(k)), *(wp.get(k)), *(vp.get(k)), *(up.get(k-1)), 
                                     *(up.get(k)), *(zp.get(k)), timeStamp_->at(k) );
    }
  } // applyAdjointHessian_22


}; // details::SerialConstraint

} // namespace details

template<typename Real>
using SerialConstraint = details::SerialConstraint<Real>;

} // namespace ROL


#endif // ROL_SERIALCONSTRAINT_HPP

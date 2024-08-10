// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKS_CONSTRAINTSERIAL_HPP
#define TANKS_CONSTRAINTSERIAL_HPP

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_TimeStamp.hpp"
#include "ROL_VectorWorkspace.hpp"

namespace Tanks {

using namespace std;

template<typename Real> 
class SerialConstraint : public ROL::Constraint_SimOpt<Real> {

  using V  = ROL::Vector<Real>;
  using PV = ROL::PartitionedVector<Real>;
  using SV = StateVector<Real>;
  using CV = ControlVector<Real>;

  using size_type = typename vector<Real>::size_type;
  template<typename T> using Ptr = ROL::Ptr<T>;

private:

  Ptr<DynamicConstraint<Real>> con_;    // Constraint for a single time step

  ROL::VectorWorkspace<Real>  workspace_;

  Ptr<const V> ui_;  // Initial condition

  size_type Nt_;            // Number of time steps

  ROL::TimeStamp<Real> ts_; // placeholder

  PV& partition ( V& x ) const { return static_cast<PV&>(x); }
  const PV& partition ( const V& x ) const { return static_cast<const PV&>(x); }

  bool skipInitialCond_;

public:

  size_type numTimeSteps() const { return Nt_; }

  SerialConstraint() {}

  SerialConstraint( ROL::ParameterList& pl ) :
    con_( DynamicConstraint<Real>::create(pl) ), 
    Nt_( static_cast<size_type>(pl.get( "Number of Time Steps", 100 ) ) ),
    skipInitialCond_(false) { 

    Ptr<SV> ui = SV::create(pl);  

    Real h0 = pl.get( "Initial Fluid Level", 2.0 );
  
    size_type rows = static_cast<size_type>( pl.get("Number of Rows", 3) );
    size_type cols = static_cast<size_type>( pl.get("Number of Columns", 3) );

    for( size_type i=0; i<rows; ++i ) 
      for( size_type j=0; j<cols; ++j ) 
        ui->h(i,j) = h0;

    ui_ = ui;
  }

  SerialConstraint( const Ptr<DynamicConstraint<Real>> & con,
                    const Ptr<const ROL::Vector<Real>> & ui,
                    size_type Nt) 
    : con_( con )
    , ui_( ui )
    , Nt_( Nt ) {
  }

  void setDynamicConstraint( const Ptr<DynamicConstraint<Real>> & con)
  { con_ = con; }

  void setInitialCondition( const Ptr<const ROL::Vector<Real>> & ui)
  { ui_ = ui; }

  void setNumberOfTimeSteps(size_type Nt)
  { Nt_ = Nt; }

  void setSkipInitialCondition(bool skip)
  { skipInitialCond_ = skip; }

  void solve( V& c, V& u, const V& z, Real& tol ) override { 

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
 
    // First interval solve uses initial condition
    if(!skipInitialCond_) 
      con_->solve( *(cp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );
    
    for( size_type k=1; k<Nt_; ++k ) 
      con_->solve( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
  }  

  void value( V& c, const V& u, const V& z, Real& tol ) override {

    auto& cp = partition(c);
    auto& up = partition(u);
    auto& zp = partition(z);
  
    // First interval value uses initial condition
    if(!skipInitialCond_) 
      con_->value( *(cp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );
    
    for( size_type k=1; k<Nt_; ++k ) 
      con_->value( *(cp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
  } 

  void applyJacobian_1( V& jv, const V& v, const V& u, 
                       const V& z, Real& tol ) override {

    auto& jvp = partition(jv);   auto& vp = partition(v);
    auto& up  = partition(u);    auto& zp = partition(z);

    auto  tmp  = workspace_.clone(jv);
    auto& tmpp = partition(*tmp); 

    if(!skipInitialCond_) 
      con_->applyJacobian_un( *(jvp.get(0)),  *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );
    
    for( size_type k=1; k<Nt_; ++k ) {
      con_->applyJacobian_uo( *(tmpp.get(k)), *(vp.get(k-1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
      con_->applyJacobian_un( *(jvp.get(k)),  *(vp.get(k)),   *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
    }

    jvp.plus(*tmp);
  }

   void applyJacobian_2( V& jv, const V& v, const V& u,
                         const V &z, Real &tol ) override { 

     auto& jvp = partition(jv);   auto& vp = partition(v);
     auto& up  = partition(u);    auto& zp = partition(z);

     if(!skipInitialCond_) 
       con_->applyJacobian_z( *(jvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );

     for( size_type k=1; k<Nt_; ++k ) 
       con_->applyJacobian_z( *(jvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
   }


   void applyInverseJacobian_1( V& ijv, const V& v, const V& u,
                                const V& z, Real& tol) override {

     auto& ijvp = partition(ijv);  auto& vp = partition(v);
     auto& up   = partition(u);    auto& zp = partition(z);

     auto  tmp  = workspace_.clone(ijv);
     auto& tmpp = partition(*tmp); 

     con_->applyInverseJacobian_un( *(ijvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );

     for( size_type k=1; k<Nt_; ++k ) {
       con_->applyJacobian_uo( *(tmpp.get(k)), *(ijvp.get(k-1)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
       tmpp.get(k)->scale(-1.0);
       tmpp.get(k)->plus( *(vp.get(k)) );
       con_->applyInverseJacobian_un( *(ijvp.get(k)), *(tmpp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
     }
     
   }


   void applyAdjointJacobian_1( V& ajv, const V& v, const V& u,
                                const V& z, const V& dualv, Real& tol) override {

     auto& ajvp  = partition(ajv);   auto& vp = partition(v);
     auto& up    = partition(u);     auto& zp = partition(z);

     auto  tmp  = workspace_.clone(ajv); 
     auto& tmpp = partition(*tmp); 

     if(!skipInitialCond_) 
       con_->applyAdjointJacobian_un( *(ajvp.get(0)),  *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );
    
     for( size_type k=1; k<Nt_; ++k ) {
       con_->applyAdjointJacobian_un( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
       con_->applyAdjointJacobian_uo( *(tmpp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );
       ajvp.get(k-1)->plus(*(tmpp.get(k)));
     }
     
   }

   void applyAdjointJacobian_2( V& ajv,  const V& v, const V& u,
                                const V& z, Real& tol ) override {

     auto& ajvp  = partition(ajv);   auto& vp = partition(v);
     auto& up    = partition(u);     auto& zp = partition(z);

     if(!skipInitialCond_) 
       con_->applyAdjointJacobian_z( *(ajvp.get(0)), *(vp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );

     for( size_type k=1; k<Nt_; ++k ) 
       con_->applyAdjointJacobian_z( *(ajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );

   }

   void applyInverseAdjointJacobian_1( V& iajv, const V& v, const V& u,
                                       const V& z, Real& tol) override {

     auto& iajvp = partition(iajv);  auto& vp = partition(v);
     auto& up    = partition(u);     auto& zp = partition(z);

     auto  tmp  = workspace_.clone(iajv);
     auto& tmpp = partition(*tmp); 

     size_type k = Nt_-1;

     con_->applyInverseAdjointJacobian_un( *(iajvp.get(k)), *(vp.get(k)), *(up.get(k-1)), *(up.get(k)), *(zp.get(k)), ts_ );

     for( size_type k=Nt_-2; k>0; --k ) {
       con_->applyAdjointJacobian_uo( *(tmpp.get(k)), *(iajvp.get(k+1)), *(up.get(k)), *(up.get(k+1)), *(zp.get(k+1)), ts_ );
       tmpp.get(k)->scale(-1.0);
       tmpp.get(k)->plus( *(vp.get(k) ) );
       con_->applyInverseAdjointJacobian_un( *(iajvp.get(k)), *(tmpp.get(k)), *(up.get(k)), *(up.get(k+1)), *(zp.get(k+1)), ts_ );             
     }

     con_->applyAdjointJacobian_uo( *(tmpp.get(0)), *(iajvp.get(1)), *(up.get(0)), *(up.get(1)), *(zp.get(1)), ts_ );
     tmpp.get(0)->scale(-1.0);
     tmpp.get(0)->plus( *(vp.get(0) ) );
     con_->applyInverseAdjointJacobian_un( *(iajvp.get(0)), *(tmpp.get(0)), *ui_, *(up.get(0)), *(zp.get(0)), ts_ );           
     
   }


}; // Tanks::ConstraintSerial


} // namespace Tanks



#endif // TANKS_CONSTRAINTSERIAL_HPP


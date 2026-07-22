// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKS_DYNAMICOBJECTIVE_HPP
#define TANKS_DYNAMICOBJECTIVE_HPP

#include "ROL_ParameterList.hpp"
#include "ROL_DynamicObjective.hpp"

#include "Tanks_StateVector.hpp"
#include "Tanks_ControlVector.hpp"
#include "LowerBandedMatrix.hpp"

#include <utility>

/** \class Tanks_DynamicObjective based on the new DynamicObjective interface
*/

namespace Tanks {

template<typename Real>
class DynamicObjective : public ROL::DynamicObjective<Real> {

  using State     = StateVector<Real>;
  using Control   = ControlVector<Real>;
  using Matrix    = LowerBandedMatrix<Real>;
  using V         = ROL::Vector<Real>;
  using TS        = ROL::TimeStamp<Real>;
  using size_type = typename State::size_type;

private:
  size_type rows_;             // Number of tank rows
  size_type cols_;             // Number of tank columns

  //---------- Time Discretization ---------------------------------------------
  Real T_;                     // Total time
  Real theta_;                 // Implicit/Explicit splitting factor
  int  Nt_;                    // Number of Time steps
  Real dt_;                    // Time step size

  size_type    Ntanks_;        // Total number of tanks 

  Real htarg_;

  //--------- Subvector addressing ---------------------------------------------
  size_type  h_, Qout_, Qin_,  z_;

public: 

  DynamicObjective( ROL::ParameterList& pl ) :
  // ----------- Begin Initializer List ----------------//
  rows_   ( pl.get( "Number of Rows",        3      ) ),
  cols_   ( pl.get( "Number of Columns",     3      ) ),
  T_      ( pl.get( "Total Time",            20.0   ) ),
  Nt_     ( pl.get( "Number of Time Steps",  100    ) ),
  //----------------------------------------------------//
  dt_     ( T_/Nt_                                    ),
  Ntanks_ (rows_*cols_                                ),
  htarg_  ( pl.get( "Target Fluid Level",    3.0    ) ),
  h_      (0                                          ),
  Qout_   (2*Ntanks_                                  ),
  Qin_    (Ntanks_                                    ),
  z_      (0                                          )
  // ------------- End Initializer List ----------------//
  {}

  Real value( const V& uo, const V& un, const V& z, const TS& timeStamp ) const {
    auto& uo_state = to_state(uo);
    auto& un_state = to_state(un);
    Real result = 0;
    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        Real hdiff = un_state.h(i,j)-htarg_;
        result += hdiff*hdiff;
      }
    }
    if( timeStamp.k > 0 ) {
      for( size_type i=0; i<rows_; ++i ) {
        for( size_type j=0; j<cols_; ++j ) {
          Real hdiff = uo_state.h(i,j)-htarg_;
          result += hdiff*hdiff;
        }
      }
    }
    return static_cast<Real>(0.5)*dt_*result;
  }

  //----------------------------------------------------------------------------
  // Gradient Terms
  void gradient_uo( V& g, const V& uo, const V& un,
                          const V& z, const TS& timeStamp ) const {
    if( timeStamp.k > 0 ) {
      auto& g_state  = to_state(g);
      auto& uo_state = to_state(uo);
      for( size_type i=0; i<rows_; ++i ) {
        for( size_type j=0; j<cols_; ++j ) {
          g_state.h(i,j)    = dt_*(uo_state.h(i,j)-htarg_);
          g_state.Qin(i,j)  = static_cast<Real>(0);
          g_state.Qout(i,j) = static_cast<Real>(0);
        }
      }
    }
    else {
      g.zero();
    }
  }

  void gradient_un( V& g, const V& uo, const V& un,
                          const V& z, const TS& timeStamp ) const {
    auto& g_state  = to_state(g);
    auto& un_state = to_state(un);
    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        g_state.h(i,j)    = dt_*(un_state.h(i,j)-htarg_);
        g_state.Qin(i,j)  = static_cast<Real>(0);
        g_state.Qout(i,j) = static_cast<Real>(0);
      }
    }
  }

  void gradient_z( V& g, const V& uo, const V& un,
                         const V& z, const TS& timeStamp ) const {
    g.zero();
  }

  //----------------------------------------------------------------------------
  // Hessian-Vector product terms
  void hessVec_uo_uo( V& hv, const V& v, const V& uo, const V& un,
                             const V& z, const TS& timeStamp ) const {
    if( timeStamp.k > 0 ) {
      auto& hv_state = to_state(hv);
      auto& v_state  = to_state(v);
      for( size_type i=0; i<rows_; ++i ) {
        for( size_type j=0; j<cols_; ++j ) {
          hv_state.h(i,j)    = dt_*v_state.h(i,j);
          hv_state.Qin(i,j)  = static_cast<Real>(0);
          hv_state.Qout(i,j) = static_cast<Real>(0);
        }
      }
    }
    else {
      hv.zero();
    }
  }

  void hessVec_un_un( V& hv, const V& v, const V& uo, const V& un,
                             const V& z, const TS& timeStamp ) const {
    auto& hv_state = to_state(hv);
    auto& v_state  = to_state(v);
    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        hv_state.h(i,j)    = dt_*v_state.h(i,j);
        hv_state.Qin(i,j)  = static_cast<Real>(0);
        hv_state.Qout(i,j) = static_cast<Real>(0);
      }
    }
  }

  void hessVec_z_z( V& hv, const V& v, const V& uo, const V& un,
                           const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  void hessVec_uo_un( V& hv, const V& v, const V& uo, const V& un,
                             const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  void hessVec_uo_z( V& hv, const V& v, const V& uo, const V& un,
                            const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  void hessVec_un_uo( V& hv, const V& v, const V& uo, const V& un,
                             const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  void hessVec_un_z( V& hv, const V& v, const V& uo, const V& un,
                            const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  void hessVec_z_uo( V& hv, const V& v, const V& uo, const V& un,
                            const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

  void hessVec_z_un( V& hv, const V& v, const V& uo, const V& un,
                            const V& z, const TS& timeStamp ) const {
    hv.zero();
  }

}; // Tanks::DynamicObjective

} // namespace Tanks

#endif // TANKS_DYNAMICOBJECTIVE_HPP

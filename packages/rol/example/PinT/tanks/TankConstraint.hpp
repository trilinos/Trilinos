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
#ifndef TANKCONSTRAINT_HPP
#define TANKCONSTRAINT_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "ROL_Constraint_TimeSimOpt.hpp"
#include "ROL_StdVector.hpp"

#include "TankState.hpp"

#include <utility>

/** \class TankConstraint
    \brief Compute time-step for the coupled tank network

*/

namespace details {

using namespace std;

using ROL::Vector;

template<typename Real>
class TankConstraint : public ROL::Constraint_TimeSimOpt<Real> {
private:

 // Tank Array Numbers
 int numRows_;          // Number of tank rows
 int numCols_;          // Number of tank columns
 vector<int> pt_rows_;  // Pass-through rows
 vector<int> pt_cols_;  // Pass-through columns

 //---------- Physical Parameters ---------------------------------------------
 Real Cv_;              // Valve constant 
 Real rho_;             // Density of fluid
 Real h0_;              // Initial fluid level
 Real H_;               // Height of tanks
 Real A_;               // Cross-sectional area of tanks
 Real g_;               // Gravity constant
 Real Qin00_;           // Corner inflow 


 //---------- Time Discretization ---------------------------------------------
 Real T_;               // Total time
 Real dt_;              // Time step size
 Real theta_;           // Implicit/Explicit splitting factor
 int  Nt_;              // Number of Time steps


 //---------- Lumped Coefficients ---------------------------------------------
 Real a1_;              // dt*(1-theta)/A
 Real a2_;              // dt*theta/A
 Real b_;               // Cv*rho*g
 
private:

  auto tank_state( const Vector<Real>& x ) -> TankState<Real,true> {
    return move( make_tank_state( x, numRows_, numCols_ ) );
  }

  auto tank_state( Vector<Real>& x ) -> TankState<Real,false> {
    return move( make_tank_state( x, numRows_, numCols_, pt_rows_, pt_cols_ ) );
  }

public:

  TankConstraint( Teuchos::ParameterList& pl ) :
    // ---------------- Initializer List -----------------//
    numRows_( pl.get( "Number of Rows",          3      ) ),
    numCols_( pl.get( "Number of Columns",       3      ) ),
    Cv_     ( pl.get( "Valve Constant",          1.0e-2 ) ),
    rho_    ( pl.get( "Density of Fluid",        1.0e3  ) ),
    h0_     ( pl.get( "Initial Fluid Level",     2.0    ) ),
    H_      ( pl.get( "Height of Tank",          10.0   ) ),
    A_      ( pl.get( "Cross-sectional Area",    10.0   ) ),
    g_      ( pl.get( "Gravity Constant",        9.8    ) ),
    Qin00_  ( pl.get( "Corner Inflow",           100.0  ) ),
    T_      ( pl.get( "Total Time",              20.0   ) ),
    theta_  ( pl.get( "Theta",                   0.5    ) ),
    Nt_     ( pl.get( "Number of Time Steps",    100    ) ) {
    //----------------------------------------------------//

    auto& ptrows = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Rows"    );
    auto& ptcols = Teuchos::getArrayFromStringParameter<int>( pl, "Pass-Through Columns" );

    pt_rows_ = ptrows.toVector();
    pt_cols_ = ptcols.toVector();

    dt_ = T_/Nt_;
    a1_ = dt_*(1.0-theta_)/A_;
    a2_ = dt_*theta_/A_;
    b_  = Cv_*rho_*g_;

  } // end Constructor


  void value(Vector<Real> &c, const Vector<Real> &u_old,
             const Vector<Real> &u_new, const Vector<Real> &z,
             Real &tol) override {

    auto c_s  = tank_state(c);
    auto uo_s = tank_state(u_old);
    auto un_s = tank_state(u_new);


   

  }

  void solve(Vector<Real> &c, const Vector<Real> &u_old, 
             Vector<Real> &u_new, const Vector<Real> &z, 
                     Real &tol) override {

    auto c_s  = tank_state(c);
    auto uo_s = tank_state(u_old);
    auto un_s = tank_state(u_new);




  }





}; // class TankConstraint

} // namespace details

using details::TankConstraint;

#endif // TANKCONSTRAINT_HPP


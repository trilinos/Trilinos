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
#ifndef TANKSTATE_HPP
#define TANKSTATE_HPP

#include <memory>

#include "Teuchos_ParameterList.hpp"
#include "LowerBandedMatrix.hpp"


namespace details {

using namespace std;



template<typename Real> 
class TankState {
public:

  using Vector     = ROL::Vector<Real>;
  using StdVector  = ROL::StdVector<Real>;
  using Matrix     = LowerBandedMatrix<Real>;
  using size_type  = typename vector<Real>::size_type;
   
  TankState( Teuchos::ParameterList& pl );

  void value( vector<Real>& c, const vector<Real>* u_old, 
              const vector<Real>& u_new, const vector<Real>& z ) const;

  // Given a vector u=(h,Qin,Qout) and f, update Qin and Qout from h
  void compute_flow( vector<Real>& u, const vector<Real>& f ) const;


  void solve_level(       vector<Real>& c,           
                          vector<Real>& u_new, 
                    const vector<Real>& u_old, 
                    const vector<Real>& f     ) const;

  // Subvector Accessor Methods
        Real& h( vector<Real>& x,       size_type r, size_type c ) const { return x[cols_*r+c]; }
  const Real& h( const vector<Real>& x, size_type r, size_type c ) const { return x[cols_*r+c]; }

        Real& Qout(       vector<Real>& x, size_type r, size_type c ) const { return x[Ntanks_+cols_*r+c]; }
  const Real& Qout( const vector<Real>& x, size_type r, size_type c ) const { return x[Ntanks_+cols_*r+c]; }
  
        Real& Qin(       vector<Real>& x, size_type r, size_type c ) const { return x[2*Ntanks_+cols_*r+c]; }
  const Real& Qin( const vector<Real>& x, size_type r, size_type c ) const { return x[2*Ntanks_+cols_*r+c]; }

private:
  size_type rows_;             // Number of tank rows
  size_type cols_;             // Number of tank columns

  //---------- Physical Parameters ---------------------------------------------
  Real Cv_;                    // Valve constant 
  Real rho_;                   // Density of fluid
  Real h0_;                    // Initial fluid level
  Real H_;                     // Height of tanks
  Real A_;                     // Cross-sectional area of tanks
  Real g_;                     // Gravity constant
  
  //---------- Time Discretization ---------------------------------------------
  Real T_;                     // Total time
  Real dt_;                    // Time step size
  Real theta_;                 // Implicit/Explicit splitting factor
  int  Nt_;                    // Number of Time steps

  size_type    Ntanks_;        // Total number of tanks 
  vector<Real> p_;             // Passthrough coefficients
 
  Real         coeff1_;        // Cv*rho*g
  Real         kappa_;         // Cv*rho*g*dt/2A
  Real         betaL_;         // (theta-1)*dt/A
  Real         betaR_;         // theta*dt/A

  shared_ptr<Matrix> L_, R_;

}; // class TankState


} // namespace details

using details::TankState;

#include "TankState_Impl.hpp"

#endif // TANKSTATE_HPP


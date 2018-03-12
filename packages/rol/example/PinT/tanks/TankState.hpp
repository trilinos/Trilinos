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

#include "ROL_Ptr.hpp"
#include <utility>

/** \class TankState 
    \brief Provides a convenient access to the blocks and elements of the
           tank system state vector where VectorImpl<Real>::getVector()
           returns a pointer to a subscriptable container
*/

namespace details {

using namespace std;

using ROL::Vector;
using ROL::StdVector;

// Casting helpers
template<typename Real, template<typename> class VectorType>
inline
// const // Certain versions of Intel compiler object to the const modifier
auto getVector( const Vector<Real>& x ) -> 
  const decltype(*declval<const VectorType<Real>>().getVector())& {
  return *(dynamic_cast<const VectorType<Real>&>(x).getVector());
}

template<typename Real, template<typename> class VectorType>
inline auto getVector( Vector<Real>& x ) -> 
  decltype(*declval<VectorType<Real>>().getVector())& {
  return *(dynamic_cast<VectorType<Real>&>(x).getVector());
}


template<typename Real, bool isConst> class TankState;

template<typename Real>
class TankState<Real,true> {
private:

  using container_type = vector<Real>;
  const container_type& state_;
  int r, c, N;

public:

  TankState( const Vector<Real>& state, int rows, int cols ) : 
    state_( getVector<StdVector>(state) ), r(rows), c(cols), N(r*c) {}

  // Row-column accessors
  const Real& h   ( int i, int j ) const { return state_[    c*j+i]; }  // Level
  const Real& Qin ( int i, int j ) const { return state_[  N+c*j+i]; }  // Inflow
  const Real& Qout( int i, int j ) const { return state_[2*N+c*j+i]; }  // Outflow

  // Get the row and column indices from the single index
  void getRowCol( const int k, int& i, int& j ) const {
    j = k % c;  i = (k - j)/c;
  }

}; // class TankState


template<typename Real> // Non-constant specialization
class TankState<Real,false> : public TankState<Real,true> {
private:
  using container_type = vector<Real>;
  container_type& state_;
  int r,c,N;

public:

  TankState( Vector<Real>& state, int rows, int cols ) :
    TankState<Real,true>::TankState(state,rows,cols), 
    state_( getVector<StdVector>(state) ), r(rows), c(cols), N(r*c) {}

  Real& h   ( int i, int j ) { return state_[    c*j+i]; }  // Level
  Real& Qin ( int i, int j ) { return state_[  N+c*j+i]; }  // Inflow
  Real& Qout( int i, int j ) { return state_[2*N+c*j+i]; }  // Outflow

}; // class TankState

template<typename Real>
auto make_tank_state( const Vector<Real>& x, int r, int c ) -> TankState<Real,true> {
  TankState<Real,true> ts( x, r, c );
  return move(ts);
}

template<typename Real>
auto make_tank_state( Vector<Real>& x, int r, int c ) -> TankState<Real,false> {
  TankState<Real,false> ts( x, r, c );
  return move(ts);
}

} // namespace details

using details::TankState;
using details::make_tank_state;

#endif // TANKSTATE_HPP


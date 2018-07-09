// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (1614) Sandia Corporation
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
#ifndef TANKS_STATEVECTOR_HPP
#define TANKS_STATEVECTOR_HPP

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"

#include <iomanip>
#include <iostream>

namespace Tanks {

using namespace std;

// Forward declaration
template<typename> class ControlVector;

template<typename Real>
class StateVector : public ROL::StdVector<Real> {
public:
  using size_type = typename vector<Real>::size_type;

private:
  ROL::Ptr<vector<Real>> vec_;
  size_type rows_,cols_, N_;
  string name_;

public:
  using ROL::StdVector<Real>::getVector;

  StateVector( const ROL::Ptr<vector<Real>>& vec, size_type rows, 
              size_type cols, const string& name="anonynous" ) : 
    ROL::StdVector<Real>(vec), vec_(vec), rows_(rows), cols_(cols), N_(rows*cols), name_(name) {
  }

  StateVector( size_type rows, size_type cols, const string& name="anonynous" ) : 
    StateVector(ROL::makePtr<vector<Real>>(3*rows*cols),rows,cols,name) { 
  }

  static ROL::Ptr<StateVector> create( ROL::ParameterList& pl, const string& name="anonymous" ) {
    auto nrows  = static_cast<size_type>( pl.get("Number of Rows",3) );
    auto ncols  = static_cast<size_type>( pl.get("Number of Columns",3) );
    return ROL::makePtr<StateVector>( nrows, ncols, name );
  }

  ROL::Ptr<ROL::Vector<Real>> clone() const { 
    return ROL::makePtr<StateVector>( rows_, cols_, "clone of " + name_ ); 
  }

  ROL::Ptr<StateVector<Real>> clone( const string& name ) const { 
    return ROL::makePtr<StateVector>( rows_, cols_, name ); 
  }

  Real& h( size_type i, size_type j ) { return vec_->at(cols_*i+j); }

  const Real& h( size_type i, size_type j ) const { return vec_->at(cols_*i+j); }

  Real& Qin( size_type i, size_type j ) { return vec_->at(cols_*i+j+N_); }
  const Real& Qin( size_type i, size_type j ) const { return vec_->at(cols_*i+j+N_); }

  Real& Qout( size_type i, size_type j ) { return vec_->at(cols_*i+j+2*N_); }
  const Real& Qout( size_type i, size_type j ) const { return vec_->at(cols_*i+j+2*N_); }

  using ROL::StdVector<Real>::axpy;
  using ROL::StdVector<Real>::set;

  void set( const StateVector& x, size_type begin, size_type xbegin );
  void set( const ControlVector<Real>& x, size_type begin );
  void axpy( Real alpha, const StateVector& x, size_type begin, size_type xbegin );
  void axpy( Real alpha, const ControlVector<Real>& x, size_type begin );
  void hadamard( const StateVector& x, size_type begin, size_type xbegin );
  void hadamard( const ControlVector<Real>& x, size_type begin );

  void print( ostream& os ) const {
    auto& u = *vec_;
    os << endl << name_ << endl;
    os << setw(5) << "row" << setw(5) << "col" << setw(7) << "   |   " << 
          setw(18) << "h"  << setw(18) << "Qin" << setw(18) << "Qout" << endl;
    os << string(13,'-') + '+' + string(60,'-') << endl;

    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        size_type k = cols_*i+j;
        os << setw(5)  << i << setw(5) << j << setw(7) << "   |   " << 
              setw(18) << u[k] << setw(18) << u[k+N_] << setw(18) << u[k+2*N_] << endl;
      }
    }
  }
}; // StateVector

//------------------------------------------------------------------------------

template<typename Real> 
inline StateVector<Real>& to_state( ROL::Vector<Real>& x ) {
  return dynamic_cast<StateVector<Real>&>(x);
}

template<typename Real> 
inline const StateVector<Real>& to_state( const ROL::Vector<Real>& x ) {
  return dynamic_cast<const StateVector<Real>&>(x);
}

template<typename Real> 
inline const ROL::Ptr<StateVector<Real>>& to_state( const ROL::Ptr<ROL::Vector<Real>>& x ) {
  return ROL::dynamicPtrCast<StateVector<Real>>(x);
}

template<typename Real> 
inline const ROL::Ptr<const StateVector<Real>>& to_state( const ROL::Ptr<const ROL::Vector<Real>>& x ) {
  return ROL::dynamicPtrCast<const StateVector<Real>>(x);
}

//------------------------------------------------------------------------------

} // namespace Tanks


#include "Tanks_StateVector_Impl.hpp"

#endif // TANKS_STATEVECTOR_HPP


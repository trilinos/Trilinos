// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKS_CONTROLVECTOR_HPP
#define TANKS_CONTROLVECTOR_HPP

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"

#include <iomanip>
#include <iostream>

namespace Tanks {

using namespace std;

// Forward declaration 
template<typename> class StateVector;

template<typename Real>
class ControlVector : public ROL::StdVector<Real> {
public:
  using size_type = typename vector<Real>::size_type;

private:
  ROL::Ptr<vector<Real>> vec_;
  size_type rows_,cols_, N_;
  string name_;
 
public:
  using ROL::StdVector<Real>::getVector;

  ControlVector( const ROL::Ptr<vector<Real>>& vec, size_type rows, 
              size_type cols, const string& name="anonymous" ) : 
    ROL::StdVector<Real>(vec), vec_(vec), rows_(rows), cols_(cols), name_(name) {}

  ControlVector( size_type rows, size_type cols, const string& name="anonymous" ) : 
    ControlVector(ROL::makePtr<vector<Real>>(rows*cols),rows,cols,name) {}

  static ROL::Ptr<ControlVector> create( ROL::ParameterList& pl, const string& name="anonymous" ) {
    auto nrows  = static_cast<size_type>( pl.get("Number of Rows",3) );
    auto ncols  = static_cast<size_type>( pl.get("Number of Columns",3) );
    return ROL::makePtr<ControlVector>( nrows, ncols, name );
  }


  ROL::Ptr<ROL::Vector<Real>> clone() const { 
    return ROL::makePtr<ControlVector>( rows_, cols_, "clone of " + name_ ); 
  }

  ROL::Ptr<ControlVector> clone( const string& name ) const { 
    return ROL::makePtr<ControlVector>( rows_, cols_, name ); 
  }
  
  Real& operator[]( size_type k ) { return vec_->at(k); }
  const Real& operator[]( size_type k ) const { return vec_->at(k); }
  Real& operator()( size_type i, size_type j ) { return vec_->at(cols_*i+j); }
  const Real& operator()( size_type i, size_type j ) const { return vec_->at(cols_*i+j); }

  using ROL::StdVector<Real>::axpy;
  using ROL::StdVector<Real>::set;

  void set( const StateVector<Real>& x, size_type xbegin );
  void set( const ControlVector& x );
  void axpy( Real alpha, const StateVector<Real>& x, size_type xbegin );
  void axpy( Real alpha, const ControlVector& x );
  void hadamard( const StateVector<Real>& x, size_type xbegin );
  void hadamard( const ControlVector& x );
  
  void print( ostream& os ) const {
    auto& z = *vec_;

    os << endl << name_ << endl;
    os << setw(5) << "row" << setw(5) << "col" << setw(7) << "   |   " << 
          setw(18) << "value"  << endl;
    os << string(13,'-') + '+' + string(60,'-') << endl;

    for( size_type i=0; i<rows_; ++i ) {
      for( size_type j=0; j<cols_; ++j ) {
        size_type k = cols_*i+j;
        os << setw(5)  << i << setw(5) << j << setw(7) << "   |   " << 
              setw(16) << z[k] << endl;
      }
    }
  }

}; // StateVector


//------------------------------------------------------------------------------

template<typename Real> 
inline ControlVector<Real>& to_control( ROL::Vector<Real>& x ) {
  return dynamic_cast<ControlVector<Real>&>(x);
}

template<typename Real> 
inline const ControlVector<Real>& to_control( const ROL::Vector<Real>& x ) {
  return dynamic_cast<const ControlVector<Real>&>(x);
}

template<typename Real> 
inline const ROL::Ptr<ControlVector<Real>>& to_control( const ROL::Ptr<ROL::Vector<Real>>& x ) {
  return ROL::dynamicPtrCast<ControlVector<Real>>(x);
}

template<typename Real> 
const ROL::Ptr<const ControlVector<Real>>& to_control( const ROL::Ptr<const ROL::Vector<Real>>& x ) {
  return ROL::dynamicPtrCast<const ControlVector<Real>>(x);
}

//------------------------------------------------------------------------------


} // namespace Tanks

#include "Tanks_ControlVector_Impl.hpp"

#endif // TANKS_CONTROLVECTOR_HPP


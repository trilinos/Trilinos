// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef CONTROLVECTORIMPL_HPP
#define CONTROLVECTORIMPL_HPP

namespace Tanks {

using namespace std;

using size_type = typename vector<double>::size_type;


template<typename Real>
void ControlVector<Real>::set( const StateVector<Real>& x, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) = xp->at(xbegin+i);
}
  
template<typename Real>
void ControlVector<Real>::set( const ControlVector<Real>& x ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) = xp->at(i);
}

template<typename Real>
void ControlVector<Real>::axpy( Real alpha, const StateVector<Real>& x, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) += alpha*xp->at(xbegin+i);
}

template<typename Real>
void ControlVector<Real>::axpy( Real alpha, const ControlVector<Real>& x ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) += alpha*xp->at(i);
}

template<typename Real>
void ControlVector<Real>::hadamard( const StateVector<Real>& x, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) *= xp->at(xbegin+i);
}
  
template<typename Real>
void ControlVector<Real>::hadamard( const ControlVector<Real>& x ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) *= xp->at(i);
}

} // namespace Tanks

#endif // CONTROLVECTOR_IMPL_HPP


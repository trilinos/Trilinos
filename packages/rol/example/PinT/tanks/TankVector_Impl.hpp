// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef TANKVECTORIMPL_HPP
#define TANKVECTORIMPL_HPP

namespace details {

using namespace std;


using size_type = typename vector<double>::size_type;


template<typename Real>
void TankStateVector<Real>::set( const TankStateVector<Real>& x, 
                                 size_type begin, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) = xp->at(xbegin+i);
}

template<typename Real>
void TankStateVector<Real>::set( const TankControlVector<Real>& x, size_type begin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) = xp->at(i);
}

template<typename Real>
void TankStateVector<Real>::axpy( Real alpha, const TankStateVector<Real>& x, 
                                  size_type begin, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i )  vec_->at(begin+i) += alpha*xp->at(xbegin+i);
}

template<typename Real>
void TankStateVector<Real>::axpy( Real alpha, const TankControlVector<Real>& x, size_type begin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) += alpha*xp->at(i);
}

template<typename Real>
void TankStateVector<Real>::hadamard( const TankStateVector<Real>& x, 
                                 size_type begin, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) *= xp->at(xbegin+i);
}

template<typename Real>
void TankStateVector<Real>::hadamard( const TankControlVector<Real>& x, size_type begin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) *= xp->at(i);
}



template<typename Real>
void TankControlVector<Real>::set( const TankStateVector<Real>& x, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) = xp->at(xbegin+i);
}
  
template<typename Real>
void TankControlVector<Real>::set( const TankControlVector<Real>& x ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) = xp->at(i);
}

template<typename Real>
void TankControlVector<Real>::axpy( Real alpha, const TankStateVector<Real>& x, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) += alpha*xp->at(xbegin+i);
}

template<typename Real>
void TankControlVector<Real>::axpy( Real alpha, const TankControlVector<Real>& x ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) += alpha*xp->at(i);
}

template<typename Real>
void TankControlVector<Real>::hadamard( const TankStateVector<Real>& x, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) *= xp->at(xbegin+i);
}
  
template<typename Real>
void TankControlVector<Real>::hadamard( const TankControlVector<Real>& x ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(i) *= xp->at(i);
}

} // namespace details
#endif // TANKVECTOR_IMPL_HPP


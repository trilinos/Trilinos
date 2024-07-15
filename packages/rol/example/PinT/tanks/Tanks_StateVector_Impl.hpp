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

namespace Tanks {

using namespace std;

using size_type = typename vector<double>::size_type;

template<typename Real>
void StateVector<Real>::set( const StateVector<Real>& x, 
                                 size_type begin, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) = xp->at(xbegin+i);
}

template<typename Real>
void StateVector<Real>::set( const ControlVector<Real>& x, size_type begin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) = xp->at(i);
}

template<typename Real>
void StateVector<Real>::axpy( Real alpha, const StateVector<Real>& x, 
                                  size_type begin, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i )  vec_->at(begin+i) += alpha*xp->at(xbegin+i);
}

template<typename Real>
void StateVector<Real>::axpy( Real alpha, const ControlVector<Real>& x, size_type begin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) += alpha*xp->at(i);
}

template<typename Real>
void StateVector<Real>::hadamard( const StateVector<Real>& x, 
                                 size_type begin, size_type xbegin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) *= xp->at(xbegin+i);
}

template<typename Real>
void StateVector<Real>::hadamard( const ControlVector<Real>& x, size_type begin ) {
 auto xp = x.getVector(); 
 for( size_type i=0; i<N_; ++i ) vec_->at(begin+i) *= xp->at(i);
}

} // Tanks

#endif // STATEVECTOR_IMPL_HPP


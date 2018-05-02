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


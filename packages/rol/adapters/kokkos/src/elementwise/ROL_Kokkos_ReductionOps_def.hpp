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
#ifndef ROL_KOKKOS_REDUCTIONOPS_DEF_HPP
#define ROL_KOKKOS_REDUCTIONOPS_DEF_HPP

namespace ROL {
namespace Elementwise {

template<typename Real, typename Device>
std::unique_ptr<KokkosReductionOp<Real,Device>>
KokkosReductionOp<Real,Device>::create( const ::ROL::Elementwise::ReductionOp<Real>& rop ) {
  KokkosReductionOp<Real,Device>::Factory factory;
  rop.accept(factory);
  return std::move(factory.rop_);
}

//--------------------------------------------------------------------------------
// Euclidean Norm Squared

template<typename Real, typename Device>
Real KokkosEuclideanNormSquared<Real,Device>::reduce( const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Real result = 0;
  Kokkos::parallel_reduce( x.dimension(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
    lresult += d(i)*d(i);
  }, result );
  return result;
}

template<typename Real, typename Device>
void KokkosReductionOp<Real,Device>::Factory::visit( const ::ROL::Elementwise::EuclideanNormSquared<Real>& ) {
  rop_ = std::unique_ptr<KokkosEuclideanNormSquared<Real,Device>>( new KokkosEuclideanNormSquared<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Logical And

//template<typename Real, typename Device>
//Real KokkosReductionAnd<Real,Device>::reduce( const ::ROL::KokkosVector<Real,Device>& x ) const {
//  LAndFunctor f;
//  f.values = x.view_device();
//  Real init{1}; 
//  Real land_scalar = init;
//  Kokkos::LAnd<Real> reducer_scalar(land_scalar);
//  Kokkos::parallel_reduce( Kokkos::RangePolicy<Device>(0,x.dimension()), f, reducer_scalar );
//  return land_scalar;
//}
//
template<typename Real, typename Device>
void KokkosReductionOp<Real,Device>::Factory::visit( const ::ROL::Elementwise::ReductionAnd<Real>& ) {
//  rop_ = makePtr<KokkosReductionAnd<Real,Device>>();
}

//--------------------------------------------------------------------------------
// Maximum Value

template<typename Real, typename Device>
Real KokkosReductionMax<Real,Device>::reduce( const ::ROL::KokkosVector<Real,Device>& x ) const {
  MaxFunctor f;
  f.values = x.view_device();
  Real init = std::numeric_limits<Real>::min();
  Real max_scalar = init;
  Kokkos::Max<Real> reducer_scalar(max_scalar);
  Kokkos::parallel_reduce( Kokkos::RangePolicy<Device>(0,x.dimension()), f, reducer_scalar );
  return max_scalar;
}

template<typename Real, typename Device>
void KokkosReductionOp<Real,Device>::Factory::visit( const ::ROL::Elementwise::ReductionMax<Real>& ) {
  rop_ = std::unique_ptr<KokkosReductionMax<Real,Device>>( new KokkosReductionMax<Real,Device>() );
}


//--------------------------------------------------------------------------------
// Minimum Value

template<typename Real, typename Device>
Real KokkosReductionMin<Real,Device>::reduce( const ::ROL::KokkosVector<Real,Device>& x ) const {
  MinFunctor f;
  f.values = x.view_device();
  Real init = std::numeric_limits<Real>::max();
  Real min_scalar = init;
  Kokkos::Min<Real> reducer_scalar(min_scalar);
  Kokkos::parallel_reduce( Kokkos::RangePolicy<Device>(0,x.dimension()), f, reducer_scalar );
  return min_scalar;
}

template<typename Real, typename Device>
void KokkosReductionOp<Real,Device>::Factory::visit( const ::ROL::Elementwise::ReductionMin<Real>& ) {
  rop_ = std::unique_ptr<KokkosReductionMin<Real,Device>>( new KokkosReductionMin<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Sum All Elements 

template<typename Real, typename Device>
Real KokkosReductionSum<Real,Device>::reduce( const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Real result = 0;
  Kokkos::parallel_reduce( x.dimension(), KOKKOS_LAMBDA( const int i, Real& lresult ) {
    lresult += d(i);
  }, result );
  return result;
}

template<typename Real, typename Device>
void KokkosReductionOp<Real,Device>::Factory::visit( const ::ROL::Elementwise::ReductionSum<Real>& ) {
  rop_ = std::unique_ptr<KokkosReductionSum<Real,Device>>( new KokkosReductionSum<Real,Device>() );
}



} // namespace Elementwise
} // namespace ROL

#endif  // ROL_KOKKOS_REDUCTIONOPS_DEF_HPP

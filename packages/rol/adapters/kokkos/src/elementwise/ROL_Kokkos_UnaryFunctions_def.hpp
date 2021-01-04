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
#ifndef ROL_KOKKOS_UNARYFUNCTIONS_DEF_HPP
#define ROL_KOKKOS_UNARYFUNCTIONS_DEF_HPP

namespace ROL {

namespace Elementwise {

template<typename Real, typename Device>
std::unique_ptr<KokkosUnaryFunction<Real,Device>>
KokkosUnaryFunction<Real,Device>::create( const ::ROL::Elementwise::UnaryFunction<Real>& uf ) {
  KokkosUnaryFunction<Real,Device>::Factory factory;
  uf.accept(factory);
  return std::move(factory.ufun_);
}

//--------------------------------------------------------------------------------
// Absolute Value Function

template<typename Real, typename Device>
void KokkosAbsoluteValue<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { d(i) = std::abs(d(i)); } );
}


template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::AbsoluteValue<Real>& ) {
  ufun_ = std::unique_ptr<KokkosAbsoluteValue<Real,Device>>(
                         new KokkosAbsoluteValue<Real,Device>());
}

//--------------------------------------------------------------------------------
// Fill Function

template<typename Real, typename Device>
void KokkosFill<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Real value = value_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { d(i) = value; } );
}


template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Fill<Real>& f ) {
  ufun_ = std::unique_ptr<KokkosFill<Real,Device>>(
                         new KokkosFill<Real,Device>(std::forward<Real>(f.get_value())));
}

//--------------------------------------------------------------------------------
// Heaviside Function

template<typename Real, typename Device>
void KokkosHeaviside<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    d(i) = d(i) < 0 ? 0 : ( d(i) > 0 ? 1.0 : 0.5 );} );
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Heaviside<Real>& ) {
  ufun_ = std::unique_ptr<KokkosHeaviside<Real,Device>>(
                         new KokkosHeaviside<Real,Device>());
}
  
//--------------------------------------------------------------------------------
// Natural Logarithm Function

template<typename Real, typename Device>
void KokkosLogarithm<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { d(i) = std::log(d(i)); } );
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Logarithm<Real>& ) {
  ufun_ = std::unique_ptr<KokkosLogarithm<Real,Device>>(new KokkosLogarithm<Real,Device>() );
}
  
//--------------------------------------------------------------------------------
// Power Function

template<typename Real, typename Device>
void KokkosPower<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  // Introduce local copy of member variable as KOKKOS_LAMBDA will not capture 
  // class variables and we do not require C++17 
  Real p = exponent_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( int i ) { d(i) = std::pow(d(i),p); } );
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Power<Real>& pf ) {
  ufun_ = std::unique_ptr<KokkosPower<Real,Device>>(new KokkosPower<Real,Device>(std::forward<Real>(pf.get_exponent())));
}

//--------------------------------------------------------------------------------
// Reciprocal Function

template<typename Real, typename Device>
void KokkosReciprocal<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const { 
  auto d = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { d(i) = 1/d(i); } );
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Reciprocal<Real>& ) {
  ufun_ = std::unique_ptr<KokkosReciprocal<Real,Device>>(new KokkosReciprocal<Real,Device>());
}

//--------------------------------------------------------------------------------
// Round Function

template<typename Real, typename Device>
void KokkosRound<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const { 
  auto d = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    d(i) = std::round(d(i)); 
  });
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Round<Real>& ) {
  ufun_ = std::unique_ptr<KokkosRound<Real,Device>>(new KokkosRound<Real,Device>());
}
  
//--------------------------------------------------------------------------------
// Scale Function

template<typename Real, typename Device>
void KokkosScale<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Real value = value_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { d(i) *= value; } );
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Scale<Real>& f ) {
  ufun_ = std::unique_ptr<KokkosScale<Real,Device>>(new KokkosScale<Real,Device>(std::forward<Real>(f.get_value())));
}

//--------------------------------------------------------------------------------
// Shift Function

template<typename Real, typename Device>
void KokkosShift<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Real value = value_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { d(i) += value; } );
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Shift<Real>& f ) {
  ufun_ = std::unique_ptr<KokkosShift<Real,Device>>(new KokkosShift<Real,Device>(std::forward<Real>(f.get_value())));
}

//--------------------------------------------------------------------------------
// Sign Function

template<typename Real, typename Device>
void KokkosSign<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    d(i) = d(i) > 0 ? Real(1) : Real(-1);
  });
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Sign<Real>& ) {
  ufun_ = std::unique_ptr<KokkosSign<Real,Device>>(new KokkosSign<Real,Device>());
}

//--------------------------------------------------------------------------------
// Square Root Function

template<typename Real, typename Device>
void KokkosSquareRoot<Real, Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { d(i) = std::sqrt(d(i)); } );
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::SquareRoot<Real>& ) {
  ufun_ = std::unique_ptr<KokkosSquareRoot<Real,Device>>(new KokkosSquareRoot<Real,Device>());
}

//--------------------------------------------------------------------------------
// Threshold Lower Function

template<typename Real, typename Device>
void KokkosThresholdLower<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const {
  auto d = x.view_device();
  Real t = threshold_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    d(i) = t < d(i) ? t : d(i);
  });
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::ThresholdLower<Real>& tl ) {
  ufun_ = std::unique_ptr<KokkosThresholdLower<Real,Device>>(new KokkosThresholdLower<Real,Device>(std::forward<Real>(tl.get_threshold())));
}

//--------------------------------------------------------------------------------
// Threshold Upper Function

template<typename Real, typename Device>
void KokkosThresholdUpper<Real,Device>::apply( ::ROL::KokkosVector<Real,Device>& x ) const { 
  auto d = x.view_device();
  Real t = threshold_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    d(i) = t > d(i) ? t : d(i);
  });
}

template<typename Real, typename Device>
void KokkosUnaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::ThresholdUpper<Real>& tu ) {
  ufun_ = std::unique_ptr<KokkosThresholdUpper<Real,Device>>(new KokkosThresholdUpper<Real,Device>(std::forward<Real>(tu.get_threshold())));
}


} // namespace Elementwise
} // namespace ROL
#endif  // ROL_KOKKOS_UNARYFUNCTIONS_DEF_HPP


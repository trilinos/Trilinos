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
#ifndef ROL_KOKKOS_BINARYFUNCTIONS_DEF_HPP
#define ROL_KOKKOS_BINARYFUNCTIONS_DEF_HPP

namespace ROL {
namespace Elementwise {

template<typename Real, typename Device>
std::unique_ptr<KokkosBinaryFunction<Real,Device>>
KokkosBinaryFunction<Real,Device>::create( const ::ROL::Elementwise::BinaryFunction<Real>& bf ) {
  KokkosBinaryFunction<Real,Device>::Factory factory;
  bf.accept(factory);
  return std::move(factory.bfun_);
}

//--------------------------------------------------------------------------------
// Axpy 

template<typename Real, typename Device>
void KokkosAxpy<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                     const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Real alpha = alpha_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { yd(i) += alpha*xd(i); } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Axpy<Real>& axpy ) {
  bfun_ = std::unique_ptr<KokkosAxpy<Real,Device>>( new KokkosAxpy<Real,Device>( std::forward<Real>(axpy.get_alpha() )));
}


//--------------------------------------------------------------------------------
// Aypx 

template<typename Real, typename Device>
void KokkosAypx<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                     const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Real alpha = alpha_;
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { yd(i) = alpha*yd(i) + xd(i); } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Aypx<Real>& aypx ) {
  bfun_ = std::unique_ptr<KokkosAypx<Real,Device>>( new KokkosAypx<Real,Device>( std::forward<Real>(aypx.get_alpha() )));
}

//--------------------------------------------------------------------------------
// Divide

template<typename Real, typename Device>
void KokkosDivide<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                       const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { yd(i) /= xd(i); } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Divide<Real>& ) {
  bfun_ = std::unique_ptr<KokkosDivide<Real,Device>>( new KokkosDivide<Real,Device>() );
}

//--------------------------------------------------------------------------------
// DivideAndInvert

template<typename Real, typename Device>
void KokkosDivideAndInvert<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                                const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { yd(i) = xd(i) / yd(i); } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::DivideAndInvert<Real>& ) {
  bfun_ = std::unique_ptr<KokkosDivideAndInvert<Real,Device>>( new KokkosDivideAndInvert<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Greater 

template<typename Real, typename Device>
void KokkosGreater<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                        const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    yd(i) = static_cast<Real>(yd(i) > xd(i));
  } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Greater<Real>& ) {
  bfun_ = std::unique_ptr<KokkosGreater<Real,Device>>( new KokkosGreater<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Lesser 

template<typename Real, typename Device>
void KokkosLesser<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                       const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    yd(i) = static_cast<Real>(yd(i) < xd(i));
  } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Lesser<Real>& ) {
  bfun_ = std::unique_ptr<KokkosLesser<Real,Device>>( new KokkosLesser<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Max

template<typename Real, typename Device>
void KokkosMax<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                   const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    yd(i) = yd(i) > xd(i) ? yd(i) : xd(i);
  } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Max<Real>& ) {
  bfun_ = std::unique_ptr<KokkosMax<Real,Device>>( new KokkosMax<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Min

template<typename Real, typename Device>
void KokkosMin<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                   const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    yd(i) = yd(i) < xd(i) ? yd(i) : xd(i);
  } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Min<Real>& ) {
  bfun_ = std::unique_ptr<KokkosMin<Real,Device>>( new KokkosMin<Real,Device>() );
}


//--------------------------------------------------------------------------------
// Multiply

template<typename Real, typename Device>
void KokkosMultiply<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                    const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    yd(i) = yd(i) * xd(i);
  } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Multiply<Real>& ) {
  bfun_ = std::unique_ptr<KokkosMultiply<Real,Device>>( new KokkosMultiply<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Plus

template<typename Real, typename Device>
void KokkosPlus<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                     const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    yd(i) = yd(i) + xd(i);
  } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Plus<Real>& ) {
  bfun_ = std::unique_ptr<KokkosPlus<Real,Device>>( new KokkosPlus<Real,Device>() );
}

//--------------------------------------------------------------------------------
// Set

template<typename Real, typename Device>
void KokkosSet<Real,Device>::apply(       ::ROL::KokkosVector<Real,Device>& y,
                                    const ::ROL::KokkosVector<Real,Device>& x ) const {
  auto yd = y.view_device();
  const auto xd = x.view_device();
  Kokkos::parallel_for( x.dimension(), KOKKOS_LAMBDA( const int i ) { 
    yd(i) = xd(i);
  } );
}

template<typename Real, typename Device>
void KokkosBinaryFunction<Real,Device>::Factory::visit( const ::ROL::Elementwise::Set<Real>& ) {
  bfun_ = std::unique_ptr<KokkosSet<Real,Device>>( new KokkosSet<Real,Device>() );
}

} // namespace Elementwise
} // namespace ROL

#endif  // ROL_KOKKOS_BINARYFUNCTIONS_DEF_HPP

// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_UQ_PCE_FWD_HPP
#define KOKKOS_VIEW_UQ_PCE_FWD_HPP

// Only include forward declarations so any overloads appear before they
// might be used inside Kokkos
#include "Kokkos_Core_fwd.hpp"
#include "Kokkos_View.hpp"

//----------------------------------------------------------------------------

namespace Sacado {
  namespace UQ {
    template <typename Storage >
    class PCE;
  }
}

namespace Kokkos {

  namespace Impl {
    template<class Space, class T, class ... P>
    struct MirrorType;
  }

}

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct ViewPCEContiguous;

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

namespace Kokkos {

// Declare overloads of create_mirror() so they are in scope
// Kokkos_Core.hpp is included below

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror(
  const Kokkos::View<T,P...> & src,
  typename std::enable_if<
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Experimental::Impl::ViewPCEContiguous >::value &&
    !std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
        Kokkos::LayoutStride >::value >::type * = 0);

template< class T , class ... P >
inline
typename Kokkos::View<T,P...>::HostMirror
create_mirror(
  const Kokkos::View<T,P...> & src,
  typename std::enable_if<
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Experimental::Impl::ViewPCEContiguous >::value &&
    std::is_same< typename Kokkos::ViewTraits<T,P...>::array_layout,
      Kokkos::LayoutStride >::value >::type * = 0);

template<class Space, class T, class ... P>
typename Impl::MirrorType<Space,T,P ...>::view_type
create_mirror(
  const Space&,
  const Kokkos::View<T,P...> & src,
  typename std::enable_if<
    std::is_same< typename ViewTraits<T,P...>::specialize ,
      Kokkos::Experimental::Impl::ViewPCEContiguous >::value >::type * = 0);

// Overload of deep_copy for UQ::PCE views intializing to a constant scalar
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
                 )>::type * = 0 );

// Overload of deep_copy for UQ::PCE views intializing to a constant UQ::PCE
template< class DT, class ... DP >
void deep_copy(
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
                 )>::type * = 0 );

// Overload of deep_copy for UQ::PCE views intializing to a constant scalar
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::array_type::value_type & value
  , typename std::enable_if<(
  Kokkos::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
                 )>::type * = 0 );

// Overload of deep_copy for UQ::PCE views intializing to a constant UQ::PCE
template< class ExecSpace , class DT, class ... DP >
void deep_copy(
  const ExecSpace &,
  const View<DT,DP...> & view ,
  const typename View<DT,DP...>::value_type & value
  , typename std::enable_if<(
  Kokkos::is_execution_space< ExecSpace >::value &&
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
                 )>::type * = 0 );

/* Specialize for deep copy of UQ::PCE */
template< class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  &&
  std::is_same< typename ViewTraits<ST,SP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
                  )>::type * = 0 );

/* Specialize for deep copy of UQ::PCE */
template< class ExecSpace, class DT , class ... DP , class ST , class ... SP >
inline
void deep_copy( const ExecSpace &,
                const View<DT,DP...> & dst ,
                const View<ST,SP...> & src
  , typename std::enable_if<(
  std::is_same< typename ViewTraits<DT,DP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
  &&
  std::is_same< typename ViewTraits<ST,SP...>::specialize
              , Kokkos::Experimental::Impl::ViewPCEContiguous >::value
                  )>::type * = 0 );

} // namespace Kokkos

#endif /* #ifndef KOKKOS_VIEW_UQ_PCE_FWD_HPP */

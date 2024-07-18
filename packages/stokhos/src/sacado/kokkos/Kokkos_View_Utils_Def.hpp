// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_UTILS_DEF_HPP
#define KOKKOS_VIEW_UTILS_DEF_HPP

#include "Kokkos_View_Utils.hpp"
//#include "Kokkos_ExecPolicy.hpp"
//#include "Kokkos_Parallel.hpp"

namespace Kokkos {

namespace Impl {

// Specialization for deep_copy( view, view::value_type ) for Cuda

template <class OutputView, typename Enabled>
struct StokhosViewFill
{
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::execution_space execution_space ;

  const OutputView output ;
  const_value_type input ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_t i0 ) const
  {
    const size_t n1 = output.extent(1);
    const size_t n2 = output.extent(2);
    const size_t n3 = output.extent(3);
    const size_t n4 = output.extent(4);
    const size_t n5 = output.extent(5);
    const size_t n6 = output.extent(6);
    const size_t n7 = output.extent(7);

    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_t i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_t i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_t i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_t i7 = 0 ; i7 < n7 ; ++i7 ) {
      output.access(i0,i1,i2,i3,i4,i5,i6,i7) = input ;
    }}}}}}}
  }

  StokhosViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      const size_t n0 = output.extent(0);
      Kokkos::RangePolicy<execution_space> policy( 0, n0 );
      Kokkos::parallel_for( policy, *this );
      execution_space().fence();
    }
};
} // namespace Impl

} // namespace Kokkos

#endif // KOKKOS_VIEW_UTILS_DEF_HPP

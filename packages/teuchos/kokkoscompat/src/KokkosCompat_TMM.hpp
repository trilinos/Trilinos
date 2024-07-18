// @HEADER
// *****************************************************************************
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//
// Copyright 2012 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOSCOMPAT_TMM_HPP
#define KOKKOSCOMPAT_TMM_HPP

/// \file KokkosCompat_TMM.hpp
/// \brief Compatibility between Kokkos Array View and Teuchos memory management classes.
///
/// \warning This file and everything in it should be considered an
///   implementation detail of Tpetra.  We make no promises of
///   backwards compatibility for any interface in this file, nor do
///   we even promise that this header file will continue to exist.

#include <TeuchosKokkosCompat_config.h>
#include <Kokkos_Core.hpp>
#include <Teuchos_ArrayView.hpp>

#if 0

namespace Kokkos {
namespace Impl {

  // template<class DestViewType>
  // struct ViewAssignment<DestViewType, Teuchos::ArrayView<typename DestViewType::value_type> {}; // ???

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  template<class ValueType>
  struct ViewAssignment<View<ValueType, Threads>, Teuchos::ArrayView<const ValueType> > {
    ViewAssignment (View<ValueType, Threads>& dst,
                    const Teuchos::ArrayView<const ValueType>& src,
                    typename enable_if< (
                      is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                               typename ViewTraits<ST,SL,SD,SM>::array_layout >::value
                      &&
                      ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) == unsigned(ViewTraits<ST,SL,SD,SM>::rank) )
                      &&
                      ( unsigned(ViewTraits<DT,DL,DD,DM>::rank_dynamic) >= unsigned(ViewTraits<ST,SL,SD,SM>::rank_dynamic) )
                    ) >::type * = 0 )



      View<ValueType, LayoutRight, Threads, MemoryUnmanaged> srcView; // (src.getRawPtr (), src.size ());
      srcView.m_shape.N0 = src.size ();


    }

    static void deep_copy (const DestViewType& dst, const Teuchos::ArrayView<const ValueType>& src);
  };


} // namespace Impl
} // namespace Kokkos

namespace Kokkos {
  namespace Compat {

    template<class ViewType>
    inline void
    deep_copy (const ViewType& dst, const Teuchos::ArrayView<DT>& src);

    template<class ValueType>
    inline void
    deep_copy (const Kokkos::View<ValueType, Kokkos::Threads>& dst,
               const Teuchos::ArrayView<ValueType>& src)
    {
      using Kokkos::LayoutRight;
      using Kokkos::Threads;
      using Kokkos::MemoryUnmanaged;
      using Kokkos::View;

      View<ValueType, LayoutRight, Threads, MemoryUnmanaged> srcView (src.getRawPtr (), src.size ());
      Kokkos::deep_copy (dst, srcView);
    }

  } // namespace Compat
} // namespace Kokkos

#endif // 0

#endif // KOKKOSCOMPAT_TMM_HPP

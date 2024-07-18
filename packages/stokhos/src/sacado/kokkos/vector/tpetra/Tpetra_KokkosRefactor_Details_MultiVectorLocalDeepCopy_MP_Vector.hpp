// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_MP_VECTOR_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_MP_VECTOR_HPP

#include "Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"

namespace Tpetra {
namespace Details {

  template<class DT, class ... DP,
           class ST, class ... SP,
           class DstWhichVecsType,
           class SrcWhichVecsType>
  typename std::enable_if<
    Kokkos::is_view_mp_vector< Kokkos::View<DT,DP...> >::value &&
    Kokkos::is_view_mp_vector< Kokkos::View<ST,SP...> >::value >::type
  localDeepCopy (
    const Kokkos::View<DT,DP...>& dst,
    const Kokkos::View<ST,SP...>& src,
    const bool dstConstStride, const bool srcConstStride,
    const DstWhichVecsType& dstWhichVecs,
    const SrcWhichVecsType& srcWhichVecs)
  {
    typedef Kokkos::View<DT,DP...> DstViewType;
    typedef Kokkos::View<ST,SP...> SrcViewType;
    typename Kokkos::FlatArrayType<DstViewType>::type dst_flat = dst;
    typename Kokkos::FlatArrayType<SrcViewType>::type src_flat = src;
    localDeepCopy( dst_flat, src_flat, dstConstStride, srcConstStride,
                   dstWhichVecs, srcWhichVecs );
  }

  template<class DT, class ... DP,
           class ST, class ... SP>
  typename std::enable_if<
    Kokkos::is_view_mp_vector< Kokkos::View<DT,DP...> >::value &&
    Kokkos::is_view_mp_vector< Kokkos::View<ST,SP...> >::value >::type
  localDeepCopyConstStride (
    const Kokkos::View<DT,DP...>& dst,
    const Kokkos::View<ST,SP...>& src)
  {
    typedef Kokkos::View<DT,DP...> DstViewType;
    typedef Kokkos::View<ST,SP...> SrcViewType;
    typename Kokkos::FlatArrayType<DstViewType>::type dst_flat = dst;
    typename Kokkos::FlatArrayType<SrcViewType>::type src_flat = src;
    localDeepCopyConstStride( dst_flat, src_flat );
  }

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_MP_VECTOR_HPP

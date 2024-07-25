// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_TEUCHOS_COMM_ADAPTERS_UQ_PCE_HPP
#define KOKKOS_TEUCHOS_COMM_ADAPTERS_UQ_PCE_HPP

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_TEUCHOSKOKKOSCOMM)

#include "Sacado_UQ_PCE.hpp"
#include "Kokkos_View_UQ_PCE.hpp"
#include "Kokkos_TeuchosCommAdapters.hpp"

//----------------------------------------------------------------------------
// Overloads of Teuchos Comm View functions for Sacado::UQ::PCE scalar type
//----------------------------------------------------------------------------

namespace Teuchos {

//! Variant of send() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename D, typename ... P>
typename std::enable_if<Kokkos::is_view_uq_pce< Kokkos::View<D,P...> >::value>::type
send (const Kokkos::View<D,P...>& sendBuffer,
      const Ordinal count,
      const int destRank,
      const int tag,
      const Comm<Ordinal>& comm)
{
  typedef Kokkos::View<D,P...> view_type;
  typedef typename Kokkos::FlatArrayType<view_type>::type flat_array_type;

  flat_array_type array = sendBuffer;
  Ordinal array_count = count * Kokkos::dimension_scalar(sendBuffer);
  send(array, array_count, destRank, tag, comm);
}

//! Variant of ssend() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename D, typename ... P>
typename std::enable_if<Kokkos::is_view_uq_pce< Kokkos::View<D,P...> >::value>::type
ssend (const Kokkos::View<D,P...>& sendBuffer,
       const Ordinal count,
       const int destRank,
       const int tag,
       const Comm<Ordinal>& comm)
{
  typedef Kokkos::View<D,P...> view_type;
  typedef typename Kokkos::FlatArrayType<view_type>::type flat_array_type;

  flat_array_type array = sendBuffer;
  Ordinal array_count = count * Kokkos::dimension_scalar(sendBuffer);
  ssend(array, array_count, destRank, tag, comm);
}

//! Variant of readySend() that accepts a message tag.
template<typename Ordinal, typename D, typename ... P>
typename std::enable_if<Kokkos::is_view_uq_pce< Kokkos::View<D,P...> >::value>::type
readySend (const Kokkos::View<D,P...>& sendBuffer,
           const Ordinal count,
           const int destRank,
           const int tag,
           const Comm<Ordinal>& comm)
{
  typedef Kokkos::View<D,P...> view_type;
  typedef typename Kokkos::FlatArrayType<view_type>::type flat_array_type;

  flat_array_type array = sendBuffer;
  Ordinal array_count = count * Kokkos::dimension_scalar(sendBuffer);
  readySend(array, array_count, destRank, tag, comm);
}

//! Variant of isend() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename D, typename ... P>
typename std::enable_if<Kokkos::is_view_uq_pce< Kokkos::View<D,P...> >::value, RCP<CommRequest<Ordinal> > >::type
isend (const Kokkos::View<D,P...>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<Ordinal>& comm)
{
  typedef Kokkos::View<D,P...> view_type;
  typedef typename Kokkos::FlatArrayType<view_type>::type flat_array_type;

  flat_array_type array = sendBuffer;
  return isend(array, destRank, tag, comm);
}

//! Variant of ireceive that takes a tag argument (and restores the correct order of arguments).
template<typename Ordinal, typename D, typename ... P>
typename std::enable_if<Kokkos::is_view_uq_pce< Kokkos::View<D,P...> >::value, RCP<CommRequest<Ordinal> > >::type
ireceive (const Kokkos::View<D,P...>& recvBuffer,
          const int sourceRank,
          const int tag,
          const Comm<Ordinal>& comm)
{
  typedef Kokkos::View<D,P...> view_type;
  typedef typename Kokkos::FlatArrayType<view_type>::type flat_array_type;

  flat_array_type array = recvBuffer;
  return ireceive(array, sourceRank, tag, comm);
}

}

#endif

#endif /* #ifndef KOKKOS_TEUCHOS_COMM_ADAPTERS_UQ_PCE_HPP */

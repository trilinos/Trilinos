//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_TEUCHOS_COMM_ADAPTERS_HPP
#define KOKKOS_TEUCHOS_COMM_ADAPTERS_HPP

#include <TeuchosKokkosComm_config.h>
#include <Teuchos_CommHelpers.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosCompat_View.hpp>

/// \file Kokkos_TeuchosCommAdapters.hpp
/// \brief Adapters for Teuchos::Comm functions for Kokkos:View.
///
/// Currently these just pass the Kokkos::View pointer along to the existing
/// Teuchos::Comm functions, but in the future could be more adapted
/// implementations.
///
/// Not all comm functions have been overloaded.

namespace Teuchos {

//! Variant of send() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, class ViewType>
typename std::enable_if<(Kokkos::is_view<ViewType>::value)>::type
send (const ViewType& sendBuffer,
      const Ordinal count,
      const int destRank,
      const int tag,
      const Comm<Ordinal>& comm)
{
  send(sendBuffer.data(), count, destRank, tag, comm);
}

//! Variant of ssend() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, class ViewType>
typename std::enable_if<(Kokkos::is_view<ViewType>::value)>::type
ssend (const ViewType& sendBuffer,
       const Ordinal count,
       const int destRank,
       const int tag,
       const Comm<Ordinal>& comm)
{
  ssend(sendBuffer.data(), count, destRank, tag, comm);
}

//! Variant of readySend() that accepts a message tag.
template<typename Ordinal, class ViewType>
typename std::enable_if<(Kokkos::is_view<ViewType>::value)>::type
readySend (const ViewType& sendBuffer,
           const Ordinal count,
           const int destRank,
           const int tag,
           const Comm<Ordinal>& comm)
{
  readySend(sendBuffer.data(), count, destRank, tag, comm);
}

//! Variant of isend() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, class ViewType>
typename std::enable_if<(Kokkos::is_view<ViewType>::value),RCP<CommRequest<Ordinal> >>::type
isend (const ViewType& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<Ordinal>& comm)
{
  using Kokkos::Compat::persistingView;
  // See Issue #1454: https://github.com/trilinos/Trilinos/issues/1454
  typename ViewType::const_type sendBuffer_const = sendBuffer;
  return isend (persistingView (sendBuffer_const), destRank, tag, comm);
}

//! Variant of ireceive that takes a tag argument (and restores the correct order of arguments).
template<typename Ordinal, class ViewType>
typename std::enable_if<(Kokkos::is_view<ViewType>::value),RCP<CommRequest<Ordinal> >>::type
ireceive (const ViewType& recvBuffer,
          const int sourceRank,
          const int tag,
          const Comm<Ordinal>& comm)
{
  using Kokkos::Compat::persistingView;
  return ireceive(persistingView(recvBuffer), sourceRank, tag, comm);
}


template<typename Ordinal, typename SendViewType, typename RecvViewType>
typename std::enable_if<(Kokkos::is_view<SendViewType>::value && Kokkos::is_view<RecvViewType>::value)>::type
reduceAll (const SendViewType& sendBuf,
           const RecvViewType& recvBuf,
           const EReductionType reductionType,
           const Comm<Ordinal>& comm)
{
  // We can't use the array of intrinsic scalar type
  // ((non_)const_array_intrinsic_type) here, because we're doing a
  // reduction.  That means we need to compute with the actual value
  // type.
  typedef typename SendViewType::value_type send_value_type;
  typedef typename RecvViewType::value_type recv_value_type;

  const bool typesDiffer =
    ! std::is_same<send_value_type, recv_value_type>::value;
  TEUCHOS_TEST_FOR_EXCEPTION(
    typesDiffer, std::invalid_argument, "Teuchos::reduceAll: Send and receive "
    "Views contain data of different types.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    SendViewType::rank > 1 || RecvViewType::rank > 1, std::invalid_argument,
    "Teuchos::reduceAll: Both send and receive Views must have rank 1.  "
    "The send View's rank is " << SendViewType::rank << " and the receive "
    "View's rank is " << RecvViewType::rank << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    sendBuf.extent (0) != recvBuf.extent (0), std::invalid_argument,
    "Send and receive buffer lengths do not match.  sendBuf.extent(0) = "
    << sendBuf.extent (0) << " != recvBuf.extent(0) = "
    << recvBuf.extent (0) << ".");

  // mfh 04 Nov 2014: Don't let Teuchos::SerialComm do a deep copy;
  // that always happens on the host, since SerialComm doesn't know
  // about Kokkos.
  if (comm.getSize () == 1) {
    Kokkos::deep_copy (recvBuf, sendBuf);
  }
  else {
    const Ordinal count = static_cast<Ordinal> (sendBuf.extent (0));
    reduceAll (comm, reductionType, count,
               sendBuf.data (),
               recvBuf.data ());
  }
}

template<typename Ordinal, typename Serializer,
         class SendViewType,
         class RecvViewType>
typename std::enable_if<(Kokkos::is_view<SendViewType>::value && Kokkos::is_view<RecvViewType>::value)>::type
reduceAll(const Comm<Ordinal>& comm,
          const Serializer& serializer,
          const EReductionType reductType,
          const Ordinal count,
          const SendViewType& sendBuffer,
          const RecvViewType& recvBuffer)
{
  // We can't use the array of intrinsic scalar type
  // ((non_)const_array_intrinsic_type) here, because we're doing a
  // reduction.  That means we need to compute with the actual value
  // type.
  typedef typename SendViewType::value_type send_value_type;
  typedef typename RecvViewType::value_type recv_value_type;

  const bool typesDiffer =
    ! std::is_same<send_value_type, recv_value_type>::value;
  TEUCHOS_TEST_FOR_EXCEPTION(
    typesDiffer, std::invalid_argument, "Teuchos::reduceAll: Send and receive "
    "Views contain data of different types.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    SendViewType::rank > 1 || RecvViewType::rank > 1, std::invalid_argument,
    "Teuchos::reduceAll: Both send and receive Views must have rank 1.  "
    "The send View's rank is " << SendViewType::rank << " and the receive "
    "View's rank is " << RecvViewType::rank << ".");

  // mfh 04 Nov 2014: Don't let Teuchos::SerialComm do a deep copy;
  // that always happens on the host, since SerialComm doesn't know
  // about Kokkos.
  if (comm.getSize () == 1) {
    Kokkos::deep_copy (recvBuffer, sendBuffer);
  }
  else {
    reduceAll (comm, serializer, reductType, count,
               sendBuffer.data (),
               recvBuffer.data ());
  }
}

template<typename Ordinal,
         class ViewType>
typename std::enable_if<(Kokkos::is_view<ViewType>::value)>::type
broadcast(const Comm<Ordinal>& comm,
               const int rootRank,
               const Ordinal count,
               const ViewType& buffer)
{
  broadcast( comm, rootRank, count, buffer.data() );
}

template<typename Ordinal,
         class ViewType,
         typename Serializer>
typename std::enable_if<(Kokkos::is_view<ViewType>::value)>::type
broadcast(const Comm<Ordinal>& comm,
               const Serializer& serializer,
               const int rootRank,
               const Ordinal count,
               const ViewType& buffer)
{
  broadcast( comm, serializer, rootRank, count, buffer.data() );
}

} // namespace Teuchos

#endif // KOKKOS_TEUCHOS_COMM_ADAPTERS_HPP

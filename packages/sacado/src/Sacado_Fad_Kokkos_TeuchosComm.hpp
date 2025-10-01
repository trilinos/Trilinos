// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_VIEW_SUPPORT_INCLUDES
#error "This file can only be included by Sacado_Fad_Kokkos_View_Support.hpp"
#endif

#include "Kokkos_TeuchosCommAdapters.hpp"

// =====================================================================
// This file includes overloads for Teuchos communication functions for
// Kokkos::View s of FadTypes
// -- reduceAll
// -- broadcast

namespace Teuchos {

template <typename Ordinal, class SD, class... SP, class RD, class... RP>
typename std::enable_if<
    Kokkos::is_view_fad<Kokkos::View<SD, SP...>>::value &&
    Kokkos::is_view_fad<Kokkos::View<RD, RP...>>::value>::type
reduceAll(const Comm<Ordinal> &comm, const EReductionType reductType,
          const Ordinal count, const Kokkos::View<SD, SP...> &sendBuffer,
          const Kokkos::View<RD, RP...> &recvBuffer) {
  // We can't implement reduceAll by extracting the underlying array (since we
  // can't reduce across the derivative dimension) and we can't just extract
  // a pointer due to ViewFad.  In principle we could handle ViewFad in the
  // serializer, but for the time being we just copy the view's into local
  // buffers (on the host).
  typedef Kokkos::View<SD, SP...> SendViewType;
  typedef Kokkos::View<RD, RP...> RecvViewType;
  typedef typename SendViewType::value_type send_value_type;
  typedef typename RecvViewType::value_type recv_value_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
      SendViewType::rank > 1 || RecvViewType::rank > 1, std::invalid_argument,
      "Teuchos::reduceAll: Both send and receive Views must have rank 1.  "
      "The send View's rank is "
          << SendViewType::rank
          << " and the receive "
             "View's rank is "
          << RecvViewType::rank << ".");

  // Copy send buffer into local array
  Teuchos::Array<send_value_type> localSendBuffer(count);
  typename SendViewType::host_mirror_type hostSendBuffer =
      Kokkos::create_mirror_view(sendBuffer);
  Kokkos::deep_copy(hostSendBuffer, sendBuffer);
  for (Ordinal i = 0; i < count; ++i)
    localSendBuffer[i] = hostSendBuffer(i);

  // Copy receive buffer into local array (necessary to initialize Fad types
  // properly)
  Teuchos::Array<recv_value_type> localRecvBuffer(count);
  typename RecvViewType::host_mirror_type hostRecvBuffer =
      Kokkos::create_mirror_view(recvBuffer);
  Kokkos::deep_copy(hostRecvBuffer, recvBuffer);
  for (Ordinal i = 0; i < count; ++i)
    localRecvBuffer[i] = hostRecvBuffer(i);

  // Do reduce-all
  reduceAll(comm, reductType, count, localSendBuffer.getRawPtr(),
            localRecvBuffer.getRawPtr());

  // Copy back into original buffer
  for (Ordinal i = 0; i < count; ++i)
    hostRecvBuffer(i) = localRecvBuffer[i];
  Kokkos::deep_copy(recvBuffer, hostRecvBuffer);
}

template <typename Ordinal, typename Serializer, class SD, class... SP,
          class RD, class... RP>
typename std::enable_if<
    Kokkos::is_view_fad<Kokkos::View<SD, SP...>>::value &&
    Kokkos::is_view_fad<Kokkos::View<RD, RP...>>::value>::type
reduceAll(const Comm<Ordinal> &comm, const Serializer &serializer,
          const EReductionType reductType, const Ordinal count,
          const Kokkos::View<SD, SP...> &sendBuffer,
          const Kokkos::View<RD, RP...> &recvBuffer) {
  // We can't implement reduceAll by extracting the underlying array (since we
  // can't reduce across the derivative dimension) and we can't just extract
  // a pointer due to ViewFad.  In principle we could handle ViewFad in the
  // serializer, but for the time being we just copy the view's into local
  // buffers (on the host).
  typedef Kokkos::View<SD, SP...> SendViewType;
  typedef Kokkos::View<RD, RP...> RecvViewType;
  typedef typename SendViewType::value_type send_value_type;
  typedef typename RecvViewType::value_type recv_value_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
      SendViewType::rank > 1 || RecvViewType::rank > 1, std::invalid_argument,
      "Teuchos::reduceAll: Both send and receive Views must have rank 1.  "
      "The send View's rank is "
          << SendViewType::rank
          << " and the receive "
             "View's rank is "
          << RecvViewType::rank << ".");

  // Copy send buffer into local array
  Teuchos::Array<send_value_type> localSendBuffer(count);
  typename SendViewType::host_mirror_type hostSendBuffer =
      Kokkos::create_mirror_view(sendBuffer);
  Kokkos::deep_copy(hostSendBuffer, sendBuffer);
  for (Ordinal i = 0; i < count; ++i)
    localSendBuffer[i] = hostSendBuffer(i);

  // Copy receive buffer into local array (necessary to initialize Fad types
  // properly)
  Teuchos::Array<recv_value_type> localRecvBuffer(count);
  typename RecvViewType::host_mirror_type hostRecvBuffer =
      Kokkos::create_mirror_view(recvBuffer);
  Kokkos::deep_copy(hostRecvBuffer, recvBuffer);
  for (Ordinal i = 0; i < count; ++i)
    localRecvBuffer[i] = hostRecvBuffer(i);

  // Do reduce-all
  reduceAll(comm, serializer, reductType, count, localSendBuffer.getRawPtr(),
            localRecvBuffer.getRawPtr());

  // Copy back into original buffer
  for (Ordinal i = 0; i < count; ++i)
    hostRecvBuffer(i) = localRecvBuffer[i];
  Kokkos::deep_copy(recvBuffer, hostRecvBuffer);
}

template <typename Ordinal, class D, class... P>
typename std::enable_if<Kokkos::is_view_fad<Kokkos::View<D, P...>>::value>::type
broadcast(const Comm<Ordinal> &comm, const int rootRank, const Ordinal count,
          const Kokkos::View<D, P...> &buffer) {
  auto array_buffer = Sacado::as_scalar_view(buffer);
  Ordinal array_count = count * Sacado::dimension_scalar(buffer);
  broadcast(comm, rootRank, array_count, array_buffer);
}

template <typename Ordinal, typename Serializer, class D, class... P>
typename std::enable_if<Kokkos::is_view_fad<Kokkos::View<D, P...>>::value>::type
broadcast(const Comm<Ordinal> &comm, const Serializer &serializer,
          const int rootRank, const Ordinal count,
          const Kokkos::View<D, P...> &buffer) {
  auto array_buffer = Sacado::as_scalar_view(buffer);
  Ordinal array_count = count * Sacado::dimension_scalar(buffer);
  broadcast(comm, *(serializer.getValueSerializer()), rootRank, array_count,
            array_buffer);
}

} // namespace Teuchos

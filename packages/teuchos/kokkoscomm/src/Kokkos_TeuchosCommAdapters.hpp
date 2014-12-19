/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
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
template<typename Ordinal, typename T, typename L, typename D, typename M, typename S>
void
send (const Kokkos::View<T,L,D,M,S>& sendBuffer,
      const Ordinal count,
      const int destRank,
      const int tag,
      const Comm<Ordinal>& comm)
{
  send(sendBuffer.ptr_on_device(), count, destRank, tag, comm);
}

//! Variant of ssend() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename T, typename L, typename D, typename M, typename S>
void
ssend (const Kokkos::View<T,L,D,M,S>& sendBuffer,
       const Ordinal count,
       const int destRank,
       const int tag,
       const Comm<Ordinal>& comm)
{
  ssend(sendBuffer.ptr_on_device(), count, destRank, tag, comm);
}

//! Variant of readySend() that accepts a message tag.
template<typename Ordinal, typename T, typename L, typename D, typename M, typename S>
void
readySend (const Kokkos::View<T,L,D,M,S>& sendBuffer,
           const Ordinal count,
           const int destRank,
           const int tag,
           const Comm<Ordinal>& comm)
{
  readySend(sendBuffer.ptr_on_device(), count, destRank, tag, comm);
}

//! Variant of isend() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename T, typename L, typename D, typename M, typename S>
RCP<CommRequest<Ordinal> >
isend (const Kokkos::View<T,L,D,M,S>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<Ordinal>& comm)
{
  using Kokkos::Compat::persistingView;
  return isend(persistingView(sendBuffer), destRank, tag, comm);
}

//! Variant of ireceive that takes a tag argument (and restores the correct order of arguments).
template<typename Ordinal, typename T, typename L, typename D, typename M, typename S>
RCP<CommRequest<Ordinal> >
ireceive (const Kokkos::View<T,L,D,M,S>& recvBuffer,
          const int sourceRank,
          const int tag,
          const Comm<Ordinal>& comm)
{
  using Kokkos::Compat::persistingView;
  return ireceive(persistingView(recvBuffer), sourceRank, tag, comm);
}


template<typename Ordinal, typename SendViewType, typename RecvViewType>
void
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
  typedef send_value_type val_type;

  const bool typesDiffer =
    ! Kokkos::Impl::is_same<send_value_type, recv_value_type>::value;
  TEUCHOS_TEST_FOR_EXCEPTION(
    typesDiffer, std::invalid_argument, "Teuchos::reduceAll: Send and receive "
    "Views contain data of different types.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    SendViewType::rank > 1 || RecvViewType::rank > 1, std::invalid_argument,
    "Teuchos::reduceAll: Both send and receive Views must have rank 1.  "
    "The send View's rank is " << SendViewType::rank << " and the receive "
    "View's rank is " << RecvViewType::rank << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    sendBuf.dimension_0 () != recvBuf.dimension_0 (), std::invalid_argument,
    "Send and receive buffer lengths do not match.  sendBuf.dimension_0() = "
    << sendBuf.dimension_0 () << " != recvBuf.dimension_0() = "
    << recvBuf.dimension_0 () << ".");

  // mfh 04 Nov 2014: Don't let Teuchos::SerialComm do a deep copy;
  // that always happens on the host, since SerialComm doesn't know
  // about Kokkos.
  if (comm.getSize () == 1) {
    Kokkos::deep_copy (recvBuf, sendBuf);
  }
  else {
    const Ordinal count = static_cast<Ordinal> (sendBuf.dimension_0 ());
    reduceAll<Ordinal, val_type> (comm, reductionType, count,
                                  sendBuf.ptr_on_device (),
                                  recvBuf.ptr_on_device ());
  }
}

template<typename Ordinal,
         typename T, typename L, typename D, typename M, typename S>
void broadcast(const Comm<Ordinal>& comm,
               const int rootRank,
               const Ordinal count,
               const Kokkos::View<T,L,D,M,S>& buffer)
{
  broadcast( comm, rootRank, count, buffer.ptr_on_device() );
}

template<typename Ordinal,
         typename T, typename L, typename D, typename M, typename S,
         typename Serializer>
void broadcast(const Comm<Ordinal>& comm,
               const Serializer& serializer,
               const int rootRank,
               const Ordinal count,
               const Kokkos::View<T,L,D,M,S>& buffer)
{
  broadcast( comm, serializer, rootRank, count, buffer.ptr_on_device() );
}

} // namespace Teuchos

#endif // KOKKOS_TEUCHOS_COMM_ADAPTERS_HPP

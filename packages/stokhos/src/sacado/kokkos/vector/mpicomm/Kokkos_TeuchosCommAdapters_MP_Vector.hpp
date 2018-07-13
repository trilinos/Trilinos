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

#ifndef KOKKOS_TEUCHOS_COMM_ADAPTERS_MP_VECTOR_HPP
#define KOKKOS_TEUCHOS_COMM_ADAPTERS_MP_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"
#if defined(HAVE_STOKHOS_TEUCHOSKOKKOSCOMM)

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_TeuchosCommAdapters.hpp"

//----------------------------------------------------------------------------
// Overloads of Teuchos Comm View functions for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Teuchos {

//! Variant of send() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename D, typename ... P>
typename std::enable_if<Kokkos::is_view_mp_vector< Kokkos::View<D,P...> >::value>::type
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
typename std::enable_if<Kokkos::is_view_mp_vector< Kokkos::View<D,P...> >::value>::type
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
typename std::enable_if<Kokkos::is_view_mp_vector< Kokkos::View<D,P...> >::value>::type
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
typename std::enable_if<Kokkos::is_view_mp_vector< Kokkos::View<D,P...> >::value, RCP<CommRequest<Ordinal> > >::type
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
typename std::enable_if<Kokkos::is_view_mp_vector< Kokkos::View<D,P...> >::value, RCP<CommRequest<Ordinal> > >::type
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

#endif /* #ifndef KOKKOS_TEUCHOS_COMM_ADAPTERS_MP_VECTOR_HPP */

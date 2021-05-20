// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_DETAILS_MPI_HPP
#define TPETRA_DETAILS_MPI_HPP

/// \file Tpetra_Details_MPI.hpp
/// \brief Declaration of Tpetra::Details::MPI
///
/// \warning This file and its contents are implementation details of
///   Tpetra.  Users must not rely on them.
///
/// Tpetra::Details::MPI::alltoallv wraps MPI_Alltoallv.  That is the
/// only thing in this file upon which Tpetra developers should rely.
/// Tpetra developers should not rely on anything else in this file.
/// <i>Users</i> may not rely on <i>anything</i> in this file!
///

#include "TpetraCore_config.h"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#  include "Tpetra_Details_MpiTypeTraits.hpp"
#endif // HAVE_TPETRACORE_MPI
#include "Tpetra_Details_temporaryViewUtils.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"
#include <memory>
#include <stdexcept>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // forward declaration of Comm
  template<class OrdinalType> class Comm;
} // namespace Teuchos
#endif // NOT DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {
namespace Details {
namespace MPI {

#ifdef HAVE_TPETRACORE_MPI
std::string getMpiErrorString (const int errCode);
#endif

/// \brief all-to-all, for rank-1 Kokkos::View objects.
///
/// \tparam SendViewType Type of the send buffer
/// \tparam RecvViewType Type of the receive buffer
/// \tparam CountsViewType Type of the counts buffers
/// \tparam OffsetViewType Type of the offset buffers
///
/// \param sendbuf [in] Send buffer; must be a rank-1 Kokkos::View.
/// \param sendcounts [in] Send count buffer; must be a rank-1 Kokkos::View.
/// \param sdispls [in] Send displacement buffer; must be a rank-1 Kokkos::View.
/// \param recvbuf [in] Receive buffer; must be a rank-1 Kokkos::View.
/// \param recvcounts [in] Receive count buffer; must be a rank-1 Kokkos::View.
/// \param rdispls [in] Receive displacement buffer; must be a rank-1 Kokkos::View.
/// \param comm [in] Communicator.
template<class SendViewType, class RecvViewType, class CountViewType, class OffsetViewType>
void
alltoallv (const SendViewType& sendbuf,
           const CountViewType& sendcounts,
           const OffsetViewType& sdispls,
           const RecvViewType& recvbuf,
           const CountViewType& recvcounts,
           const OffsetViewType& rdispls,
           const ::Teuchos::Comm<int>& comm)
{
  static_assert (Kokkos::Impl::is_view<SendViewType>::value,
                 "SendDataViewType must be a Kokkos::View specialization.");
  static_assert (Kokkos::Impl::is_view<RecvViewType>::value,
                 "RecvDataViewType must be a Kokkos::View specialization.");
  static_assert (Kokkos::Impl::is_view<CountViewType>::value,
                 "CountViewType must be a Kokkos::View specialization.");
  static_assert (Kokkos::Impl::is_view<OffsetViewType>::value,
                 "OffsetViewType must be a Kokkos::View specialization.");
  static_assert (static_cast<int> (SendViewType::rank) == 1,
                 "SendViewType must have rank 1.");
  static_assert (static_cast<int> (RecvViewType::rank) == 1,
                 "RecvViewType must have rank 1.");
  static_assert (static_cast<int> (CountViewType::rank) == 1,
                 "CountViewType must have rank 1.");
  static_assert (static_cast<int> (OffsetViewType::rank) == 1,
                 "OffsetViewType must have rank 1.");
  typedef typename RecvViewType::non_const_value_type packet_type;
  static_assert (std::is_same<typename RecvViewType::value_type,
                   packet_type>::value,
                 "RecvViewType must be a nonconst Kokkos::View.");
  //Make sure layouts are contiguous (don't accept strided 1D view)
  static_assert (!std::is_same<typename SendViewType::array_layout, Kokkos::LayoutStride>::value,
      "Input/Output views must be contiguous (not LayoutStride)");
  static_assert (!std::is_same<typename RecvViewType::array_layout, Kokkos::LayoutStride>::value,
      "Input/Output views must be contiguous (not LayoutStride)");
  static_assert (!std::is_same<typename CountViewType::array_layout, Kokkos::LayoutStride>::value,
      "Input/Output views must be contiguous (not LayoutStride)");
  static_assert (!std::is_same<typename OffsetViewType::array_layout, Kokkos::LayoutStride>::value,
      "Input/Output views must be contiguous (not LayoutStride)");

#ifdef HAVE_TPETRACORE_MPI
  MPI_Comm rawComm = ::Tpetra::Details::extractMpiCommFromTeuchos (comm);
  const int err = MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_CHAR,
                                recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_CHAR,
                                rawComm);
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::runtime_error,
     "MPI_Alltoallv failed with error \""
     << getMpiErrorString (err));
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                             "Tpetra::Details::MPI::alltoallv not implemented.");
#endif // HAVE_TPETRACORE_MPI
}

} // namespace MPI
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_MPI_HPP

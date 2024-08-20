// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_IALLREDUCE_HPP
#define TPETRA_DETAILS_IALLREDUCE_HPP

/// \file Tpetra_Details_iallreduce.hpp
/// \brief Declaration of Tpetra::Details::iallreduce
///
/// \warning This file and its contents are implementation details of
///   Tpetra.  Users must not rely on them.
///
/// Tpetra::Details::iallreduce wraps MPI_Iallreduce.  That is the
/// only thing in this file upon which Tpetra developers should rely.
/// Tpetra developers should not rely on anything else in this file.
/// <i>Users</i> may not rely on <i>anything</i> in this file!
///
/// If you want to find the only thing in this file that you are
/// supposed to use, search for "SKIP DOWN TO HERE" (no quotes).
/// "You" only refers to Tpetra developers.  Users, this file is not
/// for you!

#include "TpetraCore_config.h"
#include "Teuchos_EReductionType.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#  include "Tpetra_Details_MpiTypeTraits.hpp"
#endif // HAVE_TPETRACORE_MPI
#include "Tpetra_Details_temporaryViewUtils.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_Core.hpp"
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <functional>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // forward declaration of Comm
  template<class OrdinalType> class Comm;
} // namespace Teuchos
#endif // NOT DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_MPI
std::string getMpiErrorString (const int errCode);
#endif

/// \brief Base class for the request (more or less a future)
///   representing a pending nonblocking MPI operation.
///
/// We recommend using <tt>auto</tt> to interact with this class.  For
/// example, use <tt>auto</tt> for the type of the return value of
/// iallreduce() (see below).
class CommRequest {
public:
  //! Destructor (virtual for memory safety of derived classes).
  virtual ~CommRequest () {}

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.  This operation must also be idempotent.
  virtual void wait () {}

  /// \brief Cancel the pending communication request.
  ///
  /// This operation must be idempotent.
  virtual void cancel () {}
};

// Don't rely on anything in this namespace.
namespace Impl {

//! Return an "empty" comm request (waiting on it does nothing).
std::shared_ptr<CommRequest>
emptyCommRequest ();

#ifdef HAVE_TPETRACORE_MPI
#if MPI_VERSION >= 3
template<typename InputViewType, typename OutputViewType, typename ResultViewType>
struct MpiRequest : public CommRequest
{
  MpiRequest(const InputViewType& send, const OutputViewType& recv, const ResultViewType& result, MPI_Request req_)
    : sendBuf(send), recvBuf(recv), resultBuf(result), req(req_)
  {}

  ~MpiRequest()
  {
    //this is a no-op if wait() or cancel() have already been called
    cancel();
  }

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.  This operation must also be idempotent.
  void wait () override
  {
    if (req != MPI_REQUEST_NULL) {
      const int err = MPI_Wait (&req, MPI_STATUS_IGNORE);
      TEUCHOS_TEST_FOR_EXCEPTION
        (err != MPI_SUCCESS, std::runtime_error,
         "MpiCommRequest::wait: MPI_Wait failed with error \""
         << getMpiErrorString (err));
      // MPI_Wait should set the MPI_Request to MPI_REQUEST_NULL on
      // success.  We'll do it here just to be conservative.
      req = MPI_REQUEST_NULL;
      //Since recvBuf contains the result, copy it to the user's resultBuf.
      Kokkos::deep_copy(resultBuf, recvBuf);
    }
  }

  /// \brief Cancel the pending communication request.
  ///
  /// This operation must be idempotent.
  void cancel () override
  {
    //BMK: Per https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node126.htm,
    //MPI_Cancel cannot be used for collectives like iallreduce.
    req = MPI_REQUEST_NULL;
  }

private:
  InputViewType sendBuf;
  OutputViewType recvBuf;
  ResultViewType resultBuf;
  //This request is active if and only if req != MPI_REQUEST_NULL.
  MPI_Request req;
};

/// \brief Low-level wrapper for MPI_Iallreduce (async).
///   May only be used if MPI version >= 3.
MPI_Request
iallreduceRaw (const void* sendbuf,
               void* recvbuf,
               const int count,
               MPI_Datatype mpiDatatype,
               const Teuchos::EReductionType op,
               MPI_Comm comm);
#endif

/// \brief Low-level wrapper for MPI_Allreduce (blocking)
void
allreduceRaw  (const void* sendbuf,
               void* recvbuf,
               const int count,
               MPI_Datatype mpiDatatype,
               const Teuchos::EReductionType op,
               MPI_Comm comm);

template<class InputViewType, class OutputViewType>
std::shared_ptr<CommRequest>
iallreduceImpl (const InputViewType& sendbuf,
            const OutputViewType& recvbuf,
            const ::Teuchos::EReductionType op,
            const ::Teuchos::Comm<int>& comm)
{
  using Packet = typename InputViewType::non_const_value_type;
  if(comm.getSize() == 1)
  {
    Kokkos::deep_copy(recvbuf, sendbuf);
    return emptyCommRequest();
  }
  Packet examplePacket;
  MPI_Datatype mpiDatatype = sendbuf.extent(0) ?
    MpiTypeTraits<Packet>::getType (examplePacket) :
    MPI_BYTE;
  bool datatypeNeedsFree = MpiTypeTraits<Packet>::needsFree;
  MPI_Comm rawComm = ::Tpetra::Details::extractMpiCommFromTeuchos (comm);
  //Note BMK: Nonblocking collectives like iallreduce cannot use GPU buffers.
  //See https://www.open-mpi.org/faq/?category=runcuda#mpi-cuda-support
  auto sendMPI = Tpetra::Details::TempView::toMPISafe<InputViewType, false>(sendbuf);
  auto recvMPI = Tpetra::Details::TempView::toMPISafe<OutputViewType, false>(recvbuf);
  std::shared_ptr<CommRequest> req;
  //Next, if input/output alias and comm is an intercomm, make a deep copy of input.
  //Not possible to do in-place allreduce for intercomm.
  if(isInterComm(comm) && sendMPI.data() == recvMPI.data())
  {
    //Can't do in-place collective on an intercomm,
    //so use a separate 1D copy as the input.
    Kokkos::View<Packet*, Kokkos::HostSpace> tempInput(Kokkos::ViewAllocateWithoutInitializing("tempInput"), sendMPI.extent(0));
    for(size_t i = 0; i < sendMPI.extent(0); i++)
      tempInput(i) = sendMPI.data()[i];
#if MPI_VERSION >= 3
    //MPI 3+: use async allreduce
    MPI_Request mpiReq = iallreduceRaw((const void*) tempInput.data(), (void*) recvMPI.data(), tempInput.extent(0), mpiDatatype, op, rawComm);
    req = std::shared_ptr<CommRequest>(new MpiRequest<decltype(tempInput), decltype(recvMPI), OutputViewType>(tempInput, recvMPI, recvbuf, mpiReq));
#else
    //Older MPI: Iallreduce not available. Instead do blocking all-reduce and return empty request.
    allreduceRaw((const void*) sendMPI.data(), (void*) recvMPI.data(), sendMPI.extent(0), mpiDatatype, op, rawComm);
    Kokkos::deep_copy(recvbuf, recvMPI);
    req = emptyCommRequest();
#endif
  }
  else
  {
#if MPI_VERSION >= 3
    //MPI 3+: use async allreduce
    MPI_Request mpiReq = iallreduceRaw((const void*) sendMPI.data(), (void*) recvMPI.data(), sendMPI.extent(0), mpiDatatype, op, rawComm);
    req = std::shared_ptr<CommRequest>(new MpiRequest<decltype(sendMPI), decltype(recvMPI), OutputViewType>(sendMPI, recvMPI, recvbuf, mpiReq));
#else
    //Older MPI: Iallreduce not available. Instead do blocking all-reduce and return empty request.
    allreduceRaw((const void*) sendMPI.data(), (void*) recvMPI.data(), sendMPI.extent(0), mpiDatatype, op, rawComm);
    Kokkos::deep_copy(recvbuf, recvMPI);
    req = emptyCommRequest();
#endif
  }
  if(datatypeNeedsFree)
    MPI_Type_free(&mpiDatatype);
  return req;
}

#else

//No MPI: reduction is always the same as input.
template<class InputViewType, class OutputViewType>
std::shared_ptr<CommRequest>
iallreduceImpl (const InputViewType& sendbuf,
            const OutputViewType& recvbuf,
            const ::Teuchos::EReductionType,
            const ::Teuchos::Comm<int>&)
{
  Kokkos::deep_copy(recvbuf, sendbuf);
  return emptyCommRequest();
}

#endif // HAVE_TPETRACORE_MPI

} // namespace Impl

//
// SKIP DOWN TO HERE
//

/// \brief Nonblocking all-reduce, for either rank-1 or rank-0
///   Kokkos::View objects.
///
/// \tparam InputViewType Type of the send buffer
/// \tparam OutputViewType Type of the receive buffer
///
/// This function wraps MPI_Iallreduce.  It does a nonblocking
/// all-reduce over the input communicator \c comm, from \c sendbuf
/// into \c recvbuf, using \c op as the reduction operator.  The
/// function returns without blocking; the all-reduce only blocks for
/// completion when one calls wait() on the returned request.
///
/// \param sendbuf [in] Input buffer; must be either a rank-1 or
///   rank-0 Kokkos::View, and must have the same rank as recvbuf.
/// \param recvbuf [in] Output buffer; must be either a rank-1 or
///   rank-0 Kokkos::View, and must have the same rank as sendbuf.
/// \param op [in] Teuchos enum representing the reduction operator.
/// \param comm [in] Communicator over which to do the all-reduce.
///
/// \c sendbuf and \c recvbuf must either be disjoint, or identical
/// (point to the same array).  They may not partially overlap.
/// Furthermore, if they are identical, the input communicator must be
/// an intracommunicator.  It may not be an intercommunicator.  If you
/// don't know what an intercommunicator is, you probably just have an
/// intracommunicator, so everything is fine.
template<class InputViewType, class OutputViewType>
std::shared_ptr<CommRequest>
iallreduce (const InputViewType& sendbuf,
            const OutputViewType& recvbuf,
            const ::Teuchos::EReductionType op,
            const ::Teuchos::Comm<int>& comm)
{
  static_assert (Kokkos::is_view<InputViewType>::value,
                 "InputViewType must be a Kokkos::View specialization.");
  static_assert (Kokkos::is_view<OutputViewType>::value,
                 "OutputViewType must be a Kokkos::View specialization.");
  constexpr int rank = static_cast<int> (OutputViewType::rank);
  static_assert (static_cast<int> (InputViewType::rank) == rank,
                 "InputViewType and OutputViewType must have the same rank.");
  static_assert (rank == 0 || rank == 1,
                 "InputViewType and OutputViewType must both have "
                 "rank 0 or rank 1.");
  typedef typename OutputViewType::non_const_value_type packet_type;
  static_assert (std::is_same<typename OutputViewType::value_type,
                   packet_type>::value,
                 "OutputViewType must be a nonconst Kokkos::View.");
  static_assert (std::is_same<typename InputViewType::non_const_value_type,
                   packet_type>::value,
                 "InputViewType and OutputViewType must be Views "
                 "whose entries have the same type.");
  //Make sure layouts are contiguous (don't accept strided 1D view)
  static_assert (!std::is_same<typename InputViewType::array_layout, Kokkos::LayoutStride>::value,
      "Input/Output views must be contiguous (not LayoutStride)");
  static_assert (!std::is_same<typename OutputViewType::array_layout, Kokkos::LayoutStride>::value,
      "Input/Output views must be contiguous (not LayoutStride)");

  return Impl::iallreduceImpl<InputViewType, OutputViewType> (sendbuf, recvbuf, op, comm);
}

std::shared_ptr<CommRequest>
iallreduce (const int localValue,
            int& globalValue,
            const ::Teuchos::EReductionType op,
            const ::Teuchos::Comm<int>& comm);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_IALLREDUCE_HPP

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
  virtual void wait () = 0;

  /// \brief Cancel the pending communication request.
  ///
  /// This operation must be idempotent.
  virtual void cancel () = 0;
};

// Don't rely on anything in this namespace.  You shouldn't be relying
// on anything in Tpetra::Details anyway.  You _especially_ shouldn't
// be relying on anything in Tpetra::Details::Impl.  Consider
// yourselves warned!
namespace Impl {

//! Return an "empty" comm request (waiting on it does nothing).
std::shared_ptr<CommRequest>
emptyCommRequest ();

#ifdef HAVE_TPETRACORE_MPI

/// \class MpiCommRequest
/// \brief MPI implementation of CommRequest.
///
/// This class wraps MPI_Request, which is MPI's reification of a
/// nonblocking communication operation.
///
/// Users would not normally create an instance of this class.  Calls
/// to nonblocking communication operations may return a pointer to a
/// CommRequest, which may be an instance of this class.
///
/// Users might wish to create an MpiCommRequest directly if they want
/// to encapsulate an MPI_Request returned by an external library or
/// by their own code.
class MpiCommRequest : public CommRequest {
public:
  //! Default constructor; creates an empty ("null") request.
  MpiCommRequest ();

  //! Constructor (from a raw MPI_Request).
  MpiCommRequest (const MPI_Request req);

  //! Destructor; cancels the request if it is still pending.
  virtual ~MpiCommRequest ();

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.  (For example, a receive must have a matching
  /// send, otherwise a wait on the receive will wait forever.)
  virtual void wait ();

  /// \brief Wait on this communication request to complete, and
  ///   return the resulting MPI_Status as an output argument.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.  (For example, a receive must have a matching
  /// send, otherwise a wait on the receive will wait forever.)
  virtual void waitWithStatus (MPI_Status& status);

  //! Cancel the pending communication request.
  virtual void cancel ();

private:
  MPI_Request req_;
};
#endif // HAVE_TPETRACORE_MPI

/// \brief Part of the Work-around for not having MPI >= 3.
///
/// Substitute for MPI_Iallreduce (which requires MPI 3): defer
/// calling MPI_Allreduce until wait.  It's ugly, but it works.
class DeferredActionCommRequest : public CommRequest {
public:
  /// \brief Default constructor (take no action on wait).
  DeferredActionCommRequest ();

  /// \brief Constructor that takes an action to defer.
  ///
  /// \brief action [in] An action to take on wait().
  ///
  /// We only take that action the first time you call wait().
  /// This means that wait() is idempotent.
  DeferredActionCommRequest (std::function<void () > action);

  //! Wait on this communication request to complete.
  void wait ();

  /// \brief Cancel the pending communication request, without taking
  ///   the specified action.
  ///
  /// Subsequent calls to wait() do nothing.
  void cancel ();

private:
  //! The deferred action to perform on wait().
  std::function<void(void) > action_;
  //! Whether the action has been taken.
  bool actionTaken_;
};

/// \brief Object representing a pending ::Tpetra::Details::iallreduce
///   operation.
///
/// This subclass keeps the send and receive buffers.  Since
/// ::Kokkos::View reference-counts, this ensures that the buffers
/// will not be deallocated until the iallreduce completes.  The
/// buffer references get cleared on wait().
///
/// Tpetra developers should not use this directly; they should
/// instead create instances of this via the wrapIallreduceCommRequest
/// function (see below).
///
/// \tparam PacketType Type of each entry of the send and receive
///   buffers.
/// \tparam SendLayoutType array_layout of the send buffer.  Must be
///   Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam SendDeviceType Kokkos::Device specialization used by the
///   send buffer.
/// \tparam RecvLayoutType array_layout of the receive buffer.  Must
///   be Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam RecvDeviceType Kokkos::Device specialization used by the
///   receive buffer.  It's OK for this to differ from SendDeviceType.
///   We assume that MPI implementations can handle this.  (This is a
///   reasonable assumption with CUDA-enabled MPI implementations.)
/// \tparam rank Integer rank of the send and receive buffers.  Must
///   be either 0 or 1.
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType,
         const int rank>
class IallreduceCommRequest;

//! Partial pecialization for rank-1 send and receive buffers.
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType>
class IallreduceCommRequest<PacketType,
                            SendLayoutType,
                            SendDeviceType,
                            RecvLayoutType,
                            RecvDeviceType,
                            1> :
    public CommRequest {
public:
  //! Default constructor.
  IallreduceCommRequest ()
  {}

  //! Type of the send buffer.
  typedef ::Kokkos::View<const PacketType*, SendLayoutType,
                         SendDeviceType> send_buffer_type;
  //! Type of the receive buffer.
  typedef ::Kokkos::View<PacketType*, RecvLayoutType,
                         RecvDeviceType> recv_buffer_type;

  /// \brief Constructor that takes a wrapped request (representing
  ///   the pending MPI_Iallreduce operation itself), and saved
  ///   buffers.
  IallreduceCommRequest (const std::shared_ptr<CommRequest>& req,
                         const send_buffer_type& sendbuf,
                         const recv_buffer_type& recvbuf) :
    req_ (req),
    sendbuf_ (sendbuf),
    recvbuf_ (recvbuf)
  {
    static_assert (std::is_same<SendLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<SendLayoutType, ::Kokkos::LayoutRight>::value,
                   "SendLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
    static_assert (std::is_same<RecvLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<RecvLayoutType, ::Kokkos::LayoutRight>::value,
                   "RecvLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
  }

  //! Destructor (virtual for memory safety).
  virtual ~IallreduceCommRequest () {
    if (req_.get () != NULL) {
      // We're in a destructor, so don't throw.  We'll just try our best
      // to handle whatever happens here without throwing.
      try {
        req_->cancel ();
      }
      catch (...) {}

      try {
        req_ = std::shared_ptr<CommRequest> ();
      }
      catch (...) {}
    }
  }

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.
  void
  wait ()
  {
    if (req_.get () != NULL) {
      req_->wait ();
      // Relinquish request handle.
      req_ = std::shared_ptr<CommRequest> ();
    }
    // Relinquish references to saved buffers, if we have not already
    // done so.  This operation is idempotent.
    sendbuf_ = send_buffer_type ();
    recvbuf_ = recv_buffer_type ();
  }

  /// \brief Cancel the pending communication request.
  ///
  /// This operation must be idempotent.
  void
  cancel ()
  {
    if (req_.get () != NULL) {
      req_->cancel ();
      // Relinquish request handle.
      req_ = std::shared_ptr<CommRequest> ();
    }
    // Relinquish references to saved buffers, if we have not already
    // done so.  This operation is idempotent.
    sendbuf_ = send_buffer_type ();
    recvbuf_ = recv_buffer_type ();
  }

private:
  //! The wrapped request
  std::shared_ptr<CommRequest> req_;

  //! Saved send buffer from iallreduce call
  send_buffer_type sendbuf_;

  //! Saved recv buffer from iallreduce call
  recv_buffer_type recvbuf_;
};

//! Partial pecialization for rank-0 (single-value) send and receive buffers.
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType>
class IallreduceCommRequest<PacketType,
                            SendLayoutType,
                            SendDeviceType,
                            RecvLayoutType,
                            RecvDeviceType,
                            0> :
    public CommRequest {
public:
  //! Default constructor.
  IallreduceCommRequest ()
  {}

  //! Type of the send buffer.
  typedef ::Kokkos::View<const PacketType, SendLayoutType,
                         SendDeviceType> send_buffer_type;
  //! Type of the receive buffer.
  typedef ::Kokkos::View<PacketType, RecvLayoutType,
                         RecvDeviceType> recv_buffer_type;

  /// \brief Constructor that takes a wrapped request (representing
  ///   the pending MPI_Iallreduce operation itself), and saved
  ///   buffers.
  IallreduceCommRequest (const std::shared_ptr<CommRequest>& req,
                         const send_buffer_type& sendbuf,
                         const recv_buffer_type& recvbuf) :
    req_ (req),
    sendbuf_ (sendbuf),
    recvbuf_ (recvbuf)
  {
    static_assert (std::is_same<SendLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<SendLayoutType, ::Kokkos::LayoutRight>::value,
                   "SendLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
    static_assert (std::is_same<RecvLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<RecvLayoutType, ::Kokkos::LayoutRight>::value,
                   "RecvLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
  }

  //! Destructor (virtual for memory safety).
  virtual ~IallreduceCommRequest () {
    if (req_.get () != NULL) {
      // We're in a destructor, so don't throw.  We'll just try our best
      // to handle whatever happens here without throwing.
      try {
        req_->cancel ();
      }
      catch (...) {}

      try {
        req_ = std::shared_ptr<CommRequest> ();
      }
      catch (...) {}
    }
  }

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.
  void
  wait ()
  {
    if (req_.get () != NULL) {
      req_->wait ();
      // Relinquish request handle.
      req_ = std::shared_ptr<CommRequest> ();
    }
    // Relinquish references to saved buffers, if we have not already
    // done so.  This operation is idempotent.
    sendbuf_ = send_buffer_type ();
    recvbuf_ = recv_buffer_type ();
  }

  /// \brief Cancel the pending communication request.
  ///
  /// This operation must be idempotent.
  void
  cancel ()
  {
    if (req_.get () != NULL) {
      req_->cancel ();
      // Relinquish request handle.
      req_ = std::shared_ptr<CommRequest> ();
    }
    // Relinquish references to saved buffers, if we have not already
    // done so.  This operation is idempotent.
    sendbuf_ = send_buffer_type ();
    recvbuf_ = recv_buffer_type ();
  }

private:
  //! The wrapped request
  std::shared_ptr<CommRequest> req_;

  //! Saved send buffer from iallreduce call
  send_buffer_type sendbuf_;

  //! Saved recv buffer from iallreduce call
  recv_buffer_type recvbuf_;
};

/// \brief Function for wrapping the CommRequest to be returned from
///   ::Tpetra::Details::iallreduce; overload for rank-1 send and
///   receive buffers.
///
/// The object returned from this function keeps the send and receive
/// buffers.  Since ::Kokkos::View reference-counts, this ensures that
/// the buffers will not be deallocated until the iallreduce
/// completes.  The buffer references get cleared on wait() or
/// cancel().
///
/// \param req [in] Pointer to CommRequest from
///   ::Tpetra::Details::iallreduce.
/// \param sendbuf [in] Send buffer for ::Tpetra::Details::iallreduce.
/// \param recvbuf [in] Receive buffer for
///   ::Tpetra::Details::iallreduce.  This is an input argument,
///   because this function itself does not modify the contents of the
///   receive buffer; it just keeps a reference to it.
///
/// \return Pointer to wrapped CommRequest.
///
/// \tparam PacketType Type of each entry of the send and receive
///   buffers.
/// \tparam SendLayoutType array_layout of the send buffer.  Must be
///   Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam SendDeviceType Kokkos::Device specialization used by the
///   send buffer.
/// \tparam RecvLayoutType array_layout of the receive buffer.  Must
///   be Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam RecvDeviceType Kokkos::Device specialization used by the
///   receive buffer.  It's OK for this to differ from SendDeviceType.
///   We assume that MPI implementations can handle this.  (This is a
///   reasonable assumption with CUDA-enabled MPI implementations.)
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType>
std::shared_ptr<CommRequest>
wrapIallreduceCommRequest (const std::shared_ptr<CommRequest>& req,
                           const ::Kokkos::View<const PacketType*,
                             SendLayoutType,
                             SendDeviceType>& sendbuf,
                           const ::Kokkos::View<PacketType*,
                             RecvLayoutType,
                             RecvDeviceType>& recvbuf)
{
  static_assert (std::is_same<SendLayoutType, ::Kokkos::LayoutLeft>::value ||
                 std::is_same<SendLayoutType, ::Kokkos::LayoutRight>::value,
                 "SendLayoutType must be either Kokkos::LayoutLeft or "
                 "Kokkos::LayoutRight.");
  static_assert (std::is_same<RecvLayoutType, ::Kokkos::LayoutLeft>::value ||
                 std::is_same<RecvLayoutType, ::Kokkos::LayoutRight>::value,
                 "RecvLayoutType must be either Kokkos::LayoutLeft or "
                 "Kokkos::LayoutRight.");
  typedef IallreduceCommRequest<PacketType, SendLayoutType, SendDeviceType,
    RecvLayoutType, RecvDeviceType, 1> req_type;
  return std::shared_ptr<CommRequest> (new req_type (req, sendbuf, recvbuf));
}

/// \brief Function for wrapping the CommRequest to be returned from
///   ::Tpetra::Details::iallreduce; overload for rank-0 send and
///   receive buffers.
///
/// The object returned from this function keeps the send and receive
/// buffers.  Since ::Kokkos::View reference-counts, this ensures that
/// the buffers will not be deallocated until the iallreduce
/// completes.  The buffer references get cleared on wait() or
/// cancel().
///
/// \param req [in] Pointer to CommRequest from
///   ::Tpetra::Details::iallreduce.
/// \param sendbuf [in] Send buffer for ::Tpetra::Details::iallreduce.
/// \param recvbuf [in] Receive buffer for
///   ::Tpetra::Details::iallreduce.  This is an input argument,
///   because this function itself does not modify the contents of the
///   receive buffer; it just keeps a reference to it.
///
/// \return Pointer to wrapped CommRequest.
///
/// \tparam PacketType Type of each entry of the send and receive
///   buffers.
/// \tparam SendLayoutType array_layout of the send buffer.  Must be
///   Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam SendDeviceType Kokkos::Device specialization used by the
///   send buffer.
/// \tparam RecvLayoutType array_layout of the receive buffer.  Must
///   be Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam RecvDeviceType Kokkos::Device specialization used by the
///   receive buffer.  It's OK for this to differ from SendDeviceType.
///   We assume that MPI implementations can handle this.  (This is a
///   reasonable assumption with CUDA-enabled MPI implementations.)
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType>
std::shared_ptr<CommRequest>
wrapIallreduceCommRequest (const std::shared_ptr<CommRequest>& req,
                           const ::Kokkos::View<const PacketType,
                             SendLayoutType,
                             SendDeviceType>& sendbuf,
                           const ::Kokkos::View<PacketType,
                             RecvLayoutType,
                             RecvDeviceType>& recvbuf)
{
  static_assert (std::is_same<SendLayoutType, ::Kokkos::LayoutLeft>::value ||
                 std::is_same<SendLayoutType, ::Kokkos::LayoutRight>::value,
                 "SendLayoutType must be either Kokkos::LayoutLeft or "
                 "Kokkos::LayoutRight.");
  static_assert (std::is_same<RecvLayoutType, ::Kokkos::LayoutLeft>::value ||
                 std::is_same<RecvLayoutType, ::Kokkos::LayoutRight>::value,
                 "RecvLayoutType must be either Kokkos::LayoutLeft or "
                 "Kokkos::LayoutRight.");
  typedef IallreduceCommRequest<PacketType, SendLayoutType, SendDeviceType,
    RecvLayoutType, RecvDeviceType, 0> req_type;
  return std::shared_ptr<CommRequest> (new req_type (req, sendbuf, recvbuf));
}

#ifdef HAVE_TPETRACORE_MPI

/// \brief Bottom (lowest)-level implementation of
///   ::Tpetra::Details::iallreduce (see bottom of this header file).
///
/// This doesn't need to know about the "Packet" type, because the
/// MPI_Datatype expresses that information (how MPI should
/// communicate Packet instances) in a run-time way.
///
/// This doesn't need to know about the Kokkos "Device" type, because
/// we assume that MPI implementations can read CUDA device memory,
/// host memory, or indeed from any memory space.
std::shared_ptr<CommRequest>
iallreduceRawVoid (const void* sendbuf,
                   void* recvbuf,
                   const int count,
                   MPI_Datatype mpiDatatype,
                   const bool mpiDatatypeNeedsFree,
                   const Teuchos::EReductionType op,
                   MPI_Comm comm);

#endif // HAVE_TPETRACORE_MPI

/// \brief Second lowest-level implementation of
///   ::Tpetra::Details::iallreduce (see bottom of this header file).
///
/// \tparam Packet Type of each entry of the send and receive buffer.
///
/// This doesn't need to know about the Kokkos "Device" type, because
/// we assume that MPI implementations can read CUDA device memory,
/// host memory, or indeed from any memory space.
template<class Packet>
std::shared_ptr<CommRequest>
iallreduceRaw (const Packet sendbuf[],
               Packet recvbuf[],
               const int count,
               const ::Teuchos::EReductionType op,
               const ::Teuchos::Comm<int>& comm)
{
#ifdef HAVE_TPETRACORE_MPI
  using ::Tpetra::Details::MpiTypeTraits;

  // mfh 12 Nov 2016: It's reasonable to assume that if users call
  // this function with count > 1, then users should expect that all
  // sent and received data have the same MPI_Datatype.  If count is
  // zero, then it doesn't matter what datatype we use, as long as
  // it's not MPI_DATATYPE_NULL.  For the latter, see e.g., the
  // discussion here:
  //
  // http://lists.mpi-forum.org/pipermail/mpi-forum/2016-January/003159.html

  // FIXME (mfh 14 Nov 2016) For derived MPI_Datatype, it may make
  // sense for us to cache them somewhere for later reuse, rather than
  // constructing and freeing a new one each time.

  MPI_Datatype mpiDatatype = (count > 0) ?
    MpiTypeTraits<Packet>::getType (sendbuf[0]) :
    MPI_BYTE;
  MPI_Comm rawComm = ::Tpetra::Details::extractMpiCommFromTeuchos (comm);
  return iallreduceRawVoid (sendbuf, recvbuf, count, mpiDatatype,
                            MpiTypeTraits<Packet>::needsFree, op, rawComm);
#else // NOT HAVE_TPETRACORE_MPI
  throw std::logic_error ("Tpetra::Details::Impl::iallreduceRaw: This function "
                          "should never be called if MPI is not enabled.  "
                          "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRACORE_MPI
}

/// \brief Implementation of ::Tpetra::Details::iallreduce.
///
/// \tparam PacketType Type of each entry of the send and receive
///   buffers.
/// \tparam SendLayoutType array_layout of the send buffer.  Must be
///   Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam SendDeviceType Kokkos::Device specialization used by the
///   send buffer.
/// \tparam RecvLayoutType array_layout of the receive buffer.  Must
///   be Kokkos::LayoutLeft or Kokkos::LayoutRight.
/// \tparam RecvDeviceType Kokkos::Device specialization used by the
///   receive buffer.  It's OK for this to differ from SendDeviceType.
///   We assume that MPI implementations can handle this.  (This is a
///   reasonable assumption with CUDA-enabled MPI implementations.)
/// \tparam rank Integer rank of the send and receive buffers.  Must
///   be either 0 or 1.
///
/// The actual implementation of ::Tpetra::Details::iallreduce lives
/// in the partial specializations of this struct.  See below.  We do
/// this complicated struct thing in order to make it easier to
/// overload ::Tpetra::Details::iallreduce for different kinds of
/// output arguments (e.g., rank-0 or rank-1 Kokkos::View).
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType,
         const int rank>
struct Iallreduce {};

//! Partial specialization of Iallreduce for rank-1 send and receive buffers.
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType>
struct Iallreduce<PacketType,
                  SendLayoutType,
                  SendDeviceType,
                  RecvLayoutType,
                  RecvDeviceType,
                  1> {
  typedef ::Kokkos::View<const PacketType*,
                         SendLayoutType,
                         SendDeviceType> send_buffer_type;
  typedef ::Kokkos::View<PacketType*,
                         RecvLayoutType,
                         RecvDeviceType> recv_buffer_type;

  static std::shared_ptr<CommRequest>
  iallreduce (const send_buffer_type& sendbuf,
              const recv_buffer_type& recvbuf,
              const ::Teuchos::EReductionType op,
              const ::Teuchos::Comm<int>& comm)
  {
    static_assert (! std::is_const<PacketType>::value,
                   "PacketType must be a nonconst type.");
    static_assert (std::is_same<SendLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<SendLayoutType, ::Kokkos::LayoutRight>::value,
                   "SendLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
    static_assert (std::is_same<RecvLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<RecvLayoutType, ::Kokkos::LayoutRight>::value,
                   "RecvLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
#ifdef HAVE_TPETRACORE_MPI
    // Avoid instantiating Impl::iallreduceRaw for both T and const T,
    // by canonicalizing to the non-const version of T.
    typedef typename std::remove_const<PacketType>::type packet_type;

    std::shared_ptr<CommRequest> req =
      Impl::iallreduceRaw<packet_type> (sendbuf.data (),
                                        recvbuf.data (),
                                        static_cast<int> (sendbuf.extent (0)),
                                        op, comm);
    return Impl::wrapIallreduceCommRequest (req, sendbuf, recvbuf);
#else // NOT HAVE_TPETRACORE_MPI

    // MPI is disabled, so comm is a SerialComm.
    // Avoid needing to check the SerialComm case in Impl::iallreduce.
    // That lets Impl::iallreduce not need to know about DeviceType.
    ::Kokkos::deep_copy (recvbuf, sendbuf);
    // This request has already finished.  There's nothing more to do.
    return Impl::emptyCommRequest ();
#endif // HAVE_TPETRACORE_MPI
  }
};

//! Partial specialization of Iallreduce for rank-0 send and receive buffers.
template<class PacketType,
         class SendLayoutType,
         class SendDeviceType,
         class RecvLayoutType,
         class RecvDeviceType>
struct Iallreduce<PacketType,
                  SendLayoutType,
                  SendDeviceType,
                  RecvLayoutType,
                  RecvDeviceType,
                  0> {
  // mfh 02 Jan 2017: For rank-0 Views, LayoutLeft and LayoutRight are
  // unfortunately NOT mutually assignable, even though it would make
  // sense for this to be the case.  That's why we have to template on
  // the layout types.
  typedef ::Kokkos::View<const PacketType,
                         SendLayoutType,
                         SendDeviceType> send_buffer_type;
  typedef ::Kokkos::View<PacketType,
                         RecvLayoutType,
                         RecvDeviceType> recv_buffer_type;

  static std::shared_ptr<CommRequest>
  iallreduce (const send_buffer_type& sendbuf,
              const recv_buffer_type& recvbuf,
              const ::Teuchos::EReductionType op,
              const ::Teuchos::Comm<int>& comm)
  {
    static_assert (! std::is_const<PacketType>::value,
                   "PacketType must be a nonconst type.");
    static_assert (std::is_same<SendLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<SendLayoutType, ::Kokkos::LayoutRight>::value,
                   "SendLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
    static_assert (std::is_same<RecvLayoutType, ::Kokkos::LayoutLeft>::value ||
                   std::is_same<RecvLayoutType, ::Kokkos::LayoutRight>::value,
                   "RecvLayoutType must be either Kokkos::LayoutLeft or "
                   "Kokkos::LayoutRight.");
#ifdef HAVE_TPETRACORE_MPI
    // Avoid instantiating Impl::iallreduceRaw for both T and const T,
    // by canonicalizing to the non-const version of T.
    typedef typename std::remove_const<PacketType>::type packet_type;

    std::shared_ptr<CommRequest> req =
      Impl::iallreduceRaw<packet_type> (sendbuf.data (),
                                        recvbuf.data (),
                                        static_cast<int> (1),
                                        op, comm);
    return Impl::wrapIallreduceCommRequest (req, sendbuf, recvbuf);
#else // NOT HAVE_TPETRACORE_MPI

    // MPI is disabled, so comm is a SerialComm.
    // Avoid needing to check the SerialComm case in Impl::iallreduce.
    // That lets Impl::iallreduce not need to know about DeviceType.
    ::Kokkos::deep_copy (recvbuf, sendbuf);
    // This request has already finished.  There's nothing more to do.
    return Impl::emptyCommRequest ();
#endif // HAVE_TPETRACORE_MPI
  }
};

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
  static_assert (Kokkos::Impl::is_view<InputViewType>::value,
                 "InputViewType must be a Kokkos::View specialization.");
  static_assert (Kokkos::Impl::is_view<OutputViewType>::value,
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
  typedef typename InputViewType::array_layout send_layout_type;
  typedef typename OutputViewType::array_layout recv_layout_type;
  typedef typename InputViewType::device_type send_device_type;
  typedef typename OutputViewType::device_type recv_device_type;

  typedef Impl::Iallreduce<packet_type, send_layout_type, send_device_type,
    recv_layout_type, recv_device_type, rank> impl_type;
  return impl_type::iallreduce (sendbuf, recvbuf, op, comm);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_IALLREDUCE_HPP

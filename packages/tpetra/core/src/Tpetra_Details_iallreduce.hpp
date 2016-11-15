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

#include "Teuchos_CommHelpers.hpp" // where EReductionType enum is defined
#include "Teuchos_Details_MpiTypeTraits.hpp"
#include "Kokkos_Core.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {

// Don't rely on anything in this namespace.  You shouldn't be relying
// on anything in Tpetra::Details anyway.  You _especially_ shouldn't
// be relying on anything in Tpetra::Details::Impl.  Consider
// yourselves warned!

namespace Impl {

//! Return an "empty" comm request (waiting on it does nothing).
::Teuchos::RCP< ::Teuchos::CommRequest<int> >
emptyCommRequest ();

//! Return an "empty" comm status.
::Teuchos::RCP< ::Teuchos::CommStatus<int> >
emptyCommStatus ();

/// \brief Part of the Work-around for not having MPI >= 3.
///
/// Substitute for MPI_Iallreduce (which requires MPI 3): defer
/// calling MPI_Allreduce until wait.  It's ugly, but it works.
class DeferredActionRequest : public ::Teuchos::CommRequest<int> {
public:
  /// \brief Default constructor (take no action on wait).
  DeferredActionRequest ();

  /// \brief Constructor that takes an action to defer.
  ///
  /// \brief action [in] An action to take on wait.
  DeferredActionRequest (std::function<void () > action);

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.
  ::Teuchos::RCP< ::Teuchos::CommStatus<int> > wait ();

private:
  std::function<void(void) > action_;
  bool actionTaken_;
};

/// \brief ::Teuchos::CommRequest<int> subclass specifically to be
///   returned from ::Tpetra::Details::iallreduce.
///
/// This subclass keeps the send and receive buffers.  Since
/// ::Kokkos::View reference-counts, this ensures that the buffers
/// will not be deallocated until the iallreduce completes.  The
/// buffer references get cleared on wait().
///
/// Tpetra developers should not use this directly; they should
/// instead create instances of this via the wrapIallreduceCommRequest
/// function (see below).
template<class PacketType, class DeviceType>
class IallreduceCommRequest : public ::Teuchos::CommRequest<int> {
public:
  //! Default constructor.
  IallreduceCommRequest () :
    status_ (emptyCommStatus ())
  {}

  //! Constructor that takes a request, and saved data.
  IallreduceCommRequest (const ::Teuchos::RCP< ::Teuchos::CommRequest<int> > req,
                         const ::Kokkos::View<const PacketType*, DeviceType>& sendbuf,
                         const ::Kokkos::View<PacketType*, DeviceType>& recvbuf) :
    req_ (req),
    sendbuf_ (sendbuf),
    recvbuf_ (recvbuf),
    status_ (emptyCommStatus ())
  {}

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.
  ::Teuchos::RCP< ::Teuchos::CommStatus<int> >
  wait ()
  {
    if (! req_.is_null ()) {
      status_ = req_->wait ();
    }
    // Relinquish references to saved buffers, if we have not already
    // done so.  This operation is idempotent.
    sendbuf_ = ::Kokkos::View<const PacketType*, DeviceType> ();
    recvbuf_ = ::Kokkos::View<PacketType*, DeviceType> ();

    return status_;
  }

private:
  //! The wrapped request
  ::Teuchos::RCP< ::Teuchos::CommRequest<int> > req_;

  //! Saved send buffer from iallreduce call
  ::Kokkos::View<const PacketType*, DeviceType> sendbuf_;

  //! Saved recv buffer from iallreduce call
  ::Kokkos::View<PacketType*, DeviceType> recvbuf_;

  //! Saved status (so that multiple wait() calls are idempotent)
  ::Teuchos::RCP< ::Teuchos::CommStatus<int> > status_;
};


/// \brief Function for wrapping the ::Teuchos::CommRequest<int> to be
///   returned from ::Tpetra::Details::iallreduce.
///
/// The object returned from this function keeps the send and receive
/// buffers.  Since ::Kokkos::View reference-counts, this ensures that
/// the buffers will not be deallocated until the iallreduce
/// completes.  The buffer references get cleared on wait().
template<class PacketType, class DeviceType>
::Teuchos::RCP< ::Teuchos::CommRequest<int> >
wrapIallreduceCommRequest (const ::Teuchos::RCP< ::Teuchos::CommRequest<int> >& req,
                           const ::Kokkos::View<const PacketType*, DeviceType>& sendbuf,
                           const ::Kokkos::View<PacketType*, DeviceType>& recvbuf)
{
  return ::Teuchos::rcp (new IallreduceCommRequest<PacketType, DeviceType> (req, sendbuf, recvbuf));
}

#ifdef HAVE_TEUCHOS_MPI

/// \brief Lowest-level implementation of ::Tpetra::Details::iallreduce.
///
/// This doesn't need to know about Packet, because the MPI_Datatype
/// expresses that information (how MPI should communicate Packet
/// instances) in a run-time way.
///
/// This doesn't need to know about Device, because we assume that MPI
/// implementations can read CUDA device memory, host memory, or
/// indeed from any memory space.  If they can't, Tpetra::DistObject
/// must then first copy to a space from which they can.
Teuchos::RCP<Teuchos::CommRequest<int> >
iallreduceRawVoid (const void* sendbuf,
                   void* recvbuf,
                   const int count,
                   MPI_Datatype mpiDatatype,
                   const bool mpiDatatypeNeedsFree,
                   const Teuchos::EReductionType op,
                   const Teuchos::Comm<int>& comm);

#endif // HAVE_TEUCHOS_MPI

/// \brief Medium-level implementation of ::Tpetra::Details::iallreduce.
///
/// This doesn't need to know about Device, because we assume that MPI
/// implementations can read CUDA device memory, host memory, or
/// indeed from any memory space.  If they can't, Tpetra::DistObject
/// must then first copy to a space from which they can.
template<class Packet>
::Teuchos::RCP< ::Teuchos::CommRequest<int> >
iallreduceRaw (const Packet sendbuf[],
               Packet recvbuf[],
               const int count,
               const ::Teuchos::EReductionType op,
               const ::Teuchos::Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using ::Teuchos::Details::MpiTypeTraits;
  using ::Teuchos::CommRequest;
  using ::Teuchos::RCP;

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
  return iallreduceRawVoid (sendbuf, recvbuf, count, mpiDatatype,
                            MpiTypeTraits<Packet>::needsFree, op, comm);
#else // NOT HAVE_TEUCHOS_MPI
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Tpetra::Details::Impl::iallreduceRaw: "
     "This function should never be called if MPI is not enabled.  "
     "Please report this bug to the Tpetra developers.");
#endif // HAVE_TEUCHOS_MPI
}

} // namespace Impl


/// \brief Nonblocking all-reduce
///
/// \tparam PacketType Type of the object(s) to send and receive
/// \tparam DeviceType Kokkos::Device specialization
///
/// This function wraps MPI_Iallreduce.  It does a nonblocking
/// all-reduce over the input communicator \c comm, from \c sendbuf
/// into \c recvbuf, using \c op as the reduction operator.  The
/// function returns without blocking; the all-reduce only blocks for
/// completion when one calls wait() on the returned request.
///
/// \param sendbuf [in] Input array of PacketType.
/// \param recvbuf [out] Output array of PacketType.
/// \param op [in] Teuchos enum representing the reduction operator.
/// \param comm [in] Communicator over which to do the all-reduce.
///
/// \c sendbuf and \c recvbuf must either be disjoint, or identical
/// (point to the same array).  They may not partially overlap.
/// Furthermore, if they are identical, the input communicator must be
/// an intracommunicator.  It may not be an intercommunicator.  If you
/// don't know what an intercommunicator is, you probably just have an
/// intracommunicator, so everything is fine.
template<class PacketType, class DeviceType>
::Teuchos::RCP< ::Teuchos::CommRequest<int> >
iallreduce (const ::Kokkos::View<const PacketType*, DeviceType>& sendbuf,
            const ::Kokkos::View<PacketType*, DeviceType>& recvbuf,
            const ::Teuchos::EReductionType op,
            const ::Teuchos::Comm<int>& comm)
{
  static_assert (! std::is_const<PacketType>::value,
                 "PacketType must be a nonconst type.");

  // Avoid instantiating Impl::iallreduceRaw for both T and const T,
  // by canonicalizing to the non-const version of T.
  typedef typename std::remove_const<PacketType>::type packet_type;

  if (comm.getSize () == 1) {
    // Avoid needing to check the SerialComm case in Impl::iallreduce.
    // That lets Impl::iallreduce not need to know about DeviceType.
    ::Kokkos::deep_copy (recvbuf, sendbuf);
    // This request has already finished.  There's nothing more to do.
    return Impl::emptyCommRequest ();
  }
  else {
    using ::Teuchos::CommRequest;
    using ::Teuchos::inOutArg;
    using ::Teuchos::RCP;
    using ::Teuchos::set_extra_data;

    RCP<CommRequest<int> > req =
      Impl::iallreduceRaw<packet_type> (sendbuf.ptr_on_device (),
                                        recvbuf.ptr_on_device (),
                                        static_cast<int> (sendbuf.dimension_0 ()),
                                        op, comm);
    // NOTE (mfh 14 Nov 2016) We can't use Teuchos::set_extra_data to
    // attach sendbuf and recvbuf to the Teuchos::RCP, because that
    // weirdly requires operator<< to be defined (via Teuchos::any)
    // for Kokkos::View.  It's not.  Our work-around is to use a
    // custom CommRequest subclass that "wraps" the Kokkos::View in a
    // Teuchos::OpaqueWrapper.  This makes sense, because Kokkos::View
    // has shallow copy semantics.  It also would easily permit us to
    // replace RCP with std::shared_ptr, if we want.

    // set_extra_data (sendbuf, "sendbuf", inOutArg (req));
    // set_extra_data (recvbuf, "recvbuf", inOutArg (req));
    // return req;

    return Impl::wrapIallreduceCommRequest<packet_type, DeviceType> (req, sendbuf, recvbuf);
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_IALLREDUCE_HPP

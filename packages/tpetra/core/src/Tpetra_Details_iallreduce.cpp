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

#include "Tpetra_Details_iallreduce.hpp"

#ifdef HAVE_TPETRACORE_MPI
#  include "Teuchos_DefaultMpiComm.hpp" // only needs to be in .cpp file
#endif // HAVE_TPETRACORE_MPI
#include "Teuchos_DefaultSerialComm.hpp" // only needs to be in .cpp file

#ifdef HAVE_TPETRACORE_MPI
namespace { // (anonymous)
std::string getMpiErrorString (const int errCode) {
  // Space for storing the error string returned by MPI.
  // Leave room for null termination, since I don't know if MPI does this.
  char errString [MPI_MAX_ERROR_STRING+1];
  int errStringLen = MPI_MAX_ERROR_STRING; // output argument
  (void) MPI_Error_string (errCode, errString, &errStringLen);
  // errStringLen on output is the number of characters written.
  // I'm not sure (the MPI 3.0 Standard doesn't say) if this
  // includes the '\0', so I'll make sure.  We reserved space for
  // the extra '\0' if needed.
  if (errString[errStringLen-1] != '\0') {
    errString[errStringLen] = '\0';
  }
  return std::string (errString); // This copies the original string.
}
} // namespace (anonymous)
#endif // HAVE_TPETRACORE_MPI

namespace Tpetra {
namespace Details {
namespace Impl {

#ifdef HAVE_TPETRACORE_MPI
MpiCommRequest::
MpiCommRequest () :
  req_ (MPI_REQUEST_NULL)
{}

MpiCommRequest::
MpiCommRequest (const MPI_Request req) :
  req_ (req)
{}

void
MpiCommRequest::
waitWithStatus (MPI_Status& status)
{
  if (req_ != MPI_REQUEST_NULL) {
    MPI_Request req = req_;
    const int err = MPI_Wait (&req, &status);
    TEUCHOS_TEST_FOR_EXCEPTION
      (err != MPI_SUCCESS, std::runtime_error,
       "MpiCommRequest::waitWithStatus: MPI_Wait failed with error \""
       << getMpiErrorString (err));
    // MPI_Wait should set the MPI_Request to MPI_REQUEST_NULL on
    // success.  We'll do it here just to be conservative.
    req_ = MPI_REQUEST_NULL;
  }
}

void
MpiCommRequest::
wait ()
{
  if (req_ != MPI_REQUEST_NULL) {
    MPI_Request req = req_;
    const int err = MPI_Wait (&req, MPI_STATUS_IGNORE);
    TEUCHOS_TEST_FOR_EXCEPTION
      (err != MPI_SUCCESS, std::runtime_error,
       "MpiCommRequest::wait: MPI_Wait failed with error \""
       << getMpiErrorString (err));
    // MPI_Wait should set the MPI_Request to MPI_REQUEST_NULL on
    // success.  We'll do it here just to be conservative.
    req_ = MPI_REQUEST_NULL;
  }
}

void
MpiCommRequest::
cancel ()
{
  if (req_ != MPI_REQUEST_NULL) {
    const int err = MPI_Cancel (&req_);
    TEUCHOS_TEST_FOR_EXCEPTION
      (err != MPI_SUCCESS, std::runtime_error,
       "MpiCommRequest::cancel: MPI_Cancel failed with the following error: "
       << getMpiErrorString (err));

    // Wait on the request.  MPI requires doing this after cancel, and
    // promises that a wait after cancel is a local operation.
    this->wait ();

    // The returned status may still be useful; for example, one may
    // call MPI_Test_cancelled to test an MPI_Status from a
    // nonblocking send.  For now, we'll ignore it.
  }
}

MpiCommRequest::
~MpiCommRequest ()
{
  if (req_ != MPI_REQUEST_NULL) {
    // We're in a destructor, so don't throw errors.  However, if
    // MPI_Cancel fails, it's probably a bad idea to call MPI_Wait.
    const int err = MPI_Cancel (&req_);
    if (err == MPI_SUCCESS) {
      // The MPI_Cancel succeeded.  Now wait on the request.  Ignore
      // any reported error, since we can't do anything about those in
      // the destructor (other than kill the program).  If successful,
      // MPI_Wait will set the MPI_Request to MPI_REQUEST_NULL.  We
      // ignore the returned MPI_Status, since if users let the
      // request fall out of scope, then they must not care about the
      // status.
      //
      // mfh 21 Oct 2012: The MPI standard requires completing a
      // canceled request by calling a function like MPI_Wait,
      // MPI_Test, or MPI_Request_free.  MPI_Wait on a canceled
      // request behaves like a local operation (it does not
      // communicate or block waiting for communication).  One could
      // also call MPI_Request_free instead of MPI_Wait, but
      // MPI_Request_free is intended more for persistent requests
      // (created with functions like MPI_Recv_init).
      (void) MPI_Wait (&req_, MPI_STATUS_IGNORE);
    }
  }
}

#endif // HAVE_TPETRACORE_MPI

std::shared_ptr<CommRequest>
emptyCommRequest ()
{
  return std::shared_ptr<CommRequest> (new DeferredActionCommRequest ());
}

DeferredActionCommRequest::
DeferredActionCommRequest () :
  action_ ([] () {}), // do nothing
  actionTaken_ (false)
{}

DeferredActionCommRequest::
DeferredActionCommRequest (std::function<void () > action) :
  action_ (action),
  actionTaken_ (false)
{}

void
DeferredActionCommRequest::
wait ()
{
  if (! actionTaken_) {
    action_ ();
    actionTaken_ = true;
  }
}

void
DeferredActionCommRequest::
cancel ()
{
  actionTaken_ = true;
}

#ifdef HAVE_TPETRACORE_MPI

std::shared_ptr<CommRequest>
iallreduceRawVoid (const void* sendbuf,
                   void* recvbuf,
                   const int count,
                   MPI_Datatype mpiDatatype,
                   const bool mpiDatatypeNeedsFree,
                   const Teuchos::EReductionType op,
                   MPI_Comm comm)
{
  MPI_Op rawOp = ::Teuchos::Details::getMpiOpForEReductionType (op);

#if MPI_VERSION >= 3
  const bool useMpi3 = true;
#else
  const bool useMpi3 = false;
#endif // MPI_VERSION >= 3

  // Fix for #852: always build the fall-back (MPI_VERSION < 3)
  // implementation.
  if (useMpi3) {
#if MPI_VERSION >= 3
    MPI_Request rawRequest = MPI_REQUEST_NULL;
    int err = MPI_SUCCESS;
    if (sendbuf == recvbuf) {
      // Fix for #850.  This only works if rawComm is an
      // intracommunicator.  Intercommunicators don't have an in-place
      // option for collectives.
      err = MPI_Iallreduce (MPI_IN_PLACE, recvbuf, count, mpiDatatype,
                            rawOp, comm, &rawRequest);
    }
    else {
      err = MPI_Iallreduce (sendbuf, recvbuf, count, mpiDatatype,
                            rawOp, comm, &rawRequest);
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (err != MPI_SUCCESS, std::runtime_error,
       "MPI_Iallreduce failed with the following error: "
       << getMpiErrorString (err));
    if (mpiDatatypeNeedsFree) {
      // As long as the MPI_Datatype goes into MPI_Iallreduce, it's OK
      // to free it, even if the MPI_Iallreduce has not yet completed.
      // There's no sense in checking the error code here.
      (void) MPI_Type_free (&mpiDatatype);
    }
    return std::shared_ptr<CommRequest> (new MpiCommRequest (rawRequest));
#else
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Should never get here.  "
       "Please report this bug to the Tpetra developers.");
#endif // MPI_VERSION >= 3
  }
  else { // ! useMpi3
    // We don't have MPI_Iallreduce.  The next best thing is to defer an
    // MPI_Allreduce call until wait.  We do this by returning a
    // "DeferredActionCommRequest," which is just a wrapped
    // std::function.
    //
    // NOTE (mfh 12 Nov 2016, 14 Nov 2016) One issue with this approach
    // is that we have to make sure that the MPI_Datatype won't go away
    // before the MPI_Allreduce gets called.  We handle this for now by
    // calling MPI_Type_dup and stashing the destructor in the request.
    // (Don't use the MPI_COMM_SELF trick here, unless you first check
    // whether you've seen that MPI_Datatype before -- otherwise you'll
    // get memory growth linear in the number of iallreduce calls.)
    return std::shared_ptr<CommRequest> (new DeferredActionCommRequest ([=] () {
          // It could be that this action got deferred beyond
          // MPI_Finalize.  In that case, do nothing.
          int mpiInitialized = 0;
          (void) MPI_Initialized (&mpiInitialized);
          int mpiFinalized = 0;
          (void) MPI_Finalized (&mpiFinalized);
          if (mpiFinalized == 0 && mpiInitialized != 0) {
            // FIXME (mfh 14 Nov 2016) Unfortunately, there is no
            // MPI_Op_dup, so I can't guarantee that the input MPI_Op
            // will still exist to the point where it is actually
            // used.
            //
            // FIXME (mfh 14 Nov 2016) Best practice would be to
            // duplicate the input MPI_Comm, so that we can ensure its
            // survival to this point.  However, we can't guarantee
            // survival of the input MPI_Op, so we might as well just
            // not bother.
            if (mpiDatatypeNeedsFree) {
              // Copy the original MPI_Datatype, so that we can safely
              // defer this call past survival of the original.
              MPI_Datatype dupDatatype;
              (void) MPI_Type_dup (mpiDatatype, &dupDatatype);
#if MPI_VERSION >= 3
              if (sendbuf == recvbuf) {
                (void) MPI_Allreduce (MPI_IN_PLACE, recvbuf, count, dupDatatype,
                                      rawOp, comm);
              }
              else {
                (void) MPI_Allreduce (sendbuf, recvbuf, count, dupDatatype,
                                      rawOp, comm);
              }
#else // MPI_VERSION < 3
              if (sendbuf == recvbuf) {
                (void) MPI_Allreduce (MPI_IN_PLACE, recvbuf,
                                      count, dupDatatype, rawOp, comm);
              }
              else {
                // OpenMPI 1.6.5 insists on void*, not const void*, for sendbuf.
                (void) MPI_Allreduce (const_cast<void*> (sendbuf), recvbuf,
                                      count, dupDatatype, rawOp, comm);
              }
#endif // MPI_VERSION >= 3
              (void) MPI_Type_free (&dupDatatype);
            }
            else {
#if MPI_VERSION >= 3
              if (sendbuf == recvbuf) {
                (void) MPI_Allreduce (MPI_IN_PLACE, recvbuf, count, mpiDatatype,
                                      rawOp, comm);
              }
              else {
                (void) MPI_Allreduce (sendbuf, recvbuf, count, mpiDatatype,
                                      rawOp, comm);
              }
#else // MPI_VERSION < 3
              if (sendbuf == recvbuf) {
                (void) MPI_Allreduce (MPI_IN_PLACE, recvbuf,
                                      count, mpiDatatype, rawOp, comm);
              }
              else {
                // OpenMPI 1.6.5 insists on void*, not const void*, for sendbuf.
                (void) MPI_Allreduce (const_cast<void*> (sendbuf), recvbuf,
                                      count, mpiDatatype, rawOp, comm);
              }
#endif // MPI_VERSION >= 3
            }
          }
        }));
  } // useMpi3
}

#endif // HAVE_TPETRACORE_MPI

} // namespace Impl

std::shared_ptr<CommRequest>
iallreduce (const int localValue,
            int& globalValue,
            const ::Teuchos::EReductionType op,
            const ::Teuchos::Comm<int>& comm)
{
  using input_view_type =
    Kokkos::View<const int*, Kokkos::HostSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using output_view_type =
    Kokkos::View<int*, Kokkos::HostSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  input_view_type localView (&localValue, 1);
  output_view_type globalView (&globalValue, 1);
  return ::Tpetra::Details::iallreduce (localView, globalView, op, comm);
}

} // namespace Details
} // namespace Tpetra



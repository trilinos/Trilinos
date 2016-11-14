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

//#include "Teuchos_config.h" // HAVE_TEUCHOS_MPI macro def (already included)
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_Details_MpiTypeTraits.hpp"
#  include "Teuchos_DefaultMpiComm.hpp" // only needs to be in .cpp file
#endif // HAVE_TEUCHOS_MPI

namespace Tpetra {
namespace Details {
namespace Impl {

::Teuchos::RCP< ::Teuchos::CommRequest<int> >
emptyCommRequest ()
{
#ifdef HAVE_TEUCHOS_MPI
  typedef ::Teuchos::MpiCommRequestBase<int> req_type;
  return ::Teuchos::rcp (new req_type (MPI_REQUEST_NULL));
#else // HAVE_TEUCHOS_MPI
  return ::Teuchos::rcp (new ::Tpetra::Details::Impl::DeferredActionRequest ());
#endif // HAVE_TEUCHOS_MPI
}

::Teuchos::RCP< ::Teuchos::CommStatus<int> >
emptyCommStatus ()
{
  return ::Teuchos::rcp (new ::Teuchos::SerialCommStatus<int> ());
}

DeferredActionRequest::
DeferredActionRequest () :
  action_ ([] () {}), // do nothing
  actionTaken_ (false)
{}

DeferredActionRequest::
DeferredActionRequest (std::function<void () > action) :
  action_ (action),
  actionTaken_ (false)
{}

::Teuchos::RCP< ::Teuchos::CommStatus<int> >
DeferredActionRequest::wait ()
{
  if (! actionTaken_) {
    action_ ();
  }
  return emptyCommStatus ();
}

#ifdef HAVE_TEUCHOS_MPI

Teuchos::RCP<Teuchos::CommRequest<int> >
iallreduceRawVoid (const void* sendbuf,
                   void* recvbuf,
                   const int count,
                   MPI_Datatype mpiDatatype,
                   const bool mpiDatatypeNeedsFree,
                   const Teuchos::EReductionType op,
                   const Teuchos::Comm<int>& comm)
{
  using ::Teuchos::MpiComm;
  using ::Teuchos::rcp;

  // We've excluded the SerialComm case in the wrapper function
  // already, by checking whether comm.getSize() == 1.  Thus, we only
  // need to test the MpiComm case.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Tpetra::Details::Impl::iallreduceRawVoid: "
       "Not implemented for SerialComm, or for any other Comm subclass other "
       "than MpiComm.  We should never reach this point.  "
       "Please report this bug to the Tpetra developers.");
  }
  else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawComm = * (mpiComm->getRawMpiComm ());
    MPI_Op rawOp = ::Teuchos::Details::getMpiOpForEReductionType (op);

#if MPI_VERSION >= 3
    MPI_Request rawRequest = MPI_REQUEST_NULL;
    const int err = MPI_Iallreduce (sendbuf, recvbuf, count, mpiDatatype,
                                    rawOp, rawComm, &rawRequest);
    TEUCHOS_TEST_FOR_EXCEPTION
      (err != MPI_SUCCESS, std::runtime_error,
       "MPI_Iallreduce failed with the following error: "
       << ::Teuchos::Details::getMpiErrorString (err));
    if (mpiDatatypeNeedsFree) {
      // As long as the MPI_Datatype goes into MPI_Iallreduce, it's OK
      // to free it, even if the MPI_Iallreduce has not yet completed.
      // There's no sense in checking the error code here.
      (void) MPI_Type_free (&mpiDatatype);
    }
    return rcp (new ::Teuchos::MpiCommRequestBase<int> (rawRequest));
#else
    // We don't have MPI_Iallreduce.  The next best thing is to defer
    // an MPI_Allreduce call until wait.  We do this by returning a
    // "DeferredActionRequest," which is just a wrapped std::function.
    //
    // NOTE (mfh 12 Nov 2016, 14 Nov 2016) One issue with this
    // approach is that we have to make sure that the MPI_Datatype
    // won't go away before the MPI_Allreduce gets called.  We handle
    // this for now by calling MPI_Type_dup and stashing the
    // destructor in the request.  (Don't use the MPI_COMM_SELF trick
    // here, unless you first check whether you've seen that
    // MPI_Datatype before -- otherwise you'll get memory growth
    // linear in the number of iallreduce calls.)
    return rcp (new Impl::DeferredActionRequest ([=] () {
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
              (void) MPI_Type_dup (mpiDatatype, dupDatatype);
              (void) MPI_Allreduce (sendbuf, recvbuf, count, dupDatatype,
                                    rawOp, rawComm);
              (void) MPI_Type_free (&dupDatatype);
            }
            else {
              (void) MPI_Allreduce (sendbuf, recvbuf, count, mpiDatatype,
                                    rawOp, rawComm);
            }
          }
        }));
#endif // MPI_VERSION >= 3
  }
}

#endif // HAVE_TEUCHOS_MPI

} // namespace Impl
} // namespace Details
} // namespace Tpetra



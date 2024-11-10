// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_iallreduce.hpp"

#ifdef HAVE_TPETRACORE_MPI
#  include "Teuchos_DefaultMpiComm.hpp" // only needs to be in .cpp file
#endif // HAVE_TPETRACORE_MPI
#include "Teuchos_DefaultSerialComm.hpp" // only needs to be in .cpp file

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_MPI
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
#endif // HAVE_TPETRACORE_MPI

namespace Impl {

std::shared_ptr<CommRequest>
emptyCommRequest ()
{
  return std::shared_ptr<CommRequest> (new CommRequest());
}

#ifdef HAVE_TPETRACORE_MPI

#if MPI_VERSION >= 3
MPI_Request
iallreduceRaw (const void* sendbuf,
               void* recvbuf,
               const int count,
               MPI_Datatype mpiDatatype,
               const Teuchos::EReductionType op,
               MPI_Comm comm)
{
  MPI_Op rawOp = ::Teuchos::Details::getMpiOpForEReductionType (op);
  MPI_Request req = MPI_REQUEST_NULL;
  int err = MPI_SUCCESS;
  if (sendbuf == recvbuf) {
    // Fix for #850.  This only works if rawComm is an
    // intracommunicator.  Intercommunicators don't have an in-place
    // option for collectives.
    err = MPI_Iallreduce (MPI_IN_PLACE, recvbuf, count, mpiDatatype,
                          rawOp, comm, &req);
  }
  else {
    err = MPI_Iallreduce (sendbuf, recvbuf, count, mpiDatatype,
                          rawOp, comm, &req);
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::runtime_error,
     "MPI_Iallreduce failed with the following error: "
     << getMpiErrorString (err));
  return req;
}
#endif  //MPI >= 3

void
allreduceRaw  (const void* sendbuf,
               void* recvbuf,
               const int count,
               MPI_Datatype mpiDatatype,
               const Teuchos::EReductionType op,
               MPI_Comm comm)
{
  MPI_Op rawOp = ::Teuchos::Details::getMpiOpForEReductionType (op);
  int err = MPI_SUCCESS;
  if (sendbuf == recvbuf) {
    err = MPI_Allreduce (MPI_IN_PLACE, recvbuf,
        count, mpiDatatype, rawOp, comm);
  }
  else {
    // OpenMPI 1.6.5 insists on void*, not const void*, for sendbuf.
    (void) MPI_Allreduce (const_cast<void*> (sendbuf), recvbuf,
        count, mpiDatatype, rawOp, comm);
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (err != MPI_SUCCESS, std::runtime_error,
     "MPI_Allreduce failed with the following error: "
     << getMpiErrorString (err));
}

#endif // HAVE_TPETRACORE_MPI

} // namespace Impl

std::shared_ptr<CommRequest>
iallreduce (const int localValue,
            int& globalValue,
            const ::Teuchos::EReductionType op,
            const ::Teuchos::Comm<int>& comm)
{
  //Input: needs to be an owning view containing localValue
  Kokkos::View<int*, Kokkos::HostSpace> localView(
      Kokkos::ViewAllocateWithoutInitializing("localValue"), 1);
  localView(0) = localValue;
  Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    globalView(&globalValue, 1);
  return ::Tpetra::Details::iallreduce<decltype(localView), decltype(globalView)>(localView, globalView, op, comm);
}

} // namespace Details
} // namespace Tpetra


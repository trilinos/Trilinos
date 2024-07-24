// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_PackTriples.hpp"

#ifdef HAVE_TPETRACORE_MPI

namespace { // (anonymous)
  std::string
  mpiErrorCodeToString (const int errCode)
  {
    if (errCode == MPI_SUCCESS) {
      return "MPI_SUCCESS";
    }
    else {
      char rawErrString[MPI_MAX_ERROR_STRING];
      int len = 0;
      int err = MPI_Error_string (errCode, rawErrString, &len);
      if (err != MPI_SUCCESS) {
        // Assume that the string wasn't written.  This means it might
        // not be null terminated, so make it a valid empty string by
        // writing the null termination character to it.
        if (MPI_MAX_ERROR_STRING > 0) {
          rawErrString[0] = '\0';
        }
      }
      return std::string (rawErrString);
    }
  }
} // namespace (anonymous)

#endif // HAVE_TPETRACORE_MPI

namespace Tpetra {
namespace Details {

#ifdef HAVE_TPETRACORE_MPI

int
countPackTriplesCountMpi (MPI_Comm comm,
                          int& size,
                          std::ostream* errStrm)
{
  using std::endl;

  int curSize = 0;
  const int errCode = MPI_Pack_size (1, MPI_INT, comm, &curSize);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "countPackTripleMpi: MPI_Pack_size failed on "
               << "MPI_INT call.  MPI reports: "
               << mpiErrorCodeToString (errCode) << endl;
    }
    return errCode;
  }
  size = curSize; // "commit" the result

  return errCode;
}

int
packTriplesCountMpi (const int numEnt,
                     char outBuf[],
                     const int outBufSize,
                     int& outBufCurPos,
                     MPI_Comm comm,
                     std::ostream* errStrm)
{
  using std::endl;

  // mfh 19 Jan 2017: Some (generally older) MPI implementations want
  // the first argument of MPI_Pack to be a pointer to nonconst, even
  // though it's an input argument that MPI_Pack does not modify.
  int theNumEnt = numEnt;
  const int errCode = MPI_Pack (&theNumEnt, 1, MPI_INT, outBuf,
                                outBufSize, &outBufCurPos, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "packTriplesCountMpi: MPI_Pack failed with outBufSize="
               << outBufSize << " and outBufCurPos=" << outBufCurPos
               << ".  MPI reports: " << mpiErrorCodeToString (errCode)
               << endl;
    }
    return errCode;
  }
  return errCode;
}

int
unpackTriplesCountMpi (const char inBuf[],
                       const int inBufSize,
                       int& inBufCurPos,
                       int& numEnt,
                       MPI_Comm comm,
                       std::ostream* errStrm)
{
  using std::endl;
  int errCode = MPI_SUCCESS;

  // mfh 19 Jan 2017: Some (generally older) MPI implementations want
  // the first argument to be a nonconst pointer, even though
  // MPI_Unpack does not modify that argument.
  int theNumEnt = 0;
  errCode = MPI_Unpack (const_cast<char*> (inBuf), inBufSize, &inBufCurPos,
                        &theNumEnt, 1, MPI_INT, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "unpackTriplesCountMpi: MPI_Unpack failed on gblRowInd: "
               << "inBufSize=" << inBufSize
               << ", inBufCurPos=" << inBufCurPos
               << ".  MPI reports: " << mpiErrorCodeToString (errCode)
               << endl;
    }
    return errCode;
  }
  if (theNumEnt < 0) {
    if (errStrm != NULL) {
      *errStrm << "unpackTriplesCountMpi: The unpacked number of entries "
               << theNumEnt << " is negative." << endl;
    }
    return -1;
  }
  numEnt = theNumEnt; // "commit" the change to the output argument
  return errCode;
}

#endif // HAVE_TPETRACORE_MPI

int
#ifdef HAVE_TPETRACORE_MPI
countPackTriplesCount (const ::Teuchos::Comm<int>& comm,
                       int& size,
                       std::ostream* errStrm)
#else // NOT HAVE_TPETRACORE_MPI
countPackTriplesCount (const ::Teuchos::Comm<int>& /* comm */,
                       int& size,
                       std::ostream* errStrm)
#endif // HAVE_TPETRACORE_MPI
{
#ifdef HAVE_TPETRACORE_MPI
  using ::Tpetra::Details::extractMpiCommFromTeuchos;
  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  return countPackTriplesCountMpi (mpiComm, size, errStrm);
#else // NOT HAVE_TPETRACORE_MPI
  using std::endl;
  if (errStrm != NULL) {
    *errStrm << "countPackTriplesCount: Not implemented (no need; there's no "
      "need to pack or unpack anything if there's only one process)." << endl;
  }
  return -1;
#endif // HAVE_TPETRACORE_MPI
}

int
#ifdef HAVE_TPETRACORE_MPI
packTriplesCount (const int numEnt,
                  char outBuf[],
                  const int outBufSize,
                  int& outBufCurPos,
                  const ::Teuchos::Comm<int>& comm,
                  std::ostream* errStrm)
#else // NOT HAVE_TPETRACORE_MPI
packTriplesCount (const int /* numEnt */,
                  char /* outBuf */[],
                  const int /* outBufSize */,
                  int& /* outBufCurPos */,
                  const ::Teuchos::Comm<int>& /* comm */,
                  std::ostream* errStrm)
#endif // HAVE_TPETRACORE_MPI
{
#ifdef HAVE_TPETRACORE_MPI
  using ::Tpetra::Details::extractMpiCommFromTeuchos;
  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  return packTriplesCountMpi (numEnt, outBuf, outBufSize,
                              outBufCurPos, mpiComm, errStrm);
#else // NOT HAVE_TPETRACORE_MPI
  if (errStrm != NULL) {
    *errStrm << "packTriplesCount: Not implemented (no need; there's no need "
      "to pack or unpack anything if there's only one process)." << std::endl;
  }
  return -1;
#endif // HAVE_TPETRACORE_MPI
}

int
#ifdef HAVE_TPETRACORE_MPI
unpackTriplesCount (const char inBuf[],
                    const int inBufSize,
                    int& inBufCurPos,
                    int& numEnt, // output argument!
                    const ::Teuchos::Comm<int>& comm,
                    std::ostream* errStrm)
#else // NOT HAVE_TPETRACORE_MPI
unpackTriplesCount (const char /* inBuf */[],
                    const int /* inBufSize */,
                    int& /* inBufCurPos */,
                    int& /* numEnt */,
                    const ::Teuchos::Comm<int>& /* comm */,
                    std::ostream* errStrm)
#endif // HAVE_TPETRACORE_MPI
{
#ifdef HAVE_TPETRACORE_MPI
  using ::Tpetra::Details::extractMpiCommFromTeuchos;

  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  const int errCode =
    unpackTriplesCountMpi (inBuf, inBufSize, inBufCurPos,
                           numEnt, mpiComm, errStrm);
  return errCode;

#else // NOT HAVE_TPETRACORE_MPI
  if (errStrm != NULL) {
    *errStrm << "unpackTriplesCount: Not implemented (no need; there's no need "
      "to pack or unpack anything if there's only one process)." << std::endl;
  }
  return -1;
#endif // HAVE_TPETRACORE_MPI
}

} // namespace Details
} // namespace Tpetra


// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PACKTRIPLES_HPP
#define TPETRA_DETAILS_PACKTRIPLES_HPP

#include "TpetraCore_config.h"
#include "Teuchos_Comm.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#  include "Tpetra_Details_MpiTypeTraits.hpp"
#endif // HAVE_TPETRACORE_MPI
#include <ostream>

namespace Tpetra {
namespace Details {

//
// Search for "SKIP DOWN TO HERE" (omit quotes) for the "public"
// interface.  I put "public" in quotes because it's public only for
// Tpetra developers, NOT for Tpetra users.
//

#ifdef HAVE_TPETRACORE_MPI

// It's not useful to pack or unpack if not using MPI.  Thus, we only
// need to define these with MPI.  For a non-MPI build, we just need
// some outer stub interface that throws an exception.

int
countPackTriplesCountMpi (MPI_Comm comm,
                          int& size,
                          std::ostream* errStrm = NULL);

int
packTriplesCountMpi (const int numEnt,
                     char outBuf[],
                     const int outBufSize,
                     int& outBufCurPos,
                     MPI_Comm comm,
                     std::ostream* errStrm = NULL);

template<class ScalarType, class OrdinalType>
int
countPackTriplesMpi (MPI_Datatype ordinalDt,
                     MPI_Datatype scalarDt,
                     const int numEnt,
                     MPI_Comm comm,
                     int& size,
                     std::ostream* errStrm = NULL)
{
  using std::endl;
  int errCode = MPI_SUCCESS;

  int totalSize = 0; // return value

  // Count the global row and column indices.
  {
    int curSize = 0;
    // We're packing two ordinals, the row resp. column index.
    errCode = MPI_Pack_size (2, ordinalDt, comm, &curSize);
    if (errCode != MPI_SUCCESS) {
      if (errStrm != NULL) {
        *errStrm << "countPackTripleMpi: MPI_Pack_size failed on "
          "ordinalDt call" << endl;
      }
      return errCode;
    }
    totalSize += curSize;
  }
  // Count the matrix value.
  {
    int curSize = 0;
    errCode = MPI_Pack_size (1, scalarDt, comm, &curSize);
    if (errCode != MPI_SUCCESS) {
      if (errStrm != NULL) {
        *errStrm << "countPackTripleMpi: MPI_Pack_size failed on "
          "scalarDt call" << endl;
      }
      return errCode;
    }
    totalSize += curSize; // one Scalar
  }
  totalSize *= numEnt; // all the entries we want to pack

  size = totalSize; // "commit" the result
  return errCode;
}

template<class ScalarType, class OrdinalType>
int
packTripleMpi (const OrdinalType gblRowInd,
               const OrdinalType gblColInd,
               MPI_Datatype ordinalDt,
               const ScalarType& val,
               MPI_Datatype scalarDt,
               char outBuf[],
               const int outBufSize,
               int& outBufCurPos,
               MPI_Comm comm,
               std::ostream* errStrm = NULL)
{
  using std::endl;
  int errCode = MPI_SUCCESS;

  // Combine the two indices into a single MPI_Pack call.
  OrdinalType inds[2] = {gblRowInd, gblColInd};
  errCode = MPI_Pack (inds, 2, ordinalDt,
                      outBuf, outBufSize, &outBufCurPos, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "packTripleMpi: MPI_Pack failed for indices i=" << gblRowInd
               << ", j=" << gblColInd << ", where outBufSize=" << outBufSize
               << " and outBufCurPos=" << outBufCurPos << "." << endl;
    }
    return errCode;
  }
  // mfh 17,20 Jan 2017: Some (generally older) MPI implementations
  // want the first argument to be a pointer to nonconst, even though
  // MPI_Pack does not modify that argument.
  errCode = MPI_Pack (const_cast<ScalarType*> (&val), 1, scalarDt,
                      outBuf, outBufSize, &outBufCurPos, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "packTripleMpi: MPI_Pack failed on val = " << val
               << endl;
    }
    return errCode;
  }
  return errCode;
}

template<class ScalarType, class OrdinalType>
int
packTriplesMpi (const OrdinalType gblRowInds[],
                const OrdinalType gblColInds[],
                MPI_Datatype ordinalDt,
                const ScalarType vals[],
                MPI_Datatype scalarDt,
                const int numEnt,
                char outBuf[],
                const int outBufSize,
                int& outBufCurPos,
                MPI_Comm comm,
                std::ostream* errStrm = NULL)
{
  using std::endl;
  int errCode = MPI_SUCCESS;

  for (int k = 0; k < numEnt; ++k) {
    errCode = packTripleMpi (gblRowInds[k], gblColInds[k], ordinalDt,
                             vals[k], scalarDt, outBuf, outBufSize,
                             outBufCurPos, comm, errStrm);
    if (errCode != MPI_SUCCESS) {
      if (errStrm != NULL) {
        *errStrm << "packTriplesMpi: packTripleMpi failed at entry k=" << k
                 << endl;
      }
      return errCode;
    }
  }
  return errCode;
}

template<class ScalarType, class OrdinalType>
int
unpackTripleMpi (const char inBuf[],
                 const int inBufSize,
                 int& inBufCurPos,
                 OrdinalType& gblRowInd,
                 OrdinalType& gblColInd,
                 MPI_Datatype ordinalDt,
                 ScalarType& val,
                 MPI_Datatype scalarDt,
                 MPI_Comm comm,
                 std::ostream* errStrm = NULL)
{
  using std::endl;
  int errCode = MPI_SUCCESS;

  // mfh 17 Jan 2017: Some (generally older) MPI implementations want
  // the input buffer argument to be a pointer to nonconst, even
  // though MPI_Unpack does not modify that argument.

  // Combine the two indices into a single MPI_Pack call.
  OrdinalType inds[2] = {static_cast<OrdinalType> (0),
                         static_cast<OrdinalType> (0)};
  errCode = MPI_Unpack (const_cast<char*> (inBuf), inBufSize, &inBufCurPos,
                        inds, 2, ordinalDt, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "unpackTripleMpi: MPI_Unpack failed when unpacking indices: "
        "inBufSize=" << inBufSize << ", inBufCurPos=" << inBufCurPos << endl;
    }
    return errCode;
  }
  gblRowInd = inds[0];
  gblColInd = inds[1];

  // mfh 17 Jan 2017: Some (generally older) MPI implementations want
  // the input buffer argument to be a pointer to nonconst, even
  // though MPI_Unpack does not modify that argument.
  errCode = MPI_Unpack (const_cast<char*> (inBuf), inBufSize, &inBufCurPos,
                        &val, 1, scalarDt, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "unpackTripleMpi: MPI_Unpack failed when unpacking value: "
        "inBufSize=" << inBufSize << ", inBufCurPos=" << inBufCurPos << endl;
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
                       std::ostream* errStrm = NULL);

template<class ScalarType, class OrdinalType>
int
unpackTriplesMpi (const char inBuf[],
                  const int inBufSize,
                  int& inBufCurPos,
                  OrdinalType gblRowInds[],
                  OrdinalType gblColInds[],
                  MPI_Datatype ordinalDt,
                  ScalarType vals[],
                  MPI_Datatype scalarDt,
                  const int numEnt, // input arg, from unpackTriplesCountMpi
                  MPI_Comm comm,
                  std::ostream* errStrm = NULL)
{
  using std::endl;
  int errCode = MPI_SUCCESS;

  for (int k = 0; k < numEnt; ++k) {
    errCode = unpackTripleMpi (inBuf, inBufSize, inBufCurPos,
                               gblRowInds[k], gblColInds[k], ordinalDt,
                               vals[k], scalarDt, comm, errStrm);
    if (errCode != MPI_SUCCESS) {
      if (errStrm != NULL) {
        *errStrm << "unpackTriplesMpi: packTripleMpi failed at entry k=" << k
                 << ": inBufSize=" << inBufSize
                 << ", inBufCurPos=" << inBufCurPos
                 << "." << endl;
      }
      return errCode;
    }
  }
  return errCode;
}

#endif // HAVE_TPETRACORE_MPI

//
// SKIP DOWN TO HERE FOR "PUBLIC" INTERFACE
//

/// \brief Compute the buffer size required by packTriples for packing
///   the number of matrix entries ("triples").
///
/// countPackTriples tells me an upper bound on how much buffer space
/// I need to hold numEnt triples.  packTriplesCount actually packs
/// numEnt, the number of triples.  countPackTriplesCount tells me an
/// upper bound on how much buffer space I need to hold the number of
/// triples, not the triples themselves.
///
/// \param comm [in] Communicator used in sending and receiving the
///   packed entries.  (MPI wants this, so we have to include it.).
/// \param size [out] Pack buffer size in bytes (sizeof(char)).
/// \param errStrm [out] If nonnull, print any error messages to this
///   stream, else don't print error messages.
///
/// \return Error code.  MPI_SUCCESS (0) if successful, else nonzero.
///
/// \warning It only makes sense to call this function if using MPI.
///   If <i>not</i> building with MPI, this function is a stub that
///   returns nonzero.
int
countPackTriplesCount (const ::Teuchos::Comm<int>& comm,
                       int& size,
                       std::ostream* errStrm = NULL);

/// \brief Compute the buffer size required by packTriples for packing
///   \c numEnt number of (i,j,A(i,j)) matrix entries ("triples").
///
/// This function is NOT the same thing as packTriplesCount.
/// countPackTriples tells me an upper bound on how much buffer space
/// I need to hold numEnt triples.  packTriplesCount actually packs
/// numEnt, the number of triples.  countPackTriplesCount tells me an
/// upper bound on how much buffer space I need to hold the number of
/// triples, not the triples themselves.
///
/// \tparam ScalarType Type of each matrix entry A(i,j).
/// \tparam OrdinalType Type of each matrix index i or j.
///
/// \param numEnt [in] Number of matrix entries ("triples") to pack.
/// \param comm [in] Communicator used in sending and receiving the
///   packed entries.  (MPI wants this, so we have to include it.).
/// \param size [out] Pack buffer size in bytes (sizeof(char)).
/// \param errStrm [out] If nonnull, print any error messages to this
///   stream, else don't print error messages.
///
/// \return Error code.  MPI_SUCCESS (0) if successful, else nonzero.
///
/// \warning It only makes sense to call this function if using MPI.
///   If <i>not</i> building with MPI, this function is a stub that
///   returns nonzero.
template<class ScalarType, class OrdinalType>
int
countPackTriples (const int numEnt,
                  const ::Teuchos::Comm<int>& comm,
                  int& size, // output argument
                  std::ostream* errStrm = NULL)
{
#ifdef HAVE_TPETRACORE_MPI
  using ::Tpetra::Details::extractMpiCommFromTeuchos;
  using ::Tpetra::Details::MpiTypeTraits;

  static_assert (MpiTypeTraits<ScalarType>::isSpecialized, "countPackTriples: "
                 "ScalarType lacks an MpiTypeTraits specialization.");
  static_assert (MpiTypeTraits<OrdinalType>::isSpecialized, "countPackTriples: "
                 "OrdinalType lacks an MpiTypeTraits specialization.");

  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  MPI_Datatype ordinalDt = MpiTypeTraits<OrdinalType>::getType ();
  MPI_Datatype scalarDt = MpiTypeTraits<ScalarType>::getType ();

  const int errCode =
    countPackTriplesMpi<ScalarType, OrdinalType> (ordinalDt, scalarDt,
                                                  numEnt, mpiComm,
                                                  size, errStrm);
  if (MpiTypeTraits<ScalarType>::needsFree) {
    (void) MPI_Type_free (&scalarDt);
  }
  if (MpiTypeTraits<OrdinalType>::needsFree) {
    (void) MPI_Type_free (&ordinalDt);
  }
  return errCode;

#else // NOT HAVE_TPETRACORE_MPI
  if (errStrm != NULL) {
    *errStrm << "countPackTriples: Not implemented (no need; there's no need "
      "to pack or unpack anything if there's only one process)." << std::endl;
  }
  return -1;
#endif // HAVE_TPETRACORE_MPI
}

/// \brief Pack the count (number) of matrix triples.
///
/// This function is NOT the same thing as countPackTriples.
/// countPackTriples tells me an upper bound on how much buffer space
/// I need to hold numEnt triples.  packTriplesCount actually packs
/// numEnt, the number of triples.  countPackTriplesCount tells me an
/// upper bound on how much buffer space I need to hold the number of
/// triples, not the triples themselves.
///
/// \param numEnt [in] Number of matrix entries ("triples") to pack.
/// \param outBuf [out] Output buffer.
/// \param outBufSize [out] Total output buffer size in bytes.
/// \param outBufCurPos [in/out] Current position from which to start
///   writing to the output buffer.  This corresponds to the
///   'position' in/out argument of MPI_Pack.
/// \param comm [in] Communicator used in sending and receiving the
///   packed entries.  (MPI wants this, so we have to include it.).
/// \param errStrm [out] If nonnull, print any error messages to this
///   stream, else don't print error messages.
///
/// \return Error code.  MPI_SUCCESS (0) if successful, else nonzero.
///
/// \warning It only makes sense to call this function if using MPI.
///   If <i>not</i> building with MPI, this function is a stub that
///   returns nonzero.
int
packTriplesCount (const int numEnt,
                  char outBuf[],
                  const int outBufSize,
                  int& outBufCurPos,
                  const ::Teuchos::Comm<int>& comm,
                  std::ostream* errStrm = NULL);

/// \brief Pack matrix entries ("triples" (i, j, A(i,j))) into the
///   given output buffer.
///
/// \tparam ScalarType Type of each matrix entry A(i,j).
/// \tparam OrdinalType Type of each matrix index i or j.
///
/// \param gblRowInds [in] Row indices to pack.
/// \param gblColInds [in] Column indices to pack.
/// \param val [in] Matrix values A(i,j) to pack.
/// \param numEnt [in] Number of matrix entries ("triples") to pack.
/// \param outBuf [out] Output buffer.
/// \param outBufSize [out] Total output buffer size in bytes.
/// \param outBufCurPos [in/out] Current position from which to start
///   writing to the output buffer.  This corresponds to the
///   'position' in/out argument of MPI_Pack.
/// \param comm [in] Communicator used in sending and receiving the
///   packed entries.  (MPI wants this, so we have to include it.).
/// \param errStrm [out] If nonnull, print any error messages to this
///   stream, else don't print error messages.
///
/// \return Error code.  MPI_SUCCESS (0) if successful, else nonzero.
///
/// \warning It only makes sense to call this function if using MPI.
///   If <i>not</i> building with MPI, this function is a stub that
///   returns nonzero.
template<class ScalarType, class OrdinalType>
int
#ifdef HAVE_TPETRACORE_MPI
packTriples (const OrdinalType gblRowInds[],
             const OrdinalType gblColInds[],
             const ScalarType vals[],
             const int numEnt,
             char outBuf[],
             const int outBufSize,
             int& outBufCurPos,
             const ::Teuchos::Comm<int>& comm,
             std::ostream* errStrm = NULL)
#else // NOT HAVE_TPETRACORE_MPI
packTriples (const OrdinalType /* gblRowInds */ [],
             const OrdinalType /* gblColInds */ [],
             const ScalarType /* vals */ [],
             const int /* numEnt */,
             char /* outBuf */ [],
             const int /* outBufSize */,
             int& /* outBufCurPos */,
             const ::Teuchos::Comm<int>& /* comm */,
             std::ostream* errStrm = NULL)
#endif // HAVE_TPETRACORE_MPI
{
#ifdef HAVE_TPETRACORE_MPI
  using ::Tpetra::Details::extractMpiCommFromTeuchos;
  using ::Tpetra::Details::MpiTypeTraits;

  static_assert (MpiTypeTraits<ScalarType>::isSpecialized, "packTriples: "
                 "ScalarType lacks an MpiTypeTraits specialization.");
  static_assert (MpiTypeTraits<OrdinalType>::isSpecialized, "packTriples: "
                 "OrdinalType lacks an MpiTypeTraits specialization.");

  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  MPI_Datatype ordinalDt = MpiTypeTraits<OrdinalType>::getType ();
  MPI_Datatype scalarDt = MpiTypeTraits<ScalarType>::getType ();

  const int errCode =
    packTriplesMpi (gblRowInds, gblColInds, ordinalDt, vals, scalarDt,
                    numEnt, outBuf, outBufSize, outBufCurPos, mpiComm, errStrm);
  if (MpiTypeTraits<ScalarType>::needsFree) {
    (void) MPI_Type_free (&scalarDt);
  }
  if (MpiTypeTraits<OrdinalType>::needsFree) {
    (void) MPI_Type_free (&ordinalDt);
  }
  return errCode;

#else // NOT HAVE_TPETRACORE_MPI
  if (errStrm != NULL) {
    *errStrm << "packTriples: Not implemented (no need; there's no need to "
      "pack or unpack anything if there's only one process)." << std::endl;
  }
  return -1;
#endif // HAVE_TPETRACORE_MPI
}

/// \brief Unpack just the count of triples from the given input
///   buffer.
///
/// We store the count of triples as an \c int, because MPI buffer
/// sizes are \c int.
///
/// \param inBuf [in] Input buffer.
/// \param inBufSize [out] Total input buffer size in bytes.
/// \param inBufCurPos [in/out] Current position from which to start
///   reading from the input buffer.  This corresponds to the
///   'position' in/out argument of MPI_Unpack.
/// \param numEnt [out] Number of matrix entries ("triples") that were
///   packed.
/// \param comm [in] Communicator used in sending and receiving the
///   packed entries.  (MPI wants this, so we have to include it.).
/// \param errStrm [out] If nonnull, print any error messages to this
///   stream, else don't print error messages.
///
/// \return Error code.  MPI_SUCCESS (0) if successful, else nonzero.
///
/// \warning It only makes sense to call this function if using MPI.
///   If <i>not</i> building with MPI, this function is a stub that
///   returns nonzero.
int
unpackTriplesCount (const char inBuf[],
                    const int inBufSize,
                    int& inBufCurPos,
                    int& numEnt, // output argument!
                    const ::Teuchos::Comm<int>& comm,
                    std::ostream* errStrm = NULL);

/// \brief Unpack matrix entries ("triples" (i, j, A(i,j))) from the
///   given input buffer.
///
/// \tparam ScalarType Type of each matrix entry A(i,j).
/// \tparam OrdinalType Type of each matrix index i or j.
///
/// \param inBuf [in] Input buffer.
/// \param inBufSize [out] Total pack buffer size in bytes (sizeof(char)).
/// \param inBufCurPos [in/out] Current position from which to start
///   reading from the input buffer.  This corresponds to the
///   'position' in/out argument of MPI_Unpack.
/// \param gblRowInds [out] Row indices unpacked.
/// \param gblColInds [out] Column indices unpacked.
/// \param val [out] Matrix values A(i,j) unpacked.
/// \param numEnt [in] Number of matrix entries ("triples") to unpack.
///   If you don't know it, then you should have senders pack the
///   triples count as the first thing in the buffer, and unpack it
///   first via unpackTriplesCount().
/// \param comm [in] Communicator used in sending and receiving the
///   packed entries.  (MPI wants this, so we have to include it.).
/// \param errStrm [out] If nonnull, print any error messages to this
///   stream, else don't print error messages.
///
/// \return Error code.  MPI_SUCCESS (0) if successful, else nonzero.
///
/// \warning It only makes sense to call this function if using MPI.
///   If <i>not</i> building with MPI, this function is a stub that
///   returns nonzero.
template<class ScalarType, class OrdinalType>
int
#ifdef HAVE_TPETRACORE_MPI
unpackTriples (const char inBuf[],
               const int inBufSize,
               int& inBufCurPos,
               OrdinalType gblRowInds[],
               OrdinalType gblColInds[],
               ScalarType vals[],
               const int numEnt,
               const ::Teuchos::Comm<int>& comm,
               std::ostream* errStrm = NULL)
#else // NOT HAVE_TPETRACORE_MPI
unpackTriples (const char /* inBuf */ [],
               const int /* inBufSize */,
               int& /* inBufCurPos */,
               OrdinalType /* gblRowInds */ [],
               OrdinalType /* gblColInds */ [],
               ScalarType /* vals */ [],
               const int /* numEnt */,
               const ::Teuchos::Comm<int>& /* comm */,
               std::ostream* errStrm = NULL)
#endif // HAVE_TPETRACORE_MPI
{
#ifdef HAVE_TPETRACORE_MPI
  using ::Tpetra::Details::extractMpiCommFromTeuchos;
  using ::Tpetra::Details::MpiTypeTraits;

  static_assert (MpiTypeTraits<ScalarType>::isSpecialized, "unpackTriples: "
                 "ScalarType lacks an MpiTypeTraits specialization.");
  static_assert (MpiTypeTraits<OrdinalType>::isSpecialized, "unpackTriples: "
                 "OrdinalType lacks an MpiTypeTraits specialization.");

  MPI_Comm mpiComm = extractMpiCommFromTeuchos (comm);
  MPI_Datatype ordinalDt = MpiTypeTraits<OrdinalType>::getType ();
  MPI_Datatype scalarDt = MpiTypeTraits<ScalarType>::getType ();

  const int errCode =
    unpackTriplesMpi (inBuf, inBufSize, inBufCurPos,
                      gblRowInds, gblColInds, ordinalDt,
                      vals, scalarDt, numEnt, mpiComm, errStrm);
  if (MpiTypeTraits<ScalarType>::needsFree) {
    (void) MPI_Type_free (&scalarDt);
  }
  if (MpiTypeTraits<OrdinalType>::needsFree) {
    (void) MPI_Type_free (&ordinalDt);
  }
  return errCode;

#else // NOT HAVE_TPETRACORE_MPI
  if (errStrm != NULL) {
    *errStrm << "unpackTriples: Not implemented (no need; there's no need to "
      "pack or unpack anything if there's only one process)." << std::endl;
  }
  return -1;
#endif // HAVE_TPETRACORE_MPI
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKTRIPLES_HPP

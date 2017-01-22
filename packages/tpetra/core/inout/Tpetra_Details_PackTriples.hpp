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

#ifndef TPETRA_DETAILS_PACKTRIPLES_HPP
#define TPETRA_DETAILS_PACKTRIPLES_HPP

#include "TpetraCore_config.h"
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

  size = 0;
  {
    int curSize = 0;
    errCode = MPI_Pack_size (1, ordinalDt, comm, &curSize);
    if (errCode != MPI_SUCCESS) {
      if (errStrm != NULL) {
        *errStrm << "countPackTripleMpi: MPI_Pack_size failed on "
          "ordinalDt call" << endl;
      }
      return errCode;
    }
    size += 2 * curSize; // two Ordinals
  }
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
    size += curSize; // one Scalar
  }
  size *= numEnt;
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

  // mfh 17 Jan 2017: Some (generally older) MPI implementations want
  // this to be a void* instead of a const void*.
  errCode = MPI_Pack (const_cast<OrdinalType*> (&gblRowInd), 1, ordinalDt,
                      outBuf, outBufSize, &outBufCurPos, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "packTripleMpi: MPI_Pack failed on gblRowInd = " << gblRowInd
               << endl;
    }
    return errCode;
  }
  // mfh 17 Jan 2017: Some (generally older) MPI implementations want
  // this to be a void* instead of a const void*.
  errCode = MPI_Pack (const_cast<OrdinalType*> (&gblColInd), 1, ordinalDt,
                      outBuf, outBufSize, &outBufCurPos, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "packTripleMpi: MPI_Pack failed on gblColInd = " << gblColInd
               << endl;
    }
    return errCode;
  }
  // mfh 17 Jan 2017: Some (generally older) MPI implementations want
  // this to be a void* instead of a const void*.
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
  // this to be a void* instead of a const void*.
  errCode = MPI_Unpack (const_cast<char*> (inBuf), inBufSize, &inBufCurPos,
                        &gblRowInd, 1, ordinalDt, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "unpackTripleMpi: MPI_Unpack failed on gblRowInd" << endl;
    }
    return errCode;
  }
  // mfh 17 Jan 2017: Some (generally older) MPI implementations want
  // this to be a void* instead of a const void*.
  errCode = MPI_Unpack (const_cast<char*> (inBuf), inBufSize, &inBufCurPos,
                        &gblColInd, 1, ordinalDt, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "unpackTripleMpi: MPI_Unpack failed on gblColInd" << endl;
    }
    return errCode;
  }
  // mfh 17 Jan 2017: Some (generally older) MPI implementations want
  // this to be a void* instead of a const void*.
  errCode = MPI_Unpack (const_cast<char*> (inBuf), inBufSize, &inBufCurPos,
                        &val, 1, scalarDt, comm);
  if (errCode != MPI_SUCCESS) {
    if (errStrm != NULL) {
      *errStrm << "unpackTripleMpi: MPI_Unpack failed on val" << endl;
    }
    return errCode;
  }
  return errCode;
}

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
                  const int numEnt, // we know this from DistObject
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
        *errStrm << "packTriplesMpi: packTripleMpi failed at entry k=" << k
                 << endl;
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
             const Teuchos::Comm<int>& comm,
             std::ostream* errStrm = NULL)
#else // NOT HAVE_TPETRACORE_MPI
packTriples (const OrdinalType /* gblRowInds */ [],
             const OrdinalType /* gblColInds */ [],
             const ScalarType /* vals */ [],
             const int /* numEnt */,
             char /* outBuf */ [],
             const int /* outBufSize */,
             int& /* outBufCurPos */,
             const Teuchos::Comm<int>& /* comm */,
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

template<class ScalarType, class OrdinalType>
int
#ifdef HAVE_TPETRACORE_MPI
unpackTriples (const char inBuf[],
               const int inBufSize,
               int& inBufCurPos,
               OrdinalType gblRowInds[],
               OrdinalType gblColInds[],
               ScalarType vals[],
               const int numEnt, // we know this from DistObject
               const Teuchos::Comm<int>& comm,
               std::ostream* errStrm = NULL)
#else // NOT HAVE_TPETRACORE_MPI
unpackTriples (const char /* inBuf */ [],
               const int /* inBufSize */,
               int& /* inBufCurPos */,
               OrdinalType /* gblRowInds */ [],
               OrdinalType /* gblColInds */ [],
               ScalarType /* vals */ [],
               const int /* numEnt */,
               const Teuchos::Comm<int>& /* comm */,
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

template<class ScalarType, class OrdinalType>
int
countPackTriples (const int numEnt,
                  const Teuchos::Comm<int>& comm,
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

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKTRIPLES_HPP

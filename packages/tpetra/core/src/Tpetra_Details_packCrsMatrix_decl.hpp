// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PACKCRSMATRIX_DECL_HPP
#define TPETRA_DETAILS_PACKCRSMATRIX_DECL_HPP

#include "TpetraCore_config.h"
#include "Kokkos_DualView.hpp"
#include "Tpetra_DistObject_decl.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"

/// \file Tpetra_Details_packCrsMatrix_decl.hpp
/// \brief Functions for packing the entries of a Tpetra::CrsMatrix
///   for communication, in the case where it is valid to go to the
///   KokkosSparse::CrsMatrix (local sparse matrix data structure)
///   directly.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.
///
/// Data (bytes) describing the row of the CRS matrix are "packed"
/// (concatenated) in to a (view of) char* object in the following order:
///
///   1. number of entries (LocalOrdinal)
///   2. global column indices (GlobalOrdinal)
///   3. proces IDs (optional, int)
///   4. row values (Scalar)
///
/// The functions in this file are companions to
/// Tpetra_Details_unpackCrsMatrix.hpp, i.e., Tpetra_Details_unpackCrsMatrix.hpp
/// implements the reverse of the packing order described above to ensure proper
/// unpacking.

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
// Forward declaration of Array
template<class T> class Array;
// Forward declaration of ArrayView
template<class T> class ArrayView;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

/// \brief Pack specified entries of the given local sparse matrix for
///   communication.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param exportLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// This is the public interface to the pack machinery
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrix (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
               Teuchos::Array<char>& exports,
               const Teuchos::ArrayView<size_t>& numPacketsPerLID,
               const Teuchos::ArrayView<const LO>& exportLIDs,
               size_t& constantNumPackets);

/// \brief Pack specified entries of the given local sparse matrix for
///   communication, for "new" DistObject interface.
///
/// \tparam ST The type of the entries of the matrix.  This must be
///   the same as the Scalar template parameter of Tpetra::CrsMatrix.
/// \tparam LO The type of local indices.  This must be the same as
///   the LocalOrdinal template parameter of Tpetra::CrsMatrix.
/// \tparam GO The type of global indices.  This must be the same as
///   the GlobalOrdinal template parameter of Tpetra::CrsMatrix.
/// \tparam NT The Node type.  This must be the same as the Node
///   template parameter of Tpetra::CrsMatrix.
///
/// \param sourceMatrix [in] The "source" matrix to pack.
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] On output,
///   numPacketsPerLID.d_view[k] is the number of bytes packed for row
///   exportLIDs.d_view[k] of the local matrix.
///
/// \param exportLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Same as the constantNumPackets
///   output argument of Tpetra::DistObject::packAndPrepare (which
///   see).
///
/// This method implements CrsMatrix::packNew, and thus
/// CrsMatrix::packAndPrepare, for the case where the matrix to
/// pack has a valid KokkosSparse::CrsMatrix.
template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrixNew (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                  Kokkos::DualView<char*,
                    typename DistObject<char, LO, GO, NT>::buffer_device_type>& exports,
                  const Kokkos::DualView<size_t*,
                    typename DistObject<char, LO, GO, NT>::buffer_device_type>& numPacketsPerLID,
                  const Kokkos::DualView<const LO*,
                    typename DistObject<char, LO, GO, NT>::buffer_device_type>& exportLIDs,
                  size_t& constantNumPackets);

/// \brief Pack specified entries of the given local sparse matrix for
///   communication.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param exportLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// This is the public interface to the pack machinery
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrixWithOwningPIDs (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                             Kokkos::DualView<char*, typename DistObject<char, LO, GO, NT>::buffer_device_type>& exports_dv,
                             const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                             const Teuchos::ArrayView<const LO>& exportLIDs,
                             const Teuchos::ArrayView<const int>& sourcePIDs,
                             size_t& constantNumPackets);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKCRSMATRIX_DECL_HPP

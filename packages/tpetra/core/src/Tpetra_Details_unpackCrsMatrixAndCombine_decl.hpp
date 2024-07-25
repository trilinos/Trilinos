// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DECL_HPP
#define TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DECL_HPP

#include "TpetraCore_config.h"
#include "Tpetra_CombineMode.hpp"
#include "Kokkos_DualView.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_DistObject_decl.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_Details_unpackCrsMatrixAndCombine_decl.hpp
/// \brief Declaration of functions for unpacking the entries of a
///   Tpetra::CrsMatrix for communication, in the case where it is
///   valid to go to the KokkosSparse::CrsMatrix (local sparse matrix
///   data structure) directly.
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
/// Tpetra_Details_packCrsMatrix.hpp, i.e., Tpetra_Details_packCrsMatrix.hpp
/// implements the packing order described above to ensure proper unpacking.

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

/// \brief Unpack the imported column indices and values, and combine
///   into matrix.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param imports [in] Input pack buffer
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param importLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// \param combineMode [in] the mode to use for combining values
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
///
/// This is the public interface to the unpack and combine machinery and
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename ST, typename LO, typename GO, typename NT>
void
unpackCrsMatrixAndCombine (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                           const Teuchos::ArrayView<const char>& imports,
                           const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
                           const Teuchos::ArrayView<const LO>& importLIDs,
                           size_t constantNumPackets,
                           CombineMode combineMode);

template<typename ST, typename LO, typename GO, typename NT>
void
unpackCrsMatrixAndCombineNew(
  const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
  Kokkos::DualView<char*,
    typename DistObject<char, LO, GO, NT>::buffer_device_type> imports,
  Kokkos::DualView<size_t*,
    typename DistObject<char, LO, GO, NT>::buffer_device_type> numPacketsPerLID,
  const Kokkos::DualView<const LO*,
    typename DistObject<char, LO, GO, NT>::buffer_device_type>& importLIDs,
  const size_t constantNumPackets,
  const CombineMode combineMode);

/// \brief Special version of Tpetra::Details::unpackCrsMatrixAndCombine
///   that also unpacks owning process ranks.
///
/// Perform the count for unpacking the imported column indices pids,
/// and values, and combining them into matrix.  Return (a ceiling on)
/// the number of local stored entries ("nonzeros") in the matrix.  If
/// there are no shared rows in the sourceMatrix this count is exact.
///
/// Note: This routine also counts the copyAndPermute nonzeros in
/// addition to those that come in via import.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam Node The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param imports [in] Input pack buffer
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param importLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// \param combineMode [in] the mode to use for combining values
///
/// \param numSameIds [in]
///
/// \param permuteToLIDs [in]
///
/// \param permuteFromLIDs [in]
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
//
/// \warning This method is intended for expert developer use
///   only, and should never be called by user code.
///
/// Note: This is the public interface to the unpack and combine machinery and
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
size_t
unpackAndCombineWithOwningPIDsCount (
    const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & sourceMatrix,
    const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
    const Teuchos::ArrayView<const char> &imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    size_t constantNumPackets,
    CombineMode combineMode,
    size_t numSameIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs);

/// \brief unpackAndCombineIntoCrsArrays
///
/// Note: The SourcePids vector (on input) should contain owning PIDs
/// for each column in the (source) ColMap, as from
/// Tpetra::Import_Util::getPids, with the "-1 for local" option being
/// used.
///
/// Note: The TargetPids vector (on output) will contain owning PIDs
/// for each entry in the matrix, with the "-1 for local" for locally
/// owned entries.
///
/// Note: This method does the work previously done in unpackAndCombineWithOwningPIDsCount,
/// namely, calculating the local number of nonzeros, and allocates CRS
/// arrays of the correct sizes.

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
unpackAndCombineIntoCrsArrays (
    const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & sourceMatrix,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void>,
    const Kokkos::View<const char*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>
          ,void, void >,
    const Kokkos::View<const size_t*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>
          ,void, void >,
    const size_t numSameIDs,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void>,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void>,
    size_t TargetNumRows,
    const int MyTargetPID,
    Teuchos::ArrayRCP<size_t>& CRS_rowptr,
    Teuchos::ArrayRCP<GlobalOrdinal>& CRS_colind,
    Teuchos::ArrayRCP<Scalar>& CRS_vals,
    const Teuchos::ArrayView<const int>& SourcePids,
    Teuchos::Array<int>& TargetPids);

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
unpackAndCombineIntoCrsArrays (
    const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & sourceMatrix,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void>,
    const Kokkos::View<const char*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>
          ,void, void >,
    const Kokkos::View<const size_t*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>
          ,void, void >,
    const size_t numSameIDs,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void>,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void>,
    size_t TargetNumRows,
    const int MyTargetPID,
    Kokkos::View<size_t*,typename Node::device_type>& /*crs_rowptr_d*/,
    Kokkos::View<GlobalOrdinal*,typename Node::device_type>&     /*crs_colind_d*/,
    Kokkos::View<typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type*,typename Node::device_type>& /*crs_vals_d*/,
    const Teuchos::ArrayView<const int>& SourcePids,
    Kokkos::View<int*,typename Node::device_type>& /*TargetPids*/);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DECL_HPP

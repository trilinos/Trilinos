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

#ifndef TPETRA_DETAILS_UNPACKCRSGRAPHANDCOMBINE_DECL_HPP
#define TPETRA_DETAILS_UNPACKCRSGRAPHANDCOMBINE_DECL_HPP

#include "TpetraCore_config.h"
#include "Tpetra_CombineMode.hpp"
#include "Kokkos_DualView.hpp"
#include "Tpetra_DistObject_decl.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"

/// \file Tpetra_Details_unpackCrsGraphAndCombine_decl.hpp
/// \brief Declaration of functions for unpacking the entries of a
///   Tpetra::CrsGraph for communication, in the case where it is
///   valid to go to the KokkosSparse::CrsGraph (local sparse graph
///   data structure) directly.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.
///
/// Data (bytes) describing the row of the CRS graph are "packed"
/// (concatenated) in to a (view of) packet_type* object in the following order:
///
///   1. global column indices (GO)
///   2. process IDs (optional, int)
///
/// The functions in this file are companions to
/// Tpetra_Details_packCrsGraph.hpp, i.e., Tpetra_Details_packCrsGraph.hpp
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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Forward declaration of Distributor
class Distributor;
#endif // DOXYGEN_SHOULD_SKIP_THIS

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

/// \brief Unpack the imported column indices and combine
///   into graph.
///
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceGraph [in] the CrsGraph source
///
/// \param imports [in] Input pack buffer
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local graph.
///
/// \param importLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// \param combineMode [in] the mode to use for combining
///
/// \param atomic [in] whether or not do atomic adds/replaces in to the graph
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
///
/// This is the public interface to the unpack and combine machinery and
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsGraph migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<class LO, class GO, class NT>
void
unpackCrsGraphAndCombine(
    const CrsGraph<LO, GO, NT>& sourceGraph,
    const Teuchos::ArrayView<const typename CrsGraph<LO,GO,NT>::packet_type>& imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    const Teuchos::ArrayView<const LO>& importLIDs,
    size_t constantNumPackets,
    Distributor & distor,
    CombineMode combineMode,
    const bool atomic);

template<class LO, class GO, class NT>
void
unpackCrsGraphAndCombineNew(
    const CrsGraph<LO, GO, NT>& sourceGraph,
    const Kokkos::DualView<const typename CrsGraph<LO,GO,NT>::packet_type*,
                           typename CrsGraph<LO,GO,NT>::buffer_device_type>&
                           imports,
    const Kokkos::DualView<const size_t*,
                           typename CrsGraph<LO,GO,NT>::buffer_device_type>&
                           numPacketsPerLID,
    const Kokkos::DualView<const LO*, typename NT::device_type>& importLIDs,
    const size_t constantNumPackets,
    Distributor & distor,
    const CombineMode combineMode,
    const bool atomic);

/// \brief Special version of Tpetra::Details::unpackCrsGraphAndCombine
///   that also unpacks owning process ranks.
///
/// Perform the count for unpacking the imported column indices and pids,
/// and combining them into graph.  Return (a ceiling on)
/// the number of local stored entries ("nonzeros") in the graph.  If
/// there are no shared rows in the sourceGraph this count is exact.
///
/// Note: This routine also counts the copyAndPermute nonzeros in
/// addition to those that come in via import.
///
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceGraph [in] the CrsGraph source
///
/// \param imports [in] Input pack buffer
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local graph.
///
/// \param importLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// \param combineMode [in] the mode to use for combining
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
/// CrsGraph migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<class LO, class GO, class NT>
size_t
unpackAndCombineWithOwningPIDsCount(
    const CrsGraph<LO, GO, NT> & sourceGraph,
    const Teuchos::ArrayView<const LO> &importLIDs,
    const Teuchos::ArrayView<const typename CrsGraph<LO,GO,NT>::packet_type> &imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    size_t constantNumPackets,
    Distributor &distor,
    CombineMode combineMode,
    size_t numSameIDs,
    const Teuchos::ArrayView<const LO>& permuteToLIDs,
    const Teuchos::ArrayView<const LO>& permuteFromLIDs);

/// \brief unpackAndCombineIntoCrsArrays
///
/// \note You should call unpackAndCombineWithOwningPIDsCount first
///   and allocate all arrays accordingly, before calling this
///   function.
///
/// Note: The SourcePids vector (on input) should contain owning PIDs
/// for each column in the (source) ColMap, as from
/// Tpetra::Import_Util::getPids, with the "-1 for local" option being
/// used.
///
/// Note: The TargetPids vector (on output) will contain owning PIDs
/// for each entry in the graph, with the "-1 for local" for locally
/// owned entries.
template<class LO, class GO, class NT>
void
unpackAndCombineIntoCrsArrays(
    const CrsGraph<LO, GO, NT> & sourceGraph,
    const Teuchos::ArrayView<const LO>& importLIDs,
    const Teuchos::ArrayView<const typename CrsGraph<LO,GO,NT>::packet_type>& imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    const size_t constantNumPackets,
    Distributor& distor,
    const CombineMode combineMode,
    const size_t numSameIDs,
    const Teuchos::ArrayView<const LO>& permuteToLIDs,
    const Teuchos::ArrayView<const LO>& permuteFromLIDs,
    size_t TargetNumRows,
    size_t TargetNumNonzeros,
    const int MyTargetPID,
    const Teuchos::ArrayView<size_t>& CRS_rowptr,
    const Teuchos::ArrayView<GO>& CRS_colind,
    const Teuchos::ArrayView<const int>& SourcePids,
    Teuchos::Array<int>& TargetPids);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_UNPACKCRSGRAPHANDCOMBINE_DECL_HPP

// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PACKCRSGRAPH_DECL_HPP
#define TPETRA_DETAILS_PACKCRSGRAPH_DECL_HPP

#include "TpetraCore_config.h"
#include "Kokkos_DualView.hpp"
#include "Tpetra_DistObject_decl.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"

/// \file Tpetra_Details_packCrsGraph_decl.hpp
/// \brief Functions for packing the entries of a Tpetra::CrsGraph
///   for communication, in the case where it is valid to go to the
///   KokkosSparse::CrsGraph (local sparse graph data structure)
///   directly.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.
///
/// Data (bytes) describing the row of the CRS graph are "packed"
/// (concatenated) in to a (view of) packet_type* object in the following order:
///
///   1. number of entries (LocalOrdinal)
///   2. global column indices (GlobalOrdinal)
///   3. proces IDs (optional, int)
///
/// The functions in this file are companions to
/// Tpetra_Details_unpackCrsGraph.hpp, i.e., Tpetra_Details_unpackCrsGraph.hpp
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

/// \brief Pack specified entries of the given local sparse graph for
///   communication.
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
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local graph.
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
/// CrsGraph migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename LO, typename GO, typename NT>
void
packCrsGraph (const CrsGraph<LO, GO, NT>& sourceGraph,
               Teuchos::Array<typename CrsGraph<LO,GO,NT>::packet_type>& exports,
               const Teuchos::ArrayView<size_t>& numPacketsPerLID,
               const Teuchos::ArrayView<const LO>& exportLIDs,
               size_t& constantNumPackets);

/// \brief Pack specified entries of the given local sparse graph for
///   communication, for "new" DistObject interface.
///
/// \tparam LO The type of local indices.  This must be the same as
///   the LocalOrdinal template parameter of Tpetra::CrsGraph.
/// \tparam GO The type of global indices.  This must be the same as
///   the GlobalOrdinal template parameter of Tpetra::CrsGraph.
/// \tparam NT The Node type.  This must be the same as the Node
///   template parameter of Tpetra::CrsGraph.
///
/// \param sourceGraph [in] The "source" graph to pack.
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] On output,
///   numPacketsPerLID.d_view[k] is the number of bytes packed for row
///   exportLIDs.d_view[k] of the local graph.
///
/// \param exportLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Same as the constantNumPackets
///   output argument of Tpetra::DistObject::packAndPrepare (which
///   see).
///
/// This method implements CrsGraph::packNew, and thus
/// CrsGraph::packAndPrepare, for the case where the graph to
/// pack has a valid KokkosSparse::CrsGraph.
template<typename LO, typename GO, typename NT>
void
packCrsGraphNew (const CrsGraph<LO, GO, NT>& sourceGraph,
                 const Kokkos::DualView<
                   const LO*,
                   typename CrsGraph<LO, GO, NT>::buffer_device_type
                 >& exportLIDs,
                 const Kokkos::DualView<
                   const int*,
                   typename CrsGraph<LO, GO, NT>::buffer_device_type
                 >& exportPIDs,
                 Kokkos::DualView<
                   typename CrsGraph<LO, GO, NT>::packet_type*,
                   typename CrsGraph<LO, GO, NT>::buffer_device_type
                 >& exports,
                 Kokkos::DualView<
                   size_t*,
                   typename CrsGraph<LO, GO, NT>::buffer_device_type
                 > numPacketsPerLID,
                 size_t& constantNumPackets,
                 const bool pack_pids);

/// \brief Pack specified entries of the given local sparse graph for
///   communication.
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
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local graph.
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
/// CrsGraph migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename LO, typename GO, typename NT>
void
packCrsGraphWithOwningPIDs (const CrsGraph<LO,GO,NT>& sourceGraph,
                            Kokkos::DualView<typename CrsGraph<LO,GO,NT>::packet_type*,
                                             typename CrsGraph<LO,GO,NT>::buffer_device_type>&
                                             exports_dv,
                            const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                            const Teuchos::ArrayView<const LO>& exportLIDs,
                            const Teuchos::ArrayView<const int>& sourcePIDs,
                            size_t& constantNumPackets);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKCRSGRAPH_DECL_HPP

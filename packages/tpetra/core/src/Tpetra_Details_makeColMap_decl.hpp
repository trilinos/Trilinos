#ifndef TPETRA_DETAILS_MAKECOLMAP_DECL_HPP
#define TPETRA_DETAILS_MAKECOLMAP_DECL_HPP

// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tpetra_Details_makeColMap_decl.hpp
/// \brief Declaration of Tpetra::Details::makeColMap, a function for
///   creating the column Map of a Tpetra::CrsGraph
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra.  Users may not rely on this file or its contents.
///   They may change or disappear at any time.
///
/// This file declares the Tpetra::Details::makeColMap function, which
/// creates the column Map of a Tpetra::CrsGraph.

#include "TpetraCore_config.h"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_RowGraph_fwd.hpp"
#include "Kokkos_Core.hpp"
#include <ostream>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // forward declaration of Array
  template<class T> class Array;

  // forward declaration of RCP
  template<class T> class RCP;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {

namespace Details {

/// \brief Make the graph's column Map.
///
/// \tparam LO Local ordinal type; the type of local indices in the graph.
/// \tparam GO Global ordinal type; the type of global indices in the graph.
/// \tparam NT Node type; the third template parameter of Tpetra::CrsGraph.
///
/// \param colMap [out] On output: pointer to the column Map for the
///   given graph.  This is only valid if the returned error code is
///   zero on <i>all</i> processes in the communicator of
///   <tt>domMap</tt> (see below).  This may be the same (literally
///   the same object as) as <tt>domMap</tt>, depending on the graph.
/// \param remotePIDs [out] The process ranks corresponding to the
///   column Map's "remote" (not on the calling process in the domain
///   Map) indices.
/// \param domMap [in] The domain Map to use for creating the column
///   Map.  This need not be the same as graph.getDomainMap().  It's
///   OK for the latter to be null, in fact.  <tt>domMap</tt> needs to
///   be passed in by RCP, because it's possible for the returned
///   column Map <tt>colMap</tt> (see above) to equal <tt>domMap</tt>.
/// \param graph [in] The graph for which to make a column Map.  This
///   function does NOT modify the graph's column Map, if it happens
///   to have one already.  Thus, this function supports graph
///   modification.
/// \param sortEachProcsGids [in] Whether to sort column Map GIDs
///   associated with each remote process in ascending order.  This is
///   \c true by default.  If \c false, leave GIDs in their original
///   order as discovered in the graph by iterating in ascending
///   order through the local rows of the graph.
/// \param errStrm [out] If nonnull, print error messages to this.
///
/// \return Error code; zero if and only if successful.  This value is
///   local to the calling process.  On error, the code may be zero on
///   some processes, and nonzero on other processes.  You, the user,
///   are responsible for propagating that error state to all
///   processes.
///
/// This function <i>always</i> makes a column Map, even if the graph
/// already has one.  This makes it possible to change the graph's
/// structure, and have its column Map and corresponding Import update
/// in the same way.
///
/// The sortEachProcsGids argument corresponds to
/// sortGhostsAssociatedWithEachProcessor_ in CrsGraph.  This function
/// always groups remote GIDs by process rank, so that all remote GIDs
/// with the same owning rank occur contiguously.  The
/// sortEachProcsGids argument (see above) whether this function sorts
/// remote GIDs in increasing order within those groups.  This
/// function sorts by default.  This behavior differs from Epetra,
/// which does not sort remote GIDs with the same owning process.
/// means "sort remote GIDs."  If you don't want to sort, for
/// compatibility with Epetra, set sortEachProcsGids to false.
template <class LO, class GO, class NT>
int
makeColMap (Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& colMap,
            Teuchos::Array<int>& remotePIDs,
            const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& domMap,
            const RowGraph<LO, GO, NT>& graph,
            const bool sortEachProcsGids = true,
            std::ostream* errStrm = NULL);

/// \brief Construct a column map for the given set of gids (always sorting remote GIDs within each remote process).
/// \param colMap [out] Will be set to the new column map.
/// \param domMap [in] The domain map, used to determine which global columns are locally owned.
template <class LO, class GO, class NT>
int
makeColMap (Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>& colMap,
            const Teuchos::RCP<const Tpetra::Map<LO, GO, NT>>& domMap,
            Kokkos::View<GO*, typename NT::memory_space> gids,
            std::ostream* errStrm = NULL);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_MAKECOLMAP_DECL_HPP

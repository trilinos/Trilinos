// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_GLOBALLEXICOGRPHICINDEXMANAGER_DECL_HPP
#define MUELU_GLOBALLEXICOGRPHICINDEXMANAGER_DECL_HPP

// use for Teuchos:Comm<T>
#include "Teuchos_CommHelpers.hpp"

#include <Xpetra_Map_fwd.hpp>

#include <MueLu_ConfigDefs.hpp>
#include <MueLu_IndexManager.hpp>
#include <MueLu_GlobalLexicographicIndexManager_fwd.hpp>

/*****************************************************************************

****************************************************************************/

namespace MueLu {

/*!
    @class Aggregates
    @brief Container class for aggregation information.

    @ingroup Aggregation

    Structure holding aggregate information. Right now, nAggregates, IsRoot,
    Vertex2AggId, procWinner are populated.  This allows us to look at a node
    and determine the aggregate to which it has been assigned and the id of the
    processor that owns this aggregate. It is not so easy to determine vertices
    within the kth aggregate or the size of the kth aggregate. Thus, it might be
    useful to have a secondary structure which would be a rectangular CrsGraph
    where rows (or vertices) correspond to aggregates and colunmns (or edges)
    correspond to nodes. While not strictly necessary, it might be convenient.
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class GlobalLexicographicIndexManager : public IndexManager<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_GLOBALLEXICOGRAPHICINDEXMANAGER_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  GlobalLexicographicIndexManager();

  GlobalLexicographicIndexManager(const RCP<const Teuchos::Comm<int> > comm, const bool coupled,
                                  const int NumDimensions, const int interpolationOrder,
                                  const Array<GO> GFineNodesPerDir,
                                  const Array<LO> LFineNodesPerDir, const Array<LO> CoarseRate,
                                  const GO MinGlobalIndex);

  virtual ~GlobalLexicographicIndexManager() {}

  void computeGlobalCoarseParameters();

  void getGhostedNodesData(const RCP<const Map> fineMap,
                           Array<LO>& ghostedNodeCoarseLIDs,
                           Array<int>& ghostedNodeCoarsePIDs,
                           Array<GO>& ghostedNodeCoarseGIDs) const;

  void getCoarseNodesData(const RCP<const Map> fineCoordinatesMap,
                          Array<GO>& coarseNodeCoarseGIDs,
                          Array<GO>& coarseNodeFineGIDs) const;

  std::vector<std::vector<GO> > getCoarseMeshData() const;

  void getFineNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const;

  void getFineNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const;

  void getFineNodeGhostedTuple(const LO myLID, LO& i, LO& j, LO& k) const;

  void getFineNodeGID(const GO i, const GO j, const GO k, GO& myGID) const;

  void getFineNodeLID(const LO i, const LO j, const LO k, LO& myLID) const;

  void getCoarseNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const;

  void getCoarseNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const;

  void getCoarseNodeGID(const GO i, const GO j, const GO k, GO& myGID) const;

  void getCoarseNodeLID(const LO i, const LO j, const LO k, LO& myLID) const;

  void getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const;

  void getCoarseNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const;

  void getGhostedNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const;

  void getGhostedNodeCoarseLID(const LO i, const LO j, const LO k, LO& myLID) const;

 private:
};

}  // namespace MueLu

#define MUELU_GLOBALLEXICOGRPHICINDEXMANAGER_SHORT
#endif  // MUELU_GLOBALLEXICOGRPHICINDEXMANAGER_DECL_HPP

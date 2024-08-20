// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DECL_HPP
#define MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DECL_HPP

// use for Teuchos:Comm<T>
#include "Teuchos_CommHelpers.hpp"

#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_IndexManager.hpp"
#include "MueLu_LocalLexicographicIndexManager_fwd.hpp"

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

template <class LocalOrdinal, class GlobalOrdinal, class Node>
class LocalLexicographicIndexManager : public IndexManager<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  LocalLexicographicIndexManager() = default;

  LocalLexicographicIndexManager(const RCP<const Teuchos::Comm<int> > comm, const bool coupled,
                                 const int NumDimensions, const int interpolationOrder,
                                 const int MyRank, const int NumRanks,
                                 const Array<GO> GFineNodesPerDir,
                                 const Array<LO> LFineNodesPerDir,
                                 const Array<LO> CoarseRate, const Array<GO> MeshData);

  virtual ~LocalLexicographicIndexManager() {}

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
  const int myRank;    ///< Local rank ID.
  const int numRanks;  ///< Number of ranks used to decompose the problem.

  // Iterator delimiting the entries in meshData that correspond to the block that owns the local
  // part of the mesh.
  typename std::vector<std::vector<GO> >::iterator myBlockStart, myBlockEnd;

  int pi, pj, pk;  ///< Number of processors in each diretcion.

  int numBlocks;  ///< Number of mesh block.
  int myBlock;    ///< local mesh block ID.

  int myRankIndex;                               ///< local process index for record in meshData after sorting.
  Array<int> rankIndices;                        ///< mapping between rank ID and reordered rank ID.
  std::vector<std::vector<GO> > meshData;        ///< layout of indices accross all processes.
  std::vector<std::vector<GO> > coarseMeshData;  ///< layout of indices accross all processes after coarsening.

  void sortLocalLexicographicData();

  void computeCoarseLocalLexicographicData();

  void getGIDLocalLexicographic(const LO iGhosted, const LO jGhosted, const LO kGhosted,
                                const Array<LO> coarseNodeFineIndices, GO& myGID, LO& myPID,
                                LO& myLID) const;
};

}  // namespace MueLu

#define MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_SHORT
#endif  // MUELU_LOCALLEXICOGRAPHICINDEXMANAGER_DECL_HPP

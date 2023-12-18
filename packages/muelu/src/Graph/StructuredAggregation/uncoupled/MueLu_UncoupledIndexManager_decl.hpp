// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu        (jhu@sandia.gov)
//                    Ray Tuminaro       (rstumin@sandia.gov)
//                    Luc Berger-Vergiat (lberge@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_UNCOUPLEDINDEXMANAGER_DECL_HPP
#define MUELU_UNCOUPLEDINDEXMANAGER_DECL_HPP

// use for Teuchos:Comm<T>
#include "Teuchos_CommHelpers.hpp"

#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_IndexManager.hpp"
#include "MueLu_UncoupledIndexManager_fwd.hpp"

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
class UncoupledIndexManager : public IndexManager<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_UNCOUPLEDINDEXMANAGER_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  // LBV: I doubt that it makes sense to have
  //  this particular constructor since it is
  //  not used anywhere and parameters cannot
  //  all accessible after construction.
  UncoupledIndexManager() = default;

  UncoupledIndexManager(const RCP<const Teuchos::Comm<int> > comm, const bool coupled,
                        const int NumDimensions, const int interpolationOrder,
                        const int MyRank, const int NumRanks,
                        const Array<GO> GFineNodesPerDir,
                        const Array<LO> LFineNodesPerDir,
                        const Array<LO> CoarseRate,
                        const bool singleCoarsePoint);

  virtual ~UncoupledIndexManager() {}

  void computeGlobalCoarseParameters();

  std::vector<std::vector<GO> > getCoarseMeshData() const;

  void getGhostedNodesData(const RCP<const Map> fineMap,
                           Array<LO>& ghostedNodeCoarseLIDs,
                           Array<int>& ghostedNodeCoarsePIDs,
                           Array<GO>& ghostedNodeCoarseGIDs) const;

  void getCoarseNodesData(const RCP<const Map> fineCoordinatesMap,
                          Array<GO>& coarseNodeCoarseGIDs,
                          Array<GO>& coarseNodeFineGIDs) const;

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
};

}  // namespace MueLu

#define MUELU_UNCOUPLEDINDEXMANAGER_SHORT
#endif  // MUELU_UNCOUPLEDINDEXMANAGER_DECL_HPP

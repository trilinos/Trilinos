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
#ifndef MUELU_INDEXMANAGER_DECL_HPP
#define MUELU_INDEXMANAGER_DECL_HPP

// use for Teuchos:Comm<T>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"

#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_IndexManager_fwd.hpp"

/*****************************************************************************

****************************************************************************/

namespace MueLu {

/*!
    @class IndexManager
    @brief Container class for mesh layout and indices calculation.

    @ingroup Aggregation

    Structure holding mesh parameters for structured mesh. Based on these
    parameters the IndexManager computes indices in different index spaces
    and it also provides utilites for coarsening. Finally it also provides
    an interface to different mesh layouts. Currently global lexicographic
    and local lexicographic mesh orderings are supported.
*/

template <class LocalOrdinal, class GlobalOrdinal, class Node>
class IndexManager : public BaseClass {
#undef MUELU_INDEXMANAGER_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 private:
 protected:
  const RCP<const Teuchos::Comm<int> > comm_;  ///< Communicator used by uncoupled aggregation
  const bool coupled_;                         ///< Flag for coupled vs uncoupled aggregation mode, if true aggregation is coupled.
  const bool singleCoarsePoint_;               ///< Flag telling us if can reduce dimensions to a single layer.
  const int numDimensions;                     ///< Number of spacial dimensions in the problem
  const int interpolationOrder_;               ///< Interpolation order used by grid transfer operators using these aggregates.

  Array<int> coarseRate;  ///< coarsening rate in each direction
  Array<int> endRate;     ///< adapted coarsening rate at the edge of the mesh in each direction.

  GO gNumFineNodes;                  ///< global number of nodes.
  GO gNumFineNodes10;                ///< global number of nodes per 0-1 slice.
  const Array<GO> gFineNodesPerDir;  ///< global number of nodes per direction.

  LO lNumFineNodes;                  ///< local number of nodes.
  LO lNumFineNodes10;                ///< local number of nodes per 0-1 slice.
  const Array<LO> lFineNodesPerDir;  ///< local number of nodes per direction.

  GO gNumCoarseNodes;            ///< global number of nodes remaining after coarsening.
  GO gNumCoarseNodes10;          ///< global number of nodes per 0-1 slice remaining after coarsening.
  Array<GO> gCoarseNodesPerDir;  ///< global number of nodes per direction remaining after coarsening.

  LO lNumCoarseNodes;            ///< local number of nodes remaining after coarsening.
  LO lNumCoarseNodes10;          ///< local number of nodes per 0-1 slice remaining after coarsening.
  Array<LO> lCoarseNodesPerDir;  ///< local number of nodes per direction remaing after coarsening.

  LO numGhostNodes;              ///< local number of ghost nodes
  LO numGhostedNodes;            ///< local number of ghosted nodes (i.e. ghost + coarse nodes).
  LO numGhostedNodes10;          ///< local number of ghosted nodes (i.e. ghost + coarse nodes) per 0-1 slice.
  Array<LO> ghostedNodesPerDir;  ///< local number of ghosted nodes (i.e. ghost + coarse nodes) per direction

  GO minGlobalIndex;                 ///< lowest GID of any node in the local process
  Array<LO> offsets;                 ///< distance between lowest (resp. highest) index to the lowest (resp. highest) ghostedNodeIndex in that direction.
  Array<LO> coarseNodeOffsets;       ///< distance between lowest (resp. highest) index to the lowest (resp. highest) coarseNodeIndex in that direction.
  Array<GO> startIndices;            ///< lowest global tuple (i,j,k) of a node on the local process
  Array<GO> startGhostedCoarseNode;  ///< lowest coarse global tuple (i,j,k) of a node remaing on the local process after coarsening.

  bool meshEdge[6]       = {false};  ///< flags indicating if we run into the edge of the mesh in ilo, ihi, jlo, jhi, klo or khi.
  bool ghostInterface[6] = {false};  ///< flags indicating if ghost points are needed at ilo, ihi, jlo, jhi, klo and khi boundaries.
  bool ghostedDir[6]     = {false};  ///< flags indicating if ghost points are needed at ilo, ihi, jlo, jhi, klo and khi boundaries.

 public:
  IndexManager() = default;

  IndexManager(const RCP<const Teuchos::Comm<int> > comm, const bool coupled,
               const bool singleCoarsePoint, const int NumDimensions,
               const int interpolationOrder, const Array<GO> GFineNodesPerDir,
               const Array<LO> LFineNodesPerDir);

  virtual ~IndexManager() {}

  //! Sets basic parameters used to compute indices on the mesh.
  //! This method requires you to have set this->coarseRate and this->startIndices.
  void computeMeshParameters();

  virtual void computeGlobalCoarseParameters() = 0;

  virtual void getGhostedNodesData(const RCP<const Map> fineMap,
                                   Array<LO>& ghostedNodeCoarseLIDs,
                                   Array<int>& ghostedNodeCoarsePIDs,
                                   Array<GO>& ghostedNodeCoarseGIDs) const = 0;

  virtual void getCoarseNodesData(const RCP<const Map> fineCoordinatesMap,
                                  Array<GO>& coarseNodeCoarseGIDs,
                                  Array<GO>& coarseNodeFineGIDs) const = 0;

  bool isAggregationCoupled() const { return coupled_; }

  bool isSingleCoarsePoint() const { return singleCoarsePoint_; }

  int getNumDimensions() const { return numDimensions; }

  int getInterpolationOrder() const { return interpolationOrder_; }

  GO getNumGlobalFineNodes() const { return gNumFineNodes; }

  GO getNumGlobalCoarseNodes() const { return gNumCoarseNodes; }

  LO getNumLocalFineNodes() const { return lNumFineNodes; }

  LO getNumLocalCoarseNodes() const { return lNumCoarseNodes; }

  LO getNumLocalGhostedNodes() const { return numGhostedNodes; }

  Array<int> getCoarseningRates() const { return coarseRate; }

  int getCoarseningRate(const int dim) const { return coarseRate[dim]; }

  Array<int> getCoarseningEndRates() const { return endRate; }

  int getCoarseningEndRate(const int dim) const { return endRate[dim]; }

  bool getMeshEdge(const int dir) const { return meshEdge[dir]; }

  bool getGhostInterface(const int dir) const { return ghostInterface[dir]; }

  Array<LO> getOffsets() const { return offsets; }

  LO getOffset(int const dim) const { return offsets[dim]; }

  Array<LO> getCoarseNodeOffsets() const { return coarseNodeOffsets; }

  LO getCoarseNodeOffset(int const dim) const { return coarseNodeOffsets[dim]; }

  Array<GO> getStartIndices() const { return startIndices; }

  GO getStartIndex(int const dim) const { return startIndices[dim]; }

  Array<GO> getStartGhostedCoarseNodes() const { return startGhostedCoarseNode; }

  GO getStartGhostedCoarseNode(int const dim) const { return startGhostedCoarseNode[dim]; }

  Array<LO> getLocalFineNodesPerDir() const { return lFineNodesPerDir; }

  LO getLocalFineNodesInDir(const int dim) const { return lFineNodesPerDir[dim]; }

  Array<GO> getGlobalFineNodesPerDir() const { return gFineNodesPerDir; }

  GO getGlobalFineNodesInDir(const int dim) const { return gFineNodesPerDir[dim]; }

  Array<LO> getLocalCoarseNodesPerDir() const { return lCoarseNodesPerDir; }

  LO getLocalCoarseNodesInDir(const int dim) const { return lCoarseNodesPerDir[dim]; }

  Array<GO> getGlobalCoarseNodesPerDir() const { return gCoarseNodesPerDir; }

  GO getGlobalCoarseNodesInDir(const int dim) const { return gCoarseNodesPerDir[dim]; }

  Array<LO> getGhostedNodesPerDir() const { return ghostedNodesPerDir; }

  LO getGhostedNodesInDir(const int dim) const { return ghostedNodesPerDir[dim]; }

  virtual std::vector<std::vector<GO> > getCoarseMeshData() const = 0;

  virtual void getFineNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const = 0;

  virtual void getFineNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const = 0;

  virtual void getFineNodeGhostedTuple(const LO myLID, LO& i, LO& j, LO& k) const = 0;

  virtual void getFineNodeGID(const GO i, const GO j, const GO k, GO& myGID) const = 0;

  virtual void getFineNodeLID(const LO i, const LO j, const LO k, LO& myLID) const = 0;

  virtual void getCoarseNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const = 0;

  virtual void getCoarseNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const = 0;

  virtual void getCoarseNodeGID(const GO i, const GO j, const GO k, GO& myGID) const = 0;

  virtual void getCoarseNodeLID(const LO i, const LO j, const LO k, LO& myLID) const = 0;

  virtual void getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const = 0;

  virtual void getCoarseNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const = 0;

  virtual void getGhostedNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const = 0;

  virtual void getGhostedNodeCoarseLID(const LO i, const LO j, const LO k, LO& myLID) const = 0;
};

}  // namespace MueLu

#define MUELU_INDEXMANAGER_SHORT
#endif  // MUELU_INDEXMANAGER_DECL_HPP

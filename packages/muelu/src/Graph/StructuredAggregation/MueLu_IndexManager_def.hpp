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
#ifndef MUELU_INDEXMANAGER_DEF_HPP
#define MUELU_INDEXMANAGER_DEF_HPP

#include "Teuchos_OrdinalTraits.hpp"

#include "MueLu_ConfigDefs.hpp"
#include <MueLu_IndexManager_decl.hpp>

/*****************************************************************************

****************************************************************************/

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
IndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    IndexManager(const RCP<const Teuchos::Comm<int> > comm,
                 const bool coupled,
                 const bool singleCoarsePoint,
                 const int NumDimensions,
                 const int interpolationOrder,
                 const Array<GO> GFineNodesPerDir,
                 const Array<LO> LFineNodesPerDir)
  : comm_(comm)
  , coupled_(coupled)
  , singleCoarsePoint_(singleCoarsePoint)
  , numDimensions(NumDimensions)
  , interpolationOrder_(interpolationOrder)
  , gFineNodesPerDir(GFineNodesPerDir)
  , lFineNodesPerDir(LFineNodesPerDir) {
  coarseRate.resize(3);
  endRate.resize(3);
  gCoarseNodesPerDir.resize(3);
  lCoarseNodesPerDir.resize(3);
  ghostedNodesPerDir.resize(3);

  offsets.resize(3);
  coarseNodeOffsets.resize(3);
  startIndices.resize(6);
  startGhostedCoarseNode.resize(3);

}  // Constructor

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void IndexManager<LocalOrdinal, GlobalOrdinal, Node>::
    computeMeshParameters() {
  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_INDEXMANAGER_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  if (coupled_) {
    gNumFineNodes10 = gFineNodesPerDir[1] * gFineNodesPerDir[0];
    gNumFineNodes   = gFineNodesPerDir[2] * gNumFineNodes10;
  } else {
    gNumFineNodes10 = Teuchos::OrdinalTraits<GO>::invalid();
    gNumFineNodes   = Teuchos::OrdinalTraits<GO>::invalid();
  }
  lNumFineNodes10 = lFineNodesPerDir[1] * lFineNodesPerDir[0];
  lNumFineNodes   = lFineNodesPerDir[2] * lNumFineNodes10;
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      if (coupled_) {
        if (startIndices[dim] == 0) {
          meshEdge[2 * dim] = true;
        }
        if (startIndices[dim + 3] + 1 == gFineNodesPerDir[dim]) {
          meshEdge[2 * dim + 1] = true;
          endRate[dim]          = startIndices[dim + 3] % coarseRate[dim];
        }
      } else {  // With uncoupled problem each rank might require a different endRate
        meshEdge[2 * dim]     = true;
        meshEdge[2 * dim + 1] = true;
        endRate[dim]          = (lFineNodesPerDir[dim] - 1) % coarseRate[dim];
      }
      if (endRate[dim] == 0) {
        endRate[dim] = coarseRate[dim];
      }

      // If uncoupled aggregation is used, offsets[dim] = 0, so nothing to do.
      if (coupled_) {
        offsets[dim] = Teuchos::as<LO>(startIndices[dim]) % coarseRate[dim];
        if (offsets[dim] == 0) {
          coarseNodeOffsets[dim] = 0;
        } else if (startIndices[dim] + endRate[dim] == lFineNodesPerDir[dim]) {
          coarseNodeOffsets[dim] = endRate[dim] - offsets[dim];
        } else {
          coarseNodeOffsets[dim] = coarseRate[dim] - offsets[dim];
        }

        if (interpolationOrder_ == 0) {
          int rem = startIndices[dim] % coarseRate[dim];
          if ((rem != 0) && (rem <= Teuchos::as<double>(coarseRate[dim]) / 2.0)) {
            ghostInterface[2 * dim] = true;
          }
          rem = startIndices[dim + 3] % coarseRate[dim];
          // uncoupled by nature does not require ghosts nodes
          if (coupled_ && (startIndices[dim + 3] != gFineNodesPerDir[dim] - 1) &&
              (rem > Teuchos::as<double>(coarseRate[dim]) / 2.0)) {
            ghostInterface[2 * dim + 1] = true;
          }

        } else if (interpolationOrder_ == 1) {
          if (coupled_ && (startIndices[dim] % coarseRate[dim] != 0 ||
                           startIndices[dim] == gFineNodesPerDir[dim] - 1)) {
            ghostInterface[2 * dim] = true;
          }
          if (coupled_ && (startIndices[dim + 3] != gFineNodesPerDir[dim] - 1) &&
              ((lFineNodesPerDir[dim] == 1) || (startIndices[dim + 3] % coarseRate[dim] != 0))) {
            ghostInterface[2 * dim + 1] = true;
          }
        }
      }
    } else {  // Default value for dim >= numDimensions
      endRate[dim] = 1;
    }
  }

  *out << "singleCoarsePoint? " << singleCoarsePoint_ << std::endl;
  *out << "gFineNodesPerDir: " << gFineNodesPerDir << std::endl;
  *out << "lFineNodesPerDir: " << lFineNodesPerDir << std::endl;
  *out << "endRate: " << endRate << std::endl;
  *out << "ghostInterface: {" << ghostInterface[0] << ", " << ghostInterface[1] << ", "
       << ghostInterface[2] << ", " << ghostInterface[3] << ", " << ghostInterface[4] << ", "
       << ghostInterface[5] << "}" << std::endl;
  *out << "meshEdge: {" << meshEdge[0] << ", " << meshEdge[1] << ", "
       << meshEdge[2] << ", " << meshEdge[3] << ", " << meshEdge[4] << ", "
       << meshEdge[5] << "}" << std::endl;
  *out << "startIndices: " << startIndices << std::endl;
  *out << "offsets: " << offsets << std::endl;
  *out << "coarseNodeOffsets: " << coarseNodeOffsets << std::endl;

  // Here one element can represent either the degenerate case of one node or the more general
  // case of two nodes, i.e. x---x is a 1D element with two nodes and x is a 1D element with
  // one node. This helps generating a 3D space from tensorial products...
  // A good way to handle this would be to generalize the algorithm to take into account the
  // discretization order used in each direction, at least in the FEM sense, since a 0 degree
  // discretization will have a unique node per element. This way 1D discretization can be
  // viewed as a 3D problem with one 0 degree element in the y direction and one 0 degre
  // element in the z direction.
  // !!! Operations below are aftecting both local and global values that have two         !!!
  // different orientations. Orientations can be interchanged using mapDirG2L and mapDirL2G.
  // coarseRate, endRate and offsets are in the global basis, as well as all the variables
  // starting with a g.
  // !!! while the variables starting with an l are in the local basis.                    !!!
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Check whether the partition includes the "end" of the mesh which means that endRate
      // will apply. Also make sure that endRate is not 0 which means that the mesh does not
      // require a particular treatment at the boundaries.
      if (meshEdge[2 * dim + 1]) {
        lCoarseNodesPerDir[dim] = (lFineNodesPerDir[dim] - endRate[dim] + offsets[dim] - 1) / coarseRate[dim] + 1;
        if (offsets[dim] == 0) {
          ++lCoarseNodesPerDir[dim];
        }
        // We might want to coarsening the direction
        // into a single layer if there are not enough
        // points left to form two aggregates
        if (singleCoarsePoint_ && lFineNodesPerDir[dim] - 1 < coarseRate[dim]) {
          lCoarseNodesPerDir[dim] = 1;
        }
      } else {
        lCoarseNodesPerDir[dim] = (lFineNodesPerDir[dim] + offsets[dim] - 1) / coarseRate[dim];
        if (offsets[dim] == 0) {
          ++lCoarseNodesPerDir[dim];
        }
      }

      // The first branch of this if-statement will be used if the rank contains only one layer
      // of nodes in direction i, that layer must also coincide with the boundary of the mesh
      // and coarseRate[i] == endRate[i]...
      if (interpolationOrder_ == 0) {
        startGhostedCoarseNode[dim] = startIndices[dim] / coarseRate[dim];
        int rem                     = startIndices[dim] % coarseRate[dim];
        if (rem > (Teuchos::as<double>(coarseRate[dim]) / 2.0)) {
          ++startGhostedCoarseNode[dim];
        }
      } else {
        if ((startIndices[dim] == gFineNodesPerDir[dim] - 1) &&
            (startIndices[dim] % coarseRate[dim] == 0)) {
          startGhostedCoarseNode[dim] = startIndices[dim] / coarseRate[dim] - 1;
        } else {
          startGhostedCoarseNode[dim] = startIndices[dim] / coarseRate[dim];
        }
      }

      // This array is passed to the RAPFactory and eventually becomes gFineNodePerDir on the next
      // level.
      gCoarseNodesPerDir[dim] = (gFineNodesPerDir[dim] - 1) / coarseRate[dim];
      if ((gFineNodesPerDir[dim] - 1) % coarseRate[dim] == 0) {
        ++gCoarseNodesPerDir[dim];
      } else {
        gCoarseNodesPerDir[dim] += 2;
      }
    } else {  // Default value for dim >= numDimensions
      // endRate[dim] = 1;
      gCoarseNodesPerDir[dim] = 1;
      lCoarseNodesPerDir[dim] = 1;
    }  // if (dim < numDimensions)

    // This would happen if the rank does not own any nodes but in that case a subcommunicator
    // should be used so this should really not be a concern.
    if (lFineNodesPerDir[dim] < 1) {
      lCoarseNodesPerDir[dim] = 0;
    }
    ghostedNodesPerDir[dim] = lCoarseNodesPerDir[dim];
    // Check whether face *low needs ghost nodes
    if (ghostInterface[2 * dim]) {
      ghostedNodesPerDir[dim] += 1;
    }
    // Check whether face *hi needs ghost nodes
    if (ghostInterface[2 * dim + 1]) {
      ghostedNodesPerDir[dim] += 1;
    }
  }  // Loop for dim=0:3

  // With uncoupled aggregation we need to communicate to compute the global number of coarse points
  if (!coupled_) {
    for (int dim = 0; dim < 3; ++dim) {
      gCoarseNodesPerDir[dim] = -1;
    }
  }

  // Compute cummulative values
  lNumCoarseNodes10 = lCoarseNodesPerDir[0] * lCoarseNodesPerDir[1];
  lNumCoarseNodes   = lNumCoarseNodes10 * lCoarseNodesPerDir[2];
  numGhostedNodes10 = ghostedNodesPerDir[1] * ghostedNodesPerDir[0];
  numGhostedNodes   = numGhostedNodes10 * ghostedNodesPerDir[2];
  numGhostNodes     = numGhostedNodes - lNumCoarseNodes;

  *out << "lCoarseNodesPerDir: " << lCoarseNodesPerDir << std::endl;
  *out << "gCoarseNodesPerDir: " << gCoarseNodesPerDir << std::endl;
  *out << "ghostedNodesPerDir: " << ghostedNodesPerDir << std::endl;
  *out << "lNumCoarseNodes=" << lNumCoarseNodes << std::endl;
  *out << "numGhostedNodes=" << numGhostedNodes << std::endl;
}

}  // namespace MueLu

#define MUELU_INDEXMANAGER_SHORT
#endif  // MUELU_INDEXMANAGER_DEF_HPP

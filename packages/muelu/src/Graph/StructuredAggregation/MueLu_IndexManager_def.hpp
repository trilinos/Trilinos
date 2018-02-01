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

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include <MueLu_IndexManager_decl.hpp>

/*****************************************************************************

****************************************************************************/

namespace MueLu {

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  IndexManager<LocalOrdinal, GlobalOrdinal, Node>::IndexManager(const int NumDimensions,
                                                                const Array<GO> GFineNodesPerDir,
                                                                const Array<LO> LFineNodesPerDir) :
    numDimensions(NumDimensions), gFineNodesPerDir(GFineNodesPerDir),
    lFineNodesPerDir(LFineNodesPerDir) {

    coarseRate.resize(3);
    endRate.resize(3);
    gCoarseNodesPerDir.resize(3);
    lCoarseNodesPerDir.resize(3);
    ghostedCoarseNodesPerDir.resize(3);

    offsets.resize(6);
    startIndices.resize(6);
    startGhostedCoarseNode.resize(3);

  }

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  void IndexManager<LocalOrdinal, GlobalOrdinal, Node>::computeMeshParameters() {

    gNumFineNodes10 = gFineNodesPerDir[1]*gFineNodesPerDir[0];
    gNumFineNodes   = gFineNodesPerDir[2]*gNumFineNodes10;
    lNumFineNodes10 = lFineNodesPerDir[1]*lFineNodesPerDir[0];
    lNumFineNodes   = lFineNodesPerDir[2]*lNumFineNodes10;
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        offsets[dim]     = Teuchos::as<LO>(startIndices[dim]) % coarseRate[dim];
        offsets[dim + 3] = Teuchos::as<LO>(startIndices[dim]) % coarseRate[dim];

        if(startIndices[dim] % coarseRate[dim] != 0 ||
           startIndices[dim] == gFineNodesPerDir[dim]-1) {
          ghostInterface[2*dim] = true;
        }
        if((startIndices[dim + 3] != gFineNodesPerDir[dim] - 1) &&
           ((lFineNodesPerDir[dim] == 1) || (startIndices[dim + 3] % coarseRate[dim] != 0))) {
          ghostInterface[2*dim+1] = true;
        }
      }
    }

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
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        // This array is passed to the RAPFactory and eventually becomes gFineNodePerDir on the next
        // level.
        gCoarseNodesPerDir[dim] = (gFineNodesPerDir[dim] - 1) / coarseRate[dim];
        endRate[dim] = Teuchos::as<LO>((gFineNodesPerDir[dim] - 1) % coarseRate[dim]);
        if(endRate[dim] == 0) {
          endRate[dim] = coarseRate[dim];
          ++gCoarseNodesPerDir[dim];
        } else {
          gCoarseNodesPerDir[dim] += 2;
        }

        // Check whether the partition includes the "end" of the mesh which means that endRate
        // will apply. Also make sure that endRate is not 0 which means that the mesh does not
        // require a particular treatment at the boundaries.
        if( (startIndices[dim] + lFineNodesPerDir[dim]) == gFineNodesPerDir[dim] ) {
          lCoarseNodesPerDir[dim] = (lFineNodesPerDir[dim] - endRate[dim] + offsets[dim] - 1)
            / coarseRate[dim] + 1;
          if(offsets[dim] == 0) {++lCoarseNodesPerDir[dim];}
        } else {
          lCoarseNodesPerDir[dim] = (lFineNodesPerDir[dim] + offsets[dim] - 1) / coarseRate[dim];
          if(offsets[dim] == 0) {++lCoarseNodesPerDir[dim];}
        }

        // The first branch of this if-statement will be used if the rank contains only one layer
        // of nodes in direction i, that layer must also coincide with the boundary of the mesh
        // and coarseRate[i] == endRate[i]...
        if((startIndices[dim] == gFineNodesPerDir[dim] - 1) &&
           (startIndices[dim] % coarseRate[dim] == 0)) {
          startGhostedCoarseNode[dim] = startIndices[dim] / coarseRate[dim] - 1;
        } else {
          startGhostedCoarseNode[dim] = startIndices[dim] / coarseRate[dim];
        }
      } else { // Default value for dim >= numDimensions
        endRate[dim] = 1;
        gCoarseNodesPerDir[dim] = 1;
        lCoarseNodesPerDir[dim] = 1;
      }

      // This would happen if the rank does not own any nodes but in that case a subcommunicator
      // should be used so this should really not be a concern.
      if(lFineNodesPerDir[dim] < 1) {lCoarseNodesPerDir[dim] = 0;}
      ghostedCoarseNodesPerDir[dim] = lCoarseNodesPerDir[dim];
      // Check whether face *low needs ghost nodes
      if(ghostInterface[2*dim]) {ghostedCoarseNodesPerDir[dim] += 1;}
      // Check whether face *hi needs ghost nodes
      if(ghostInterface[2*dim + 1]) {ghostedCoarseNodesPerDir[dim] += 1;}
    }


    // Compute cummulative values
    gNumCoarseNodes10       = gCoarseNodesPerDir[0]*gCoarseNodesPerDir[1];
    gNumCoarseNodes         = gNumCoarseNodes10*gCoarseNodesPerDir[2];
    lNumCoarseNodes10       = lCoarseNodesPerDir[0]*lCoarseNodesPerDir[1];
    lNumCoarseNodes         = lNumCoarseNodes10*lCoarseNodesPerDir[2];
    numGhostedCoarseNodes10 = ghostedCoarseNodesPerDir[1]*ghostedCoarseNodesPerDir[0];
    numGhostedCoarseNodes   = numGhostedCoarseNodes10*ghostedCoarseNodesPerDir[2];
    lNumGhostNodes          = numGhostedCoarseNodes - lNumCoarseNodes;

  }

} //namespace MueLu

#define MUELU_INDEXMANAGER_SHORT
#endif // MUELU_INDEXMANAGER_DEF_HPP

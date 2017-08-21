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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_GENERALGEOMETRICPFACTORY_DEF_HPP
#define MUELU_GENERALGEOMETRICPFACTORY_DEF_HPP

#include <stdlib.h>
#include <iomanip>


// #include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_GeneralGeometricPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    // Coarsen can come in two forms, either a single char that will be interpreted as an integer which is used as the coarsening
    // rate in every spatial dimentions, or it can be a longer string that will then be interpreted as an array of integers.
    // By default coarsen is set as "{2}", hence a coarsening rate of 2 in every spatial dimension is the default setting!
    validParamList->set<std::string >           ("Coarsen",                 "{2}", "Coarsening rate in all spatial dimensions");
    validParamList->set<int>                    ("order",                   1, "Order of the interpolation scheme used");
    validParamList->set<RCP<const FactoryBase> >("A",                         Teuchos::null, "Generating factory of the matrix A");
    validParamList->set<RCP<const FactoryBase> >("Nullspace",                 Teuchos::null, "Generating factory of the nullspace");
    validParamList->set<RCP<const FactoryBase> >("Coordinates",               Teuchos::null, "Generating factory for coorindates");
    validParamList->set<RCP<const FactoryBase> >("gNodesPerDim",              Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",              Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<std::string >("axisPermutation",                    "{0,1,2}", "Assuming a global (x,y,z) orientation, local might be (z,y,x). This vector gives a permutation from global to local orientation.");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "Coordinates");
    // Request the global number of nodes per dimensions
    if(fineLevel.GetLevelID() == 0) {
      if(fineLevel.IsAvailable("gNodesPerDim", NoFactory::get())) {
        fineLevel.DeclareInput("gNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("gNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "gNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(fineLevel, "gNodesPerDim");
    }

    // Request the local number of nodes per dimensions
    if(fineLevel.GetLevelID() == 0) {
      if(fineLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
        fineLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "lNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(fineLevel, "lNodesPerDim");
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    // Obtain general variables
    using xdMV = Xpetra::MultiVector<double,LO,GO,NO>;
    RCP<Matrix>      A             = Get< RCP<Matrix> >      (fineLevel, "A");
    RCP<MultiVector> fineNullspace = Get< RCP<MultiVector> > (fineLevel, "Nullspace");
    RCP<xdMV>        fineCoords    = Get< RCP<xdMV> >(fineLevel, "Coordinates");
    RCP<xdMV>        coarseCoords;

    // Get user-provided coarsening rate parameter (constant over all levels)
    const ParameterList& pL = GetParameterList();

    // collect general input data
    LO blkSize                 = A->GetFixedBlockSize();
    RCP<const Map> rowMap      = A->getRowMap();
    LO numDimensions           = 0;                             // Number of spatial dimensions
    Array<GO> gFineNodesPerDir(3);                              // Global number of fine points per direction
    Array<GO> gCoarseNodesPerDir(3);                            // Global number of coarse points per direction
    LO lNumFinePoints;                                          // Local number of fine points
    Array<LO> lFineNodesPerDir(3);                              // Local number of fine points per direction
    Array<LO> lCoarseNodesPerDir(3);                            // Local number of coarse points per direction
    Array<LO> mapDirL2G(3);
    Array<LO> mapDirG2L(3);

    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords==Teuchos::null, Exceptions::RuntimeError, "Coordinates cannot be accessed from fine level!");
    numDimensions  = fineCoords->getNumVectors();
    lNumFinePoints = fineCoords->getLocalLength();

    // Get the number of points in each direction
    if(fineLevel.GetLevelID() == 0) {
      gFineNodesPerDir = fineLevel.Get<Array<GO> >("gNodesPerDim", NoFactory::get());
      lFineNodesPerDir = fineLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
    } else {
      // Loading global number of nodes per diretions
      gFineNodesPerDir = Get<Array<GO> >(fineLevel, "gNodesPerDim");

      // Loading local number of nodes per diretions
      lFineNodesPerDir = Get<Array<LO> >(fineLevel, "lNodesPerDim");
    }

    // Get the coarsening rate
    std::string coarsenRate = pL.get<std::string>("Coarsen");
    ArrayRCP<LO> coarseRate = arcp<LO>(3);
    Teuchos::Array<LO> crates;
    try {
      crates = Teuchos::fromStringToArray<LO>(coarsenRate);
    } catch(const Teuchos::InvalidArrayStringRepresentation e) {
      GetOStream(Errors,-1) << " *** Coarsen must be a string convertible into an array! *** " << std::endl;
      throw e;
    }
    TEUCHOS_TEST_FOR_EXCEPTION((crates.size()>1) && (crates.size()<numDimensions),
                               Exceptions::RuntimeError,
                               "Coarsen must have at least as many components as the number of spatial dimensions in the problem.");
    for(LO i = 0; i < 3; ++i) {
      if(i < numDimensions) {
        if( crates.size()==1 ) {
          coarseRate[i] = crates[0];
        } else if( crates.size()==numDimensions ) {
          coarseRate[i] = crates[i];
        }
      } else {
        coarseRate[i] = 1;
      }
    }

    int interpolationOrder = pL.get<int>("order");
    TEUCHOS_TEST_FOR_EXCEPTION((interpolationOrder < 0) || (interpolationOrder > 1),
                               Exceptions::RuntimeError,
                               "The interpolation order can only be set to 0 or 1.");

    // Get the axis permutation from Global axis to Local axis
    std::string axisPermutation = pL.get<std::string>("axisPermutation");
    try {
      mapDirG2L = Teuchos::fromStringToArray<LO>(axisPermutation);
    } catch(const Teuchos::InvalidArrayStringRepresentation e) {
      GetOStream(Errors,-1) << " *** axisPermutation must be a string convertible into an array! *** " << std::endl;
      throw e;
    }
    for(LO i = 0; i < numDimensions; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(mapDirG2L[i] > numDimensions,
                                 Exceptions::RuntimeError,
                                 "axis permutation values must all be less than the number of spatial dimensions.");
      mapDirL2G[mapDirG2L[i]] = i;
    }

    // Find the offsets for the nodes on the current processor, this tells us the position of the first node on the partition compare to the
    // closest lower coarse point. This information is needed to know who are the coarse nodes that need to be used for interpolation and by
    // extension this tells us what are the ghostnodes we need for the computation.
    //
    //     x - o - o - x
    //     |   |   |   |
    //     o - o - o - o
    //     |   |   |   |
    //     o - o - * - o
    //     |   |   |   |
    //     x - o - o - x
    //
    //     x are coarse points, o are fine points and * is the first point on the local processor.
    //     the offsets of * are 2 and 1
    GO minGlobalIndex, maxGlobalIndex, rem;
    GO gIndices[6] = {0, 0, 0, 0, 0, 0};
    LO offsets[6] = {0, 0, 0, 0, 0, 0};
    RCP<const Map> fineCoordsMap = fineCoords->getMap();
    minGlobalIndex = fineCoordsMap->getMinGlobalIndex();
    maxGlobalIndex = fineCoordsMap->getMaxGlobalIndex();
    if(numDimensions == 1) {
      gIndices[0] = minGlobalIndex;
      offsets[0]  = Teuchos::as<LO>(gIndices[0]) % coarseRate[0];

      gIndices[3] = maxGlobalIndex;
      offsets[3]  = Teuchos::as<LO>(gIndices[3]) % coarseRate[0];
    } else if(numDimensions == 2) {
      gIndices[1] = minGlobalIndex / Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[1]  = Teuchos::as<LO>(gIndices[1]) % coarseRate[1];
      gIndices[0] = minGlobalIndex % Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[0]  = Teuchos::as<LO>(gIndices[0]) % coarseRate[0];

      gIndices[4] = maxGlobalIndex / Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[4]  = Teuchos::as<LO>(gIndices[4]) % coarseRate[1];
      gIndices[3] = maxGlobalIndex % Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[3]  = Teuchos::as<LO>(gIndices[3]) % coarseRate[3];
    } else if(numDimensions == 3) {
      gIndices[2] = minGlobalIndex / Teuchos::as<GO>(gFineNodesPerDir[0]*gFineNodesPerDir[1]);
      rem         = minGlobalIndex % Teuchos::as<GO>(gFineNodesPerDir[0]*gFineNodesPerDir[1]);
      offsets[2]  = Teuchos::as<LO>(gIndices[2]) % coarseRate[2];
      gIndices[1] = rem / Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[1]  = Teuchos::as<LO>(gIndices[1]) % coarseRate[1];
      gIndices[0] = rem % Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[0]  = Teuchos::as<LO>(gIndices[0]) % coarseRate[0];

      gIndices[5] = maxGlobalIndex / Teuchos::as<GO>(gFineNodesPerDir[0]*gFineNodesPerDir[1]);
      rem         = maxGlobalIndex % Teuchos::as<GO>(gFineNodesPerDir[0]*gFineNodesPerDir[1]);
      offsets[5]  = Teuchos::as<LO>(gIndices[5]) % coarseRate[2];
      gIndices[4] = rem / Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[4]  = Teuchos::as<LO>(gIndices[4]) % coarseRate[1];
      gIndices[3] = rem % Teuchos::as<GO>(gFineNodesPerDir[0]);
      offsets[3]  = Teuchos::as<LO>(gIndices[3]) % coarseRate[0];
    }
    /* At this point the local geometrical discovery is over */

    // Check if the partition contains nodes on a boundary, if so that boundary (face, line or point) will not require ghost nodes.
    // Always include a layer if ghost nodes for inner interface since even for faces with coarse nodes, an outward curvature might
    // require ghost nodes to compute a proper interpolation operator. Curvature could be check localy to avoid including extra
    // ghost nodes...
    bool ghostInterface[6] = {false, false, false, false, false, false};
    for(LO i=0; i < numDimensions; ++i) {
      if(  gIndices[i] != 0) {
        ghostInterface[2*i]=true;
      }
      if( (gIndices[i]+lFineNodesPerDir[mapDirL2G[i]]) != gFineNodesPerDir[i] ) {
        ghostInterface[2*i+1]=true;
      }
    }

    /* Here one element can represent either the degenerate case of one node or the more general case of two nodes.
       i.e. x---x is a 1D element with two nodes and x is a 1D element with one node. This helps generating a 3D space from tensorial products...
       A good way to handle this would be to generalize the algorithm to take into account the discretization order used in each direction,
       at least in the FEM sense, since a 0 degree discretization will have a unique node per element. This way a 1D discretization can be viewed
       as a 3D problem with one 0 degree element in the y direction and one 0 degre element in the z direction. */
    LO endRate[3] = {0, 0, 0};
    /* /!\ Operations below are aftecting both local and global values that have two different orientations.   /!\
           orientations can be interchanged using mapDirG2L and mapDirL2G. coarseRate, endRate and offsets
           are in the global basis, as well as all the variables starting with a g, while the variables
       /!\ starting with an l are in the local basis.                                                          /!\ */
    for(LO i = 0; i < 3; ++i) {
      if(i < numDimensions) {
        // This array is passed to the RAPFactory and eventually becomes gFineNodePerDir on the next level.
        gCoarseNodesPerDir[i] = (gFineNodesPerDir[i] - 1) / coarseRate[i];
        endRate[i] = Teuchos::as<LO>((gFineNodesPerDir[i] - 1) % coarseRate[i]);
        if(endRate[i] == 0) {
          endRate[i] = coarseRate[i];
          ++gCoarseNodesPerDir[i];
        } else {
          gCoarseNodesPerDir[i] += 2;
        }
      } else {
        endRate[i] = 1;
        gCoarseNodesPerDir[i] = 1;
      }
    }

    for(LO i = 0; i < 3; ++i) {
      if(i < numDimensions) {
        // Check whether the partition includes the "end" of the mesh which means that endRate will apply.
        // Also make sure that endRate is not 0 which means that the mesh does not require a particular treatment
        // at the boundaries.
        if( (gIndices[mapDirG2L[i]] + lFineNodesPerDir[i]) == gFineNodesPerDir[mapDirG2L[i]] ) {
          lCoarseNodesPerDir[i] = (lFineNodesPerDir[i] - endRate[mapDirG2L[i]] + offsets[mapDirG2L[i]] - 1) / coarseRate[mapDirG2L[i]] + 1;
          if(offsets[mapDirG2L[i]] == 0) {++lCoarseNodesPerDir[i];}
        } else {
          lCoarseNodesPerDir[i] = (lFineNodesPerDir[i] + offsets[mapDirG2L[i]] - 1) / coarseRate[mapDirG2L[i]];
          if(offsets[mapDirG2L[i]] == 0) {++lCoarseNodesPerDir[i];}
        }
      } else {
        lCoarseNodesPerDir[i] = 1;
      }
      if(lFineNodesPerDir[i] < 1) {lCoarseNodesPerDir[i] = 0;}
    }

    // Assuming linear interpolation, each fine point has contribution from 8 coarse points
    // and each coarse point value gets injected.
    // For systems of PDEs we assume that all dofs have the same P operator.
    LO lNumCoarsePoints = lCoarseNodesPerDir[0]*lCoarseNodesPerDir[1]*lCoarseNodesPerDir[2];
    LO nTerms = 8*lNumFinePoints - 7*lNumCoarsePoints;
    nTerms=nTerms*blkSize;

    // For each direction, determine how many ghost points are required.
    LO numGhosts = 0;
    const LO complementaryIndices[3][2] = {{1,2}, {0,2}, {0,1}};
    for(LO i = 0; i < 3; ++i) {
      LO tmp = 0;
      // Check whether a face in direction i needs ghost nodes
      if(ghostInterface[2*i] || ghostInterface[2*i+1]) {
        if(i == 0) {tmp = lCoarseNodesPerDir[mapDirL2G[1]]*lCoarseNodesPerDir[mapDirL2G[2]];}
        if(i == 1) {tmp = lCoarseNodesPerDir[mapDirL2G[0]]*lCoarseNodesPerDir[mapDirL2G[2]];}
        if(i == 2) {tmp = lCoarseNodesPerDir[mapDirL2G[0]]*lCoarseNodesPerDir[mapDirL2G[1]];}
      }
      // If both faces in direction i need nodes, double the number of ghost nodes
      if(ghostInterface[2*i] && ghostInterface[2*i+1]) {tmp = 2*tmp;}
      numGhosts += tmp;

      // The corners and edges need to be checked in 2D / 3D to add more ghosts nodes
      for(LO j = 0; j < 2; ++j) {
        for(LO k = 0; k < 2; ++k) {
          // Check if two adjoining faces need ghost nodes and then add their common edge
          if(ghostInterface[2*complementaryIndices[i][0]+j] && ghostInterface[2*complementaryIndices[i][1]+k]) {
            numGhosts += lCoarseNodesPerDir[mapDirL2G[i]];
            // Add corners if three adjoining faces need ghost nodes,
            // but add them only once! Hence when i == 0.
            if(ghostInterface[2*i] && (i == 0)) { numGhosts += 1; }
            if(ghostInterface[2*i+1] && (i == 0)) { numGhosts += 1; }
          }
        }
      }
    }

    // Build a list of GIDs to import the required ghost nodes.
    // The ordering of the ghosts nodes will be as natural as possible,
    // i.e. it should follow the ordering of the mesh.
    //
    // Saddly we have to more or less redo what was just done to figure out the value of numGhosts,
    // there might be some optimization possibility here...
    Array<GO> ghostsGIDs(numGhosts);
    LO countGhosts = 0;
    // Get the GID of the first point on the current partition.
    GO startingGID = minGlobalIndex;
    Array<GO> startingIndices(3);
    // We still want ghost nodes even if have with a 0 offset,
    // except when we are on a boundary
    if(ghostInterface[4] && (offsets[2] == 0)) {
      if(gIndices[2] + coarseRate[2] > gFineNodesPerDir[2]) {
        startingGID -= endRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
      } else {
        startingGID -= coarseRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
      }
    } else {
      startingGID -= offsets[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
    }
    if(ghostInterface[2] && (offsets[1] == 0)) {
      if(gIndices[1] + coarseRate[1] > gFineNodesPerDir[1]) {
        startingGID -= endRate[1]*gFineNodesPerDir[0];
      } else {
        startingGID -= coarseRate[1]*gFineNodesPerDir[0];
      }
    } else {
      startingGID -= offsets[1]*gFineNodesPerDir[0];
    }
    if(ghostInterface[0] && (offsets[0] == 0)) {
      if(gIndices[0] + coarseRate[0] > gFineNodesPerDir[0]) {
        startingGID -= endRate[0];
      } else {
        startingGID -= coarseRate[0];
      }
    } else {
      startingGID -= offsets[0];
    }

    { // scope for tmp
      GO tmp;
      startingIndices[2] = startingGID / (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      tmp = startingGID % (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      startingIndices[1] = tmp / gFineNodesPerDir[0];
      startingIndices[0] = tmp % gFineNodesPerDir[0];
    }

    GO ghostOffset[3] = {0, 0, 0};
    LO lengthZero  = lCoarseNodesPerDir[mapDirL2G[0]];
    LO lengthOne   = lCoarseNodesPerDir[mapDirL2G[1]];
    LO lengthTwo   = lCoarseNodesPerDir[mapDirL2G[2]];
    if(ghostInterface[0]) {++lengthZero;}
    if(ghostInterface[1]) {++lengthZero;}
    if(ghostInterface[2]) {++lengthOne;}
    if(ghostInterface[3]) {++lengthOne;}
    if(ghostInterface[4]) {++lengthTwo;}
    if(ghostInterface[5]) {++lengthTwo;}

    // First check the bottom face as it will have the lowest GIDs
    if(ghostInterface[4]) {
      ghostOffset[2] = startingGID;
      for(LO j = 0; j < lengthOne; ++j) {
        if( (j == lengthOne-1) && (startingIndices[1] + j*coarseRate[1] + 1 > gFineNodesPerDir[1]) ) {
          ghostOffset[1] = ((j-1)*coarseRate[1] + endRate[1])*gFineNodesPerDir[0];
        } else {
          ghostOffset[1] = j*coarseRate[1]*gFineNodesPerDir[0];
        }
        for(LO k = 0; k < lengthZero; ++k) {
          if( (k == lengthZero-1) && (startingIndices[0] + k*coarseRate[0] + 1 > gFineNodesPerDir[0]) ) {
            ghostOffset[0] = (k-1)*coarseRate[0] + endRate[0];
          } else {
            ghostOffset[0] = k*coarseRate[0];
          }
          // If the partition includes a changed rate at the edge, ghost nodes need to be picked carefully.
          // This if statement is repeated each time a ghost node is selected.
          ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
          ++countGhosts;
        }
      }
    }

    // Sweep over the lCoarseNodesPerDir[2] coarse layers in direction 2 and gather necessary ghost nodes
    // located on these layers.
    for(LO i = 0; i < lengthTwo; ++i) {
      // Exclude the cases where ghost nodes exists on the faces in directions 2, these faces are swept
      // seperatly for efficiency.
      if( !((i == lengthTwo-1) && ghostInterface[5]) && !((i == 0) && ghostInterface[4]) ) {
        // Set the ghostOffset in direction 2 taking into account a possible endRate different
        // from the regular coarseRate.
        if( (i == lengthTwo-1) && !ghostInterface[5] ) {
          ghostOffset[2] = startingGID + ((i-1)*coarseRate[2] + endRate[2])*gFineNodesPerDir[1]*gFineNodesPerDir[0];
        } else {
          ghostOffset[2] = startingGID + i*coarseRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
        }
        for(LO j = 0; j < lengthOne; ++j) {
          if( (j == 0) && ghostInterface[2] ) {
            for(LO k = 0; k < lengthZero; ++k) {
              if( (k == lengthZero-1) && (startingIndices[0] + k*coarseRate[0] + 1 > gFineNodesPerDir[0]) ) {
                if(k == 0) {
                  ghostOffset[0] = endRate[0];
                } else {
                  ghostOffset[0] = (k-1)*coarseRate[0] + endRate[0];
                }
              } else {
                ghostOffset[0] = k*coarseRate[0];
              }
              // In this case j == 0 so ghostOffset[1] is zero and can be ignored
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[0];
              ++countGhosts;
            }
          } else if( (j == lengthOne-1) && ghostInterface[3] ) {
            // Set the ghostOffset in direction 1 taking into account a possible endRate different
            // from the regular coarseRate.
            if( (j == lengthOne-1) && (startingIndices[1] + j*coarseRate[1] + 1 > gFineNodesPerDir[1]) ) {
              ghostOffset[1] = ((j-1)*coarseRate[1] + endRate[1])*gFineNodesPerDir[0];
            } else {
              ghostOffset[1] = j*coarseRate[1]*gFineNodesPerDir[0];
            }
            for(LO k = 0; k < lengthZero; ++k) {
              if( (k == lengthZero-1) && (startingIndices[0] + k*coarseRate[0] + 1 > gFineNodesPerDir[0]) ) {
                ghostOffset[0] = (k-1)*coarseRate[0] + endRate[0];
              } else {
                ghostOffset[0] = k*coarseRate[0];
              }
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
              ++countGhosts;
            }
          } else {
            // Set ghostOffset[1] for side faces sweep
            if( (j == lengthOne-1) && (startingIndices[1] + j*coarseRate[1] + 1 > gFineNodesPerDir[1]) ) {
              ghostOffset[1] = ( (j-1)*coarseRate[1] + endRate[1] )*gFineNodesPerDir[0];
            } else {
              ghostOffset[1] = j*coarseRate[1]*gFineNodesPerDir[0];
            }

            // Set ghostOffset[0], ghostsGIDs and countGhosts
            if(ghostInterface[0]) { // In that case ghostOffset[0]==0, so we can ignore it
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1];
              ++countGhosts;
            }
            if(ghostInterface[1]) { // Grab ghost point at the end of direction 0.
              if( (startingIndices[0] + (lengthZero-1)*coarseRate[0]) > gFineNodesPerDir[0] - 1 ) {
                if(lengthZero > 2) {
                  ghostOffset[0] = (lengthZero-2)*coarseRate[0] + endRate[0];
                } else {
                  ghostOffset[0] = endRate[0];
                }
              } else {
                ghostOffset[0] = (lengthZero-1)*coarseRate[0];
              }
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
              ++countGhosts;
            }
          }
        }
      }
    }

    // Finally check the top face
    if(ghostInterface[5]) {
      if( startingIndices[2] + (lengthTwo-1)*coarseRate[2] + 1 > gFineNodesPerDir[2] ) {
        ghostOffset[2] = startingGID + ((lengthTwo-2)*coarseRate[2] + endRate[2])*gFineNodesPerDir[1]*gFineNodesPerDir[0];
      } else {
        ghostOffset[2] = startingGID + (lengthTwo-1)*coarseRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
      }
      for(LO j = 0; j < lengthOne; ++j) {
        if( (j == lengthOne-1) && (startingIndices[1] + j*coarseRate[1] + 1 > gFineNodesPerDir[1]) ) { // && !ghostInterface[3] ) {
          ghostOffset[1] = ( (j-1)*coarseRate[1] + endRate[1] )*gFineNodesPerDir[0];
        } else {
          ghostOffset[1] = j*coarseRate[1]*gFineNodesPerDir[0];
        }
        for(LO k = 0; k < lengthZero; ++k) {
          if( (k == lengthZero-1) && (startingIndices[0] + k*coarseRate[0] + 1 > gFineNodesPerDir[0]) ) {// && !ghostInterface[1] ) {
            ghostOffset[0] = (k-1)*coarseRate[0] + endRate[0];
          } else {
            ghostOffset[0] = k*coarseRate[0];
          }
          ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
          ++countGhosts;
        }
      }
    }

    // All that is left to do is loop over NCpts and:
    //   - extract coarse points coordiate for coarseCoords
    //   - get coordinates for current stencil computation
    //   - compute current stencil
    //   - compute row and column indices for stencil entries
    RCP<const Map> stridedDomainMapP;
    RCP<Matrix>    P;
    MakeGeneralGeometricP(numDimensions, mapDirL2G, mapDirG2L, lFineNodesPerDir, lCoarseNodesPerDir, gCoarseNodesPerDir,
                          gFineNodesPerDir, coarseRate, endRate, offsets, ghostInterface, fineCoords, nTerms, blkSize,
                          stridedDomainMapP, A, P, coarseCoords, ghostsGIDs, interpolationOrder);

    // set StridingInformation of P
    if (A->IsView("stridedMaps") == true) {
      P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), stridedDomainMapP);
    } else {
      P->CreateView("stridedMaps", P->getRangeMap(), stridedDomainMapP);
    }

    // store the transfer operator and node coordinates on coarse level
    Set(coarseLevel, "P", P);
    Set(coarseLevel, "coarseCoordinates", coarseCoords);
    Set<Array<GO> >(coarseLevel, "gCoarseNodesPerDim", gCoarseNodesPerDir);
    Set<Array<LO> >(coarseLevel, "lCoarseNodesPerDim", lCoarseNodesPerDir);

    // rst: null space might get scaled here ... do we care. We could just inject at the cpoints, but I don't
    //  feel that this is needed.
    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(P->getDomainMap(), fineNullspace->getNumVectors());
    P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(), Teuchos::ScalarTraits<SC>::zero());
    Set(coarseLevel, "Nullspace", coarseNullspace);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MakeGeneralGeometricP(LO const ndm, const Array<LO> mapDirL2G, const Array<LO> mapDirG2L, const Array<LO> lFineNodesPerDir, const Array<LO> lCoarseNodesPerDir, Array<GO> gCoarseNodesPerDir, Array<GO> gFineNodesPerDir, ArrayRCP<LO> const coarseRate, LO const endRate[3], LO const offsets[6], bool const ghostInterface[6], const RCP<Xpetra::MultiVector<double,LO,GO,Node> >& fineCoords, LO const nnzP, LO const dofsPerNode, RCP<const Map>& stridedDomainMapP, RCP<Matrix> & Amat, RCP<Matrix>& P, RCP<Xpetra::MultiVector<double,LO,GO,Node> >& coarseCoords, Array<GO> ghostsGIDs, int interpolationOrder) const {

    /*
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     *
     * On termination, return the number of local prolongator columns owned by
     * this processor.
     *
     * Input
     * =====
     *    nNodes       Number of fine level Blk Rows owned by this processor
     *    coarseRate   Rate of coarsening in each spatial direction.
     *    endRate      Rate of coarsening in each spatial direction for the last
     *                 nodes in the mesh where an adaptive coarsening rate is
     *                 required.
     *    nTerms       Number of nonzero entries in the prolongation matrix.
     *    dofsPerNode  Number of degrees-of-freedom per mesh node.
     *
     * Output
     * =====
     *    So far nothing...
     */

    using xdMV = Xpetra::MultiVector<double,LO,GO,NO>;

    LO lNumFineNodes    = lFineNodesPerDir[0]*lFineNodesPerDir[1]*lFineNodesPerDir[2];
    LO lNumCoarseNodes  = lCoarseNodesPerDir[0]*lCoarseNodesPerDir[1]*lCoarseNodesPerDir[2];
    GO gNumCoarseNodes  = gCoarseNodesPerDir[0]*gCoarseNodesPerDir[1]*gCoarseNodesPerDir[2];
    LO lNumGhostNodes   = ghostsGIDs.size();
    GO numGloCols       = dofsPerNode*gNumCoarseNodes;

    // Build the required column map for the prolongator operator.
    // This requies to find the GIDs of the coarse nodes on the coarse mesh,
    // including the ghost nodes, on the local partition.
    // Note: the ordering of the coarse nodes is assumed
    // to be the same as the ordering of the fine nodes.

    GO gStartIndices[3];
    RCP<const Map> fineCoordsMap = fineCoords->getMap();
    // Compute the global indices of the first node on the partition.
    { // Scope for dummy
      gStartIndices[2] = fineCoordsMap->getMinGlobalIndex() / (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      GO dummy         = fineCoordsMap->getMinGlobalIndex() % (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      gStartIndices[1] = dummy / gFineNodesPerDir[0];
      gStartIndices[0] = dummy % gFineNodesPerDir[0];
    }

    Array<GO> gCoarseNodesGIDs(lNumCoarseNodes);
    LO currentNode, offset2, offset1, offset0;
    // Find the GIDs of the coarse nodes on the partition.
    for(LO ind2 = 0; ind2 < lCoarseNodesPerDir[mapDirL2G[2]]; ++ind2) {
      if(offsets[2] == 0) {
        offset2 = gStartIndices[2];
      } else {
        if(gStartIndices[2] + endRate[2] - offsets[2] == gFineNodesPerDir[2] - 1) {
          offset2 = gStartIndices[2] + endRate[2] - offsets[2];
        } else {
          offset2 = gStartIndices[2] + coarseRate[2] - offsets[2];
        }
      }
      if(offset2 + ind2*coarseRate[2] > gFineNodesPerDir[2] - 1) {
        offset2 += (ind2 - 1)*coarseRate[2] + endRate[2];
      } else {
        offset2 += ind2*coarseRate[2];
      }
      offset2 = offset2*gFineNodesPerDir[1]*gFineNodesPerDir[0];

      for(LO ind1 = 0; ind1 < lCoarseNodesPerDir[mapDirL2G[1]]; ++ind1) {
        if(offsets[1] == 0) {
          offset1 = gStartIndices[1];
        } else {
          if(gStartIndices[1] + endRate[1] - offsets[1] == gFineNodesPerDir[1] - 1) {
            offset1 = gStartIndices[1] + endRate[1] - offsets[1];
          } else {
            offset1 = gStartIndices[1] + coarseRate[1] - offsets[1];
          }
        }
        if(offset1 + ind1*coarseRate[1] > gFineNodesPerDir[1] - 1) {
          offset1 += (ind1 - 1)*coarseRate[1] + endRate[1];
        } else {
          offset1 += ind1*coarseRate[1];
        }
        offset1 = offset1*gFineNodesPerDir[0];
        for(LO ind0 = 0; ind0 < lCoarseNodesPerDir[mapDirL2G[0]]; ++ind0) {
          offset0 = gStartIndices[0] - offsets[0];
          if(offsets[0] == 0) {
            offset0 += ind0*coarseRate[0];
          } else {
            offset0 += (ind0 + 1)*coarseRate[0];
          }
          if(offset0 > gFineNodesPerDir[0] - 1) {offset0 += endRate[0] - coarseRate[0];}

          currentNode = ind2*lCoarseNodesPerDir[mapDirL2G[1]]*lCoarseNodesPerDir[mapDirL2G[0]]
                      + ind1*lCoarseNodesPerDir[mapDirL2G[0]]
                      + ind0;
          gCoarseNodesGIDs[currentNode] = offset2 + offset1 + offset0;
        }
      }
    }

    // Actual loop over all the coarse/ghost nodes to find their index on the coarse mesh
    // and the corresponding dofs that will need to be added to colMapP.
    Array<GO> colGIDs(dofsPerNode*(lNumCoarseNodes+lNumGhostNodes));
    Array<GO> coarseNodesGIDs(lNumCoarseNodes);
    GO fineNodesPerCoarseSlab    = coarseRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
    GO fineNodesEndCoarseSlab    = endRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
    GO fineNodesPerCoarsePlane   = coarseRate[1]*gFineNodesPerDir[0];
    GO fineNodesEndCoarsePlane   = endRate[1]*gFineNodesPerDir[0];
    GO coarseNodesPerCoarseLayer = gCoarseNodesPerDir[1]*gCoarseNodesPerDir[0];
    GO gCoarseNodeOnCoarseGridGID;
    LO gInd[3], lCol;
    Array<int> ghostsPIDs  (lNumGhostNodes);
    Array<LO>  ghostsLIDs  (lNumGhostNodes);
    Array<LO>  ghostsPermut(lNumGhostNodes);
    for(LO k = 0; k < lNumGhostNodes; ++k) {ghostsPermut[k] = k;}
    fineCoordsMap->getRemoteIndexList(ghostsGIDs, ghostsPIDs, ghostsLIDs);
    sh_sort_permute(ghostsPIDs.begin(),ghostsPIDs.end(), ghostsPermut.begin(),ghostsPermut.end());

    { // scope for tmpInds, tmpVars and tmp.
      GO tmpInds[3], tmpVars[2];
      LO tmp;
      // Loop over the coarse nodes of the partition and add them to colGIDs
      // that will be used to construct the column and domain maps of P as well
      // as to construct the coarse coordinates map.
      for(LO col = 0; col < lNumCoarseNodes; ++col) {
        if((endRate[2] != coarseRate[2]) && (gCoarseNodesGIDs[col] > (gCoarseNodesPerDir[2] - 2)*fineNodesPerCoarseSlab + fineNodesEndCoarseSlab - 1)) {
          tmpInds[2] = gCoarseNodesGIDs[col] / fineNodesPerCoarseSlab + 1;
          tmpVars[0] = gCoarseNodesGIDs[col] - (tmpInds[2] - 1)*fineNodesPerCoarseSlab - fineNodesEndCoarseSlab;
        } else {
          tmpInds[2] = gCoarseNodesGIDs[col] / fineNodesPerCoarseSlab;
          tmpVars[0] = gCoarseNodesGIDs[col] % fineNodesPerCoarseSlab;
        }
        if((endRate[1] != coarseRate[1]) && (tmpVars[0] > (gCoarseNodesPerDir[1] - 2)*fineNodesPerCoarsePlane + fineNodesEndCoarsePlane - 1)) {
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane + 1;
          tmpVars[1] = tmpVars[0] - (tmpInds[1] - 1)*fineNodesPerCoarsePlane - fineNodesEndCoarsePlane;
        } else {
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane;
          tmpVars[1] = tmpVars[0] % fineNodesPerCoarsePlane;
        }
        if(tmpVars[1] == gFineNodesPerDir[0] - 1) {
          tmpInds[0] = gCoarseNodesPerDir[0] - 1;
        } else {
          tmpInds[0] = tmpVars[1] / coarseRate[0];
        }
        gInd[2] = col / (lCoarseNodesPerDir[mapDirG2L[1]]*lCoarseNodesPerDir[mapDirG2L[0]]);
        tmp     = col % (lCoarseNodesPerDir[mapDirG2L[1]]*lCoarseNodesPerDir[mapDirG2L[0]]);
        gInd[1] = tmp / lCoarseNodesPerDir[mapDirG2L[0]];
        gInd[0] = tmp % lCoarseNodesPerDir[mapDirG2L[0]];
        lCol = gInd[mapDirG2L[2]]*(lCoarseNodesPerDir[1]*lCoarseNodesPerDir[0]) + gInd[mapDirG2L[1]]*lCoarseNodesPerDir[0] + gInd[mapDirG2L[0]];
        gCoarseNodeOnCoarseGridGID = tmpInds[2]*coarseNodesPerCoarseLayer + tmpInds[1]*gCoarseNodesPerDir[0] + tmpInds[0];
        coarseNodesGIDs[lCol] = gCoarseNodeOnCoarseGridGID;
        for(LO dof = 0; dof < dofsPerNode; ++dof) {
          colGIDs[dofsPerNode*lCol + dof] = dofsPerNode*gCoarseNodeOnCoarseGridGID + dof;
        }
      }
      // Now loop over the ghost nodes of the partition to add them to colGIDs
      // since they will need to be included in the column map of P
      for(LO col = lNumCoarseNodes; col < lNumCoarseNodes + lNumGhostNodes; ++col) {
        if((endRate[2] != coarseRate[2]) && (ghostsGIDs[ghostsPermut[col - lNumCoarseNodes]] > (gCoarseNodesPerDir[2] - 2)*fineNodesPerCoarseSlab + fineNodesEndCoarseSlab - 1)) {
          tmpInds[2] = ghostsGIDs[ghostsPermut[col - lNumCoarseNodes]] / fineNodesPerCoarseSlab + 1;
          tmpVars[0] = ghostsGIDs[ghostsPermut[col - lNumCoarseNodes]] - (tmpInds[2] - 1)*fineNodesPerCoarseSlab - fineNodesEndCoarseSlab;
        } else {
          tmpInds[2] = ghostsGIDs[ghostsPermut[col - lNumCoarseNodes]] / fineNodesPerCoarseSlab;
          tmpVars[0] = ghostsGIDs[ghostsPermut[col - lNumCoarseNodes]] % fineNodesPerCoarseSlab;
        }
        if((endRate[1] != coarseRate[1]) && (tmpVars[0] > (gCoarseNodesPerDir[1] - 2)*fineNodesPerCoarsePlane + fineNodesEndCoarsePlane - 1)) {
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane + 1;
          tmpVars[1] = tmpVars[0] - (tmpInds[1] - 1)*fineNodesPerCoarsePlane - fineNodesEndCoarsePlane;
        } else {
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane;
          tmpVars[1] = tmpVars[0] % fineNodesPerCoarsePlane;
        }
        if(tmpVars[1] == gFineNodesPerDir[0] - 1) {
          tmpInds[0] = gCoarseNodesPerDir[0] - 1;
        } else {
          tmpInds[0] = tmpVars[1] / coarseRate[0];
        }
        gCoarseNodeOnCoarseGridGID = tmpInds[2]*coarseNodesPerCoarseLayer + tmpInds[1]*gCoarseNodesPerDir[0] + tmpInds[0];
        for(LO dof = 0; dof < dofsPerNode; ++dof) {
          colGIDs[dofsPerNode*col + dof] = dofsPerNode*gCoarseNodeOnCoarseGridGID + dof;
        }
      }
    } // End of scope for tmpInds, tmpVars and tmp

    // Build maps necessary to create and fill complete the prolongator
    // note: rowMapP == rangeMapP and colMapP != domainMapP
    RCP<const Map> rowMapP = Amat->getDomainMap();

    RCP<const Map> domainMapP = Xpetra::MapFactory<LO,GO,NO>::Build(rowMapP->lib(),
                                                                    numGloCols,
                                                                    colGIDs.view(0, dofsPerNode*lNumCoarseNodes),
                                                                    rowMapP->getIndexBase(),
                                                                    rowMapP->getComm());

    RCP<const Map> colMapP = Xpetra::MapFactory<LO,GO,NO>::Build(rowMapP->lib(),
                                                                 Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                 colGIDs.view(0, colGIDs.size()),
                                                                 rowMapP->getIndexBase(),
                                                                 rowMapP->getComm());

    std::vector<size_t> strideInfo(1);
    strideInfo[0] = dofsPerNode;
    stridedDomainMapP = Xpetra::StridedMapFactory<LO,GO,NO>::Build(domainMapP, strideInfo);

    // Build the map for the coarse level coordinates
    RCP<const Map> coarseCoordsMap = MapFactory::Build (fineCoordsMap->lib(),
                                                        gNumCoarseNodes,
                                                        coarseNodesGIDs.view(0, lNumCoarseNodes),
                                                        fineCoordsMap->getIndexBase(),
                                                        rowMapP->getComm());

    // Do the actual import using the fineCoordsMap
    RCP<const Map> ghostMap = Xpetra::MapFactory<LO,GO,NO>::Build(fineCoordsMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), ghostsGIDs.view(0, lNumGhostNodes),
                                                               fineCoordsMap->getIndexBase(), rowMapP->getComm());
    RCP<const Import> importer = ImportFactory::Build(fineCoordsMap, ghostMap);
    RCP<xdMV> ghostCoords = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(ghostMap, ndm);
    ghostCoords->doImport(*fineCoords, *importer, Xpetra::INSERT);

    P = rcp(new CrsMatrixWrap(rowMapP, colMapP, 0, Xpetra::StaticProfile));
    RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

    ArrayRCP<size_t>  iaP;
    ArrayRCP<LO>      jaP;
    ArrayRCP<SC>     valP;

    PCrs->allocateAllValues(nnzP, iaP, jaP, valP);

    ArrayView<size_t> ia  = iaP();
    ArrayView<LO>     ja  = jaP();
    ArrayView<SC>     val = valP();

    LO nStencil = 0; LO tStencil = 0;
    LO firstCoarseNodeIndex;
    LO indices[4][3];
    ia[0] = 0;

    // Declaration and assignment of fineCoords which holds the coordinates of the fine nodes in 3D.
    // To do so we pull the nD coordinates from fineCoords and pad the rest with zero vectors...
    // We also create an nD array that will store the coordinates of the coarse grid nodes.
    // That array will eventually be used during the construction of coarseCoords, the MultiVector of
    // coarse nodes coordinates that we store on the coarse level.
    ArrayRCP< ArrayRCP<double> > lFineCoords(3);
    ArrayRCP< ArrayRCP<double> > lGhostCoords(3);
    RCP<const Map> ghostCoordsMap = ghostCoords->getMap();
    RCP<Xpetra::Vector<double, LO, GO, NO> > zeros = Xpetra::VectorFactory<double, LO, GO, NO>::Build(fineCoordsMap, true);
    RCP<Xpetra::Vector<double, LO, GO, NO> > ghostZeros = Xpetra::VectorFactory<double, LO, GO, NO>::Build(ghostCoordsMap, true);

    // Build the MultiVector holding the coarse grid points coordinates.
    coarseCoords = Xpetra::MultiVectorFactory<double,LO,GO,Node>::Build(coarseCoordsMap, Teuchos::as<size_t>(ndm));
    ArrayRCP<double> xCoarseNodes; ArrayRCP<double> yCoarseNodes; ArrayRCP<double> zCoarseNodes;

    if(ndm==1) {
      lFineCoords[0] = fineCoords->getDataNonConst(0);
      lFineCoords[1] = zeros->getDataNonConst(0);
      lFineCoords[2] = zeros->getDataNonConst(0);
      lGhostCoords[0] = ghostCoords->getDataNonConst(0);
      lGhostCoords[1] = ghostZeros->getDataNonConst(0);
      lGhostCoords[2] = ghostZeros->getDataNonConst(0);

      xCoarseNodes = coarseCoords->getDataNonConst(0);
    } else if(ndm==2) {
      lFineCoords[0] = fineCoords->getDataNonConst(0);
      lFineCoords[1] = fineCoords->getDataNonConst(1);
      lFineCoords[2] = zeros->getDataNonConst(0);
      lGhostCoords[0] = ghostCoords->getDataNonConst(0);
      lGhostCoords[1] = ghostCoords->getDataNonConst(1);
      lGhostCoords[2] = ghostZeros->getDataNonConst(0);

      xCoarseNodes= coarseCoords->getDataNonConst(0);
      yCoarseNodes= coarseCoords->getDataNonConst(1);
    } else if(ndm==3) {
      lFineCoords[0] = fineCoords->getDataNonConst(0);
      lFineCoords[1] = fineCoords->getDataNonConst(1);
      lFineCoords[2] = fineCoords->getDataNonConst(2);
      lGhostCoords[0] = ghostCoords->getDataNonConst(0);
      lGhostCoords[1] = ghostCoords->getDataNonConst(1);
      lGhostCoords[2] = ghostCoords->getDataNonConst(2);

      xCoarseNodes = coarseCoords->getDataNonConst(0);
      yCoarseNodes = coarseCoords->getDataNonConst(1);
      zCoarseNodes = coarseCoords->getDataNonConst(2);
    }

    // Finally let's load the PIDs of the coarse nodes on the coarse grid
    // to later perform permutation of stencil entries to conform with
    // Tpetra matrix ordering when calling setAllValues.

    // Loop over the fine nodes and compute their interpolation stencils
    // Get the rank of the local process to order entries by PIDs in the
    // prolongator matrix rows
    GO currentCoarseNode = 0;
    for(LO i = 0; i < lNumFineNodes; ++i) {

      // Get point indices on fine grid
      {
        std::div_t tmp;
        tmp=std::div(i,lFineNodesPerDir[0]*lFineNodesPerDir[1]);
        indices[0][2] = tmp.quot;
        tmp=std::div(tmp.rem,lFineNodesPerDir[0]);
        indices[0][1] = tmp.quot;
        indices[0][0] = tmp.rem;

        // Get ref point indices on coarse grid
        tmp=std::div(indices[0][0] - offsets[0],coarseRate[0]);
        indices[1][0] = tmp.quot;
        tmp=std::div(indices[0][1] - offsets[1],coarseRate[1]);
        indices[1][1] = tmp.quot;
        tmp=std::div(indices[0][2] - offsets[2],coarseRate[2]);
        indices[1][2] = tmp.quot;
      }

      // location "flags" indicate if the current node is on a coarse
      // face, edge or node.
      indices[2][0] = (indices[0][0] + offsets[0]) % coarseRate[0];
      indices[2][1] = (indices[0][1] + offsets[1]) % coarseRate[1];
      indices[2][2] = (indices[0][2] + offsets[2]) % coarseRate[2];

      // Get indices of ref point on fine grid
      indices[3][0] = indices[1][0]*coarseRate[0];
      indices[3][1] = indices[1][1]*coarseRate[1];
      indices[3][2] = indices[1][2]*coarseRate[2];

      if( (indices[0][0] == lFineNodesPerDir[0]-1) && !ghostInterface[1] ) {
        indices[1][0] = lCoarseNodesPerDir[0]-1;
        indices[2][0] = 0;
        indices[3][0] = lFineNodesPerDir[0]-1;
      }
      if( (indices[0][1] == lFineNodesPerDir[1]-1) && !ghostInterface[3] ) {
        indices[1][1] = lCoarseNodesPerDir[1]-1;
        indices[2][1] = 0;
        indices[3][1] = lFineNodesPerDir[1]-1;
      }
      if( (indices[0][2] == lFineNodesPerDir[2]-1) && !ghostInterface[5] ) {
        indices[1][2] = lCoarseNodesPerDir[2]-1;
        indices[2][2] = 0;
        indices[3][2] = lFineNodesPerDir[2]-1;
      }

      Array<GO> currentNodeIndices(3);
      currentNodeIndices[0] = gStartIndices[0] + Teuchos::as<GO>(indices[0][mapDirL2G[0]]);
      currentNodeIndices[1] = gStartIndices[1] + Teuchos::as<GO>(indices[0][mapDirL2G[1]]);
      currentNodeIndices[2] = gStartIndices[2] + Teuchos::as<GO>(indices[0][mapDirL2G[2]]);

      // Assuming lexicographic indexing the coarse nodes forming a prism
      // around fine node "i" are selected and store them in connec.
      // Two tricky things to be careful about:
      //    - are we using coarseRate or endRate?
      //      --> check indices and set rate correctly
      //    - are we on the east, north or top face?
      //      --> if so fix firstCoarseNodeIndex to make sure there is
      //          a node to the east, north and top of firstCoarseNodeIndex
      LO rate[3];
      if(currentNodeIndices[0] >= gFineNodesPerDir[0] - endRate[0] - 1) {
        rate[0] = endRate[0];
      } else {
        rate[0] = coarseRate[0];
      }
      if(currentNodeIndices[1] >= gFineNodesPerDir[1] - endRate[1] - 1) {
        rate[1] = endRate[1];
      } else {
        rate[1] = coarseRate[1];
      }
      if(currentNodeIndices[2] >= gFineNodesPerDir[2] - endRate[2] - 1) {
        rate[2] = endRate[2];
      } else {
        rate[2] = coarseRate[2];
      }
      if(ndm < 3) { rate[2] = 0;}
      if(ndm < 2) { rate[1] = 0;}

      // We need to check whether we are on the edge of the mesh in which case we need to adjust the coarse nodes
      firstCoarseNodeIndex = 0;
      Array<GO> firstCoarseNodeIndices(3); // These are fine grid indices, divide by coarseRate[i] to get coarse grid indices
      if((currentNodeIndices[2] == gFineNodesPerDir[2] -1) && (endRate[2] == coarseRate[2])) {
        // If we are on the last node and have a endRate == coarseRate we need to take the coarse node below the current node
        firstCoarseNodeIndices[2] = ((currentNodeIndices[2] / coarseRate[2]) - 1) * coarseRate[2];
      } else {
        firstCoarseNodeIndices[2] = (currentNodeIndices[2] / coarseRate[2]) * coarseRate[2];
      }
      if((currentNodeIndices[1] == gFineNodesPerDir[1] -1) && (endRate[1] == coarseRate[1])) {
        firstCoarseNodeIndices[1] = ((currentNodeIndices[1] / coarseRate[1]) - 1) * coarseRate[1];
      } else {
        firstCoarseNodeIndices[1] = (currentNodeIndices[1] / coarseRate[1]) * coarseRate[1];
      }
      if((currentNodeIndices[0] == gFineNodesPerDir[0] -1) && (endRate[0] == coarseRate[0])) {
        firstCoarseNodeIndices[0] = ((currentNodeIndices[0] / coarseRate[0]) - 1) * coarseRate[0];
      } else {
        firstCoarseNodeIndices[0] = (currentNodeIndices[0] / coarseRate[0]) * coarseRate[0];
      }
      firstCoarseNodeIndex += firstCoarseNodeIndices[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0]
        + firstCoarseNodeIndices[1]*gFineNodesPerDir[0] + firstCoarseNodeIndices[0];

      GO firstCoarseNodeOnCoarseGridIndex;
      {
        GO tmpInds[4];
        tmpInds[2] = firstCoarseNodeIndex / fineNodesPerCoarseSlab;
        tmpInds[3] = firstCoarseNodeIndex % fineNodesPerCoarseSlab;
        tmpInds[1] = tmpInds[3] / fineNodesPerCoarsePlane;
        tmpInds[0] = (tmpInds[3] % fineNodesPerCoarsePlane) / coarseRate[0];
        firstCoarseNodeOnCoarseGridIndex = tmpInds[2]*coarseNodesPerCoarseLayer + tmpInds[1]*gCoarseNodesPerDir[0] + tmpInds[0];
      }

      LO coarseGhosts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      if( ghostInterface[0] ) {
        if( ((indices[0][mapDirG2L[0]] < rate[0] - offsets[0]) && offsets[0] != 0)
            || (((currentNodeIndices[0] == gFineNodesPerDir[0] -1) && (indices[0][mapDirG2L[0]] < rate[0] - offsets[0] + 1)) && offsets[0] != 0)
            || ((currentNodeIndices[0] == gFineNodesPerDir[0] -1) && lFineNodesPerDir[mapDirG2L[0]] == 1) ) {
          coarseGhosts[0] = 1;
          coarseGhosts[2] = 1;
          coarseGhosts[4] = 1;
          coarseGhosts[6] = 1;
        }
      }
      if( ghostInterface[1] && (indices[0][mapDirG2L[0]] > lFineNodesPerDir[mapDirG2L[0]] - offsets[3] - 2) ) {
        coarseGhosts[1] = 1;
        coarseGhosts[3] = 1;
        coarseGhosts[5] = 1;
        coarseGhosts[7] = 1;
      }
      if( ghostInterface[2] ) {
        if( ((indices[0][mapDirG2L[1]] < rate[1] - offsets[1]) && offsets[1] != 0)
            || (((currentNodeIndices[1] == gFineNodesPerDir[1] -1) && (indices[0][mapDirG2L[1]] < rate[1] - offsets[1] + 1)) && offsets[1] != 0)
            || ((currentNodeIndices[1] == gFineNodesPerDir[1] -1) && lFineNodesPerDir[mapDirG2L[1]] == 1) ) {
          coarseGhosts[0] = 1;
          coarseGhosts[1] = 1;
          coarseGhosts[4] = 1;
          coarseGhosts[5] = 1;
        }
      }
      if( ghostInterface[3] && (indices[0][mapDirG2L[1]] > lFineNodesPerDir[mapDirG2L[1]] - offsets[4] - 2) ) {
        coarseGhosts[2] = 1;
        coarseGhosts[3] = 1;
        coarseGhosts[6] = 1;
        coarseGhosts[7] = 1;
      }
      if( ghostInterface[4] ) {
        if( ((indices[0][mapDirG2L[2]] < rate[2] - offsets[2]) && offsets[2] != 0)
            || (((currentNodeIndices[2] == gFineNodesPerDir[2] -1) && (indices[0][mapDirG2L[2]] < rate[2] - offsets[2] + 1)) && offsets[2] != 0)
            || ((currentNodeIndices[2] == gFineNodesPerDir[2] -1) && lFineNodesPerDir[mapDirG2L[2]] == 1) ) {
          coarseGhosts[0] = 1;
          coarseGhosts[1] = 1;
          coarseGhosts[2] = 1;
          coarseGhosts[3] = 1;
        }
      }
      if( ghostInterface[5] && (indices[0][mapDirG2L[2]] > lFineNodesPerDir[mapDirG2L[2]] - offsets[5] - 2) ) {
        coarseGhosts[4] = 1;
        coarseGhosts[5] = 1;
        coarseGhosts[6] = 1;
        coarseGhosts[7] = 1;
      }

      GO firstGhostNodeIndices[3], firstGhostNodeIndex;
      if(currentNodeIndices[0] == gFineNodesPerDir[0] - 1) {
        firstGhostNodeIndices[0] = (currentNodeIndices[0]-rate[0]) - (currentNodeIndices[0]-rate[0])%coarseRate[0];
      } else {
        firstGhostNodeIndices[0] = currentNodeIndices[0] - currentNodeIndices[0]%coarseRate[0];
      }
      if(currentNodeIndices[1] == gFineNodesPerDir[1] - 1) {
        firstGhostNodeIndices[1] = (currentNodeIndices[1]-rate[1]) - (currentNodeIndices[1]-rate[1])%coarseRate[1];
      } else {
        firstGhostNodeIndices[1] = currentNodeIndices[1] - currentNodeIndices[1]%coarseRate[1];
      }
      if(currentNodeIndices[2] == gFineNodesPerDir[2] - 1) {
        firstGhostNodeIndices[2] = (currentNodeIndices[2]-rate[2]) - (currentNodeIndices[2]-rate[2])%coarseRate[2];
      } else {
        firstGhostNodeIndices[2] = currentNodeIndices[2] - currentNodeIndices[2]%coarseRate[2];
      }
      firstGhostNodeIndex = firstGhostNodeIndices[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0] + firstGhostNodeIndices[1]*gFineNodesPerDir[0] + firstGhostNodeIndices[0];

      Array<GO> gCoarseNodesOnCoarseGridIndices(8);
      GO ghostNodeIndex, ghostNodeOnCoarseGridIndex;
      GO coarseNodeIndex, coarseNodeOnCoarseGridIndex;
      LO ind;
      double connec[9][3];
      LO currentGhostLID;
      Array<LO> connecPIDs(8); // A PID of -1 indicate a local element.
      Array<GO> connecLGIDs(8);
      // First load the coordinates of the node of interest
      for(LO dim = 0; dim < 3; ++dim) {connec[0][dim] = lFineCoords[dim][i];}
      // Then loop over the three spatial dimension and load the 8 coarse nodes
      // required for the linear interpolation.
      for(LO ind2 = 0; ind2 < 2; ++ind2) {
        for(LO ind1 = 0; ind1 < 2; ++ind1) {
          for(LO ind0 = 0; ind0 < 2; ++ind0) {
            ind = 4*ind2 + 2*ind1 + ind0;
            // Check whether ghost nodes are needed for the current fine point.
            if(coarseGhosts[ind] == 1) {
              // Get the global ghost node index and load its coordinates
              ghostNodeIndex = firstGhostNodeIndex + Teuchos::as<GO>(ind2*rate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0]
                                                                     + ind1*rate[1]*gFineNodesPerDir[0]
                                                                     + ind0*rate[0]);
              ghostNodeOnCoarseGridIndex = firstCoarseNodeOnCoarseGridIndex + Teuchos::as<GO>(ind2*gCoarseNodesPerDir[1]*gCoarseNodesPerDir[0]
                                                                                              + ind1*gCoarseNodesPerDir[0]
                                                                                              + ind0);
              currentGhostLID  = ghostMap->getLocalElement(ghostNodeIndex);
              connecPIDs[ind]  = ghostsPIDs[currentGhostLID];
              connecLGIDs[ind] = ghostNodeIndex;
              for(LO dim = 0; dim < 3; ++dim) {connec[ind + 1][dim] = lGhostCoords[dim][currentGhostLID];}
              gCoarseNodesOnCoarseGridIndices[4*ind2 + 2*ind1 + ind0] = ghostNodeOnCoarseGridIndex;
            } else {
              // Get the local coarse node index and load its coordinates
              coarseNodeIndex = firstCoarseNodeIndex + Teuchos::as<GO>(ind2*rate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0]
                                                                     + ind1*rate[1]*gFineNodesPerDir[0]
                                                                     + ind0*rate[0]);
              coarseNodeOnCoarseGridIndex = firstCoarseNodeOnCoarseGridIndex + Teuchos::as<GO>(ind2*gCoarseNodesPerDir[1]*gCoarseNodesPerDir[0]
                                                                                               + ind1*gCoarseNodesPerDir[0]
                                                                                               + ind0);
              for(LO dim = 0; dim < 3; ++dim) {connec[ind + 1][dim] = lFineCoords[dim][fineCoordsMap->getLocalElement(coarseNodeIndex)];}
              gCoarseNodesOnCoarseGridIndices[4*ind2 + 2*ind1 + ind0] = coarseNodeOnCoarseGridIndex;
              connecPIDs[ind]  = -1;
              connecLGIDs[ind] = coarseCoordsMap->getLocalElement(coarseNodeOnCoarseGridIndex);
            }
          }
        }
      }

      // Compute the actual geometric interpolation stencil
      SC stencil[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      ComputeStencil(ndm, currentNodeIndices, firstCoarseNodeIndices, rate, connec, interpolationOrder, stencil);

      // Finally check whether the fine node is on a coarse: node, edge or face
      // and select accordingly the non-zero values from the stencil and the
      // corresponding column indices
      LO nzIndStencil[8] = {0,0,0,0,0,0,0,0};
      if( ((currentNodeIndices[0] % coarseRate[0] == 0) || currentNodeIndices[0] == gFineNodesPerDir[0] - 1)
          && ((currentNodeIndices[1] % coarseRate[1] == 0) || currentNodeIndices[1] == gFineNodesPerDir[1] - 1)
          && ((currentNodeIndices[2] % coarseRate[2] == 0) || currentNodeIndices[2] == gFineNodesPerDir[2] - 1) ) {
        if(ndm==1) {
          xCoarseNodes[currentCoarseNode] = lFineCoords[0][i];
        } else if(ndm==2) {
          xCoarseNodes[currentCoarseNode] = lFineCoords[0][i];
          yCoarseNodes[currentCoarseNode] = lFineCoords[1][i];
        } else if(ndm==3) {
          xCoarseNodes[currentCoarseNode] = lFineCoords[0][i];
          yCoarseNodes[currentCoarseNode] = lFineCoords[1][i];
          zCoarseNodes[currentCoarseNode] = lFineCoords[2][i];
        }
        ++currentCoarseNode;

        if(currentNodeIndices[0] == gFineNodesPerDir[0] - 1) {
          nzIndStencil[0] += 1;
        }
        if((currentNodeIndices[1] == gFineNodesPerDir[1] - 1) && (ndm > 1)) {
          nzIndStencil[0] += 2;
        }
        if((currentNodeIndices[2] == gFineNodesPerDir[2] - 1) && (ndm > 2)) {
          nzIndStencil[0] += 4;
        }
        nStencil = 1;
        for(LO k = 0; k < 8; ++k) {
          nzIndStencil[k] = nzIndStencil[0];
        }
      } else {
        nStencil = 8;
        for(LO k = 0; k < 8; ++k) {
          nzIndStencil[k] = k;
        }
      }

      // Here the values are filled in the Crs matrix arrays
      // This is basically the only place these variables are modified
      // Hopefully this makes handling system of PDEs easy!

      // Loop on dofsPerNode and process each row for the current Node

      Array<LO> permutation(8);
      for(LO k = 0; k < 8; ++k) {permutation[k] = k;}

      // Sort nodes by PIDs using stable sort to keep sublist ordered by LIDs and GIDs
      sh_sort2(connecPIDs.begin(),connecPIDs.end(), permutation.begin(), permutation.end());

      for(LO j = 0; j < dofsPerNode; ++j) {
        ia[i*dofsPerNode + j + 1] = ia[i*dofsPerNode + j] + nStencil;
        for(LO k = 0; k < nStencil; ++k) {
          ja [ia[i*dofsPerNode + j] + k] = colMapP->getLocalElement(gCoarseNodesOnCoarseGridIndices[nzIndStencil[permutation[k]]]*dofsPerNode + j);
          val[ia[i*dofsPerNode + j] + k] = stencil[nzIndStencil[permutation[k]]];
        }
        // Add the stencil for each degree of freedom.
        tStencil += nStencil;
      }
    } // End loop over fine nodes

    if (rowMapP->lib() == Xpetra::UseTpetra) {
      // - Cannot resize for Epetra, as it checks for same pointers
      // - Need to resize for Tpetra, as it check ().size() == ia[numRows]
      // NOTE: these invalidate ja and val views
      jaP .resize(tStencil);
      valP.resize(tStencil);
    }

    // Set the values of the prolongation operators into the CrsMatrix P and call FillComplete
    PCrs->setAllValues(iaP, jaP, valP);
    PCrs->expertStaticFillComplete(domainMapP,rowMapP);

  } // MakeGeneralGeometricP

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeStencil(
                                const LO numDimensions, const Array<GO> currentNodeIndices,
                                const Array<GO> coarseNodeIndices, const LO rate[3],
                                const double coord[9][3], const int interpolationOrder,
                                SC stencil[8]) const {

    TEUCHOS_TEST_FOR_EXCEPTION((interpolationOrder > 1) || (interpolationOrder < 0),
                               Exceptions::RuntimeError,
                               "The interpolation order can be set to 0 or 1 only.");

    if(interpolationOrder == 0) {
      ComputeConstantInterpolationStencil(numDimensions, currentNodeIndices, coarseNodeIndices, rate, stencil);
    } else if(interpolationOrder == 1) {
      ComputeLinearInterpolationStencil(numDimensions, coord, stencil);
    }

  } // End ComputeStencil

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeConstantInterpolationStencil(
                                const LO numDimensions, const Array<GO> currentNodeIndices,
                                const Array<GO> coarseNodeIndices, const LO rate[3], SC stencil[8]) const {

    LO coarseNode = 0;
    if(numDimensions > 2) {
      if((currentNodeIndices[2] - coarseNodeIndices[2]) > (rate[2] / 2)) {
        coarseNode += 4;
      }
    }
    if(numDimensions > 1) {
      if((currentNodeIndices[1] - coarseNodeIndices[1]) > (rate[1] / 2)) {
        coarseNode += 2;
      }
    }
    if((currentNodeIndices[0] - coarseNodeIndices[0]) > (rate[0] / 2)) {
      coarseNode += 1;
    }
    stencil[coarseNode] = 1.0;

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeLinearInterpolationStencil(
                                const LO numDimensions, const double coord[9][3], SC stencil[8]) const {

    //                7         8                Find xi, eta and zeta such that
    //                x---------x
    //               /|        /|             Rx = x_p - sum N_i(xi,eta,zeta)x_i = 0
    //             5/ |      6/ |             Ry = y_p - sum N_i(xi,eta,zeta)y_i = 0
    //             x---------x  |             Rz = z_p - sum N_i(xi,eta,zeta)z_i = 0
    //             |  | *P   |  |
    //             |  x------|--x             We can do this with a Newton solver:
    //             | /3      | /4             We will start with initial guess (xi,eta,zeta) = (0,0,0)
    //             |/        |/               We compute the Jacobian and iterate until convergence...
    //  z  y       x---------x
    //  | /        1         2                Once we have (xi,eta,zeta), we can evaluate all N_i which
    //  |/                                    give us the weights for the interpolation stencil!
    //  o---x
    //

    Teuchos::SerialDenseMatrix<LO,double> Jacobian(numDimensions, numDimensions);
    Teuchos::SerialDenseVector<LO,double> residual(numDimensions);
    Teuchos::SerialDenseVector<LO,double> solutionDirection(numDimensions);
    Teuchos::SerialDenseVector<LO,double> paramCoords(numDimensions);
    Teuchos::SerialDenseSolver<LO,double> problem;
    LO numTerms = std::pow(2,numDimensions), iter = 0, max_iter = 5;
    double functions[4][8], norm_ref = 1, norm2 = 1, tol = 1e-5;
    paramCoords.size(numDimensions);

    while( (iter < max_iter) && (norm2 > tol*norm_ref) ) {
      ++iter;
      norm2 = 0;
      solutionDirection.size(numDimensions);
      residual.size(numDimensions);
      Jacobian = 0.0;

      // Compute Jacobian and Residual
      GetInterpolationFunctions(numDimensions, paramCoords, functions);
      for(LO i = 0; i < numDimensions; ++i) {
        residual(i) = coord[0][i];                 // Add coordinates from point of interest
        for(LO k = 0; k < numTerms; ++k) {
          residual(i) -= functions[0][k]*coord[k+1][i];  // Remove contribution from all coarse points
        }
        if(iter == 1) {
          norm_ref += residual(i)*residual(i);
          if(i == numDimensions - 1) {
            norm_ref = std::sqrt(norm_ref);
          }
        }

        for(LO j = 0; j < numDimensions; ++j) {
          for(LO k = 0; k < numTerms; ++k) {
            Jacobian(i,j) += functions[j+1][k]*coord[k+1][i];
          }
        }
      }

      // Set Jacobian, Vectors and solve problem
      problem.setMatrix(Teuchos::rcp(&Jacobian, false));
      problem.setVectors(Teuchos::rcp(&solutionDirection, false), Teuchos::rcp(&residual, false));
      problem.factorWithEquilibration(true);
      problem.solve();
      problem.unequilibrateLHS();

      for(LO i = 0; i < numDimensions; ++i) {
        paramCoords(i) = paramCoords(i) + solutionDirection(i);
      }

      // Recompute Residual norm
      GetInterpolationFunctions(numDimensions, paramCoords, functions); // These need to be recomputed with new paramCoords!
      for(LO i = 0; i < numDimensions; ++i) {
        double tmp = coord[0][i];
        for(LO k = 0; k < numTerms; ++k) {
          tmp -= functions[0][k]*coord[k+1][i];
        }
        norm2 += tmp*tmp;
        tmp = 0;
      }
      norm2 = std::sqrt(norm2);
    }

    // Load the interpolation values onto the stencil.
    for(LO i = 0; i < 8; ++i) {
      stencil[i] = functions[0][i];
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetInterpolationFunctions(const LO numDimensions, const Teuchos::SerialDenseVector<LO,double> parameters, double functions[4][8]) const {
    double xi = 0.0, eta = 0.0, zeta = 0.0, denominator = 0.0;
    if(numDimensions == 1) {
      xi = parameters[0];
      denominator = 2.0;
    } else if(numDimensions == 2) {
      xi  = parameters[0];
      eta = parameters[1];
      denominator = 4.0;
    } else if(numDimensions == 3) {
      xi   = parameters[0];
      eta  = parameters[1];
      zeta = parameters[2];
      denominator = 8.0;
    }

    functions[0][0] = (1.0 - xi)*(1.0 - eta)*(1.0 - zeta) / denominator;
    functions[0][1] = (1.0 + xi)*(1.0 - eta)*(1.0 - zeta) / denominator;
    functions[0][2] = (1.0 - xi)*(1.0 + eta)*(1.0 - zeta) / denominator;
    functions[0][3] = (1.0 + xi)*(1.0 + eta)*(1.0 - zeta) / denominator;
    functions[0][4] = (1.0 - xi)*(1.0 - eta)*(1.0 + zeta) / denominator;
    functions[0][5] = (1.0 + xi)*(1.0 - eta)*(1.0 + zeta) / denominator;
    functions[0][6] = (1.0 - xi)*(1.0 + eta)*(1.0 + zeta) / denominator;
    functions[0][7] = (1.0 + xi)*(1.0 + eta)*(1.0 + zeta) / denominator;

    functions[1][0] = -(1.0 - eta)*(1.0 - zeta) / denominator;
    functions[1][1] =  (1.0 - eta)*(1.0 - zeta) / denominator;
    functions[1][2] = -(1.0 + eta)*(1.0 - zeta) / denominator;
    functions[1][3] =  (1.0 + eta)*(1.0 - zeta) / denominator;
    functions[1][4] = -(1.0 - eta)*(1.0 + zeta) / denominator;
    functions[1][5] =  (1.0 - eta)*(1.0 + zeta) / denominator;
    functions[1][6] = -(1.0 + eta)*(1.0 + zeta) / denominator;
    functions[1][7] =  (1.0 + eta)*(1.0 + zeta) / denominator;

    functions[2][0] = -(1.0 - xi)*(1.0 - zeta) / denominator;
    functions[2][1] = -(1.0 + xi)*(1.0 - zeta) / denominator;
    functions[2][2] =  (1.0 - xi)*(1.0 - zeta) / denominator;
    functions[2][3] =  (1.0 + xi)*(1.0 - zeta) / denominator;
    functions[2][4] = -(1.0 - xi)*(1.0 + zeta) / denominator;
    functions[2][5] = -(1.0 + xi)*(1.0 + zeta) / denominator;
    functions[2][6] =  (1.0 - xi)*(1.0 + zeta) / denominator;
    functions[2][7] =  (1.0 + xi)*(1.0 + zeta) / denominator;

    functions[3][0] = -(1.0 - xi)*(1.0 - eta) / denominator;
    functions[3][1] = -(1.0 + xi)*(1.0 - eta) / denominator;
    functions[3][2] = -(1.0 - xi)*(1.0 + eta) / denominator;
    functions[3][3] = -(1.0 + xi)*(1.0 + eta) / denominator;
    functions[3][4] =  (1.0 - xi)*(1.0 - eta) / denominator;
    functions[3][5] =  (1.0 + xi)*(1.0 - eta) / denominator;
    functions[3][6] =  (1.0 - xi)*(1.0 + eta) / denominator;
    functions[3][7] =  (1.0 + xi)*(1.0 + eta) / denominator;

  } // End GetInterpolationFunctions

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sh_sort_permute(
                const typename Teuchos::Array<LocalOrdinal>::iterator& first1,
                const typename Teuchos::Array<LocalOrdinal>::iterator& last1,
                const typename Teuchos::Array<LocalOrdinal>::iterator& first2,
                const typename Teuchos::Array<LocalOrdinal>::iterator& last2) const
  {
    typedef typename std::iterator_traits<typename Teuchos::Array<LocalOrdinal>::iterator>::difference_type DT;
    DT n = last1 - first1;
    DT m = n / 2;
    DT z = Teuchos::OrdinalTraits<DT>::zero();
    while (m > z)
      {
        DT max = n - m;
        for (DT j = 0; j < max; j++)
          {
            for (DT k = j; k >= 0; k-=m)
              {
                if (first1[first2[k+m]] >= first1[first2[k]])
                  break;
                std::swap(first2[k+m], first2[k]);
              }
          }
        m = m/2;
      }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sh_sort2(
                const typename Teuchos::Array<LocalOrdinal>::iterator& first1,
                const typename Teuchos::Array<LocalOrdinal>::iterator& last1,
                const typename Teuchos::Array<LocalOrdinal>::iterator& first2,
                const typename Teuchos::Array<LocalOrdinal>::iterator& last2) const
  {
    typedef typename std::iterator_traits<typename Teuchos::Array<LocalOrdinal>::iterator>::difference_type DT;
    DT n = last1 - first1;
    DT m = n / 2;
    DT z = Teuchos::OrdinalTraits<DT>::zero();
    while (m > z)
      {
        DT max = n - m;
        for (DT j = 0; j < max; j++)
          {
            for (DT k = j; k >= 0; k-=m)
              {
                if (first1[k+m] >= first1[k])
                  break;
                std::swap(first1[k+m], first1[k]);
                std::swap(first2[k+m], first2[k]);
              }
          }
        m = m/2;
      }
  }

} //namespace MueLu

#define MUELU_GENERALGEOMETRICPFACTORY_SHORT
#endif // MUELU_GENERALGEOMETRICPFACTORY_DEF_HPP

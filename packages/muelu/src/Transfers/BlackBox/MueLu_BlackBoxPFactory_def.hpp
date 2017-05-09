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
#ifndef MUELU_BLACKBOXPFACTORY_DEF_HPP
#define MUELU_BLACKBOXPFACTORY_DEF_HPP

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
#include "MueLu_BlackBoxPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    // Coarsen can come in two forms, either a single char that will be interpreted as an integer which is used as the coarsening
    // rate in every spatial dimentions, or it can be a longer string that will then be interpreted as an array of integers.
    // By default coarsen is set as "{2}", hence a coarsening rate of 2 in every spatial dimension is the default setting!
    validParamList->set<std::string >           ("Coarsen",                 "{2}", "Coarsening rate in all spatial dimensions");
    validParamList->set<RCP<const FactoryBase> >("A",                         Teuchos::null, "Generating factory of the matrix A");
    validParamList->set<RCP<const FactoryBase> >("Nullspace",                 Teuchos::null, "Generating factory of the nullspace");
    validParamList->set<RCP<const FactoryBase> >("Coordinates",               Teuchos::null, "Generating factory for coorindates");
    validParamList->set<RCP<const FactoryBase> >("gNodesPerDim",              Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",              Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<std::string >           ("axisPermutation",         "", "Assuming a global (x,y,z) orientation, local might be (z,y,x). This vector gives a permutation from global to local orientation.");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
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
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    // // Print from all processes
    // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    // fancy->setShowAllFrontMatter(false).setShowProcRank(true);
    // Teuchos::FancyOStream& out  = *fancy;
    // // Print from a particular rank
    // const int procRank = Teuchos::GlobalMPISession::getRank();
    // Teuchos::oblackholestream blackhole;
    // std::ostream& out = ( procRank == 0 ? std::cout : blackhole );
    // // Do not print anything
    // Teuchos::oblackholestream blackhole;
    // std::ostream& out = blackhole;

    // Get parameter list
    const ParameterList& pL = GetParameterList();

    // obtain general variables
    RCP<Matrix>      A             = Get< RCP<Matrix> >      (fineLevel, "A");
    RCP<MultiVector> fineNullspace = Get< RCP<MultiVector> > (fineLevel, "Nullspace");
    RCP<Xpetra::MultiVector<double,LO,GO,NO> > coordinates = Get< RCP<Xpetra::MultiVector<double,LO,GO,NO> > >(fineLevel, "Coordinates");
    LO numDimensions  = coordinates->getNumVectors();
    LO BlkSize = A->GetFixedBlockSize();

    //  Get fine level geometric data: g(l)FineNodesPerDir and g(l)NumFineNodes
    Array<GO> gFineNodesPerDir(3);
    Array<LO> lFineNodesPerDir(3);
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
    for(LO i = 0; i < 3; ++i) {
      if(gFineNodesPerDir[i] == 0) {
        GetOStream(Runtime0) << "gNodesPerDim in direction " << i << " is set to 1 from 0" << std::endl;
        gFineNodesPerDir[i] = 1;
      }
      if(lFineNodesPerDir[i] == 0) {
        GetOStream(Runtime0) << "lNodesPerDim in direction " << i << " is set to 1 from 0" << std::endl;
        lFineNodesPerDir[i] = 1;
      }
    }
    GO gNumFineNodes = gFineNodesPerDir[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
    LO lNumFineNodes = lFineNodesPerDir[2]*lFineNodesPerDir[1]*lFineNodesPerDir[0];

    // Get the coarsening rate
    std::string coarsenRate = pL.get<std::string>("Coarsen");
    Array<LO> coarseRate(3);
    {
      Teuchos::Array<LO> crates;
      try {
        crates = Teuchos::fromStringToArray<LO>(coarsenRate);
      } catch(const Teuchos::InvalidArrayStringRepresentation e) {
        GetOStream(Errors,-1) << " *** Coarsen must be a string convertible into an array! *** " << std::endl;
        throw e;
      }
      TEUCHOS_TEST_FOR_EXCEPTION((crates.size() > 1) && (crates.size() < numDimensions),
                                 Exceptions::RuntimeError,
                                 "Coarsen must have at least as many components as the number of spatial dimensions in the problem.");
      for(LO i = 0; i < 3; ++i) {
        if(i < numDimensions) {
          if(crates.size() == 1) {
            coarseRate[i] = crates[0];
          } else if(i < crates.size()) {
            coarseRate[i] = crates[i];
          } else {
            GetOStream(Errors,-1) << " *** Coarsen must be at least as long as the number of spatial dimensions! *** " << std::endl;
            throw Exceptions::RuntimeError(" *** Coarsen must be at least as long as the number of spatial dimensions! *** \n");
          }
        } else {
          coarseRate[i] = 1;
        }
      }
    } // End of scope for crates

    GO gNumCoarseNodes = 0;
    LO lNumCoarseNodes = 0;
    Array<GO> gIndices(3), gCoarseNodesPerDir(3), ghostGIDs, coarseNodesGIDs, colGIDs;
    Array<LO> myOffset(3), lCoarseNodesPerDir(3), endRate(3);
    Array<bool> ghostInterface(6);
    ArrayRCP<Array<double> > coarseNodes(numDimensions);
    Array<ArrayView<const double> > fineNodes(numDimensions);
    for(LO dim = 0; dim < numDimensions; ++dim) {
      fineNodes[dim] = coordinates->getData(dim)();
    }
    GetGeometricData(coordinates, coarseRate, gFineNodesPerDir, lFineNodesPerDir, BlkSize,   // inputs (const)
                     gIndices, myOffset, ghostInterface, endRate, gCoarseNodesPerDir,        // outputs
                     lCoarseNodesPerDir, ghostGIDs, coarseNodesGIDs, colGIDs,
                     gNumCoarseNodes, lNumCoarseNodes, coarseNodes);

    // Create the MultiVector of coarse coordinates
    Xpetra::UnderlyingLib lib = coordinates->getMap()->lib();
    RCP<const Map> coarseCoordsMap = MapFactory::Build (lib,
                                                        gNumCoarseNodes,
                                                        coarseNodesGIDs.view(0, lNumCoarseNodes),
                                                        coordinates->getMap()->getIndexBase(),
                                                        coordinates->getMap()->getComm());
    Array<ArrayView<const double> > coarseCoords(numDimensions);
    for(LO dim = 0; dim < numDimensions; ++dim) {
      coarseCoords[dim] = coarseNodes[dim]();
    }
    RCP<Xpetra::MultiVector<double,LO,GO,NO> > coarseCoordinates = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(coarseCoordsMap, coarseCoords(), numDimensions);

    GO nTerms = gNumCoarseNodes                                                                               // Coarse nodes
      + (coarseRate[0] - 1)*(gCoarseNodesPerDir[0] - 1)*gCoarseNodesPerDir[1]*2                               // Edge nodes along direction 0
      + (coarseRate[1] - 1)*(gCoarseNodesPerDir[1] - 1)*gCoarseNodesPerDir[0]*2                               // Edge nodes along direction 1
      + (coarseRate[0] - 1)*(coarseRate[1] - 1)*(gCoarseNodesPerDir[0] - 1)*(gCoarseNodesPerDir[1] - 1)*4;    // Interior nodes
    nTerms = nTerms*BlkSize;

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetGeometricData(
                        RCP<Xpetra::MultiVector<double,LO,GO,NO> >& coordinates, const Array<LO> coarseRate,
                        const Array<GO> gFineNodesPerDir, const Array<LO> lFineNodesPerDir, const LO BlkSize,
                        Array<GO>& gIndices, Array<LO>& myOffset, Array<bool>& ghostInterface, Array<LO>& endRate,
                        Array<GO>& gCoarseNodesPerDir, Array<LO>& lCoarseNodesPerDir, Array<GO>& ghostGIDs,
                        Array<GO>& coarseNodesGIDs, Array<GO>& colGIDs, GO& gNumCoarseNodes, LO& lNumCoarseNodes,
                        ArrayRCP<Array<double> > coarseNodes) const {

    RCP<const Map> coordinatesMap = coordinates->getMap();
    LO numDimensions  = coordinates->getNumVectors();

    // Using the coarsening rate and the fine level data,
    // compute coarse level data

    //                              Phase 1                               //
    // ------------------------------------------------------------------ //
    // We first start by finding small informations on the mesh such as   //
    // the number of coarse nodes (local and global) and the number of    //
    // ghost nodes / the end rate of coarsening.                          //
    // ------------------------------------------------------------------ //
    GO minGlobalIndex = coordinatesMap->getMinGlobalIndex();
    {
      GO tmp;
      gIndices[2] = minGlobalIndex / (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      tmp         = minGlobalIndex % (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      gIndices[1] = tmp / gFineNodesPerDir[0];
      gIndices[0] = tmp % gFineNodesPerDir[0];

      myOffset[2] = gIndices[2] % coarseRate[2];
      myOffset[1] = gIndices[1] % coarseRate[1];
      myOffset[0] = gIndices[0] % coarseRate[0];
    }

    // Check whether ghost nodes are needed in each direction
    for(LO i=0; i < numDimensions; ++i) {
      if((gIndices[i] != 0) && (gIndices[i] % coarseRate[i] > 0)) {
        ghostInterface[2*i] = true;
      }
      if(((gIndices[i] + lFineNodesPerDir[i]) != gFineNodesPerDir[i]) && ((gIndices[i] + lFineNodesPerDir[i] - 1) % coarseRate[i] > 0)) {
        ghostInterface[2*i + 1] = true;
      }
    }

    for(LO i = 0; i < 3; ++i) {
      if(i < numDimensions) {
        lCoarseNodesPerDir[i] = (lFineNodesPerDir[i] + myOffset[i] - 1) / coarseRate[i];
        if(myOffset[i] == 0) { ++lCoarseNodesPerDir[i]; }
        gCoarseNodesPerDir[i] = (gFineNodesPerDir[i] - 1) / coarseRate[i];
        endRate[i]            = (gFineNodesPerDir[i] - 1) % coarseRate[i];
        if(endRate[i] == 0) {
          ++gCoarseNodesPerDir[i];
          endRate[i] = coarseRate[i];
        }
      } else {
        // Most quantities need to be set to 1 for extra dimensions
        // this is rather logical, an x-y plane is like a single layer
        // of nodes in the z direction...
        gCoarseNodesPerDir[i] = 1;
        lCoarseNodesPerDir[i] = 1;
        endRate[i]            = 1;
      }
    }

    gNumCoarseNodes = gCoarseNodesPerDir[0]*gCoarseNodesPerDir[1]*gCoarseNodesPerDir[2];
    lNumCoarseNodes = lCoarseNodesPerDir[0]*lCoarseNodesPerDir[1]*lCoarseNodesPerDir[2];

    // For each direction, determine how many ghost points are required.
    LO lNumGhostNodes = 0;
    {
      const int complementaryIndices[3][2] = {{1,2}, {0,2}, {0,1}};
      LO tmp = 0;
      for(LO i = 0; i < 3; ++i) {
        // Check whether a face in direction i needs ghost nodes
        if(ghostInterface[2*i] || ghostInterface[2*i+1]) {
          if(i == 0) {tmp = lCoarseNodesPerDir[1]*lCoarseNodesPerDir[2];}
          if(i == 1) {tmp = lCoarseNodesPerDir[0]*lCoarseNodesPerDir[2];}
          if(i == 2) {tmp = lCoarseNodesPerDir[0]*lCoarseNodesPerDir[1];}
        }
        // If both faces in direction i need nodes, double the number of ghost nodes
        if(ghostInterface[2*i] && ghostInterface[2*i+1]) {tmp = 2*tmp;}
        lNumGhostNodes += tmp;

        // The corners and edges need to be checked in 2D / 3D to add more ghosts nodes
        for(LO j = 0; j < 2; ++j) {
          for(LO k = 0; k < 2; ++k) {
            // Check if two adjoining faces need ghost nodes and then add their common edge
            if(ghostInterface[2*complementaryIndices[i][0]+j] && ghostInterface[2*complementaryIndices[i][1]+k]) {
              lNumGhostNodes += lCoarseNodesPerDir[i];
              // Add corners if three adjoining faces need ghost nodes,
              // but add them only once! Hence when i == 0.
              if(ghostInterface[2*i] && (i == 0)) { lNumGhostNodes += 1; }
              if(ghostInterface[2*i+1] && (i == 0)) { lNumGhostNodes += 1; }
            }
          }
        }
        tmp = 0;
      }
    } // end of scope for tmp and complementaryIndices

    //                              Phase 2                               //
    // ------------------------------------------------------------------ //
    // Build a list of GIDs to import the required ghost nodes.           //
    // The ordering of the ghosts nodes will be as natural as possible,   //
    // i.e. it should follow the GID ordering of the mesh.                //
    // ------------------------------------------------------------------ //

    // Saddly we have to more or less redo what was just done to figure out the value of lNumGhostNodes,
    // there might be some optimization possibility here...
    ghostGIDs.resize(lNumGhostNodes);
    LO countGhosts = 0;
    // Get the GID of the first point on the current partition.
    GO startingGID = minGlobalIndex;
    Array<GO> startingIndices(3);
    // We still want ghost nodes even if have with a 0 offset,
    // except when we are on a boundary
    if(ghostInterface[4] && (myOffset[2] == 0)) {
      if(gIndices[2] + coarseRate[2] > gFineNodesPerDir[2]) {
        startingGID -= endRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
      } else {
        startingGID -= coarseRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
      }
    } else {
      startingGID -= myOffset[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
    }
    if(ghostInterface[2] && (myOffset[1] == 0)) {
      if(gIndices[1] + coarseRate[1] > gFineNodesPerDir[1]) {
        startingGID -= endRate[1]*gFineNodesPerDir[0];
      } else {
        startingGID -= coarseRate[1]*gFineNodesPerDir[0];
      }
    } else {
      startingGID -= myOffset[1]*gFineNodesPerDir[0];
    }
    if(ghostInterface[0] && (myOffset[0] == 0)) {
      if(gIndices[0] + coarseRate[0] > gFineNodesPerDir[0]) {
        startingGID -= endRate[0];
      } else {
        startingGID -= coarseRate[0];
      }
    } else {
      startingGID -= myOffset[0];
    }

    { // scope for tmp
      GO tmp;
      startingIndices[2] = startingGID / (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      tmp = startingGID % (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
      startingIndices[1] = tmp / gFineNodesPerDir[0];
      startingIndices[0] = tmp % gFineNodesPerDir[0];
    }

    GO ghostOffset[3] = {0, 0, 0};
    LO lengthZero  = lCoarseNodesPerDir[0];
    LO lengthOne   = lCoarseNodesPerDir[1];
    LO lengthTwo   = lCoarseNodesPerDir[2];
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
          ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
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
              ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[0];
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
              ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
              ++countGhosts;
            }
          } else {
            // Set ghostOffset[1] for side faces sweep
            if( (j == lengthOne-1) && (startingIndices[1] + j*coarseRate[1] + 1 > gFineNodesPerDir[1]) ) {
              ghostOffset[1] = ( (j-1)*coarseRate[1] + endRate[1] )*gFineNodesPerDir[0];
            } else {
              ghostOffset[1] = j*coarseRate[1]*gFineNodesPerDir[0];
            }

            // Set ghostOffset[0], ghostGIDs and countGhosts
            if(ghostInterface[0]) { // In that case ghostOffset[0]==0, so we can ignore it
              ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1];
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
              ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
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
          if( (k == lengthZero-1) && (startingIndices[0] + k*coarseRate[0] + 1 > gFineNodesPerDir[0]) ) {
            ghostOffset[0] = (k-1)*coarseRate[0] + endRate[0];
          } else {
            ghostOffset[0] = k*coarseRate[0];
          }
          ghostGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
          ++countGhosts;
        }
      }
    }

    //                              Phase 3                               //
    // ------------------------------------------------------------------ //
    // Final phase of this function, lists are being built to create the  //
    // column and domain maps of the projection as well as the map for    //
    // coarse coordinates multivector.                                    //
    // ------------------------------------------------------------------ //

    Array<GO> gCoarseNodesGIDs(lNumCoarseNodes);
    LO currentNode, offset2, offset1, offset0;
    // Find the GIDs of the coarse nodes on the partition.
    for(LO ind2 = 0; ind2 < lCoarseNodesPerDir[2]; ++ind2) {
      if(myOffset[2] == 0) {
        offset2 = startingIndices[2] + myOffset[2];
      } else {
        if(startingIndices[2] + endRate[2] == gFineNodesPerDir[2] - 1) {
          offset2 = startingIndices[2] + endRate[2];
        } else {
          offset2 = startingIndices[2] + coarseRate[2];
        }
      }
      if(offset2 + ind2*coarseRate[2] > gFineNodesPerDir[2] - 1) {
        offset2 += (ind2 - 1)*coarseRate[2] + endRate[2];
      } else {
        offset2 += ind2*coarseRate[2];
      }
      offset2 = offset2*gFineNodesPerDir[1]*gFineNodesPerDir[0];

      for(LO ind1 = 0; ind1 < lCoarseNodesPerDir[1]; ++ind1) {
        if(myOffset[1] == 0) {
          offset1 = startingIndices[1] + myOffset[1];
        } else {
          if(startingIndices[1] + endRate[1] == gFineNodesPerDir[1] - 1) {
            offset1 = startingIndices[1] + endRate[1];
          } else {
            offset1 = startingIndices[1] + coarseRate[1];
          }
        }
        if(offset1 + ind1*coarseRate[1] > gFineNodesPerDir[1] - 1) {
          offset1 += (ind1 - 1)*coarseRate[1] + endRate[1];
        } else {
          offset1 += ind1*coarseRate[1];
        }
        offset1 = offset1*gFineNodesPerDir[0];
        for(LO ind0 = 0; ind0 < lCoarseNodesPerDir[0]; ++ind0) {
          offset0 = startingIndices[0];
          if(myOffset[0] == 0) {
            offset0 += ind0*coarseRate[0];
          } else {
            offset0 += (ind0 + 1)*coarseRate[0];
          }
          if(offset0 > gFineNodesPerDir[0] - 1) {offset0 += endRate[0] - coarseRate[0];}

          currentNode = ind2*lCoarseNodesPerDir[1]*lCoarseNodesPerDir[0]
                      + ind1*lCoarseNodesPerDir[0]
                      + ind0;
          gCoarseNodesGIDs[currentNode] = offset2 + offset1 + offset0;
        }
      }
    }

    // Actual loop over all the coarse/ghost nodes to find their index on the coarse mesh
    // and the corresponding dofs that will need to be added to colMapP.
    colGIDs.resize(BlkSize*(lNumCoarseNodes+lNumGhostNodes));
    coarseNodesGIDs.resize(lNumCoarseNodes);
    for(LO i = 0; i < numDimensions; ++i) {coarseNodes[i].resize(lNumCoarseNodes);}
    GO fineNodesPerCoarseSlab    = coarseRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
    GO fineNodesEndCoarseSlab    = endRate[2]*gFineNodesPerDir[1]*gFineNodesPerDir[0];
    GO fineNodesPerCoarsePlane   = coarseRate[1]*gFineNodesPerDir[0];
    GO fineNodesEndCoarsePlane   = endRate[1]*gFineNodesPerDir[0];
    GO coarseNodesPerCoarseLayer = gCoarseNodesPerDir[1]*gCoarseNodesPerDir[0];
    GO gCoarseNodeOnCoarseGridGID;
    LO gInd[3], lCol;
    Array<int> ghostPIDs  (lNumGhostNodes);
    Array<LO>  ghostLIDs  (lNumGhostNodes);
    Array<LO>  ghostPermut(lNumGhostNodes);
    for(LO k = 0; k < lNumGhostNodes; ++k) {ghostPermut[k] = k;}
    coordinatesMap->getRemoteIndexList(ghostGIDs, ghostPIDs, ghostLIDs);
    sh_sort_permute(ghostPIDs.begin(),ghostPIDs.end(), ghostPermut.begin(),ghostPermut.end());

    { // scope for tmpInds, tmpVars and tmp.
      GO tmpInds[3], tmpVars[2];
      LO tmp;
      // Loop over the coarse nodes of the partition and add them to colGIDs
      // that will be used to construct the column and domain maps of P as well
      // as to construct the coarse coordinates map.
      // for(LO col = 0; col < lNumCoarseNodes; ++col) { // This should most likely be replaced by loops of lCoarseNodesPerDir[] to simplify arithmetics
      LO col = 0;
      LO firstCoarseNodeInds[3], currentCoarseNode;
      for(LO dim = 0; dim < 3; ++dim) {
        if(myOffset[dim] == 0) {
          firstCoarseNodeInds[dim] = 0;
        } else {
          firstCoarseNodeInds[dim] = coarseRate[dim] - myOffset[dim];
        }
      }
      Array<ArrayRCP<const double> > fineNodes(numDimensions);
      for(LO dim = 0; dim < numDimensions; ++dim) {fineNodes[dim] = coordinates->getData(dim);}
      for(LO k = 0; k < lCoarseNodesPerDir[2]; ++k) {
        for(LO j = 0; j < lCoarseNodesPerDir[1]; ++j) {
          for(LO i = 0; i < lCoarseNodesPerDir[0]; ++i) {
            col = k*lCoarseNodesPerDir[1]*lCoarseNodesPerDir[0] + j*lCoarseNodesPerDir[0] + i;

            // Check for endRate
            currentCoarseNode = 0;
            if(firstCoarseNodeInds[0] + i*coarseRate[0] > lFineNodesPerDir[0] - 1) {
              currentCoarseNode += firstCoarseNodeInds[0] + (i-1)*coarseRate[0] + endRate[0];
            } else {
              currentCoarseNode += firstCoarseNodeInds[0] + i*coarseRate[0];
            }
            if(firstCoarseNodeInds[1] + j*coarseRate[1] > lFineNodesPerDir[1] - 1) {
              currentCoarseNode += (firstCoarseNodeInds[1] + (j-1)*coarseRate[1] + endRate[1])*lFineNodesPerDir[0];
            } else {
              currentCoarseNode += (firstCoarseNodeInds[1] + j*coarseRate[1])*lFineNodesPerDir[0];
            }
            if(firstCoarseNodeInds[2] + k*coarseRate[2] > lFineNodesPerDir[2] - 1) {
              currentCoarseNode += (firstCoarseNodeInds[2] + (k-1)*coarseRate[2] + endRate[2])*lFineNodesPerDir[1]*lFineNodesPerDir[0];
            } else {
              currentCoarseNode += (firstCoarseNodeInds[2] + k*coarseRate[2])*lFineNodesPerDir[1]*lFineNodesPerDir[0];
            }
            // Load coordinates
            for(LO dim = 0; dim < numDimensions; ++dim) {
              coarseNodes[dim][col] = fineNodes[dim][currentCoarseNode];
            }

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
            gInd[2] = col / (lCoarseNodesPerDir[1]*lCoarseNodesPerDir[0]);
            tmp     = col % (lCoarseNodesPerDir[1]*lCoarseNodesPerDir[0]);
            gInd[1] = tmp / lCoarseNodesPerDir[0];
            gInd[0] = tmp % lCoarseNodesPerDir[0];
            lCol = gInd[2]*(lCoarseNodesPerDir[1]*lCoarseNodesPerDir[0]) + gInd[1]*lCoarseNodesPerDir[0] + gInd[0];
            gCoarseNodeOnCoarseGridGID = tmpInds[2]*coarseNodesPerCoarseLayer + tmpInds[1]*gCoarseNodesPerDir[0] + tmpInds[0];
            coarseNodesGIDs[lCol] = gCoarseNodeOnCoarseGridGID;
            for(LO dof = 0; dof < BlkSize; ++dof) {
              colGIDs[BlkSize*lCol + dof] = BlkSize*gCoarseNodeOnCoarseGridGID + dof;
            }
          }
        }
      }
      // Now loop over the ghost nodes of the partition to add them to colGIDs
      // since they will need to be included in the column map of P
      for(col = lNumCoarseNodes; col < lNumCoarseNodes + lNumGhostNodes; ++col) {
        if((endRate[2] != coarseRate[2]) && (ghostGIDs[ghostPermut[col - lNumCoarseNodes]] > (gCoarseNodesPerDir[2] - 2)*fineNodesPerCoarseSlab + fineNodesEndCoarseSlab - 1)) {
          tmpInds[2] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] / fineNodesPerCoarseSlab + 1;
          tmpVars[0] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] - (tmpInds[2] - 1)*fineNodesPerCoarseSlab - fineNodesEndCoarseSlab;
        } else {
          tmpInds[2] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] / fineNodesPerCoarseSlab;
          tmpVars[0] = ghostGIDs[ghostPermut[col - lNumCoarseNodes]] % fineNodesPerCoarseSlab;
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
        for(LO dof = 0; dof < BlkSize; ++dof) {
          colGIDs[BlkSize*col + dof] = BlkSize*gCoarseNodeOnCoarseGridGID + dof;
        }
      }
    } // End of scope for tmpInds, tmpVars and tmp

  } // GetGeometricData()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sh_sort_permute(
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

} //namespace MueLu

#define MUELU_BLACKBOXPFACTORY_SHORT
#endif // MUELU_BLACKBOXPFACTORY_DEF_HPP

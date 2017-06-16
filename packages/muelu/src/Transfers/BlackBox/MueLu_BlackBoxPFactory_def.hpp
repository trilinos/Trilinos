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

    // Print from all processes
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fancy->setShowAllFrontMatter(false).setShowProcRank(true);
    Teuchos::FancyOStream& out  = *fancy;
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
    for(LO dim = 0; dim < numDimensions; ++dim) {fineNodes[dim] = coordinates->getData(dim)();}

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

    // Now create a new matrix: Aghost that contains all the data
    // locally needed to compute the prolongator on the node.
    // Here assuming that all the coarse nodes o and fine nodes +
    // are local then all the data associated with the coarse
    // nodes O and the fine nodes * needs to be imported.
    //
    //                  o--+--+--o--+--+--O
    //                  |  |  |  |  |  |  |
    //                  +--+--+--+--+--+--*
    //                  |  |  |  |  |  |  |
    //                  +--+--+--+--+--+--*
    //                  |  |  |  |  |  |  |
    //                  o--+--+--o--+--+--O
    //                  |  |  |  |  |  |  |
    //                  +--+--+--+--+--+--*
    //                  |  |  |  |  |  |  |
    //                  *--*--*--*--*--*--*
    //                  |  |  |  |  |  |  |
    //                  O--*--*--O--*--*--O
    //
    // Creating that local matrix is easy enough using proper range
    // and domain maps to import data from A.
    // As usual we need to be careful about any coarsening rate
    // change at the boundary!

    // The ingredients needed are an importer, a range map and a domain map
    Array<GO> ghostRowGIDs, nodeSteps(3);
    nodeSteps[0] = 1;
    nodeSteps[1] = gFineNodesPerDir[0];
    nodeSteps[2] = gFineNodesPerDir[0]*gFineNodesPerDir[1];
    Array<LO> range(3);
    GO startingGID = A->getRowMap()->getMinGlobalIndex();
    for(LO dim = 0; dim < 3; ++dim) {
      LO numCoarseNodes = 0;
      if(dim < numDimensions) {
        startingGID -= myOffset[dim]*nodeSteps[dim];
        numCoarseNodes = lCoarseNodesPerDir[dim];
        if(ghostInterface[2*dim]) {++numCoarseNodes;}
        if(ghostInterface[2*dim+1]) {++numCoarseNodes;}
        if(gIndices[dim] + lFineNodesPerDir[dim] == gFineNodesPerDir[dim]) {
          range[dim] = (numCoarseNodes - 2)*coarseRate[dim] + endRate[dim] + 1;
        } else {
          range[dim] = (numCoarseNodes - 1)*coarseRate[dim] + 1;
        }
      } else {
        range[dim] = 1;
      }
    }
    ghostRowGIDs.resize(range[0]*range[1]*range[2]*BlkSize);
    for(LO k = 0; k < range[2]; ++k) {
      for(LO j = 0; j < range[1]; ++j) {
        for(LO i = 0; i < range[0]; ++i) {
          for(LO l = 0; l < BlkSize; ++l) {
            ghostRowGIDs[(k*range[1]*range[0] + j*range[0] + i)*BlkSize + l] = startingGID
              + (k*gFineNodesPerDir[1]*gFineNodesPerDir[0] + j*gFineNodesPerDir[0] + i)*BlkSize + l;
          }
        }
      }
    }
    RCP<const Map> ghostedRowMap = Xpetra::MapFactory<LO,GO,NO>::Build(A->getRowMap()->lib(),
                                                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                       ghostRowGIDs(),
                                                                       A->getRowMap()->getIndexBase(),
                                                                       A->getRowMap()->getComm());
    RCP<const Import> ghostImporter = Xpetra::ImportFactory<LO,GO,NO>::Build(A->getRowMap(), ghostedRowMap);
    RCP<const Matrix> Aghost        = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(A, *ghostImporter, ghostedRowMap, ghostedRowMap);

    // Create the maps and data structures for the projection matrix
    RCP<const Map> rowMapP    = A->getDomainMap();
    RCP<const Map> domainMapP = Xpetra::MapFactory<LO,GO,NO>::Build(rowMapP->lib(),
                                                                    BlkSize*gNumCoarseNodes,
                                                                    colGIDs.view(0, BlkSize*lNumCoarseNodes),
                                                                    rowMapP->getIndexBase(),
                                                                    rowMapP->getComm());

    RCP<const Map> colMapP = Xpetra::MapFactory<LO,GO,NO>::Build(rowMapP->lib(),
                                                                 Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                 colGIDs.view(0, colGIDs.size()),
                                                                 rowMapP->getIndexBase(),
                                                                 rowMapP->getComm());

    std::vector<size_t> strideInfo(1);
    strideInfo[0] = BlkSize;
    RCP<const Map> stridedDomainMapP = Xpetra::StridedMapFactory<LO,GO,NO>::Build(domainMapP, strideInfo);

    GO gNumEdgeNodes = 0, gNumInteriorNodes = 1;
    LO complementaryIndices[3][2] = {{1,2}, {0,2}, {0,1}};
    for(LO dim = 0; dim < numDimensions; ++dim) {
      gNumEdgeNodes += (coarseRate[dim] - 1)*(gCoarseNodesPerDir[dim] - 1)
        *gCoarseNodesPerDir[complementaryIndices[dim][0]]*gCoarseNodesPerDir[complementaryIndices[dim][1]];
      gNumInteriorNodes = gNumInteriorNodes*(coarseRate[dim] - 1)*(gCoarseNodesPerDir[dim] - 1);
    }

    GO nnzP = gNumCoarseNodes + gNumEdgeNodes*2 + gNumInteriorNodes*std::pow(2.0,numDimensions);
    nnzP = nnzP*BlkSize;

    // Create the matrix itself using the above maps
    RCP<Matrix> P;
    P = rcp(new CrsMatrixWrap(rowMapP, colMapP, 0, Xpetra::StaticProfile));
    RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

    ArrayRCP<size_t>  iaP;
    ArrayRCP<LO>      jaP;
    ArrayRCP<SC>     valP;

    PCrs->allocateAllValues(nnzP, iaP, jaP, valP);

    ArrayView<size_t> ia  = iaP();
    ArrayView<LO>     ja  = jaP();
    ArrayView<SC>     val = valP();
    ia[0] = 0;


    LO numCoarseElements = 1;
    Array<LO> lCoarseElementsPerDir(3);
    for(LO dim = 0; dim < numDimensions; ++dim) {
      lCoarseElementsPerDir[dim] = lCoarseNodesPerDir[dim];
      if(ghostInterface[2*dim]) {++lCoarseElementsPerDir[dim];}
      if(!ghostInterface[2*dim+1]) {--lCoarseElementsPerDir[dim];}
      numCoarseElements = numCoarseElements*lCoarseElementsPerDir[dim];
    }

    for(LO dim = numDimensions; dim < 3; ++dim) {
      lCoarseElementsPerDir[dim] = 1;
    }

    out << "numCoarseElements=" << numCoarseElements << std::endl;
    out << "lCoarseElementsPerDir=(" << lCoarseElementsPerDir[0] << ", " << lCoarseElementsPerDir[1] << ", " << lCoarseElementsPerDir[2] << ")" << std::endl;

    // Loop over the coarse elements
    Array<LO> elemInds(3);
    for(elemInds[2] = 0; elemInds[2] < lCoarseElementsPerDir[2]; ++elemInds[2]) {
      for(elemInds[1] = 0; elemInds[1] < lCoarseElementsPerDir[1]; ++elemInds[1]) {
        for(elemInds[0] = 0; elemInds[0] < lCoarseElementsPerDir[0]; ++elemInds[0]) {
          ComputeLocalEntries(Aghost, coarseRate, endRate, BlkSize, elemInds, lCoarseElementsPerDir, range,
                              numDimensions, lFineNodesPerDir, gFineNodesPerDir, gIndices, lCoarseNodesPerDir);
        }
      }
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetGeometricData(
                        RCP<Xpetra::MultiVector<double,LO,GO,NO> >& coordinates, const Array<LO> coarseRate,
                        const Array<GO> gFineNodesPerDir, const Array<LO> lFineNodesPerDir, const LO BlkSize,
                        Array<GO>& gIndices, Array<LO>& myOffset, Array<bool>& ghostInterface, Array<LO>& endRate,
                        Array<GO>& gCoarseNodesPerDir, Array<LO>& lCoarseNodesPerDir, Array<GO>& ghostGIDs,
                        Array<GO>& coarseNodesGIDs, Array<GO>& colGIDs, GO& gNumCoarseNodes, LO& lNumCoarseNodes,
                        ArrayRCP<Array<double> > coarseNodes) const {
    // This function is extracting the geometric information from the coordinates
    // and creates the necessary data/formatting to perform locally the calculation
    // of the pronlongator.
    //
    // inputs:

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
    // column and domain maps of the projection as well as the map and    //
    // arrays for the coarse coordinates multivector.                     //
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
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeLocalEntries(
                        const RCP<const Matrix>& Aghost, const Array<LO> coarseRate,
                        const Array<LO> endRate, const LO BlkSize, const Array<LO> elemInds,
                        const Array<LO> lCoarseElementsPerDir, const Array<LO> range,
                        const LO numDimensions, const Array<LO> lFineNodesPerDir,
                        const Array<GO> gFineNodesPerDir, const Array<GO> gIndices,
                        const Array<LO> lCoarseNodesPerDir) const {

    // Print from all processes
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fancy->setShowAllFrontMatter(false).setShowProcRank(true);
    Teuchos::FancyOStream& out  = *fancy;
    // // Print from a particular rank
    // const int procRank = Teuchos::GlobalMPISession::getRank();
    // Teuchos::oblackholestream blackhole;
    // std::ostream& out = ( procRank == 0 ? std::cout : blackhole );
    // // Do not print anything
    // Teuchos::oblackholestream blackhole;
    // std::ostream& out = blackhole;

    // First extract data from Aghost and move it to the corresponding dense matrix
    // This step requires to compute the number of nodes (resp dofs) in the current
    // coarse element, then allocate storage for local dense matrices, finally loop
    // over the dofs and extract the corresponding row in Aghost before dispactching
    // its data to the dense matrices.
    Array<LO> elementNodesPerDir(3);
    for(LO dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        if((elemInds[dim] == lCoarseElementsPerDir[dim]) && (gIndices[dim] + lFineNodesPerDir[dim] == gFineNodesPerDir[dim])) {
          elementNodesPerDir[dim] = endRate[dim] + 1;
        } else {
          elementNodesPerDir[dim] = coarseRate[dim] + 1;
        }
      } else {
        elementNodesPerDir[dim] = 1;
      }
    }
    LO numNodesInElement = elementNodesPerDir[0]*elementNodesPerDir[1]*elementNodesPerDir[2];
    LO countInterior=0, countFace=0, countEdge=0, countCorner =0;
    Array<LO> dofType(numNodesInElement*BlkSize), lDofInd(numNodesInElement*BlkSize);
    Array<LO> collapseDir(numNodesInElement*BlkSize);
    for(LO ke = 0; ke < elementNodesPerDir[2]; ++ke) {
      for(LO je = 0; je < elementNodesPerDir[1]; ++je) {
        for(LO ie = 0; ie < elementNodesPerDir[0]; ++ie) {
          // Check for corner node
          if((ke == 0 || ke == elementNodesPerDir[2]-1)
             && (je == 0 || je == elementNodesPerDir[1]-1)
             && (ie == 0 || ie == elementNodesPerDir[0]-1)) {
            for(LO dof = 0; dof < BlkSize; ++dof) {
              dofType[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 0;
              lDofInd[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = BlkSize*countCorner + dof;
            }
            ++countCorner;

          // Check for edge node
          } else if (((ke == 0 || ke == elementNodesPerDir[2]-1) && (je == 0 || je == elementNodesPerDir[1]-1))
                     || ((ke == 0 || ke == elementNodesPerDir[2]-1) && (ie == 0 || ie == elementNodesPerDir[0]-1))
                     || ((je == 0 || je == elementNodesPerDir[1]-1) && (ie == 0 || ie == elementNodesPerDir[0]-1))) {
            for(LO dof = 0; dof < BlkSize; ++dof) {
              dofType[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 1;
              lDofInd[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = BlkSize*countEdge + dof;
              if((ke == 0 || ke == elementNodesPerDir[2]-1) && (je == 0 || je == elementNodesPerDir[1]-1)) {
                collapseDir[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 0;
              } else if((ke == 0 || ke == elementNodesPerDir[2]-1) && (ie == 0 || ie == elementNodesPerDir[0]-1)) {
                collapseDir[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 1;
              } else if((je == 0 || je == elementNodesPerDir[1]-1) && (ie == 0 || ie == elementNodesPerDir[0]-1)) {
                collapseDir[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 2;
              }
            }
            ++countEdge;

          // Check for face node
          } else if ((ke == 0 || ke == elementNodesPerDir[2]-1)
                     || (je == 0 || je == elementNodesPerDir[1]-1)
                     || (ie == 0 || ie == elementNodesPerDir[0]-1)) {
            for(LO dof = 0; dof < BlkSize; ++dof) {
              dofType[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 2;
              lDofInd[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = BlkSize*countFace + dof;
              if(ke == 0 || ke == elementNodesPerDir[2]-1) {
                collapseDir[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 2;
              } else if(je == 0 || je == elementNodesPerDir[1]-1) {
                collapseDir[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 1;
              } else if(ie == 0 || ie == elementNodesPerDir[0]-1) {
                collapseDir[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 0;
              }
            }
            ++countFace;

          // Otherwise it is an interior node
          } else {
            for(LO dof = 0; dof < BlkSize; ++dof) {
              dofType[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = 3;
              lDofInd[BlkSize*(ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie) + dof] = BlkSize*countInterior + dof;
            }
            ++countInterior;
          }
        }
      }
    }

    // Diagonal blocks
    Teuchos::SerialDenseMatrix<LO,double> Aii(BlkSize*elementNodesPerDir[0]*elementNodesPerDir[1]*elementNodesPerDir[2],
                                              BlkSize*elementNodesPerDir[0]*elementNodesPerDir[1]*elementNodesPerDir[2]);
    Teuchos::SerialDenseMatrix<LO,double> Aff(2*BlkSize*(elementNodesPerDir[0]*elementNodesPerDir[1]
                                                         + elementNodesPerDir[1]*elementNodesPerDir[2]
                                                         + elementNodesPerDir[0]*elementNodesPerDir[2]),
                                              2*BlkSize*(elementNodesPerDir[0]*elementNodesPerDir[1]
                                                         + elementNodesPerDir[1]*elementNodesPerDir[2]
                                                         + elementNodesPerDir[0]*elementNodesPerDir[2]));
    Teuchos::SerialDenseMatrix<LO,double> Aee(4*BlkSize*(elementNodesPerDir[0] + elementNodesPerDir[1] + elementNodesPerDir[2]),
                                              4*BlkSize*(elementNodesPerDir[0] + elementNodesPerDir[1] + elementNodesPerDir[2]));
    // Upper triangular blocks
    Teuchos::SerialDenseMatrix<LO,double> Aif(BlkSize*elementNodesPerDir[0]*elementNodesPerDir[1]*elementNodesPerDir[2],
                                              2*BlkSize*(elementNodesPerDir[0]*elementNodesPerDir[1]
                                                         + elementNodesPerDir[1]*elementNodesPerDir[2]
                                                         + elementNodesPerDir[0]*elementNodesPerDir[2]));
    Teuchos::SerialDenseMatrix<LO,double> Aie(BlkSize*elementNodesPerDir[0]*elementNodesPerDir[1]*elementNodesPerDir[2],
                                              4*BlkSize*(elementNodesPerDir[0] + elementNodesPerDir[1] + elementNodesPerDir[2]));
    Teuchos::SerialDenseMatrix<LO,double> Afe(2*BlkSize*(elementNodesPerDir[0]*elementNodesPerDir[1]
                                                         + elementNodesPerDir[1]*elementNodesPerDir[2]
                                                         + elementNodesPerDir[0]*elementNodesPerDir[2]),
                                              4*BlkSize*(elementNodesPerDir[0] + elementNodesPerDir[1] + elementNodesPerDir[2]));
    // Coarse nodes "right hand sides"
    Teuchos::SerialDenseMatrix<LO,double> Aic(BlkSize*elementNodesPerDir[0]*elementNodesPerDir[1]*elementNodesPerDir[2], 8);
    Teuchos::SerialDenseMatrix<LO,double> Afc(2*BlkSize*(elementNodesPerDir[0]*elementNodesPerDir[1]
                                                         + elementNodesPerDir[1]*elementNodesPerDir[2]
                                                         + elementNodesPerDir[0]*elementNodesPerDir[2]), 8);
    Teuchos::SerialDenseMatrix<LO,double> Aec(4*BlkSize*(elementNodesPerDir[0] + elementNodesPerDir[1] + elementNodesPerDir[2]), 8);

    ArrayView<const LO> rowIndices;
    ArrayView<const SC> rowValues;
    LO idof, jdof, iInd, jInd;
    LO indi, indf, inde, indc;
    int itype = 0, jtype = 0;
    Array<SC> stencil(std::pow(3,numDimensions));
    // LBV Note: we could skip the extraction of rows corresponding to coarse nodes
    for(LO ke = 0; ke < elementNodesPerDir[2]; ++ke) {
      for(LO je = 0; je < elementNodesPerDir[1]; ++je) {
        for(LO ie = 0; ie < elementNodesPerDir[0]; ++ie) {
          GetNodeInfo(ie, je, ke, elementNodesPerDir, &itype, iInd);
          for(LO dof0 = 0; dof0 < BlkSize; ++dof0) {
            idof = (ke*elementNodesPerDir[1]*elementNodesPerDir[0] + je*elementNodesPerDir[0] + ie)*BlkSize + dof0;
            Aghost->getLocalRowView(idof, rowIndices, rowValues);
            // To make life easy, let's reshape and reorder rowValues to
            // process it in a standardized way for both 7 points stencils
            // and 27 points stencil as well as handle ghost data.
            ReorderStencil(ie, je, ke, rowValues, elementNodesPerDir, stencil);
            // We are handling an interior node, we won't have to handle
            // ghost data but we need to check if entries go to Aif or Aii
            if(itype == 3) {
              indi = (ke - 1)*(elementNodesPerDir[1] - 2)*(elementNodesPerDir[0] - 2)
                + (je - 1)*(elementNodesPerDir[0] - 2) + ie - 1;
              Aii(indi, indi) = rowValues[3];
              if(ke == 1) {// The entrie goes to Aif
                indf = BlkSize*((je - 1)*(elementNodesPerDir[0] - 2) + ie - 1) + dof0;
                Aif(indi, indf) = rowValues[0];
              } else {
                Aii(indi, indi - elementNodesPerDir[1]*elementNodesPerDir[0]) = rowValues[0];
              }
              if(ke == elementNodesPerDir[2] - 2) {// The entrie goes to Aif
                indf = BlkSize*((elementNodesPerDir[1] - 2)*(elementNodesPerDir[0] - 2)
                                   + 2*(elementNodesPerDir[2] - 2)*(elementNodesPerDir[1] - 2 + elementNodesPerDir[0] - 2)
                                   + (je - 1)*(elementNodesPerDir[0] - 2) + ie - 1) + dof0;
                Aif(indi, indf) = rowValues[6];
              } else {
                Aii(indi, indi + elementNodesPerDir[1]*elementNodesPerDir[0]) = rowValues[6];
              }
              if(je == 1) {// The entrie goes to Aif
                indf = BlkSize*((elementNodesPerDir[1] - 2)*(elementNodesPerDir[0] - 2)
                                   + 2*(ke - 1)*(elementNodesPerDir[1] - 2 + elementNodesPerDir[0] - 2) + ie -1) + dof0;
                Aif(indi, indf) = rowValues[1];
              } else {
                Aii(indi, indi - elementNodesPerDir[0]) = rowValues[1];
              }
              if(je == elementNodesPerDir[1] - 2) {// The entrie goes to Aif
                indf = BlkSize*((elementNodesPerDir[1] - 2)*(elementNodesPerDir[0] - 2)
                                   + 2*(ke - 1)*(elementNodesPerDir[1] - 2 + elementNodesPerDir[0] - 2)
                                   + 2*(elementNodesPerDir[1] - 2) + elementNodesPerDir[0] - 2 + ie -1) + dof0;
                Aif(indi, indf) = rowValues[5];
              } else {
                Aii(indi, indi + elementNodesPerDir[0]) = rowValues[5];
              }
              if(ie == 1) {// The entrie goes to Aif
                indf = BlkSize*((elementNodesPerDir[1] - 2)*(elementNodesPerDir[0] - 2)
                                   + 2*(ke - 1)*(elementNodesPerDir[1] - 2 + elementNodesPerDir[0] - 2)
                                   + elementNodesPerDir[0] - 2 + 2*(je -1)) + dof0;
                Aif(indi, indf) = rowValues[2];
              } else {
                Aii(indi, indi - 1) = rowValues[2];
              }
              if(ie == elementNodesPerDir[0] - 2) {// The entrie goes to Aif
                indf = BlkSize*((je -1)*(elementNodesPerDir[0] - 2) + ie - 1) + dof0;
                Aif(indi, indf) = rowValues[4];
              } else {
                Aii(indi, indi + 1) = rowValues[4];
              }
              // We are handling a face node,
              // we can skip connections to interior nodes
              // we need to collapse 3D stencils to 2D stencils
            } else if(itype == 2) {
              // First check the orientation of the face
              if(ke == 0) {
                indf = (je - 1)*(elementNodesPerDir[0] - 2) + ie - 1;
                if(rowValues.size() == 5) {// We are on a mesh boundary
                  Aff(indf, indf) = rowValues[2] + rowValues[4];
                } else if(rowValues.size() == 6) {
                  Aff(indf, indf) = rowValues[0] + rowValues[3] + rowValues[5];
                }
                if(je == 1) {
                  inde = BlkSize*(ie - 1) + dof0;
                  if(rowValues.size() == 5) {// We are on a mesh boundary
                    Afe(indf, inde) = rowValues[0];
                  } else if(rowValues.size() == 6) {
                    Afe(indf, inde) = rowValues[1];
                  }
                } else if(je == elementNodesPerDir[1] - 2) {
                  inde = BlkSize*(2*(elementNodesPerDir[1] - 2) + elementNodesPerDir[0] - 2 + (ie - 1)) + dof0;
                  if(rowValues.size() == 5) {// We are on a mesh boundary
                    Afe(indf, inde) = rowValues[4];
                  } else if(rowValues.size() == 6) {
                    Afe(indf, inde) = rowValues[5];
                  }
                }
              }
            //   // we are handling an edge node,
            //   // we can skip connections to interior and face nodes
            //   // we need to collapse 3D stencils to 1D stencils
            // } else if(itype == 1) {
            //   if(jtype == 1) {
            //     Aee(lDofInd[idof],lDofInd[jdof]) = rowValues[dof1];
            //   } else if(jtype == 0) {
            //     Aec(lDofInd[idof],lDofInd[jdof]) = rowValues[dof1];
            //   }
            //   // We are handling a corner node
            } else {
              break;
            }
          }
        }
      }
    }
    // Compute the projection operators for edge and interior nodes
    //
    //         P_i = A_{ii}^{-1}A_{ic} - A_{ii}^{-1}A_{if}\hat{A}_{ff}^{-1}\hat{A}_{fc} - A_{ii}^{-1}A_{if}\hat{A}_{ff}^{-1}\hat{A}_{fe}\hat{A}_{ee}^{-1}\hat{A}_{ec}
    //         P_f = \hat{A}_{ff}^{-1}\hat{A}_{fc} -\hat{A}_{ff}^{-1}\hat{A}_{fe}\hat{A}_{ee}^{-1}\hat{A}_{ec}
    //         P_e = \hat{A}_{ee}^{-1}\hat{A}_{ec}
    //         P_c = I
    //

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReorderStencil(
                        const LO ie, const LO je, const LO ke, const ArrayView<const SC> rowValues,
                        const Array<LO> elementNodesPerDir, Array<SC>& stencil) const {
    if((ie == 0) && (je == 0) && (ke == 0)) {// corner 1
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[4];
        stencil[10] = rowValues[5];
        stencil[12] = rowValues[6];
        stencil[13] = rowValues[0];
        stencil[14] = rowValues[1];
        stencil[16] = rowValues[2];
        stencil[22] = rowValues[3];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[13 + k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        for(LO i = 0; i < 13; ++i) {
          stencil[i] = rowValues[8 + i];
        }
        stencil[15] = rowValues[21];
        for(LO i = 0; i < 4; ++i) {
          stencil[18 + i] = rowValues[22 + i];
        }
        stencil[24] = rowValues[26];
      }
    } else if((ie == elementNodesPerDir[0] - 1) && (je == 0) && (ke == 0)) {// corner 2
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[4];
        stencil[10] = rowValues[5];
        stencil[12] = rowValues[0];
        stencil[13] = rowValues[1];
        stencil[14] = rowValues[6];
        stencil[16] = rowValues[2];
        stencil[22] = rowValues[3];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[12 + k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        for(LO i = 0; i < 12; ++i) {
          stencil[i] = rowValues[8 + i];
        }
        stencil[14] = rowValues[20];
        for(LO i = 0; i < 4; ++i) {
          stencil[17 + i] = rowValues[21 + i];
        }
        stencil[23] = rowValues[25];
        stencil[26] = rowValues[26];
      }
    } else if((ie == 0) && (je == elementNodesPerDir[1] - 1) && (ke == 0)) {// corner 3
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[4];
        stencil[10] = rowValues[0];
        stencil[12] = rowValues[5];
        stencil[13] = rowValues[1];
        stencil[14] = rowValues[2];
        stencil[16] = rowValues[6];
        stencil[22] = rowValues[3];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[10 + k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        for(LO i = 0; i < 10; ++i) {
          stencil[i] = rowValues[8 + i];
        }
        stencil[12] = rowValues[18];
        for(LO i = 0; i < 4; ++i) {
          stencil[15 + i] = rowValues[19 + i];
        }
        stencil[21] = rowValues[23];
        for(LO i = 0; i < 3; ++i) {
          stencil[24 + i] = rowValues[24 + i];
        }
      }
    } else if((ie == elementNodesPerDir[0] - 1) && (je == elementNodesPerDir[1] - 1) && (ke == 0)) {// corner 4
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[4];
        stencil[10] = rowValues[0];
        stencil[12] = rowValues[1];
        stencil[13] = rowValues[2];
        stencil[14] = rowValues[5];
        stencil[16] = rowValues[6];
        stencil[22] = rowValues[3];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[9 + k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        for(LO i = 0; i < 9; ++i) {
          stencil[i] = rowValues[8 + i];
        }
        stencil[11] = rowValues[17];
        for(LO i = 0; i < 4; ++i) {
          stencil[14 + i] = rowValues[18 + i];
        }
        stencil[20] = rowValues[22];
        for(LO i = 0; i < 4; ++i) {
          stencil[23 + i] = rowValues[23 + i];
        }
      }
    } else if((ie == 0) && (je == 0) && (ke == elementNodesPerDir[2] - 1)) {// corner 5
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[0];
        stencil[10] = rowValues[4];
        stencil[12] = rowValues[5];
        stencil[13] = rowValues[1];
        stencil[14] = rowValues[2];
        stencil[16] = rowValues[3];
        stencil[22] = rowValues[6];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[4 + k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        for(LO i = 0; i < 4; ++i) {
          stencil[i] = rowValues[8 + i];
        }
        stencil[6] = rowValues[12];
        for(LO i = 0; i < 4; ++i) {
          stencil[9 + i] = rowValues[13 + i];
        }
        stencil[15] = rowValues[17];
        for(LO i = 0; i < 9; ++i) {
          stencil[18 + i] = rowValues[18 + i];
        }
      }
    } else if((ie == elementNodesPerDir[0] - 1) && (je == 0) && (ke == elementNodesPerDir[2] - 1)) {// corner 6
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[0];
        stencil[10] = rowValues[4];
        stencil[12] = rowValues[1];
        stencil[13] = rowValues[2];
        stencil[14] = rowValues[5];
        stencil[16] = rowValues[3];
        stencil[22] = rowValues[6];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[3 + k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        for(LO i = 0; i < 3; ++i) {
          stencil[i] = rowValues[8 + i];
        }
        stencil[5] = rowValues[11];
        for(LO i = 0; i < 4; ++i) {
          stencil[8 + i] = rowValues[12 + i];
        }
        stencil[14] = rowValues[16];
        for(LO i = 0; i < 10; ++i) {
          stencil[17 + i] = rowValues[17 + i];
        }
      }
    } else if((ie == 0) && (je == elementNodesPerDir[1] - 1) && (ke == elementNodesPerDir[2] - 1)) {// corner 7
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[0];
        stencil[10] = rowValues[1];
        stencil[12] = rowValues[4];
        stencil[13] = rowValues[2];
        stencil[14] = rowValues[3];
        stencil[16] = rowValues[5];
        stencil[22] = rowValues[6];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[1 + k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        stencil[0] = rowValues[8];
        stencil[3] = rowValues[9];
        for(LO i = 0; i < 4; ++i) {
          stencil[6 + i] = rowValues[10 + i];
        }
        stencil[12] = rowValues[14];
        for(LO i = 0; i < 12; ++i) {
          stencil[15 + i] = rowValues[15 + i];
        }
      }
    } else if((ie == elementNodesPerDir[0] - 1) && (je == elementNodesPerDir[1] - 1) && (ke == elementNodesPerDir[2] - 1)) {// corner 8
      if(rowValues.size() == 7) {
        stencil[ 4] = rowValues[0];
        stencil[10] = rowValues[1];
        stencil[12] = rowValues[2];
        stencil[13] = rowValues[3];
        stencil[14] = rowValues[4];
        stencil[16] = rowValues[5];
        stencil[22] = rowValues[6];
      } else if (rowValues.size() == 27) {
        for(LO k = 0; k < 2; ++k) {
          for(LO j = 0; j < 2; ++j) {
            for(LO i = 0; i < 2; ++i) {
              stencil[k*9 + j*3 + i] = rowValues[k*4 + j*2 + i];
            }
          }
        }
        stencil[2] = rowValues[8];
        for(LO i = 0; i < 4; ++i) {
          stencil[5 + i] = rowValues[9 + i];
        }
        stencil[11] = rowValues[13];
        for(LO i = 0; i < 13; ++i) {
          stencil[14 + i] = rowValues[14 + i];
        }
      }
    } else if((je == 0 || je == elementNodesPerDir[1] - 1) && (ke == 0 || ke == elementNodesPerDir[2] - 1)) {// i-axis edge
      Array<LO> iEdgeStencilOffset(3), iEdgeRowOffset(3);
      if(ke == 0) {
        iEdgeStencilOffset[1] =  0;
        iEdgeRowOffset[1] = 12;
        if(je == 0) {
          iEdgeStencilOffset[0] = 12;
          iEdgeStencilOffset[2] =  9;
          iEdgeRowOffset[2] = 21;
        } else {
          iEdgeStencilOffset[0] = 9;
          iEdgeStencilOffset[2] =  15;
          iEdgeRowOffset[2] = 21;
        }
      } else {
        iEdgeStencilOffset[1] =  18;
        iEdgeRowOffset[1] = 18;
        if(je == 0) {
          iEdgeStencilOffset[0] =  3;
          iEdgeStencilOffset[2] =  0;
          iEdgeRowOffset[2] = 12;
        } else {
          iEdgeStencilOffset[0] = 0;
          iEdgeStencilOffset[2] = 6;
          iEdgeRowOffset[2] = 12;
        }
      }
      for(LO k = 0; k < 2; ++k) {
        for(LO j = 0; j < 2; ++j) {
          for(LO i = 0; i < 3; ++i) {
            stencil[iEdgeStencilOffset[0] + k*9 + j*3 + i] = rowValues[6*k + 3*j + i];
          }
        }
      }
      for(LO i = 0; i < 9; ++i) {
        stencil[iEdgeStencilOffset[1] + i] = rowValues[iEdgeRowOffset[1] + i];
      }
      for(LO i = 0; i < 3; ++i) {
        stencil[iEdgeStencilOffset[2] + i]     = rowValues[iEdgeRowOffset[2] + i];
        stencil[iEdgeStencilOffset[2] + 9 + i] = rowValues[iEdgeRowOffset[2] + 3 + i];
      }
    } else if((ie == 0 || ie == elementNodesPerDir[0] - 1) && (ke == 0 || ke == elementNodesPerDir[2] - 1)) {// j-axis edge
      Array<LO> jEdgeStencilOffset(3), jEdgeRowOffset(3);
      if(ke == 0) {
        jEdgeStencilOffset[1] = 0;
        jEdgeRowOffset[1] = 12;
        if(ie == 0) {
          jEdgeStencilOffset[0] = 10;
          jEdgeStencilOffset[2] = 9;
          jEdgeRowOffset[2] = 21;
        } else {
          jEdgeStencilOffset[0] = 9;
          jEdgeStencilOffset[2] = 11;
          jEdgeRowOffset[2] = 21;
        }
      } else {
        jEdgeStencilOffset[1] = 18;
        jEdgeRowOffset[1] = 18;
        if(ie == 0) {
          jEdgeStencilOffset[0] = 1;
          jEdgeStencilOffset[2] = 0;
          jEdgeRowOffset[2] = 12;
        } else {
          jEdgeStencilOffset[0] = 0;
          jEdgeStencilOffset[2] = 2;
          jEdgeRowOffset[2] = 12;
        }
      }
      for(LO k = 0; k < 2; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 2; ++i) {
            stencil[jEdgeStencilOffset[0] + k*9 + j*3 + i] = rowValues[6*k + 2*j + i];
          }
        }
      }
      for(LO i = 0; i < 9; ++i) {
        stencil[jEdgeStencilOffset[1] + i] = rowValues[jEdgeRowOffset[1] + i];
      }
      for(LO i = 0; i < 3; ++i) {
        stencil[jEdgeStencilOffset[2] + 3*i]     = rowValues[jEdgeRowOffset[2] + i];
        stencil[jEdgeStencilOffset[2] + 9 + 3*i] = rowValues[jEdgeRowOffset[2] + 3 + i];
      }
    } else if((ie == 0 || ie == elementNodesPerDir[0] - 1) && (je == 0 || je == elementNodesPerDir[1] - 1)) {// k-axis edge
      Array<LO> kEdgeStencilOffset(3), kEdgeRowOffset(3);
      if(je == 0) {
        kEdgeStencilOffset[1] = 0;
        kEdgeRowOffset[1] = 12;
        if(ie == 0) {
          kEdgeStencilOffset[0] = 4;
          kEdgeStencilOffset[2] = 3;
          kEdgeRowOffset[2] = 15;
        } else {
          kEdgeStencilOffset[0] = 3;
          kEdgeStencilOffset[2] = 5;
          kEdgeRowOffset[2] = 15;
        }
      } else {
        kEdgeStencilOffset[1] = 6;
        kEdgeRowOffset[1] = 14;
        if(ie == 0) {
          kEdgeStencilOffset[0] = 1;
          kEdgeStencilOffset[2] = 0;
          kEdgeRowOffset[2] = 12;
        } else {
          kEdgeStencilOffset[0] = 0;
          kEdgeStencilOffset[2] = 2;
          kEdgeRowOffset[2] = 12;
        }
      }
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 2; ++j) {
          for(LO i = 0; i < 2; ++i) {
            stencil[kEdgeStencilOffset[0] + k*9 + j*3 + i] = rowValues[4*k + 2*j + i];
          }
        }
      }
      for(LO k = 0; k < 3; ++k) {
        for(LO i = 0; i < 3; ++i) {
          stencil[kEdgeStencilOffset[1] + 9*k + i] = rowValues[kEdgeRowOffset[1] + 5*k + i];
        }
      }
      for(LO i = 0; i < 3; ++i) {
        stencil[kEdgeStencilOffset[2] + 9*i]     = rowValues[kEdgeRowOffset[2] + 5*i];
        stencil[kEdgeStencilOffset[2] + 3 + 9*i] = rowValues[kEdgeRowOffset[2] + 1 + 5*i];
      }
    } else if(ie == 0 || ie == elementNodesPerDir[0] - 1) {// i-axis face
      Array<LO> iFaceStencilOffset(2);
      if(ie == 0) {
        iFaceStencilOffset[0] = 1;
        iFaceStencilOffset[1] = 0;
      } else {
        iFaceStencilOffset[0] = 0;
        iFaceStencilOffset[1] = 2;
      }
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 2; ++i) {
            stencil[iFaceStencilOffset[0] + k*9 + j*3 + i] = rowValues[6*k + 2*j + i];
          }
        }
      }
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          stencil[iFaceStencilOffset[1] + k*9 + j*3] = rowValues[18 + 3*k + j];
        }
      }
    } else if(je == 0 || je == elementNodesPerDir[1] - 1) {// j-axis face
      Array<LO> jFaceStencilOffset(2);
      if(je == 0) {
        jFaceStencilOffset[0] = 3;
        jFaceStencilOffset[1] = 0;
      } else {
        jFaceStencilOffset[0] = 0;
        jFaceStencilOffset[1] = 6;
      }
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 2; ++j) {
          for(LO i = 0; i < 3; ++i) {
            stencil[jFaceStencilOffset[0] + k*9 + j*3 + i] = rowValues[6*k + 3*j + i];
          }
        }
      }
      for(LO k = 0; k < 3; ++k) {
        for(LO i = 0; i < 3; ++i) {
          stencil[jFaceStencilOffset[1] + k*9 + i] = rowValues[18 + 3*k + i];
        }
      }
    } else if(ke == 0 || ke == elementNodesPerDir[2] - 1) {// k-axis face
      Array<LO> kFaceStencilOffset(2);
      if(ke == 0) {
        kFaceStencilOffset[0] =  9;
        kFaceStencilOffset[1] =  0;
      } else {
        kFaceStencilOffset[0] =  0;
        kFaceStencilOffset[1] = 18;
      }
      for(LO k = 0; k < 2; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            stencil[kFaceStencilOffset[0] + k*9 + j*3 + i] = rowValues[9*k + 3*j + i];
          }
        }
      }
      for(LO j = 0; j < 3; ++j) {
        for(LO i = 0; i < 3; ++i) {
          stencil[kFaceStencilOffset[1] + 3*j + i] = rowValues[18 + 3*j + i];
        }
      }
    } else {// the node is "interior" and nothing needs to be reordered
      for(LO i = 0; i < 27; ++i) {
        stencil[i] = rowValues[i];
      }
    }
  } // reorderStencil()

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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetNodeInfo(
                        const LO ie, const LO je, const LO ke,
                        const Array<LO> elementNodesPerDir, int* type, LO& ind) const {
    *type = 0, ind = 0;
    if((ke == 0 || ke == elementNodesPerDir[2]-1)
             && (je == 0 || je == elementNodesPerDir[1]-1)
             && (ie == 0 || ie == elementNodesPerDir[0]-1)) {
      // Corner node
      *type = 0;
      ind  = 4*ke / (elementNodesPerDir[2]-1) + 2*je / (elementNodesPerDir[1]-1) + ie / (elementNodesPerDir[0]-1);
    } else if(((ke == 0 || ke == elementNodesPerDir[2]-1) && (je == 0 || je == elementNodesPerDir[1]-1))
                     || ((ke == 0 || ke == elementNodesPerDir[2]-1) && (ie == 0 || ie == elementNodesPerDir[0]-1))
                     || ((je == 0 || je == elementNodesPerDir[1]-1) && (ie == 0 || ie == elementNodesPerDir[0]-1))) {
      // Edge node
      *type = 1;
      if(ke > 0) {ind += 2*(elementNodesPerDir[0] - 2 + elementNodesPerDir[1] - 2);}
      if(ke == elementNodesPerDir[2] - 1) {ind += 4*(elementNodesPerDir[2] - 2);}
      if((ke == 0) || (ke == elementNodesPerDir[2] - 1)) {
        if(je == 0) {
          ind += ie - 1;
        } else if(je == elementNodesPerDir[1] - 1) {
          ind += 2*(elementNodesPerDir[1] - 2) + elementNodesPerDir[0] - 2 + ie - 1;
        } else {
          ind += elementNodesPerDir[0] - 2 + 2*(je - 1) + ie / (elementNodesPerDir[0] - 1);
        }
      } else {
        ind += 4*(ke - 1) + 2*(je/(elementNodesPerDir[1] - 1)) + ie / (elementNodesPerDir[0] - 1);
      }
    } else if ((ke == 0 || ke == elementNodesPerDir[2]-1)
                     || (je == 0 || je == elementNodesPerDir[1]-1)
                     || (ie == 0 || ie == elementNodesPerDir[0]-1)) {
      // Face node
      *type = 2;
      if(ke == 0) {// current node is on "bottom" face
        ind = (je - 1)*(elementNodesPerDir[0] - 2) + ie - 1;
      } else {
        // add nodes from "bottom" face
        ind += (elementNodesPerDir[1] - 2)*(elementNodesPerDir[0] - 2);
        // add nodes from side faces
        ind += 2*(ke - 1)*(elementNodesPerDir[1] - 2 + elementNodesPerDir[0] - 2);
        if(ke == elementNodesPerDir[2]-1) {// current node is on "top" face
          ind += (je - 1)*(elementNodesPerDir[0] - 2) + ie - 1;
        } else {// current node is on a side face
          if(je == 0) {
            ind += ie - 1;
          } else if(je == elementNodesPerDir[1] - 1) {
            ind += 2*(elementNodesPerDir[1] - 2) + elementNodesPerDir[0] - 2 + ie - 1;
          } else {
            ind += elementNodesPerDir[0] - 2 + 2*(je - 1) + ie / (elementNodesPerDir[0]-1);
          }
        }
      }
    } else {
      // Interior node
      *type = 3;
      ind  = (ke - 1)*(elementNodesPerDir[1] - 2)*(elementNodesPerDir[0] - 2)
        + (je - 1)*(elementNodesPerDir[0] - 2) + ie - 1;
    }
  }

} //namespace MueLu

#define MUELU_BLACKBOXPFACTORY_SHORT
#endif // MUELU_BLACKBOXPFACTORY_DEF_HPP

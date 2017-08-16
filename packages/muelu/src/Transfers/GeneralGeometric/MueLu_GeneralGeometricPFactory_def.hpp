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
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DeclareInput(Level& fineLevel, Level& coarseLevel) const {
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
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel,
                                Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel,
                                Level& coarseLevel) const {
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
    const LO blkSize      = A->GetFixedBlockSize();
    RCP<const Map> rowMap = A->getRowMap();
    GeometricData myGeometry{};

    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords==Teuchos::null, Exceptions::RuntimeError,
                               "Coordinates cannot be accessed from fine level!");
    myGeometry.numDimensions = fineCoords->getNumVectors();
    myGeometry.lNumFineNodes = fineCoords->getLocalLength();

    // Get the number of points in each direction
    if(fineLevel.GetLevelID() == 0) {
      myGeometry.gFineNodesPerDir = fineLevel.Get<Array<GO> >("gNodesPerDim", NoFactory::get());
      myGeometry.lFineNodesPerDir = fineLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
    } else {
      // Loading global number of nodes per diretions
      myGeometry.gFineNodesPerDir = Get<Array<GO> >(fineLevel, "gNodesPerDim");

      // Loading local number of nodes per diretions
      myGeometry.lFineNodesPerDir = Get<Array<LO> >(fineLevel, "lNodesPerDim");
    }
    myGeometry.gNumFineNodes10 = myGeometry.gFineNodesPerDir[1]*myGeometry.gFineNodesPerDir[0];
    myGeometry.gNumFineNodes   = myGeometry.gFineNodesPerDir[2]*myGeometry.gNumFineNodes10;
    myGeometry.lNumFineNodes10 = myGeometry.lFineNodesPerDir[1]*myGeometry.lFineNodesPerDir[0];

    // Get the coarsening rate
    std::string coarsenRate = pL.get<std::string>("Coarsen");
    Teuchos::Array<LO> crates;
    try {
      crates = Teuchos::fromStringToArray<LO>(coarsenRate);
    } catch(const Teuchos::InvalidArrayStringRepresentation e) {
      GetOStream(Errors,-1) << " *** Coarsen must be a string convertible into an array! *** "
                            << std::endl;
      throw e;
    }
    TEUCHOS_TEST_FOR_EXCEPTION((crates.size() > 1) && (crates.size() < myGeometry.numDimensions),
                               Exceptions::RuntimeError,
                               "Coarsen must have at least as many components as the number of"
                               " spatial dimensions in the problem.");
    for(LO i = 0; i < 3; ++i) {
      if(i < myGeometry.numDimensions) {
        if(crates.size()==1) {
          myGeometry.coarseRate[i] = crates[0];
        } else if(crates.size() == myGeometry.numDimensions) {
          myGeometry.coarseRate[i] = crates[i];
        }
      } else {
        myGeometry.coarseRate[i] = 1;
      }
    }

    int interpolationOrder = pL.get<int>("order");
    TEUCHOS_TEST_FOR_EXCEPTION((interpolationOrder < 0) || (interpolationOrder > 1),
                               Exceptions::RuntimeError,
                               "The interpolation order can only be set to 0 or 1.");

    // Get the axis permutation from Global axis to Local axis
    std::string axisPermutation = pL.get<std::string>("axisPermutation");
    Array<LO> mapDirG2L(3), mapDirL2G(3);
    try {
      mapDirG2L = Teuchos::fromStringToArray<LO>(axisPermutation);
    } catch(const Teuchos::InvalidArrayStringRepresentation e) {
      GetOStream(Errors,-1) << " *** axisPermutation must be a string convertible into array! *** "
                            << std::endl;
      throw e;
    }
    for(LO i = 0; i < myGeometry.numDimensions; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(mapDirG2L[i] > myGeometry.numDimensions,
                                 Exceptions::RuntimeError,
                                 "axis permutation values must all be less than"
                                 " the number of spatial dimensions.");
      mapDirL2G[mapDirG2L[i]] = i;
    }

    // Find the offsets for the nodes on the current processor, this tells us the position of the
    // first node on the partition compare to the closest lower coarse point. This information is
    // needed to know who are the coarse nodes that need to be used for interpolation and by
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
    GO maxGlobalIndex, rem;
    RCP<const Map> fineCoordsMap = fineCoords->getMap();
    myGeometry.minGlobalIndex = fineCoordsMap->getMinGlobalIndex();
    maxGlobalIndex = fineCoordsMap->getMaxGlobalIndex();
    if(myGeometry.numDimensions == 1) {
      myGeometry.startIndices[0] = myGeometry.minGlobalIndex;
      myGeometry.offsets[0]= Teuchos::as<LO>(myGeometry.startIndices[0]) % myGeometry.coarseRate[0];

      myGeometry.startIndices[3] = maxGlobalIndex;
      myGeometry.offsets[3]= Teuchos::as<LO>(myGeometry.startIndices[3]) % myGeometry.coarseRate[0];
    } else if(myGeometry.numDimensions == 2) {
      myGeometry.startIndices[1] = myGeometry.minGlobalIndex / myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[1]= Teuchos::as<LO>(myGeometry.startIndices[1]) % myGeometry.coarseRate[1];
      myGeometry.startIndices[0] = myGeometry.minGlobalIndex % myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[0]= Teuchos::as<LO>(myGeometry.startIndices[0]) % myGeometry.coarseRate[0];

      myGeometry.startIndices[4] = maxGlobalIndex / myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[4]= Teuchos::as<LO>(myGeometry.startIndices[4]) % myGeometry.coarseRate[1];
      myGeometry.startIndices[3] = maxGlobalIndex % myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[3]= Teuchos::as<LO>(myGeometry.startIndices[3]) % myGeometry.coarseRate[3];
    } else if(myGeometry.numDimensions == 3) {
      myGeometry.startIndices[2] = myGeometry.minGlobalIndex / myGeometry.gNumFineNodes10;
      rem = myGeometry.minGlobalIndex % myGeometry.gNumFineNodes10;
      myGeometry.offsets[2]= Teuchos::as<LO>(myGeometry.startIndices[2]) % myGeometry.coarseRate[2];
      myGeometry.startIndices[1] = rem / myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[1]= Teuchos::as<LO>(myGeometry.startIndices[1]) % myGeometry.coarseRate[1];
      myGeometry.startIndices[0] = rem % myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[0]= Teuchos::as<LO>(myGeometry.startIndices[0]) % myGeometry.coarseRate[0];

      myGeometry.startIndices[5] = maxGlobalIndex / myGeometry.gNumFineNodes10;
      rem = maxGlobalIndex % myGeometry.gNumFineNodes10;
      myGeometry.offsets[5]= Teuchos::as<LO>(myGeometry.startIndices[5]) % myGeometry.coarseRate[2];
      myGeometry.startIndices[4] = rem / myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[4]= Teuchos::as<LO>(myGeometry.startIndices[4]) % myGeometry.coarseRate[1];
      myGeometry.startIndices[3] = rem % myGeometry.gFineNodesPerDir[0];
      myGeometry.offsets[3]= Teuchos::as<LO>(myGeometry.startIndices[3]) % myGeometry.coarseRate[0];
    }
    /* At this point the local geometrical discovery is over */

    Array<GO> ghostsGIDs;
    GetCoarsePoints(&myGeometry, ghostsGIDs, blkSize);

    // All that is left to do is loop over NCpts and:
    //   - extract coarse points coordiate for coarseCoords
    //   - get coordinates for current stencil computation
    //   - compute current stencil
    //   - compute row and column indices for stencil entries
    RCP<const Map> stridedDomainMapP;
    RCP<Matrix>    P;
    GO nTerms = 8*(myGeometry.lNumFineNodes - myGeometry.lNumCoarseNodes)
      + myGeometry.lNumCoarseNodes;

    GetOStream(Runtime1) << "P size = " << blkSize*myGeometry.gNumFineNodes
                         << " x " << blkSize*myGeometry.gNumCoarseNodes << std::endl;
    GetOStream(Runtime1) << "P Fine   grid : " << myGeometry.gFineNodesPerDir[0] << " -- "
                         << myGeometry.gFineNodesPerDir[1] << " -- "
                         << myGeometry.gFineNodesPerDir[2] << std::endl;
    GetOStream(Runtime1) << "P Coarse grid : " << myGeometry.gCoarseNodesPerDir[0] << " -- "
                         << myGeometry.gCoarseNodesPerDir[1] << " -- "
                         << myGeometry.gCoarseNodesPerDir[2] << std::endl;
    GetOStream(Runtime1) << "P nnz estimate: " << nTerms << std::endl;

    MakeGeneralGeometricP(myGeometry, fineCoords, nTerms, blkSize, stridedDomainMapP,
                          A, P, coarseCoords, ghostsGIDs, interpolationOrder);

    // set StridingInformation of P
    if (A->IsView("stridedMaps") == true) {
      P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), stridedDomainMapP);
    } else {
      P->CreateView("stridedMaps", P->getRangeMap(), stridedDomainMapP);
    }

    // store the transfer operator and node coordinates on coarse level
    Set(coarseLevel, "P", P);
    Set(coarseLevel, "coarseCoordinates", coarseCoords);
    Set<Array<GO> >(coarseLevel, "gCoarseNodesPerDim", myGeometry.gCoarseNodesPerDir);
    Set<Array<LO> >(coarseLevel, "lCoarseNodesPerDim", myGeometry.lCoarseNodesPerDir);

    // rst: null space might get scaled here ... do we care. We could just inject at the cpoints,
    // but I don't feel that this is needed.
    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(P->getDomainMap(),
                                                                 fineNullspace->getNumVectors());
    P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(),
             Teuchos::ScalarTraits<SC>::zero());
    Set(coarseLevel, "Nullspace", coarseNullspace);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DataInterface() const{
    // The goal here is to produce maps that globally labels the mesh lexicographically.
    // These maps will replace the current maps of A, the coordinate vector and the nullspace.
    // Ideally if the user provides the necessary maps then nothing needs to be done, otherwise
    // it could be advantageous to allow the user to register a re-labeling function. Ultimately
    // for common ordering schemes, some re-labeling can be directly implemented here.

  } // DataInterface

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetCoarsePoints(
                                GeometricData* myGeo, Array<GO>& ghostsGIDs,
                                const LO blkSize) const{
    // Check if the partition contains nodes on a boundary, if so that boundary (face, line or
    // point) will not require ghost nodes.
    // Always include a layer if ghost nodes for inner interface since even for faces with coarse
    // nodes, an outward curvature might require ghost nodes to compute a proper interpolation
    // operator. Curvature could be check localy to avoid including extra ghost nodes...
    for(LO i=0; i < 3; ++i) {
      if(i < myGeo->numDimensions &&  myGeo->startIndices[i] != 0) {
        myGeo->ghostInterface[2*i]=true;
        }
      if(i < myGeo->numDimensions
         && (myGeo->startIndices[i]+myGeo->lFineNodesPerDir[i]) != myGeo->gFineNodesPerDir[i]) {
        myGeo->ghostInterface[2*i+1]=true;
      }
    }

    // Here one element can represent either the degenerate case of one node or the more general
    // case of two nodes, i.e. x---x is a 1D element with two nodes and x is a 1D element with one
    // node. This helps generating a 3D space from tensorial products...
    // A good way to handle this would be to generalize the algorithm to take into account the
    // discretization order used in each direction, at least in the FEM sense, since a 0 degree
    // discretization will have a unique node per element. This way 1D discretization can be viewed
    // as a 3D problem with one 0 degree element in the y direction and one 0 degre element in the z
    // direction.
    // !!! Operations below are aftecting both local and global values that have two different   !!!
    // orientations. Orientations can be interchanged using mapDirG2L and mapDirL2G. coarseRate,
    // endRate and offsets are in the global basis, as well as all the variables starting with a g,
    // !!! while the variables starting with an l are in the local basis.                        !!!
    for(LO i = 0; i < 3; ++i) {
      if(i < myGeo->numDimensions) {
        // This array is passed to the RAPFactory and eventually becomes gFineNodePerDir on the next
        // level.
        myGeo->gCoarseNodesPerDir[i] = (myGeo->gFineNodesPerDir[i] - 1) / myGeo->coarseRate[i];
        myGeo->endRate[i] = Teuchos::as<LO>((myGeo->gFineNodesPerDir[i] - 1) % myGeo->coarseRate[i]);
        if(myGeo->endRate[i] == 0) {
          myGeo->endRate[i] = myGeo->coarseRate[i];
          ++myGeo->gCoarseNodesPerDir[i];
        } else {
          myGeo->gCoarseNodesPerDir[i] += 2;
        }
      } else {
        myGeo->endRate[i] = 1;
        myGeo->gCoarseNodesPerDir[i] = 1;
      }
    }

    myGeo->gNumCoarseNodes = myGeo->gCoarseNodesPerDir[0]*myGeo->gCoarseNodesPerDir[1]
      *myGeo->gCoarseNodesPerDir[2];

    for(LO i = 0; i < 3; ++i) {
      if(i < myGeo->numDimensions) {
        // Check whether the partition includes the "end" of the mesh which means that endRate will
        // apply. Also make sure that endRate is not 0 which means that the mesh does not require a
        // particular treatment at the boundaries.
        if( (myGeo->startIndices[i] + myGeo->lFineNodesPerDir[i]) == myGeo->gFineNodesPerDir[i] ) {
          myGeo->lCoarseNodesPerDir[i] = (myGeo->lFineNodesPerDir[i] - myGeo->endRate[i]
                                   + myGeo->offsets[i] - 1) / myGeo->coarseRate[i] + 1;
          if(myGeo->offsets[i] == 0) {++myGeo->lCoarseNodesPerDir[i];}
        } else {
          myGeo->lCoarseNodesPerDir[i] = (myGeo->lFineNodesPerDir[i] + myGeo->offsets[i] - 1)
            / myGeo->coarseRate[i];
          if(myGeo->offsets[i] == 0) {++myGeo->lCoarseNodesPerDir[i];}
        }
      } else {
        myGeo->lCoarseNodesPerDir[i] = 1;
      }
      if(myGeo->lFineNodesPerDir[i] < 1) {myGeo->lCoarseNodesPerDir[i] = 0;}
    }

    // Assuming linear interpolation, each fine point has contribution from 8 coarse points
    // and each coarse point value gets injected.
    // For systems of PDEs we assume that all dofs have the same P operator.
    myGeo->lNumCoarseNodes = myGeo->lCoarseNodesPerDir[0]*myGeo->lCoarseNodesPerDir[1]
      *myGeo->lCoarseNodesPerDir[2];
    LO nTerms = 8*myGeo->lNumFineNodes - 7*myGeo->lNumCoarseNodes;
    nTerms=nTerms*blkSize;

    // For each direction, determine how many ghost points are required.
    LO numGhosts = 0;
    const LO complementaryIndices[3][2] = {{1,2}, {0,2}, {0,1}};
    for(LO i = 0; i < 3; ++i) {
      LO tmp = 0;
      // Check whether a face in direction i needs ghost nodes
      if(myGeo->ghostInterface[2*i] || myGeo->ghostInterface[2*i+1]) {
        if(i == 0) {tmp = myGeo->lCoarseNodesPerDir[1]*myGeo->lCoarseNodesPerDir[2];}
        if(i == 1) {tmp = myGeo->lCoarseNodesPerDir[0]*myGeo->lCoarseNodesPerDir[2];}
        if(i == 2) {tmp = myGeo->lCoarseNodesPerDir[0]*myGeo->lCoarseNodesPerDir[1];}
      }
      // If both faces in direction i need nodes, double the number of ghost nodes
      if(myGeo->ghostInterface[2*i] && myGeo->ghostInterface[2*i+1]) {tmp = 2*tmp;}
      numGhosts += tmp;

      // The corners and edges need to be checked in 2D / 3D to add more ghosts nodes
      for(LO j = 0; j < 2; ++j) {
        for(LO k = 0; k < 2; ++k) {
          // Check if two adjoining faces need ghost nodes and then add their common edge
          if(myGeo->ghostInterface[2*complementaryIndices[i][0]+j]
             && myGeo->ghostInterface[2*complementaryIndices[i][1]+k]) {
            numGhosts += myGeo->lCoarseNodesPerDir[i];
            // Add corners if three adjoining faces need ghost nodes,
            // but add them only once! Hence when i == 0.
            if(myGeo->ghostInterface[2*i] && (i == 0)) { numGhosts += 1; }
            if(myGeo->ghostInterface[2*i+1] && (i == 0)) { numGhosts += 1; }
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
    ghostsGIDs.resize(numGhosts);
    LO countGhosts = 0;
    // Get the GID of the first point on the current partition.
    GO startingGID = myGeo->minGlobalIndex;
    Array<GO> startingIndices(3);
    // We still want ghost nodes even if have with a 0 offset,
    // except when we are on a boundary
    if(myGeo->ghostInterface[4] && (myGeo->offsets[2] == 0)) {
      if(myGeo->startIndices[2] + myGeo->coarseRate[2] > myGeo->gFineNodesPerDir[2]) {
        startingGID -= myGeo->endRate[2]*myGeo->gNumFineNodes10;
      } else {
        startingGID -= myGeo->coarseRate[2]*myGeo->gNumFineNodes10;
      }
    } else {
      startingGID -= myGeo->offsets[2]*myGeo->gNumFineNodes10;
    }
    if(myGeo->ghostInterface[2] && (myGeo->offsets[1] == 0)) {
      if(myGeo->startIndices[1] + myGeo->coarseRate[1] > myGeo->gFineNodesPerDir[1]) {
        startingGID -= myGeo->endRate[1]*myGeo->gFineNodesPerDir[0];
      } else {
        startingGID -= myGeo->coarseRate[1]*myGeo->gFineNodesPerDir[0];
      }
    } else {
      startingGID -= myGeo->offsets[1]*myGeo->gFineNodesPerDir[0];
    }
    if(myGeo->ghostInterface[0] && (myGeo->offsets[0] == 0)) {
      if(myGeo->startIndices[0] + myGeo->coarseRate[0] > myGeo->gFineNodesPerDir[0]) {
        startingGID -= myGeo->endRate[0];
      } else {
        startingGID -= myGeo->coarseRate[0];
      }
    } else {
      startingGID -= myGeo->offsets[0];
    }

    { // scope for tmp
      GO tmp;
      startingIndices[2] = startingGID / myGeo->gNumFineNodes10;
      tmp = startingGID % myGeo->gNumFineNodes10;
      startingIndices[1] = tmp / myGeo->gFineNodesPerDir[0];
      startingIndices[0] = tmp % myGeo->gFineNodesPerDir[0];
    }

    GO ghostOffset[3] = {0, 0, 0};
    LO lengthZero  = myGeo->lCoarseNodesPerDir[0];
    LO lengthOne   = myGeo->lCoarseNodesPerDir[1];
    LO lengthTwo   = myGeo->lCoarseNodesPerDir[2];
    if(myGeo->ghostInterface[0]) {++lengthZero;}
    if(myGeo->ghostInterface[1]) {++lengthZero;}
    if(myGeo->ghostInterface[2]) {++lengthOne;}
    if(myGeo->ghostInterface[3]) {++lengthOne;}
    if(myGeo->ghostInterface[4]) {++lengthTwo;}
    if(myGeo->ghostInterface[5]) {++lengthTwo;}

    // First check the bottom face as it will have the lowest GIDs
    if(myGeo->ghostInterface[4]) {
      ghostOffset[2] = startingGID;
      for(LO j = 0; j < lengthOne; ++j) {
        if((j == lengthOne-1)
           && (startingIndices[1] + j*myGeo->coarseRate[1] + 1 > myGeo->gFineNodesPerDir[1])) {
          ghostOffset[1] = ((j-1)*myGeo->coarseRate[1] + myGeo->endRate[1])*myGeo->gFineNodesPerDir[0];
        } else {
          ghostOffset[1] = j*myGeo->coarseRate[1]*myGeo->gFineNodesPerDir[0];
        }
        for(LO k = 0; k < lengthZero; ++k) {
          if((k == lengthZero-1)
             && (startingIndices[0] + k*myGeo->coarseRate[0] + 1 > myGeo->gFineNodesPerDir[0])) {
            ghostOffset[0] = (k-1)*myGeo->coarseRate[0] + myGeo->endRate[0];
          } else {
            ghostOffset[0] = k*myGeo->coarseRate[0];
          }
          // If the partition includes a changed rate at the edge, ghost nodes need to be picked
          // carefully. This if statement is repeated each time a ghost node is selected.
          ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
          ++countGhosts;
        }
      }
    }

    // Sweep over the myGeo->lCoarseNodesPerDir[2] coarse layers in direction 2 and gather necessary
    // ghost nodes located on these layers.
    for(LO i = 0; i < lengthTwo; ++i) {
      // Exclude the cases where ghost nodes exists on the faces in directions 2, these faces are
      // swept seperatly for efficiency.
      if( !((i == lengthTwo-1) && myGeo->ghostInterface[5])
          && !((i == 0) && myGeo->ghostInterface[4]) ) {
        // Set the ghostOffset in direction 2 taking into account a possible endRate different
        // from the regular coarseRate.
        if( (i == lengthTwo-1) && !myGeo->ghostInterface[5] ) {
          ghostOffset[2] = startingGID
            + ((i-1)*myGeo->coarseRate[2] + myGeo->endRate[2])*myGeo->gNumFineNodes10;
        } else {
          ghostOffset[2] = startingGID + i*myGeo->coarseRate[2]*myGeo->gNumFineNodes10;
        }
        for(LO j = 0; j < lengthOne; ++j) {
          if( (j == 0) && myGeo->ghostInterface[2] ) {
            for(LO k = 0; k < lengthZero; ++k) {
              if((k == lengthZero-1)
                 && (startingIndices[0] + k*myGeo->coarseRate[0] + 1 > myGeo->gFineNodesPerDir[0])) {
                if(k == 0) {
                  ghostOffset[0] = myGeo->endRate[0];
                } else {
                  ghostOffset[0] = (k-1)*myGeo->coarseRate[0] + myGeo->endRate[0];
                }
              } else {
                ghostOffset[0] = k*myGeo->coarseRate[0];
              }
              // In this case j == 0 so ghostOffset[1] is zero and can be ignored
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[0];
              ++countGhosts;
            }
          } else if( (j == lengthOne-1) && myGeo->ghostInterface[3] ) {
            // Set the ghostOffset in direction 1 taking into account a possible endRate different
            // from the regular myGeo->coarseRate.
            if((j == lengthOne-1)
               && (startingIndices[1] + j*myGeo->coarseRate[1] + 1 > myGeo->gFineNodesPerDir[1])) {
              ghostOffset[1] = ((j-1)*myGeo->coarseRate[1] + myGeo->endRate[1])
                *myGeo->gFineNodesPerDir[0];
            } else {
              ghostOffset[1] = j*myGeo->coarseRate[1]*myGeo->gFineNodesPerDir[0];
            }
            for(LO k = 0; k < lengthZero; ++k) {
              if((k == lengthZero-1)
                 && (startingIndices[0] + k*myGeo->coarseRate[0] + 1 > myGeo->gFineNodesPerDir[0])) {
                ghostOffset[0] = (k-1)*myGeo->coarseRate[0] + myGeo->endRate[0];
              } else {
                ghostOffset[0] = k*myGeo->coarseRate[0];
              }
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
              ++countGhosts;
            }
          } else {
            // Set ghostOffset[1] for side faces sweep
            if((j == lengthOne-1)
               && (startingIndices[1] + j*myGeo->coarseRate[1] + 1 > myGeo->gFineNodesPerDir[1])) {
              ghostOffset[1] = ((j-1)*myGeo->coarseRate[1] + myGeo->endRate[1])
                *myGeo->gFineNodesPerDir[0];
            } else {
              ghostOffset[1] = j*myGeo->coarseRate[1]*myGeo->gFineNodesPerDir[0];
            }

            // Set ghostOffset[0], ghostsGIDs and countGhosts
            if(myGeo->ghostInterface[0]) { // In that case ghostOffset[0]==0, so we can ignore it
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1];
              ++countGhosts;
            }
            if(myGeo->ghostInterface[1]) { // Grab ghost point at the end of direction 0.
              if((startingIndices[0] + (lengthZero - 1)*myGeo->coarseRate[0])
                 > myGeo->gFineNodesPerDir[0] - 1) {
                if(lengthZero > 2) {
                  ghostOffset[0] = (lengthZero-2)*myGeo->coarseRate[0] + myGeo->endRate[0];
                } else {
                  ghostOffset[0] = myGeo->endRate[0];
                }
              } else {
                ghostOffset[0] = (lengthZero-1)*myGeo->coarseRate[0];
              }
              ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
              ++countGhosts;
            }
          }
        }
      }
    }

    // Finally check the top face
    if(myGeo->ghostInterface[5]) {
      if( startingIndices[2] + (lengthTwo-1)*myGeo->coarseRate[2] + 1 > myGeo->gFineNodesPerDir[2] ) {
        ghostOffset[2] = startingGID
          + ((lengthTwo-2)*myGeo->coarseRate[2] + myGeo->endRate[2])*myGeo->gNumFineNodes10;
      } else {
        ghostOffset[2] = startingGID + (lengthTwo-1)*myGeo->coarseRate[2]*myGeo->gNumFineNodes10;
      }
      for(LO j = 0; j < lengthOne; ++j) {
        if((j == lengthOne-1)
           && (startingIndices[1] + j*myGeo->coarseRate[1] + 1 > myGeo->gFineNodesPerDir[1])) {
          ghostOffset[1] = ((j-1)*myGeo->coarseRate[1] + myGeo->endRate[1])*myGeo->gFineNodesPerDir[0];
        } else {
          ghostOffset[1] = j*myGeo->coarseRate[1]*myGeo->gFineNodesPerDir[0];
        }
        for(LO k = 0; k < lengthZero; ++k) {
          if((k == lengthZero-1)
             && (startingIndices[0] + k*myGeo->coarseRate[0] + 1 > myGeo->gFineNodesPerDir[0])) {
            ghostOffset[0] = (k-1)*myGeo->coarseRate[0] + myGeo->endRate[0];
          } else {
            ghostOffset[0] = k*myGeo->coarseRate[0];
          }
          ghostsGIDs[countGhosts] = ghostOffset[2] + ghostOffset[1] + ghostOffset[0];
          ++countGhosts;
        }
      }
    }
  } // GetCoarsePoint

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MakeGeneralGeometricP(
                                GeometricData myGeo,
                                const RCP<Xpetra::MultiVector<double,LO,GO,Node> >& fineCoords,
                                const LO nnzP, const LO dofsPerNode,
                                RCP<const Map>& stridedDomainMapP,RCP<Matrix> & Amat,RCP<Matrix>& P,
                                RCP<Xpetra::MultiVector<double,LO,GO,Node> >& coarseCoords,
                                Array<GO> ghostsGIDs, int interpolationOrder) const {

    /* On termination, return the number of local prolongator columns owned by
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

    using xdMV                 = Xpetra::MultiVector<double,LO,GO,NO>;
    Xpetra::global_size_t OTI  = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    LO lNumGhostNodes   = ghostsGIDs.size();
    GO numGloCols       = dofsPerNode*myGeo.gNumCoarseNodes;

    // Build the required column map for the prolongator operator.
    // This requies to find the GIDs of the coarse nodes on the coarse mesh,
    // including the ghost nodes, on the local partition.
    // Note: the ordering of the coarse nodes is assumed
    // to be the same as the ordering of the fine nodes.

    GO gStartIndices[3];
    RCP<const Map> fineCoordsMap = fineCoords->getMap();
    // Compute the global indices of the first node on the partition.
    { // Scope for dummy
      gStartIndices[2] = fineCoordsMap->getMinGlobalIndex()
        / (myGeo.gNumFineNodes10);
      GO dummy         = fineCoordsMap->getMinGlobalIndex()
        % (myGeo.gNumFineNodes10);
      gStartIndices[1] = dummy / myGeo.gFineNodesPerDir[0];
      gStartIndices[0] = dummy % myGeo.gFineNodesPerDir[0];
    }

    Array<GO> gCoarseNodesGIDs(myGeo.lNumCoarseNodes);
    LO currentNode, offset2, offset1, offset0;
    // Find the GIDs of the coarse nodes on the partition.
    for(LO ind2 = 0; ind2 < myGeo.lCoarseNodesPerDir[2]; ++ind2) {
      if(myGeo.offsets[2] == 0) {
        offset2 = gStartIndices[2];
      } else {
        if(gStartIndices[2] + myGeo.endRate[2] - myGeo.offsets[2] == myGeo.gFineNodesPerDir[2] - 1){
          offset2 = gStartIndices[2] + myGeo.endRate[2] - myGeo.offsets[2];
        } else {
          offset2 = gStartIndices[2] + myGeo.coarseRate[2] - myGeo.offsets[2];
        }
      }
      if(offset2 + ind2*myGeo.coarseRate[2] > myGeo.gFineNodesPerDir[2] - 1) {
        offset2 += (ind2 - 1)*myGeo.coarseRate[2] + myGeo.endRate[2];
      } else {
        offset2 += ind2*myGeo.coarseRate[2];
      }
      offset2 = offset2*myGeo.gNumFineNodes10;

      for(LO ind1 = 0; ind1 < myGeo.lCoarseNodesPerDir[1]; ++ind1) {
        if(myGeo.offsets[1] == 0) {
          offset1 = gStartIndices[1];
        } else {
          if(gStartIndices[1] + myGeo.endRate[1] - myGeo.offsets[1] == myGeo.gFineNodesPerDir[1]-1){
            offset1 = gStartIndices[1] + myGeo.endRate[1] - myGeo.offsets[1];
          } else {
            offset1 = gStartIndices[1] + myGeo.coarseRate[1] - myGeo.offsets[1];
          }
        }
        if(offset1 + ind1*myGeo.coarseRate[1] > myGeo.gFineNodesPerDir[1] - 1) {
          offset1 += (ind1 - 1)*myGeo.coarseRate[1] + myGeo.endRate[1];
        } else {
          offset1 += ind1*myGeo.coarseRate[1];
        }
        offset1 = offset1*myGeo.gFineNodesPerDir[0];
        for(LO ind0 = 0; ind0 < myGeo.lCoarseNodesPerDir[0]; ++ind0) {
          offset0 = gStartIndices[0] - myGeo.offsets[0];
          if(myGeo.offsets[0] == 0) {
            offset0 += ind0*myGeo.coarseRate[0];
          } else {
            offset0 += (ind0 + 1)*myGeo.coarseRate[0];
          }
          if(offset0 > myGeo.gFineNodesPerDir[0] - 1) {
            offset0 += myGeo.endRate[0] - myGeo.coarseRate[0];
          }

          currentNode = ind2*myGeo.lCoarseNodesPerDir[1]*myGeo.lCoarseNodesPerDir[0]
                      + ind1*myGeo.lCoarseNodesPerDir[0]
                      + ind0;
          gCoarseNodesGIDs[currentNode] = offset2 + offset1 + offset0;
        }
      }
    }

    // Actual loop over all the coarse/ghost nodes to find their index on the coarse mesh
    // and the corresponding dofs that will need to be added to colMapP.
    Array<GO> colGIDs(dofsPerNode*(myGeo.lNumCoarseNodes+lNumGhostNodes));
    Array<GO> coarseNodesGIDs(myGeo.lNumCoarseNodes);
    GO fineNodesPerCoarseSlab    = myGeo.coarseRate[2]*myGeo.gNumFineNodes10;
    GO fineNodesEndCoarseSlab    = myGeo.endRate[2]*myGeo.gNumFineNodes10;
    GO fineNodesPerCoarsePlane   = myGeo.coarseRate[1]*myGeo.gFineNodesPerDir[0];
    GO fineNodesEndCoarsePlane   = myGeo.endRate[1]*myGeo.gFineNodesPerDir[0];
    GO coarseNodesPerCoarseLayer = myGeo.gCoarseNodesPerDir[1]*myGeo.gCoarseNodesPerDir[0];
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
      for(LO col = 0; col < myGeo.lNumCoarseNodes; ++col) {
        if((myGeo.endRate[2] != myGeo.coarseRate[2])
           && (gCoarseNodesGIDs[col] > (myGeo.gCoarseNodesPerDir[2] - 2)
               *fineNodesPerCoarseSlab + fineNodesEndCoarseSlab - 1)){
          tmpInds[2] = gCoarseNodesGIDs[col] / fineNodesPerCoarseSlab + 1;
          tmpVars[0] = gCoarseNodesGIDs[col] - (tmpInds[2] - 1)*fineNodesPerCoarseSlab
            - fineNodesEndCoarseSlab;
        } else {
          tmpInds[2] = gCoarseNodesGIDs[col] / fineNodesPerCoarseSlab;
          tmpVars[0] = gCoarseNodesGIDs[col] % fineNodesPerCoarseSlab;
        }
        if((myGeo.endRate[1] != myGeo.coarseRate[1])
           && (tmpVars[0] > (myGeo.gCoarseNodesPerDir[1] - 2)
               *fineNodesPerCoarsePlane + fineNodesEndCoarsePlane-1)){
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane + 1;
          tmpVars[1] = tmpVars[0] - (tmpInds[1] - 1)*fineNodesPerCoarsePlane
            - fineNodesEndCoarsePlane;
        } else {
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane;
          tmpVars[1] = tmpVars[0] % fineNodesPerCoarsePlane;
        }
        if(tmpVars[1] == myGeo.gFineNodesPerDir[0] - 1) {
          tmpInds[0] = myGeo.gCoarseNodesPerDir[0] - 1;
        } else {
          tmpInds[0] = tmpVars[1] / myGeo.coarseRate[0];
        }
        gInd[2] = col / (myGeo.lCoarseNodesPerDir[1]*myGeo.lCoarseNodesPerDir[0]);
        tmp     = col % (myGeo.lCoarseNodesPerDir[1]*myGeo.lCoarseNodesPerDir[0]);
        gInd[1] = tmp / myGeo.lCoarseNodesPerDir[0];
        gInd[0] = tmp % myGeo.lCoarseNodesPerDir[0];
        lCol = gInd[2]*(myGeo.lCoarseNodesPerDir[1]*myGeo.lCoarseNodesPerDir[0])
          + gInd[1]*myGeo.lCoarseNodesPerDir[0] + gInd[0];
        gCoarseNodeOnCoarseGridGID = tmpInds[2]*coarseNodesPerCoarseLayer
          + tmpInds[1]*myGeo.gCoarseNodesPerDir[0] + tmpInds[0];
        coarseNodesGIDs[lCol] = gCoarseNodeOnCoarseGridGID;
        for(LO dof = 0; dof < dofsPerNode; ++dof) {
          colGIDs[dofsPerNode*lCol + dof] = dofsPerNode*gCoarseNodeOnCoarseGridGID + dof;
        }
      }
      // Now loop over the ghost nodes of the partition to add them to colGIDs
      // since they will need to be included in the column map of P
      for(LO col = myGeo.lNumCoarseNodes; col < myGeo.lNumCoarseNodes + lNumGhostNodes; ++col) {
        if((myGeo.endRate[2] != myGeo.coarseRate[2])
           && (ghostsGIDs[ghostsPermut[col - myGeo.lNumCoarseNodes]] >
               (myGeo.gCoarseNodesPerDir[2] - 2)*fineNodesPerCoarseSlab
               + fineNodesEndCoarseSlab - 1)) {
          tmpInds[2] = ghostsGIDs[ghostsPermut[col - myGeo.lNumCoarseNodes]]
            / fineNodesPerCoarseSlab + 1;
          tmpVars[0] = ghostsGIDs[ghostsPermut[col - myGeo.lNumCoarseNodes]]
            - (tmpInds[2] - 1)*fineNodesPerCoarseSlab - fineNodesEndCoarseSlab;
        } else {
          tmpInds[2] = ghostsGIDs[ghostsPermut[col - myGeo.lNumCoarseNodes]]
            / fineNodesPerCoarseSlab;
          tmpVars[0] = ghostsGIDs[ghostsPermut[col - myGeo.lNumCoarseNodes]]
            % fineNodesPerCoarseSlab;
        }
        if((myGeo.endRate[1] != myGeo.coarseRate[1])
           && (tmpVars[0] > (myGeo.gCoarseNodesPerDir[1] - 2)
               *fineNodesPerCoarsePlane + fineNodesEndCoarsePlane-1)){
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane + 1;
          tmpVars[1] = tmpVars[0] - (tmpInds[1] - 1)*fineNodesPerCoarsePlane
            - fineNodesEndCoarsePlane;
        } else {
          tmpInds[1] = tmpVars[0] / fineNodesPerCoarsePlane;
          tmpVars[1] = tmpVars[0] % fineNodesPerCoarsePlane;
        }
        if(tmpVars[1] == myGeo.gFineNodesPerDir[0] - 1) {
          tmpInds[0] = myGeo.gCoarseNodesPerDir[0] - 1;
        } else {
          tmpInds[0] = tmpVars[1] / myGeo.coarseRate[0];
        }
        gCoarseNodeOnCoarseGridGID = tmpInds[2]*coarseNodesPerCoarseLayer
          + tmpInds[1]*myGeo.gCoarseNodesPerDir[0] + tmpInds[0];
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
                                                                    colGIDs.view(0, dofsPerNode
                                                                                 *myGeo.lNumCoarseNodes),
                                                                    rowMapP->getIndexBase(),
                                                                    rowMapP->getComm());

    RCP<const Map> colMapP = Xpetra::MapFactory<LO,GO,NO>::Build(rowMapP->lib(),OTI,
                                                                 colGIDs.view(0, colGIDs.size()),
                                                                 rowMapP->getIndexBase(),
                                                                 rowMapP->getComm());

    std::vector<size_t> strideInfo(1);
    strideInfo[0] = dofsPerNode;
    stridedDomainMapP = Xpetra::StridedMapFactory<LO,GO,NO>::Build(domainMapP, strideInfo);

    // Build the map for the coarse level coordinates
    RCP<const Map> coarseCoordsMap = MapFactory::Build (fineCoordsMap->lib(),
                                                        myGeo.gNumCoarseNodes,
                                                        coarseNodesGIDs.view(0, myGeo.lNumCoarseNodes),
                                                        fineCoordsMap->getIndexBase(),
                                                        rowMapP->getComm());

    // Do the actual import using the fineCoordsMap
    RCP<const Map> ghostMap = Xpetra::MapFactory<LO,GO,NO>::Build(fineCoordsMap->lib(), OTI,
                                                                  ghostsGIDs.view(0, lNumGhostNodes),
                                                                  fineCoordsMap->getIndexBase(),
                                                                  rowMapP->getComm());
    RCP<const Import> importer = ImportFactory::Build(fineCoordsMap, ghostMap);
    RCP<xdMV> ghostCoords = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(ghostMap,
                                                                               myGeo.numDimensions);
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
    // That array will eventually be used during the construction of coarseCoords, the MultiVector
    // of coarse nodes coordinates that we store on the coarse level.
    ArrayRCP< ArrayRCP<double> > lFineCoords(3);
    ArrayRCP< ArrayRCP<double> > lGhostCoords(3);
    RCP<const Map> ghostCoordsMap = ghostCoords->getMap();
    RCP<Xpetra::Vector<double,LO,GO,NO> > zeros
      = Xpetra::VectorFactory<double,LO,GO,NO>::Build(fineCoordsMap, true);
    RCP<Xpetra::Vector<double,LO,GO,NO> > ghostZeros
      = Xpetra::VectorFactory<double,LO,GO,NO>::Build(ghostCoordsMap, true);

    // Build the MultiVector holding the coarse grid points coordinates.
    coarseCoords = Xpetra::MultiVectorFactory<double,LO,GO,Node>::Build(coarseCoordsMap,
                                              Teuchos::as<size_t>(myGeo.numDimensions));
    ArrayRCP<double> xCoarseNodes; ArrayRCP<double> yCoarseNodes; ArrayRCP<double> zCoarseNodes;

    if(myGeo.numDimensions==1) {
      lFineCoords[0] = fineCoords->getDataNonConst(0);
      lFineCoords[1] = zeros->getDataNonConst(0);
      lFineCoords[2] = zeros->getDataNonConst(0);
      lGhostCoords[0] = ghostCoords->getDataNonConst(0);
      lGhostCoords[1] = ghostZeros->getDataNonConst(0);
      lGhostCoords[2] = ghostZeros->getDataNonConst(0);

      xCoarseNodes = coarseCoords->getDataNonConst(0);
    } else if(myGeo.numDimensions==2) {
      lFineCoords[0] = fineCoords->getDataNonConst(0);
      lFineCoords[1] = fineCoords->getDataNonConst(1);
      lFineCoords[2] = zeros->getDataNonConst(0);
      lGhostCoords[0] = ghostCoords->getDataNonConst(0);
      lGhostCoords[1] = ghostCoords->getDataNonConst(1);
      lGhostCoords[2] = ghostZeros->getDataNonConst(0);

      xCoarseNodes= coarseCoords->getDataNonConst(0);
      yCoarseNodes= coarseCoords->getDataNonConst(1);
    } else if(myGeo.numDimensions==3) {
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
    for(LO i = 0; i < myGeo.lNumFineNodes; ++i) {

      // Get point indices on fine grid
      {
        std::div_t tmp;
        tmp=std::div(i,myGeo.lFineNodesPerDir[0]*myGeo.lFineNodesPerDir[1]);
        indices[0][2] = tmp.quot;
        tmp=std::div(tmp.rem,myGeo.lFineNodesPerDir[0]);
        indices[0][1] = tmp.quot;
        indices[0][0] = tmp.rem;

        // Get ref point indices on coarse grid
        tmp=std::div(indices[0][0] - myGeo.offsets[0],myGeo.coarseRate[0]);
        indices[1][0] = tmp.quot;
        tmp=std::div(indices[0][1] - myGeo.offsets[1],myGeo.coarseRate[1]);
        indices[1][1] = tmp.quot;
        tmp=std::div(indices[0][2] - myGeo.offsets[2],myGeo.coarseRate[2]);
        indices[1][2] = tmp.quot;
      }

      // location "flags" indicate if the current node is on a coarse
      // face, edge or node.
      indices[2][0] = (indices[0][0] + myGeo.offsets[0]) % myGeo.coarseRate[0];
      indices[2][1] = (indices[0][1] + myGeo.offsets[1]) % myGeo.coarseRate[1];
      indices[2][2] = (indices[0][2] + myGeo.offsets[2]) % myGeo.coarseRate[2];

      // Get indices of ref point on fine grid
      indices[3][0] = indices[1][0]*myGeo.coarseRate[0];
      indices[3][1] = indices[1][1]*myGeo.coarseRate[1];
      indices[3][2] = indices[1][2]*myGeo.coarseRate[2];

      if( (indices[0][0] == myGeo.lFineNodesPerDir[0]-1) && !myGeo.ghostInterface[1] ) {
        indices[1][0] = myGeo.lCoarseNodesPerDir[0]-1;
        indices[2][0] = 0;
        indices[3][0] = myGeo.lFineNodesPerDir[0]-1;
      }
      if( (indices[0][1] == myGeo.lFineNodesPerDir[1]-1) && !myGeo.ghostInterface[3] ) {
        indices[1][1] = myGeo.lCoarseNodesPerDir[1]-1;
        indices[2][1] = 0;
        indices[3][1] = myGeo.lFineNodesPerDir[1]-1;
      }
      if( (indices[0][2] == myGeo.lFineNodesPerDir[2]-1) && !myGeo.ghostInterface[5] ) {
        indices[1][2] = myGeo.lCoarseNodesPerDir[2]-1;
        indices[2][2] = 0;
        indices[3][2] = myGeo.lFineNodesPerDir[2]-1;
      }

      Array<GO> currentNodeIndices(3);
      currentNodeIndices[0] = gStartIndices[0] + Teuchos::as<GO>(indices[0][0]);
      currentNodeIndices[1] = gStartIndices[1] + Teuchos::as<GO>(indices[0][1]);
      currentNodeIndices[2] = gStartIndices[2] + Teuchos::as<GO>(indices[0][2]);

      // Assuming lexicographic indexing the coarse nodes forming a prism
      // around fine node "i" are selected and store them in connec.
      // Two tricky things to be careful about:
      //    - are we using coarseRate or endRate?
      //      --> check indices and set rate correctly
      //    - are we on the east, north or top face?
      //      --> if so fix firstCoarseNodeIndex to make sure there is
      //          a node to the east, north and top of firstCoarseNodeIndex
      LO rate[3];
      if(currentNodeIndices[0] >= myGeo.gFineNodesPerDir[0] - myGeo.endRate[0] - 1) {
        rate[0] = myGeo.endRate[0];
      } else {
        rate[0] = myGeo.coarseRate[0];
      }
      if(currentNodeIndices[1] >= myGeo.gFineNodesPerDir[1] - myGeo.endRate[1] - 1) {
        rate[1] = myGeo.endRate[1];
      } else {
        rate[1] = myGeo.coarseRate[1];
      }
      if(currentNodeIndices[2] >= myGeo.gFineNodesPerDir[2] - myGeo.endRate[2] - 1) {
        rate[2] = myGeo.endRate[2];
      } else {
        rate[2] = myGeo.coarseRate[2];
      }
      if(myGeo.numDimensions < 3) { rate[2] = 0;}
      if(myGeo.numDimensions < 2) { rate[1] = 0;}

      // We need to check whether we are on the edge of the mesh in which case we need to adjust
      // the coarse nodes
      firstCoarseNodeIndex = 0;
      Array<GO> firstCoarseNodeIndices(3); // These are fine grid indices
      if((currentNodeIndices[2] == myGeo.gFineNodesPerDir[2] -1)
         && (myGeo.endRate[2] == myGeo.coarseRate[2])) {
        // If we are on the last node and have a endRate == coarseRate we need to take the coarse
        // node below the current node
        firstCoarseNodeIndices[2] = ((currentNodeIndices[2] / myGeo.coarseRate[2]) - 1)
          * myGeo.coarseRate[2];
      } else {
        firstCoarseNodeIndices[2] = (currentNodeIndices[2] / myGeo.coarseRate[2])
          * myGeo.coarseRate[2];
      }
      if((currentNodeIndices[1] == myGeo.gFineNodesPerDir[1] -1)
         && (myGeo.endRate[1] == myGeo.coarseRate[1])) {
        firstCoarseNodeIndices[1] = ((currentNodeIndices[1] / myGeo.coarseRate[1]) - 1)
          * myGeo.coarseRate[1];
      } else {
        firstCoarseNodeIndices[1] = (currentNodeIndices[1] / myGeo.coarseRate[1])
          * myGeo.coarseRate[1];
      }
      if((currentNodeIndices[0] == myGeo.gFineNodesPerDir[0] -1)
         && (myGeo.endRate[0] == myGeo.coarseRate[0])) {
        firstCoarseNodeIndices[0] = ((currentNodeIndices[0] / myGeo.coarseRate[0]) - 1)
          * myGeo.coarseRate[0];
      } else {
        firstCoarseNodeIndices[0] = (currentNodeIndices[0] / myGeo.coarseRate[0])
          * myGeo.coarseRate[0];
      }
      firstCoarseNodeIndex += firstCoarseNodeIndices[2]*myGeo.gNumFineNodes10
        + firstCoarseNodeIndices[1]*myGeo.gFineNodesPerDir[0] + firstCoarseNodeIndices[0];

      GO firstCoarseNodeOnCoarseGridIndex;
      {
        GO tmpInds[4];
        tmpInds[2] = firstCoarseNodeIndex / fineNodesPerCoarseSlab;
        tmpInds[3] = firstCoarseNodeIndex % fineNodesPerCoarseSlab;
        tmpInds[1] = tmpInds[3] / fineNodesPerCoarsePlane;
        tmpInds[0] = (tmpInds[3] % fineNodesPerCoarsePlane) / myGeo.coarseRate[0];
        firstCoarseNodeOnCoarseGridIndex = tmpInds[2]*coarseNodesPerCoarseLayer
          + tmpInds[1]*myGeo.gCoarseNodesPerDir[0] + tmpInds[0];
      }

      LO coarseGhosts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      if( myGeo.ghostInterface[0] ) {
        if( ((indices[0][0] < rate[0] - myGeo.offsets[0]) && myGeo.offsets[0] != 0)
            || (((currentNodeIndices[0] == myGeo.gFineNodesPerDir[0] -1)
                 && (indices[0][0] < rate[0] - myGeo.offsets[0] + 1)) && myGeo.offsets[0] != 0)
            || ((currentNodeIndices[0] == myGeo.gFineNodesPerDir[0] -1)
                && myGeo.lFineNodesPerDir[0] == 1) ) {
          coarseGhosts[0] = 1;
          coarseGhosts[2] = 1;
          coarseGhosts[4] = 1;
          coarseGhosts[6] = 1;
        }
      }
      if(myGeo.ghostInterface[1]
         && (indices[0][0] > myGeo.lFineNodesPerDir[0] - myGeo.offsets[3] - 2)) {
        coarseGhosts[1] = 1;
        coarseGhosts[3] = 1;
        coarseGhosts[5] = 1;
        coarseGhosts[7] = 1;
      }
      if( myGeo.ghostInterface[2] ) {
        if( ((indices[0][1] < rate[1] - myGeo.offsets[1]) && myGeo.offsets[1] != 0)
            || (((currentNodeIndices[1] == myGeo.gFineNodesPerDir[1] -1)
                 && (indices[0][1] < rate[1] - myGeo.offsets[1] + 1)) && myGeo.offsets[1] != 0)
            || ((currentNodeIndices[1] == myGeo.gFineNodesPerDir[1] -1)
                && myGeo.lFineNodesPerDir[1] == 1) ) {
          coarseGhosts[0] = 1;
          coarseGhosts[1] = 1;
          coarseGhosts[4] = 1;
          coarseGhosts[5] = 1;
        }
      }
      if( myGeo.ghostInterface[3] && (indices[0][1] > myGeo.lFineNodesPerDir[1]
                                - myGeo.offsets[4] - 2) ) {
        coarseGhosts[2] = 1;
        coarseGhosts[3] = 1;
        coarseGhosts[6] = 1;
        coarseGhosts[7] = 1;
      }
      if( myGeo.ghostInterface[4] ) {
        if( ((indices[0][2] < rate[2] - myGeo.offsets[2]) && myGeo.offsets[2] != 0)
            || (((currentNodeIndices[2] == myGeo.gFineNodesPerDir[2] -1)
                 && (indices[0][2] < rate[2] - myGeo.offsets[2] + 1)) && myGeo.offsets[2] != 0)
            || ((currentNodeIndices[2] == myGeo.gFineNodesPerDir[2] -1)
                && myGeo.lFineNodesPerDir[2] == 1) ) {
          coarseGhosts[0] = 1;
          coarseGhosts[1] = 1;
          coarseGhosts[2] = 1;
          coarseGhosts[3] = 1;
        }
      }
      if( myGeo.ghostInterface[5]
          && (indices[0][2] > myGeo.lFineNodesPerDir[2] - myGeo.offsets[5] - 2) ) {
        coarseGhosts[4] = 1;
        coarseGhosts[5] = 1;
        coarseGhosts[6] = 1;
        coarseGhosts[7] = 1;
      }

      GO firstGhostNodeIndices[3], firstGhostNodeIndex;
      if(currentNodeIndices[0] == myGeo.gFineNodesPerDir[0] - 1) {
        firstGhostNodeIndices[0] = (currentNodeIndices[0]-rate[0]) - (currentNodeIndices[0]-rate[0])
          % myGeo.coarseRate[0];
      } else {
        firstGhostNodeIndices[0] = currentNodeIndices[0] -currentNodeIndices[0]%myGeo.coarseRate[0];
      }
      if(currentNodeIndices[1] == myGeo.gFineNodesPerDir[1] - 1) {
        firstGhostNodeIndices[1] = (currentNodeIndices[1]-rate[1]) - (currentNodeIndices[1]-rate[1])
          % myGeo.coarseRate[1];
      } else {
        firstGhostNodeIndices[1] = currentNodeIndices[1] -currentNodeIndices[1]%myGeo.coarseRate[1];
      }
      if(currentNodeIndices[2] == myGeo.gFineNodesPerDir[2] - 1) {
        firstGhostNodeIndices[2] = (currentNodeIndices[2]-rate[2]) - (currentNodeIndices[2]-rate[2])
          % myGeo.coarseRate[2];
      } else {
        firstGhostNodeIndices[2] = currentNodeIndices[2] -currentNodeIndices[2]%myGeo.coarseRate[2];
      }
      firstGhostNodeIndex = firstGhostNodeIndices[2]*myGeo.gNumFineNodes10
        + firstGhostNodeIndices[1]*myGeo.gFineNodesPerDir[0] + firstGhostNodeIndices[0];

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
              ghostNodeIndex = firstGhostNodeIndex
                + Teuchos::as<GO>(ind2*rate[2]*myGeo.gNumFineNodes10
                                  + ind1*rate[1]*myGeo.gFineNodesPerDir[0] + ind0*rate[0]);
              ghostNodeOnCoarseGridIndex = firstCoarseNodeOnCoarseGridIndex
                + Teuchos::as<GO>(ind2*myGeo.gCoarseNodesPerDir[1]*myGeo.gCoarseNodesPerDir[0]
                                  + ind1*myGeo.gCoarseNodesPerDir[0] + ind0);
              currentGhostLID  = ghostMap->getLocalElement(ghostNodeIndex);
              connecPIDs[ind]  = ghostsPIDs[currentGhostLID];
              connecLGIDs[ind] = ghostNodeIndex;
              for(LO dim = 0; dim < 3; ++dim) {
                connec[ind + 1][dim] = lGhostCoords[dim][currentGhostLID];
              }
              gCoarseNodesOnCoarseGridIndices[4*ind2 + 2*ind1 + ind0] = ghostNodeOnCoarseGridIndex;
            } else {
              // Get the local coarse node index and load its coordinates
              coarseNodeIndex = firstCoarseNodeIndex
                + Teuchos::as<GO>(ind2*rate[2]*myGeo.gNumFineNodes10
                                  + ind1*rate[1]*myGeo.gFineNodesPerDir[0] + ind0*rate[0]);
              coarseNodeOnCoarseGridIndex = firstCoarseNodeOnCoarseGridIndex
                + Teuchos::as<GO>(ind2*myGeo.gCoarseNodesPerDir[1]*myGeo.gCoarseNodesPerDir[0]
                                  + ind1*myGeo.gCoarseNodesPerDir[0] + ind0);
              for(LO dim = 0; dim < 3; ++dim) {
                connec[ind + 1][dim] =
                  lFineCoords[dim][fineCoordsMap->getLocalElement(coarseNodeIndex)];
              }
              gCoarseNodesOnCoarseGridIndices[4*ind2 + 2*ind1 + ind0] = coarseNodeOnCoarseGridIndex;
              connecPIDs[ind]  = -1;
              connecLGIDs[ind] = coarseCoordsMap->getLocalElement(coarseNodeOnCoarseGridIndex);
            }
          }
        }
      }

      // Compute the actual geometric interpolation stencil
      SC stencil[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      ComputeStencil(myGeo.numDimensions, currentNodeIndices, firstCoarseNodeIndices, rate, connec,
                     interpolationOrder, stencil);

      // Finally check whether the fine node is on a coarse: node, edge or face
      // and select accordingly the non-zero values from the stencil and the
      // corresponding column indices
      Array<LO> nzIndStencil(8);
      if( ((currentNodeIndices[0] % myGeo.coarseRate[0] == 0)
           || currentNodeIndices[0] == myGeo.gFineNodesPerDir[0] - 1)
          && ((currentNodeIndices[1] % myGeo.coarseRate[1] == 0)
              || currentNodeIndices[1] == myGeo.gFineNodesPerDir[1] - 1)
          && ((currentNodeIndices[2] % myGeo.coarseRate[2] == 0)
              || currentNodeIndices[2] == myGeo.gFineNodesPerDir[2] - 1) ) {
        if(myGeo.numDimensions==1) {
          xCoarseNodes[currentCoarseNode] = lFineCoords[0][i];
        } else if(myGeo.numDimensions==2) {
          xCoarseNodes[currentCoarseNode] = lFineCoords[0][i];
          yCoarseNodes[currentCoarseNode] = lFineCoords[1][i];
        } else if(myGeo.numDimensions==3) {
          xCoarseNodes[currentCoarseNode] = lFineCoords[0][i];
          yCoarseNodes[currentCoarseNode] = lFineCoords[1][i];
          zCoarseNodes[currentCoarseNode] = lFineCoords[2][i];
        }
        ++currentCoarseNode;

        if(currentNodeIndices[0] == myGeo.gFineNodesPerDir[0] - 1) {
          nzIndStencil[0] += 1;
        }
        if((currentNodeIndices[1] == myGeo.gFineNodesPerDir[1] - 1) && (myGeo.numDimensions > 1)) {
          nzIndStencil[0] += 2;
        }
        if((currentNodeIndices[2] == myGeo.gFineNodesPerDir[2] - 1) && (myGeo.numDimensions > 2)) {
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
          ja [ia[i*dofsPerNode + j] + k] =
            colMapP->getLocalElement(gCoarseNodesOnCoarseGridIndices[nzIndStencil[permutation[k]]]
                                     *dofsPerNode + j);
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
      ComputeConstantInterpolationStencil(numDimensions, currentNodeIndices, coarseNodeIndices,
                                          rate, stencil);
    } else if(interpolationOrder == 1) {
      ComputeLinearInterpolationStencil(numDimensions, coord, stencil);
    }

  } // End ComputeStencil

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ComputeConstantInterpolationStencil(const LO numDimensions, const Array<GO> currentNodeIndices,
                                      const Array<GO> coarseNodeIndices, const LO rate[3],
                                      SC stencil[8]) const {

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
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ComputeLinearInterpolationStencil(const LO numDimensions, const double coord[9][3], SC stencil[8])
    const {

    //                7         8                Find xi, eta and zeta such that
    //                x---------x
    //               /|        /|          Rx = x_p - sum N_i(xi,eta,zeta)x_i = 0
    //             5/ |      6/ |          Ry = y_p - sum N_i(xi,eta,zeta)y_i = 0
    //             x---------x  |          Rz = z_p - sum N_i(xi,eta,zeta)z_i = 0
    //             |  | *P   |  |
    //             |  x------|--x          We can do this with a Newton solver:
    //             | /3      | /4          We will start with initial guess (xi,eta,zeta) = (0,0,0)
    //             |/        |/            We compute the Jacobian and iterate until convergence...
    //  z  y       x---------x
    //  | /        1         2             Once we have (xi,eta,zeta), we can evaluate all N_i which
    //  |/                                 give us the weights for the interpolation stencil!
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
          residual(i) -= functions[0][k]*coord[k+1][i]; //Remove contribution from all coarse points
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
      GetInterpolationFunctions(numDimensions, paramCoords, functions);
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
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetInterpolationFunctions(const LO numDimensions,
                            const Teuchos::SerialDenseVector<LO,double> parameters,
                            double functions[4][8]) const {
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

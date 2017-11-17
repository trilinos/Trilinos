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

    // Coarsen can come in two forms, either a single char that will be interpreted as an integer
    // which is used as the coarsening rate in every spatial dimentions, or it can be a longer
    // string that will then be interpreted as an array of integers.
    // By default coarsen is set as "{2}", hence a coarsening rate of 2 in every spatial dimension
    // is the default setting!
    validParamList->set<std::string >           ("Coarsen",                 "{2}",
                                                 "Coarsening rate in all spatial dimensions");
    validParamList->set<int>                    ("order",                   1,
                                                 "Order of the interpolation scheme used");
    validParamList->set<RCP<const FactoryBase> >("A",                       Teuchos::null,
                                                 "Generating factory of the matrix A");
    validParamList->set<RCP<const FactoryBase> >("Nullspace",               Teuchos::null,
                                                 "Generating factory of the nullspace");
    validParamList->set<RCP<const FactoryBase> >("Coordinates",             Teuchos::null,
                                                 "Generating factory for coorindates");
    validParamList->set<RCP<const FactoryBase> >("gNodesPerDim",            Teuchos::null,
                                                 "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",            Teuchos::null,
                                                 "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<std::string >           ("meshLayout",              "Global Lexicographic",
                                                 "Type of mesh ordering");
    validParamList->set<RCP<const FactoryBase> >("meshData",                Teuchos::null,
                                                 "Mesh ordering associated data");

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
    RCP<GeometricData> myGeometry = rcp(new GeometricData{});

    // Load the mesh layout type and the associated mesh data
    myGeometry->meshLayout = pL.get<std::string>("meshLayout");
    if(fineLevel.GetLevelID() == 0) {
      if(myGeometry->meshLayout == "Local Lexicographic") {
        Array<GO> meshData;
        meshData = fineLevel.Get<Array<GO> >("meshData", NoFactory::get());
        TEUCHOS_TEST_FOR_EXCEPTION(meshData.empty() == true, Exceptions::RuntimeError,
                                   "The meshData array is empty, somehow the input for geometric"
                                   " multigrid are not captured correctly.");
        myGeometry->meshData.resize(rowMap->getComm()->getSize());
        for(int i = 0; i < rowMap->getComm()->getSize(); ++i) {
          myGeometry->meshData[i].resize(10);
          for(int j = 0; j < 10; ++j) {
            myGeometry->meshData[i][j] = meshData[10*i + j];
          }
        }
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords == Teuchos::null, Exceptions::RuntimeError,
                               "Coordinates cannot be accessed from fine level!");
    myGeometry->numDimensions = fineCoords->getNumVectors();

    // Get the number of points in each direction
    if(fineLevel.GetLevelID() == 0) {
      myGeometry->gFineNodesPerDir = fineLevel.Get<Array<GO> >("gNodesPerDim", NoFactory::get());
      myGeometry->lFineNodesPerDir = fineLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
    } else {
      // Loading global number of nodes per diretions
      myGeometry->gFineNodesPerDir = Get<Array<GO> >(fineLevel, "gNodesPerDim");

      // Loading local number of nodes per diretions
      myGeometry->lFineNodesPerDir = Get<Array<LO> >(fineLevel, "lNodesPerDim");
    }
    myGeometry->gNumFineNodes10 = myGeometry->gFineNodesPerDir[1]*myGeometry->gFineNodesPerDir[0];
    myGeometry->gNumFineNodes   = myGeometry->gFineNodesPerDir[2]*myGeometry->gNumFineNodes10;
    myGeometry->lNumFineNodes10 = myGeometry->lFineNodesPerDir[1]*myGeometry->lFineNodesPerDir[0];
    myGeometry->lNumFineNodes   = myGeometry->lFineNodesPerDir[2]*myGeometry->lNumFineNodes10;

    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getLocalLength()
                               != static_cast<size_t>(myGeometry->lNumFineNodes),
                               Exceptions::RuntimeError,
                               "The local number of elements in Coordinates is not equal to the"
                               " number of nodes given by: lNodesPerDim!");
    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getGlobalLength()
                               != static_cast<size_t>(myGeometry->gNumFineNodes),
                               Exceptions::RuntimeError,
                               "The global number of elements in Coordinates is not equal to the"
                               " number of nodes given by: gNodesPerDim!");

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
    TEUCHOS_TEST_FOR_EXCEPTION((crates.size() > 1) && (crates.size() < myGeometry->numDimensions),
                               Exceptions::RuntimeError,
                               "Coarsen must have at least as many components as the number of"
                               " spatial dimensions in the problem.");

    for(LO i = 0; i < 3; ++i) {
      if(i < myGeometry->numDimensions) {
        if(crates.size()==1) {
          myGeometry->coarseRate[i] = crates[0];
        } else if(crates.size() == myGeometry->numDimensions) {
          myGeometry->coarseRate[i] = crates[i];
        }
      } else {
        myGeometry->coarseRate[i] = 1;
      }
    }

    int interpolationOrder = pL.get<int>("order");
    TEUCHOS_TEST_FOR_EXCEPTION((interpolationOrder < 0) || (interpolationOrder > 1),
                               Exceptions::RuntimeError,
                               "The interpolation order can only be set to 0 or 1.");

    // Get the axis permutation from Global axis to Local axis
    Array<LO> mapDirG2L(3), mapDirL2G(3);
    for(LO i = 0; i < myGeometry->numDimensions; ++i) {
      mapDirG2L[i] = i;
    }
    for(LO i = 0; i < myGeometry->numDimensions; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(mapDirG2L[i] > myGeometry->numDimensions,
                                 Exceptions::RuntimeError,
                                 "axis permutation values must all be less than"
                                 " the number of spatial dimensions.");
      mapDirL2G[mapDirG2L[i]] = i;
    }
    RCP<const Map> fineCoordsMap = fineCoords->getMap();

    // This struct stores PIDs, LIDs and GIDs on the fine mesh and GIDs on the coarse mesh.
    RCP<NodesIDs> ghostedCoarseNodes = rcp(new NodesIDs{});
    Array<Array<GO> > lCoarseNodesGIDs(2);
    if((fineLevel.GetLevelID() == 0) && (myGeometry->meshLayout != "Global Lexicographic")) {
      MeshLayoutInterface(interpolationOrder, blkSize, fineCoordsMap, myGeometry,
                          ghostedCoarseNodes, lCoarseNodesGIDs);
    } else {
      // This function expects perfect global lexicographic ordering of nodes and will not process
      // data correctly otherwise. These restrictions allow for the simplest and most efficient
      // processing of the levels (hopefully at least).
      GetCoarsePoints(interpolationOrder, blkSize, fineCoordsMap, myGeometry, ghostedCoarseNodes,
                      lCoarseNodesGIDs);
    }

    // All that is left to do is loop over NCpts and:
    //   - extract coarse points coordiate for coarseCoords
    //   - get coordinates for current stencil computation
    //   - compute current stencil
    //   - compute row and column indices for stencil entries
    RCP<const Map> stridedDomainMapP;
    RCP<Matrix>    P;
    // Fancy formula for the number of non-zero terms
    // All coarse points are injected, other points are using polynomial interpolation
    // and have contribution from (interpolationOrder + 1)^numDimensions
    // Noticebly this leads to 1 when the order is zero, hence fancy MatMatMatMult can be used.
    GO lnnzP = Teuchos::as<LO>(std::pow(interpolationOrder + 1, myGeometry->numDimensions))
      *(myGeometry->lNumFineNodes - myGeometry->lNumCoarseNodes) + myGeometry->lNumCoarseNodes;
    lnnzP = lnnzP*blkSize;
    GO gnnzP = Teuchos::as<LO>(std::pow(interpolationOrder + 1, myGeometry->numDimensions))
      *(myGeometry->gNumFineNodes - myGeometry->gNumCoarseNodes) + myGeometry->gNumCoarseNodes;
    gnnzP = gnnzP*blkSize;

    GetOStream(Runtime1) << "P size = " << blkSize*myGeometry->gNumFineNodes
                         << " x " << blkSize*myGeometry->gNumCoarseNodes << std::endl;
    GetOStream(Runtime1) << "P Fine   grid : " << myGeometry->gFineNodesPerDir[0] << " -- "
                         << myGeometry->gFineNodesPerDir[1] << " -- "
                         << myGeometry->gFineNodesPerDir[2] << std::endl;
    GetOStream(Runtime1) << "P Coarse grid : " << myGeometry->gCoarseNodesPerDir[0] << " -- "
                         << myGeometry->gCoarseNodesPerDir[1] << " -- "
                         << myGeometry->gCoarseNodesPerDir[2] << std::endl;
    GetOStream(Runtime1) << "P nnz estimate: " << gnnzP << std::endl;

    MakeGeneralGeometricP(myGeometry, fineCoords, lnnzP, blkSize, stridedDomainMapP,
                           A, P, coarseCoords, ghostedCoarseNodes, lCoarseNodesGIDs,
                           interpolationOrder);

    // set StridingInformation of P
    if (A->IsView("stridedMaps") == true) {
      P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), stridedDomainMapP);
    } else {
      P->CreateView("stridedMaps", P->getRangeMap(), stridedDomainMapP);
    }

    // store the transfer operator and node coordinates on coarse level
    Set(coarseLevel, "P", P);
    Set(coarseLevel, "coarseCoordinates", coarseCoords);
    Set<Array<GO> >(coarseLevel, "gCoarseNodesPerDim", myGeometry->gCoarseNodesPerDir);
    Set<Array<LO> >(coarseLevel, "lCoarseNodesPerDim", myGeometry->lCoarseNodesPerDir);

    // rst: null space might get scaled here ... do we care. We could just inject at the cpoints,
    // but I don't feel that this is needed.
    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(P->getDomainMap(),
                                                                 fineNullspace->getNumVectors());
    P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(),
             Teuchos::ScalarTraits<SC>::zero());
    Set(coarseLevel, "Nullspace", coarseNullspace);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MeshLayoutInterface(const int interpolationOrder, const LO blkSize, RCP<const Map> fineCoordsMap,
                      RCP<GeometricData> myGeo, RCP<NodesIDs> ghostedCoarseNodes,
                      Array<Array<GO> >& lCoarseNodesGIDs) const{
    // The goal here is to produce maps that globally labels the mesh lexicographically.
    // These maps will replace the current maps of A, the coordinate vector and the nullspace.
    // Ideally if the user provides the necessary maps then nothing needs to be done, otherwise
    // it could be advantageous to allow the user to register a re-labeling function. Ultimately
    // for common ordering schemes, some re-labeling can be directly implemented here.

    int numRanks = fineCoordsMap->getComm()->getSize();
    int myRank   = fineCoordsMap->getComm()->getRank();

    myGeo->myBlock = myGeo->meshData[myRank][2];
    myGeo->startIndices[0] = myGeo->meshData[myRank][3];
    myGeo->startIndices[1] = myGeo->meshData[myRank][5];
    myGeo->startIndices[2] = myGeo->meshData[myRank][7];
    myGeo->startIndices[3] = myGeo->meshData[myRank][4];
    myGeo->startIndices[4] = myGeo->meshData[myRank][6];
    myGeo->startIndices[5] = myGeo->meshData[myRank][8];
    std::sort(myGeo->meshData.begin(), myGeo->meshData.end(),
              [](const std::vector<GO>& a, const std::vector<GO>& b)->bool {
                // The below function sorts ranks by blockID, kmin, jmin and imin
                if(a[2] < b[2]) {
                  return true;
                } else if(a[2] == b[2]) {
                  if(a[7] < b[7]) {
                    return true;
                  } else if(a[7] == b[7]) {
                    if(a[5] < b[5]) {
                      return true;
                    } else if(a[5] == b[5]) {
                      if(a[3] < b[3]) {return true;}
                    }
                  }
                }
                return false;
              });

    myGeo->numBlocks = myGeo->meshData[numRanks - 1][2] + 1;
    // Find the range of the current block
    auto myBlockStart = std::lower_bound(myGeo->meshData.begin(), myGeo->meshData.end(),
                                         myGeo->myBlock - 1,
                                         [](const std::vector<GO>& vec, const GO val)->bool{
                                           return (vec[2] < val) ? true : false;
                                         });
    auto myBlockEnd = std::upper_bound(myGeo->meshData.begin(), myGeo->meshData.end(),
                                       myGeo->myBlock,
                                       [](const GO val, const std::vector<GO>& vec)->bool{
                                         return (val < vec[2]) ? true : false;
                                       });
    // Assuming that i,j,k and ranges are split in pi, pj and pk processors
    // we search for these numbers as they will allow us to find quickly the PID of processors
    // owning ghost nodes.
    auto myKEnd = std::upper_bound(myBlockStart, myBlockEnd, (*myBlockStart)[3],
                                   [](const GO val, const std::vector<GO>& vec)->bool{
                                     return (val < vec[7]) ? true : false;
                                   });
    auto myJEnd = std::upper_bound(myBlockStart, myKEnd, (*myBlockStart)[3],
                                   [](const GO val, const std::vector<GO>& vec)->bool{
                                     return (val < vec[5]) ? true : false;
                                   });
    LO pi = std::distance(myBlockStart, myJEnd);
    LO pj = std::distance(myBlockStart, myKEnd) / pi;
    LO pk = std::distance(myBlockStart, myBlockEnd) / (pj*pi);

    // We also look for the index of the local rank in the current block.
    LO myRankIndex = std::distance(myGeo->meshData.begin(),
                                   std::find_if(myBlockStart, myBlockEnd,
                                                [myRank](const std::vector<GO>& vec)->bool{
                                                  return (vec[0] == myRank) ? true : false;
                                                })
                                   );

    for(int dim = 0; dim < 3; ++dim) {
      if(dim < myGeo->numDimensions) {
        myGeo->offsets[dim]= Teuchos::as<LO>(myGeo->startIndices[dim]) % myGeo->coarseRate[dim];
        myGeo->offsets[dim + 3]= Teuchos::as<LO>(myGeo->startIndices[dim]) % myGeo->coarseRate[dim];
      }
    }

    // Check if the partition contains nodes on a boundary, if so that boundary (face, line or
    // point) will not require ghost nodes.
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < myGeo->numDimensions &&  (myGeo->startIndices[dim] % myGeo->coarseRate[dim] != 0 ||
                                         myGeo->startIndices[dim] == myGeo->gFineNodesPerDir[dim]-1)) {
        myGeo->ghostInterface[2*dim] = true;
      }
      if(dim < myGeo->numDimensions
         && myGeo->startIndices[dim + 3] != myGeo->gFineNodesPerDir[dim] - 1
         && (myGeo->lFineNodesPerDir[dim] == 1 ||
             myGeo->startIndices[dim + 3] % myGeo->coarseRate[dim] != 0)) {
        myGeo->ghostInterface[2*dim+1] = true;
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
    for(int i = 0; i < 3; ++i) {
      if(i < myGeo->numDimensions) {
        // This array is passed to the RAPFactory and eventually becomes gFineNodePerDir on the next
        // level.
        myGeo->gCoarseNodesPerDir[i] = (myGeo->gFineNodesPerDir[i] - 1) / myGeo->coarseRate[i];
        myGeo->endRate[i] = Teuchos::as<LO>((myGeo->gFineNodesPerDir[i] - 1) %myGeo->coarseRate[i]);
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
      // This would happen if the rank does not own any nodes but in that case a subcommunicator
      // should be used so this should really not be a concern.
      if(myGeo->lFineNodesPerDir[i] < 1) {myGeo->lCoarseNodesPerDir[i] = 0;}
    }

    // Assuming linear interpolation, each fine point has contribution from 8 coarse points
    // and each coarse point value gets injected.
    // For systems of PDEs we assume that all dofs have the same P operator.
    myGeo->lNumCoarseNodes = myGeo->lCoarseNodesPerDir[0]*myGeo->lCoarseNodesPerDir[1]
      *myGeo->lCoarseNodesPerDir[2];

    // For each direction, determine how many points (including ghosts) are required.
    for(int dim = 0; dim < 3; ++dim) {
      // The first branch of this if-statement will be used if the rank contains only one layer
      // of nodes in direction i, that layer must also coincide with the boundary of the mesh
      // and coarseRate[i] == endRate[i]...
      if(myGeo->startIndices[dim] == myGeo->gFineNodesPerDir[dim] - 1 &&
         myGeo->startIndices[dim] % myGeo->coarseRate[dim] == 0) {
        myGeo->startGhostedCoarseNode[dim] = myGeo->startIndices[dim] / myGeo->coarseRate[dim] - 1;
      } else {
        myGeo->startGhostedCoarseNode[dim] = myGeo->startIndices[dim] / myGeo->coarseRate[dim];
      }
      myGeo->ghostedCoarseNodesPerDir[dim] = myGeo->lCoarseNodesPerDir[dim];
      // Check whether face *low needs ghost nodes
      if(myGeo->ghostInterface[2*dim]) {myGeo->ghostedCoarseNodesPerDir[dim] += 1;}
      // Check whether face *hi needs ghost nodes
      if(myGeo->ghostInterface[2*dim + 1]) {myGeo->ghostedCoarseNodesPerDir[dim] += 1;}
    }
    myGeo->lNumGhostedNodes = myGeo->ghostedCoarseNodesPerDir[2]*myGeo->ghostedCoarseNodesPerDir[1]
      *myGeo->ghostedCoarseNodesPerDir[0];
    myGeo->lNumGhostNodes = myGeo->lNumGhostedNodes - myGeo->lNumCoarseNodes;
    ghostedCoarseNodes->PIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->LIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->GIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->coarseGIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->colInds.resize(myGeo->lNumGhostedNodes);
    lCoarseNodesGIDs[0].resize(myGeo->lNumCoarseNodes);
    lCoarseNodesGIDs[1].resize(myGeo->lNumCoarseNodes);

    // Now the tricky part starts, the coarse nodes / ghosted coarse nodes need to be imported.
    // This requires finding what their GID on the fine mesh is. They need to be ordered
    // lexicographically to allow for fast sweeps through the mesh.

    // We loop over all ghosted coarse nodes by increasing global lexicographic order
    Array<LO> coarseNodeCoarseIndices(3), coarseNodeFineIndices(3);
    LO currentIndex = -1, countCoarseNodes = 0;
    for(int k = 0; k < myGeo->ghostedCoarseNodesPerDir[2]; ++k) {
      for(int j = 0; j < myGeo->ghostedCoarseNodesPerDir[1]; ++j) {
        for(int i = 0; i < myGeo->ghostedCoarseNodesPerDir[0]; ++i) {
          currentIndex = k*myGeo->ghostedCoarseNodesPerDir[1]*myGeo->ghostedCoarseNodesPerDir[0]
            + j*myGeo->ghostedCoarseNodesPerDir[0] + i;
          coarseNodeCoarseIndices[0] = myGeo->startGhostedCoarseNode[0] + i;
          coarseNodeFineIndices[0] = coarseNodeCoarseIndices[0]*myGeo->coarseRate[0];
          if(coarseNodeFineIndices[0] > myGeo->gFineNodesPerDir[0] - 1) {
            coarseNodeFineIndices[0] = myGeo->gFineNodesPerDir[0] - 1;
          }
          coarseNodeCoarseIndices[1] = myGeo->startGhostedCoarseNode[1] + j;
          coarseNodeFineIndices[1] = coarseNodeCoarseIndices[1]*myGeo->coarseRate[1];
          if(coarseNodeFineIndices[1] > myGeo->gFineNodesPerDir[1] - 1) {
            coarseNodeFineIndices[1] = myGeo->gFineNodesPerDir[1] - 1;
          }
          coarseNodeCoarseIndices[2] = myGeo->startGhostedCoarseNode[2] + k;
          coarseNodeFineIndices[2] = coarseNodeCoarseIndices[2]*myGeo->coarseRate[2];
          if(coarseNodeFineIndices[2] > myGeo->gFineNodesPerDir[2] - 1) {
            coarseNodeFineIndices[2] = myGeo->gFineNodesPerDir[2] - 1;
          }
          GO myGID = -1, myCoarseGID = -1;
          LO myLID = -1, myPID = -1;
          GetGIDLocalLexicographic(i, j, k, coarseNodeFineIndices, myGeo, myRankIndex, pi, pj, pk,
                                   myBlockStart, myBlockEnd, myGID, myPID, myLID);
          myCoarseGID = coarseNodeCoarseIndices[0]
            + coarseNodeCoarseIndices[1]*myGeo->gCoarseNodesPerDir[0]
            + coarseNodeCoarseIndices[2]*myGeo->gCoarseNodesPerDir[1]*myGeo->gCoarseNodesPerDir[0];
          ghostedCoarseNodes->PIDs[currentIndex] = myPID;
          ghostedCoarseNodes->LIDs[currentIndex] = myLID;
          ghostedCoarseNodes->GIDs[currentIndex] = myGID;
          ghostedCoarseNodes->coarseGIDs[currentIndex] = myCoarseGID;
          if(myPID == myRank){
            lCoarseNodesGIDs[0][countCoarseNodes] = myCoarseGID;
            lCoarseNodesGIDs[1][countCoarseNodes] = myGID;
            ++countCoarseNodes;
          }
        }
      }
    }

  } // End MeshLayoutInterface

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetCoarsePoints(const int interpolationOrder, const LO blkSize, RCP<const Map> fineCoordsMap,
                  RCP<GeometricData> myGeo, RCP<NodesIDs> ghostedCoarseNodes,
                  Array<Array<GO> >& lCoarseNodesGIDs) const{
    // Assuming perfect global lexicographic ordering of the mesh, produce two arrays:
    //   1) lGhostNodesIDs that stores PID, LID, GID and coarseGID associated with the coarse nodes
    //      need to compute the local part of the prolongator.
    //   2) lCoarseNodesGIDs that stores the GIDs associated with the local nodes needed to create
    //      the map of the MultiVector of coarse node coordinates.

    {
      GO tmp = 0;
      myGeo->startIndices[2] = fineCoordsMap->getMinGlobalIndex()
        / (myGeo->gFineNodesPerDir[1]*myGeo->gFineNodesPerDir[0]);
      tmp = fineCoordsMap->getMinGlobalIndex()
        % (myGeo->gFineNodesPerDir[1]*myGeo->gFineNodesPerDir[0]);
      myGeo->startIndices[1] = tmp / myGeo->gFineNodesPerDir[0];
      myGeo->startIndices[0] = tmp % myGeo->gFineNodesPerDir[0];
    } // End of scope for tmp
    for(int dim = 0; dim < 3; ++dim) {
      myGeo->startIndices[dim + 3] = myGeo->startIndices[dim] + myGeo->lFineNodesPerDir[dim] - 1;
    }

    for(int dim = 0; dim < 3; ++dim) {
      if(dim < myGeo->numDimensions) {
        myGeo->offsets[dim]= Teuchos::as<LO>(myGeo->startIndices[dim]) % myGeo->coarseRate[dim];
        myGeo->offsets[dim + 3]= Teuchos::as<LO>(myGeo->startIndices[dim]) % myGeo->coarseRate[dim];
      }
    }

    // Check if the partition contains nodes on a boundary, if so that boundary (face, line or
    // point) will not require ghost nodes, unless there is only one node in that direction.
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < myGeo->numDimensions && (myGeo->startIndices[dim] % myGeo->coarseRate[dim] != 0 ||
                                        myGeo->startIndices[dim] == myGeo->gFineNodesPerDir[dim]-1)) {
        myGeo->ghostInterface[2*dim] = true;
      }
      if(dim < myGeo->numDimensions
         && myGeo->startIndices[dim + 3] != myGeo->gFineNodesPerDir[dim] - 1
         && (myGeo->lFineNodesPerDir[dim] == 1 ||
             myGeo->startIndices[dim + 3] % myGeo->coarseRate[dim] != 0)) {
        myGeo->ghostInterface[2*dim + 1] = true;
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
    for(int i = 0; i < 3; ++i) {
      if(i < myGeo->numDimensions) {
        // This array is passed to the RAPFactory and eventually becomes gFineNodePerDir on the next
        // level.
        myGeo->gCoarseNodesPerDir[i] = (myGeo->gFineNodesPerDir[i] - 1) / myGeo->coarseRate[i];
        myGeo->endRate[i] = Teuchos::as<LO>((myGeo->gFineNodesPerDir[i] - 1) %myGeo->coarseRate[i]);
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
      // This would happen if the rank does not own any nodes but in that case a subcommunicator
      // should be used so this should really not be a concern.
      if(myGeo->lFineNodesPerDir[i] < 1) {myGeo->lCoarseNodesPerDir[i] = 0;}
    }

    // Assuming linear interpolation, each fine point has contribution from 8 coarse points
    // and each coarse point value gets injected.
    // For systems of PDEs we assume that all dofs have the same P operator.
    myGeo->lNumCoarseNodes = myGeo->lCoarseNodesPerDir[0]*myGeo->lCoarseNodesPerDir[1]
      *myGeo->lCoarseNodesPerDir[2];

    // For each direction, determine how many points (including ghosts) are required.
    bool ghostedDir[6] = {false};
    for(int dim = 0; dim < 3; ++dim) {
      // The first branch of this if-statement will be used if the rank contains only one layer
      // of nodes in direction i, that layer must also coincide with the boundary of the mesh
      // and coarseRate[i] == endRate[i]...
      if(myGeo->startIndices[dim] == myGeo->gFineNodesPerDir[dim] - 1 &&
         myGeo->startIndices[dim] % myGeo->coarseRate[dim] == 0) {
        myGeo->startGhostedCoarseNode[dim] = myGeo->startIndices[dim] / myGeo->coarseRate[dim] - 1;
      } else {
        myGeo->startGhostedCoarseNode[dim] = myGeo->startIndices[dim] / myGeo->coarseRate[dim];
      }
      myGeo->ghostedCoarseNodesPerDir[dim] = myGeo->lCoarseNodesPerDir[dim];
      // Check whether face *low needs ghost nodes
      if(myGeo->ghostInterface[2*dim]) {
        myGeo->ghostedCoarseNodesPerDir[dim] += 1;
        ghostedDir[2*dim] = true;
      }
      // Check whether face *hi needs ghost nodes
      if(myGeo->ghostInterface[2*dim + 1]) {
        myGeo->ghostedCoarseNodesPerDir[dim] += 1;
        ghostedDir[2*dim + 1] = true;
      }
    }
    myGeo->lNumGhostedNodes = myGeo->ghostedCoarseNodesPerDir[2]*myGeo->ghostedCoarseNodesPerDir[1]
      *myGeo->ghostedCoarseNodesPerDir[0];
    myGeo->lNumGhostNodes = myGeo->lNumGhostedNodes - myGeo->lNumCoarseNodes;
    ghostedCoarseNodes->PIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->LIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->GIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->coarseGIDs.resize(myGeo->lNumGhostedNodes);
    ghostedCoarseNodes->colInds.resize(myGeo->lNumGhostedNodes);
    lCoarseNodesGIDs[0].resize(myGeo->lNumCoarseNodes);
    lCoarseNodesGIDs[1].resize(myGeo->lNumCoarseNodes);

    // Now the tricky part starts, the coarse nodes / ghosted coarse nodes need to be imported.
    // This requires finding what their GID on the fine mesh is. They need to be ordered
    // lexicographically to allow for fast sweeps through the mesh.

    // We loop over all ghosted coarse nodes by increasing global lexicographic order
    Array<LO> coarseNodeCoarseIndices(3), coarseNodeFineIndices(3), ijk(3);
    LO currentIndex = -1, countCoarseNodes = 0;
    for(ijk[2] = 0; ijk[2] < myGeo->ghostedCoarseNodesPerDir[2]; ++ijk[2]) {
      for(ijk[1] = 0; ijk[1] < myGeo->ghostedCoarseNodesPerDir[1]; ++ijk[1]) {
        for(ijk[0] = 0; ijk[0] < myGeo->ghostedCoarseNodesPerDir[0]; ++ijk[0]) {
          currentIndex = ijk[2]*myGeo->ghostedCoarseNodesPerDir[1]*myGeo->ghostedCoarseNodesPerDir[0]
            + ijk[1]*myGeo->ghostedCoarseNodesPerDir[0] + ijk[0];
          coarseNodeCoarseIndices[0] = myGeo->startGhostedCoarseNode[0] + ijk[0];
          coarseNodeFineIndices[0] = coarseNodeCoarseIndices[0]*myGeo->coarseRate[0];
          if(coarseNodeFineIndices[0] > myGeo->gFineNodesPerDir[0] - 1) {
            coarseNodeFineIndices[0] = myGeo->gFineNodesPerDir[0] - 1;
          }
          coarseNodeCoarseIndices[1] = myGeo->startGhostedCoarseNode[1] + ijk[1];
          coarseNodeFineIndices[1] = coarseNodeCoarseIndices[1]*myGeo->coarseRate[1];
          if(coarseNodeFineIndices[1] > myGeo->gFineNodesPerDir[1] - 1) {
            coarseNodeFineIndices[1] = myGeo->gFineNodesPerDir[1] - 1;
          }
          coarseNodeCoarseIndices[2] = myGeo->startGhostedCoarseNode[2] + ijk[2];
          coarseNodeFineIndices[2] = coarseNodeCoarseIndices[2]*myGeo->coarseRate[2];
          if(coarseNodeFineIndices[2] > myGeo->gFineNodesPerDir[2] - 1) {
            coarseNodeFineIndices[2] = myGeo->gFineNodesPerDir[2] - 1;
          }
          GO myGID = 0, myCoarseGID = -1;
          GO factor[3] = {};
          factor[2] = myGeo->gNumFineNodes10;
          factor[1] = myGeo->gFineNodesPerDir[0];
          factor[0] = 1;
          for(int dim = 0; dim < 3; ++dim) {
            if(dim < myGeo->numDimensions) {
              if(myGeo->startIndices[dim] - myGeo->offsets[dim] + ijk[dim]*myGeo->coarseRate[dim]
                 < myGeo->gFineNodesPerDir[dim] - 1) {
                myGID += (myGeo->startIndices[dim] - myGeo->offsets[dim]
                          + ijk[dim]*myGeo->coarseRate[dim])*factor[dim];
              } else {
                myGID += (myGeo->startIndices[dim] - myGeo->offsets[dim]
                          + (ijk[dim] - 1)*myGeo->coarseRate[dim] + myGeo->endRate[dim])*factor[dim];
              }
            }
          }
          myCoarseGID = coarseNodeCoarseIndices[0]
            + coarseNodeCoarseIndices[1]*myGeo->gCoarseNodesPerDir[0]
            + coarseNodeCoarseIndices[2]*myGeo->gCoarseNodesPerDir[1]*myGeo->gCoarseNodesPerDir[0];
          ghostedCoarseNodes->GIDs[currentIndex] = myGID;
          ghostedCoarseNodes->coarseGIDs[currentIndex] = myCoarseGID;
          if((!ghostedDir[0] || ijk[0] != 0)
             && (!ghostedDir[2] || ijk[1] != 0)
             && (!ghostedDir[4] || ijk[2] != 0)
             && (!ghostedDir[1] || ijk[0] != myGeo->ghostedCoarseNodesPerDir[0] - 1)
             && (!ghostedDir[3] || ijk[1] != myGeo->ghostedCoarseNodesPerDir[1] - 1)
             && (!ghostedDir[5] || ijk[2] != myGeo->ghostedCoarseNodesPerDir[2] - 1)){
            lCoarseNodesGIDs[0][countCoarseNodes] = myCoarseGID;
            lCoarseNodesGIDs[1][countCoarseNodes] = myGID;
            ++countCoarseNodes;
          }
        }
      }
    }
    Array<int> ghostsPIDs(myGeo->lNumGhostedNodes);
    Array<LO>  ghostsLIDs(myGeo->lNumGhostedNodes);
    fineCoordsMap->getRemoteIndexList(ghostedCoarseNodes->GIDs(),
                                      ghostedCoarseNodes->PIDs(),
                                      ghostedCoarseNodes->LIDs());
  } // End GetCoarsePoint

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  MakeGeneralGeometricP(RCP<GeometricData> myGeo,
                        const RCP<Xpetra::MultiVector<double,LO,GO,Node> >& fineCoords,
                        const LO nnzP, const LO dofsPerNode,
                        RCP<const Map>& stridedDomainMapP,RCP<Matrix> & Amat, RCP<Matrix>& P,
                        RCP<Xpetra::MultiVector<double,LO,GO,Node> >& coarseCoords,
                        RCP<NodesIDs> ghostedCoarseNodes, Array<Array<GO> > coarseNodesGIDs,
                        int interpolationOrder) const {

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

    LO myRank          = Amat->getRowMap()->getComm()->getRank();
    GO numGloCols      = dofsPerNode*myGeo->gNumCoarseNodes;

    // Build maps necessary to create and fill complete the prolongator
    // note: rowMapP == rangeMapP and colMapP != domainMapP
    RCP<const Map> rowMapP = Amat->getDomainMap();

    RCP<const Map> domainMapP;

    RCP<const Map> colMapP;

    // At this point we need to create the column map which is a delicate operation.
    // The entries in that map need to be ordered as follows:
    //         1) first owned entries ordered by LID
    //         2) second order the remaining entries by PID
    //         3) entries with the same remote PID are ordered by GID.
    // One piece of good news: myGeo->lNumCoarseNodes is the number of ownedNodes and
    // myGeo->lNumGhostNodes is the number of remote nodes. The sorting can be limited to remote
    // nodes as the owned ones are alreaded ordered by LID!

    {
      std::vector<NodeID> colMapOrdering(myGeo->lNumGhostedNodes);
      for(LO ind = 0; ind < myGeo->lNumGhostedNodes; ++ind) {
        colMapOrdering[ind].GID = ghostedCoarseNodes->GIDs[ind];
        if(ghostedCoarseNodes->PIDs[ind] == myRank) {
          colMapOrdering[ind].PID = -1;
        } else {
          colMapOrdering[ind].PID = ghostedCoarseNodes->PIDs[ind];
        }
        colMapOrdering[ind].LID = ghostedCoarseNodes->LIDs[ind];
        colMapOrdering[ind].lexiInd = ind;
      }
      std::sort(colMapOrdering.begin(), colMapOrdering.end(),
                [](NodeID a, NodeID b)->bool{
                  return ( (a.PID < b.PID) || ((a.PID == b.PID) && (a.GID < b.GID)) );
                });

      Array<GO> colGIDs(dofsPerNode*myGeo->lNumGhostedNodes);
      for(LO ind = 0; ind < myGeo->lNumGhostedNodes; ++ind) {
        // Save the permutation calculated to go from Lexicographic indexing to column map indexing
        ghostedCoarseNodes->colInds[colMapOrdering[ind].lexiInd] = ind;
        for(LO dof = 0; dof < dofsPerNode; ++dof) {
          colGIDs[dofsPerNode*ind + dof] = dofsPerNode*colMapOrdering[ind].GID + dof;
        }
      }
      domainMapP = Xpetra::MapFactory<LO,GO,NO>::Build(rowMapP->lib(),
                                                       numGloCols,
                                                       colGIDs.view(0,dofsPerNode*
                                                                    myGeo->lNumCoarseNodes),
                                                       rowMapP->getIndexBase(),
                                                       rowMapP->getComm());
      colMapP = Xpetra::MapFactory<LO,GO,NO>::Build(rowMapP->lib(),
                                                    OTI,
                                                    colGIDs.view(0, colGIDs.size()),
                                                    rowMapP->getIndexBase(),
                                                    rowMapP->getComm());
    } // End of scope for colMapOrdering and colGIDs

    std::vector<size_t> strideInfo(1);
    strideInfo[0] = dofsPerNode;
    stridedDomainMapP = Xpetra::StridedMapFactory<LO,GO,NO>::Build(domainMapP, strideInfo);

    // Build the map for the coarse level coordinates, create the associated MultiVector and fill it
    // with an import from the fine coordinates MultiVector. As data is local this should not create
    // communications during the importer creation.
    RCP<const Map> coarseCoordsMap = MapFactory::Build (fineCoords->getMap()->lib(),
                                                        myGeo->gNumCoarseNodes,
                                                        coarseNodesGIDs[0](),
                                                        fineCoords->getMap()->getIndexBase(),
                                                        fineCoords->getMap()->getComm());
    RCP<const Map> coarseCoordsFineMap = MapFactory::Build (fineCoords->getMap()->lib(),
                                                            myGeo->gNumCoarseNodes,
                                                            coarseNodesGIDs[1](),
                                                            fineCoords->getMap()->getIndexBase(),
                                                            fineCoords->getMap()->getComm());

    RCP<const Import> coarseImporter = ImportFactory::Build(fineCoords->getMap(),
                                                            coarseCoordsFineMap);
    coarseCoords = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(coarseCoordsFineMap,
                                                                      myGeo->numDimensions);
    coarseCoords->doImport(*fineCoords, *coarseImporter, Xpetra::INSERT);
    coarseCoords->replaceMap(coarseCoordsMap);

    // Do the actual import using the fineCoords->getMap()
    RCP<const Map> ghostMap = Xpetra::MapFactory<LO,GO,NO>::Build(fineCoords->getMap()->lib(),
                                                               OTI,
                                                               ghostedCoarseNodes->GIDs(),
                                                               fineCoords->getMap()->getIndexBase(),
                                                               fineCoords->getMap()->getComm());
    RCP<const Import> ghostImporter = ImportFactory::Build(fineCoords->getMap(), ghostMap);
    RCP<xdMV> ghostCoords = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(ghostMap,
                                                                              myGeo->numDimensions);
    ghostCoords->doImport(*fineCoords, *ghostImporter, Xpetra::INSERT);

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

    Array<ArrayRCP<double> > ghostedCoords(3);
    {
      ArrayRCP<double> tmp(ghostCoords->getLocalLength(),  0.0);
      for(int dim = 0; dim < 3; ++dim) {
        if(dim < myGeo->numDimensions) {
          ghostedCoords[dim] = ghostCoords->getDataNonConst(dim);
        } else {
          ghostedCoords[dim] = tmp;
        }
      }
    }

    // Declaration and assignment of fineCoords which holds the coordinates of the fine nodes in 3D.
    // To do so we pull the nD coordinates from fineCoords and pad the rest with zero vectors...
    RCP<Xpetra::Vector<double,LO,GO,NO> > zeros
      = Xpetra::VectorFactory<double,LO,GO,NO>::Build(fineCoords->getMap(), true);
    ArrayRCP< ArrayRCP<double> > lFineCoords(3);
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < myGeo->numDimensions) {
        lFineCoords[dim] = fineCoords->getDataNonConst(dim);
      } else {
        lFineCoords[dim] = zeros->getDataNonConst(0);
      }
    }

    GO tStencil = 0;
    for(int currentIndex = 0; currentIndex < myGeo->lNumFineNodes; ++currentIndex) {
      Array<GO> ghostedIndices(3), firstInterpolationIndices(3);
      Array<LO> interpolationPIDs(8), interpolationLIDs(8), interpolationGIDs(8);
      Array<Array<double> > interpolationCoords(9);
      interpolationCoords[0].resize(3);
      GO firstInterpolationNodeIndex;
      int nStencil = 0;
      for(int dim = 0; dim < 3; ++dim) {
        interpolationCoords[0][dim] = lFineCoords[dim][currentIndex];
      }

      // Compute the ghosted (i,j,k) of the current node, that assumes (I,J,K)=(0,0,0) to be
      // associated with the first node in ghostCoords
      {// Scope for tmp
        ghostedIndices[2] = currentIndex / (myGeo->lFineNodesPerDir[1]*myGeo->lFineNodesPerDir[0]);
        LO tmp = currentIndex % (myGeo->lFineNodesPerDir[1]*myGeo->lFineNodesPerDir[0]);
        ghostedIndices[1] = tmp / myGeo->lFineNodesPerDir[0];
        ghostedIndices[0] = tmp % myGeo->lFineNodesPerDir[0];
        for(int dim = 0; dim < 3; ++dim) {
          ghostedIndices[dim] += myGeo->offsets[dim];
        }
        // A special case appears when the mesh is really coarse: it is possible for a rank to
        // have a single coarse node in a given direction. If this happens on the highest i, j or k
        // we need to "grab" a coarse node with a lower i, j, or k which leads us to add to the
        // value of ghostedIndices
      }
      // No we find the indices of the coarse nodes used for interpolation simply by integer
      // division.
      for(int dim = 0; dim < 3; ++dim) {
        firstInterpolationIndices[dim] = ghostedIndices[dim] / myGeo->coarseRate[dim];
        // If you are on the edge of the local domain go back one coarse node, unless there is only
        // one node on the local domain...
        if(firstInterpolationIndices[dim] == myGeo->ghostedCoarseNodesPerDir[dim] - 1
           && myGeo->ghostedCoarseNodesPerDir[dim] > 1) {
          firstInterpolationIndices[dim] -= 1;
        }
      }
      firstInterpolationNodeIndex =
        firstInterpolationIndices[2]*myGeo->ghostedCoarseNodesPerDir[1]*myGeo->ghostedCoarseNodesPerDir[0]
        + firstInterpolationIndices[1]*myGeo->ghostedCoarseNodesPerDir[0]
        + firstInterpolationIndices[0];

      // We extract the coordinates and PIDs associated with each coarse node used during
      // inteprolation in order to fill the prolongator correctly
      {
        LO ind = -1;
        for(int k = 0; k < 2; ++k) {
          for(int j = 0; j < 2; ++j) {
            for(int i = 0; i < 2; ++i) {
              ind = k*4 + j*2 + i;
              interpolationLIDs[ind] = firstInterpolationNodeIndex
                + k*myGeo->ghostedCoarseNodesPerDir[1]*myGeo->ghostedCoarseNodesPerDir[0]
                + j*myGeo->ghostedCoarseNodesPerDir[0] + i;
              if(ghostedCoarseNodes->PIDs[interpolationLIDs[ind]] == rowMapP->getComm()->getRank()){
                interpolationPIDs[ind] = -1;
              } else {
                interpolationPIDs[ind] = ghostedCoarseNodes->PIDs[interpolationLIDs[ind]];
              }
              interpolationGIDs[ind] = ghostedCoarseNodes->coarseGIDs[interpolationLIDs[ind]];

              interpolationCoords[ind + 1].resize(3);
              for(int dim = 0; dim < 3; ++dim) {
                interpolationCoords[ind + 1][dim]
                  = ghostedCoords[dim][interpolationLIDs[ind]];
              }
            }
          }
        }
      } // End of ind scope

      // Compute the actual geometric interpolation stencil
      // LO stencilLength = static_cast<LO>(std::pow(interpolationOrder + 1, 3));
      std::vector<double> stencil(8);
      Array<GO> firstCoarseNodeFineIndices(3);
      int rate[3] = {};
      for(int dim = 0; dim < 3; ++dim) {
        firstCoarseNodeFineIndices[dim] = firstInterpolationIndices[dim]*myGeo->coarseRate[dim];
        if((myGeo->startIndices[dim + 3] == myGeo->gFineNodesPerDir[dim] - 1)
           && (ghostedIndices[dim] >=
               (myGeo->ghostedCoarseNodesPerDir[dim] - 2)*myGeo->coarseRate[dim])) {
          rate[dim] = myGeo->endRate[dim];
        } else {
          rate[dim] = myGeo->coarseRate[dim];
        }
      }
      Array<GO> trueGhostedIndices(3);
      // This handles the case of a rank having a single node that also happens to be the last node
      // in any direction. It might be more efficient to re-write the algo so that this is
      // incorporated in the definition of ghostedIndices directly...
      for(int dim = 0; dim < 3; ++dim) {
        if (myGeo->startIndices[dim] == myGeo->gFineNodesPerDir[dim] - 1) {
          trueGhostedIndices[dim] = ghostedIndices[dim] + rate[dim];
        } else {
          trueGhostedIndices[dim] = ghostedIndices[dim];
        }
      }
      ComputeStencil(myGeo->numDimensions, trueGhostedIndices, firstCoarseNodeFineIndices, rate,
                     interpolationCoords, interpolationOrder, stencil);

      // Finally check whether the fine node is on a coarse: node, edge or face
      // and select accordingly the non-zero values from the stencil and the
      // corresponding column indices
      Array<LO> nzIndStencil(8), permutation(8);
      for(LO k = 0; k < 8; ++k) {permutation[k] = k;}
      if(interpolationOrder == 0) {
        nStencil = 1;
        for(LO k = 0; k < 8; ++k) {
          nzIndStencil[k] = static_cast<LO>(stencil[0]);
        }
        stencil[0] = 0.0;
        stencil[nzIndStencil[0]] = 1.0;
      } else if(interpolationOrder == 1) {
        Array<GO> currentNodeGlobalFineIndices(3);
        for(int dim = 0; dim < 3; ++dim) {
          currentNodeGlobalFineIndices[dim] = ghostedIndices[dim] - myGeo->offsets[dim]
            + myGeo->startIndices[dim];
      }
        if( ((ghostedIndices[0] % myGeo->coarseRate[0] == 0)
             || currentNodeGlobalFineIndices[0] == myGeo->gFineNodesPerDir[0] - 1)
            && ((ghostedIndices[1] % myGeo->coarseRate[1] == 0)
                || currentNodeGlobalFineIndices[1] == myGeo->gFineNodesPerDir[1] - 1)
            && ((ghostedIndices[2] % myGeo->coarseRate[2] == 0)
                || currentNodeGlobalFineIndices[2] == myGeo->gFineNodesPerDir[2] - 1) ) {
          if((currentNodeGlobalFineIndices[0] == myGeo->gFineNodesPerDir[0] - 1) ||
             (ghostedIndices[0] / myGeo->coarseRate[0] == myGeo->ghostedCoarseNodesPerDir[0] - 1)) {
            nzIndStencil[0] += 1;
          }
          if(((currentNodeGlobalFineIndices[1] == myGeo->gFineNodesPerDir[1] - 1) ||
              (ghostedIndices[1] / myGeo->coarseRate[1] == myGeo->ghostedCoarseNodesPerDir[1] - 1))
             && (myGeo->numDimensions > 1)){
            nzIndStencil[0] += 2;
          }
          if(((currentNodeGlobalFineIndices[2] == myGeo->gFineNodesPerDir[2] - 1) ||
              (ghostedIndices[2] / myGeo->coarseRate[2] == myGeo->ghostedCoarseNodesPerDir[2] - 1))
             && myGeo->numDimensions > 2) {
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
      }

      // Here the values are filled in the Crs matrix arrays
      // This is basically the only place these variables are modified
      // Hopefully this makes handling system of PDEs easy!

      // Loop on dofsPerNode and process each row for the current Node


      // Sort nodes by PIDs using stable sort to keep sublist ordered by LIDs and GIDs
      sh_sort2(interpolationPIDs.begin(),interpolationPIDs.end(),
               permutation.begin(), permutation.end());

      GO colInd;
      for(LO j = 0; j < dofsPerNode; ++j) {
        ia[currentIndex*dofsPerNode + j + 1] = ia[currentIndex*dofsPerNode + j] + nStencil;
        for(LO k = 0; k < nStencil; ++k) {
          colInd = ghostedCoarseNodes->colInds[interpolationLIDs[nzIndStencil[permutation[k]]]];
          ja[ia[currentIndex*dofsPerNode + j] + k] = colInd*dofsPerNode + j;
          val[ia[currentIndex*dofsPerNode + j] + k] = stencil[nzIndStencil[permutation[k]]];
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
  } // End MakeGeneralGeometricP

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeStencil(
                                const LO numDimensions, const Array<GO> currentNodeIndices,
                                const Array<GO> coarseNodeIndices, const LO rate[3],
                                const Array<Array<double> > coord, const int interpolationOrder,
                                std::vector<double>& stencil) const {

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
                                      std::vector<double>& stencil) const {

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
    stencil[0] = coarseNode;

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  ComputeLinearInterpolationStencil(const LO numDimensions, const Array<Array<double> > coord,
                                    std::vector<double>& stencil)
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
  } // End sh_sort_permute

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
  } // End sh_sort2

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeneralGeometricPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetGIDLocalLexicographic(const GO i, const GO j, const GO k,
                           const Array<LO> coarseNodeFineIndices, const RCP<GeometricData> myGeo,
                           const LO myRankIndex, const LO pi, const LO pj, const LO pk,
                           const typename std::vector<std::vector<GO> >::iterator blockStart,
                           const typename std::vector<std::vector<GO> >::iterator blockEnd,
                           GO& myGID, LO& myPID, LO& myLID) const {

    LO ni = -1, nj = -1, li = -1, lj = -1, lk = -1;
    LO myRankGuess = myRankIndex;
    // We try to make a logical guess as to which PID owns the current coarse node
    if(i == 0 && myGeo->ghostInterface[0]) {
      --myRankGuess;
    } else if((i == myGeo->ghostedCoarseNodesPerDir[0] - 1) && myGeo->ghostInterface[1]) {
      ++myRankGuess;
    }
    if(j == 0 && myGeo->ghostInterface[2]) {
      myRankGuess -= pi;
    } else if((j == myGeo->ghostedCoarseNodesPerDir[1] - 1) && myGeo->ghostInterface[3]) {
      myRankGuess += pi;
    }
    if(k == 0 && myGeo->ghostInterface[4]) {
      myRankGuess -= pj*pi;
    } else if((k == myGeo->ghostedCoarseNodesPerDir[2] - 1) && myGeo->ghostInterface[5]) {
      myRankGuess += pj*pi;
    }
    if(coarseNodeFineIndices[0] >= myGeo->meshData[myRankGuess][3]
       && coarseNodeFineIndices[0] <= myGeo->meshData[myRankGuess][4]
       && coarseNodeFineIndices[1] >= myGeo->meshData[myRankGuess][5]
       && coarseNodeFineIndices[1] <= myGeo->meshData[myRankGuess][6]
       && coarseNodeFineIndices[2] >= myGeo->meshData[myRankGuess][7]
       && coarseNodeFineIndices[2] <= myGeo->meshData[myRankGuess][8]) {
      myPID = myGeo->meshData[myRankGuess][0];
      ni = myGeo->meshData[myRankGuess][4] - myGeo->meshData[myRankGuess][3] + 1;
      nj = myGeo->meshData[myRankGuess][6] - myGeo->meshData[myRankGuess][5] + 1;
      li = coarseNodeFineIndices[0] - myGeo->meshData[myRankGuess][3];
      lj = coarseNodeFineIndices[1] - myGeo->meshData[myRankGuess][5];
      lk = coarseNodeFineIndices[2] - myGeo->meshData[myRankGuess][7];
      myLID = lk*nj*ni + lj*ni + li;
      myGID = myGeo->meshData[myRankGuess][9] + myLID;
    } else { // The guess failed, let us use the heavy artilery: std::find_if()
      // It could be interesting to monitor how many times this branch of the code gets
      // used as it is far more expensive than the above one...
      auto nodeRank = std::find_if(blockStart, blockEnd,
                                   [coarseNodeFineIndices](const std::vector<GO>& vec){
                                     if(coarseNodeFineIndices[0] >= vec[3]
                                        && coarseNodeFineIndices[0] <= vec[4]
                                        && coarseNodeFineIndices[1] >= vec[5]
                                        && coarseNodeFineIndices[1] <= vec[6]
                                        && coarseNodeFineIndices[2] >= vec[7]
                                        && coarseNodeFineIndices[2] <= vec[8]) {
                                       return true;
                                     } else {
                                       return false;
                                     }
                                   });
      myPID = (*nodeRank)[0];
      ni = (*nodeRank)[4] - (*nodeRank)[3] + 1;
      nj = (*nodeRank)[6] - (*nodeRank)[5] + 1;
      li = coarseNodeFineIndices[0] - (*nodeRank)[3];
      lj = coarseNodeFineIndices[1] - (*nodeRank)[5];
      lk = coarseNodeFineIndices[2] - (*nodeRank)[7];
      myLID = lk*nj*ni + lj*ni + li;
      myGID = (*nodeRank)[9] + myLID;
    }
  } // End GetGIDLocalLexicographic

} //namespace MueLu

#define MUELU_GENERALGEOMETRICPFACTORY_SHORT
#endif // MUELU_GENERALGEOMETRICPFACTORY_DEF_HPP

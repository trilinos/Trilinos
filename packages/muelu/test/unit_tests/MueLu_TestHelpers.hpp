// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TEST_HELPERS_H
#define MUELU_TEST_HELPERS_H
#include <stdio.h>  //DEBUG
#include <string>
#include <set>
#ifndef _MSC_VER
#include <dirent.h>
#endif

// Teuchos
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifdef HAVE_MUELU_EPETRA
#include "Epetra_config.h"
#endif

// Xpetra
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Parameters.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_CrsGraph.hpp"

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Level.hpp"

// Galeri
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"

namespace Galeri {
namespace Xpetra {
template <class LocalOrdinal, class GlobalOrdinal, class Map>
RCP<Map> CreateMap(const std::string& mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList& list);

#ifdef HAVE_GALERI_XPETRA
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(::Xpetra::UnderlyingLib lib, const std::string& mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList& list);
#endif
}  // namespace Xpetra
}  // namespace Galeri

#include "MueLu_NoFactory.hpp"

// Conditional Tpetra stuff
#include "TpetraCore_config.h"
#include "Xpetra_TpetraCrsGraph.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"
#include "Xpetra_TpetraBlockCrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"

#include <MueLu_TestHelpers_Common.hpp>

namespace MueLuTests {
using Teuchos::arcp;
using Teuchos::arcp_reinterpret_cast;
using Teuchos::arcpFromArrayView;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_implicit_cast;
using Teuchos::rcpFromRef;

namespace TestHelpers {

using Xpetra::global_size_t;

class Parameters {
 private:
  Parameters() {}  // static class

 public:
  static Xpetra::Parameters xpetraParameters;

  inline static RCP<const Teuchos::Comm<int> > getDefaultComm() {
    return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }

  inline static Xpetra::UnderlyingLib getLib() {
    return TestHelpers::Parameters::xpetraParameters.GetLib();
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TestFactory {
#include "MueLu_UseShortNames.hpp"

 private:
  TestFactory() {}  // static class

 public:
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
  typedef Xpetra::MultiVectorFactory<real_type, LO, GO, NO> RealValuedMultiVectorFactory;

  //
  // Method that creates a map containing a specified number of local elements per process.
  //
  static const RCP<const Map> BuildMap(LO numElementsPerProc) {
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();

    return MapFactory::Build(TestHelpers::Parameters::getLib(), INVALID, numElementsPerProc, 0, comm);

  }  // BuildMap()

  // Create a matrix as specified by parameter list options
  static RCP<Matrix> BuildMatrix(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) {
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    if (lib == Xpetra::NotSpecified)
      lib = TestHelpers::Parameters::getLib();

    GO nx, ny, nz;
    nx = ny = nz = 5;
    nx           = matrixList.get("nx", nx);
    ny           = matrixList.get("ny", ny);
    nz           = matrixList.get("nz", nz);

    std::string matrixType = matrixList.get("matrixType", "Laplace1D");
    GO numGlobalElements;  // global_size_t
    if (matrixType == "Laplace1D")
      numGlobalElements = nx;
    else if (matrixType == "Laplace2D" || matrixType == "Star2D" || matrixType == "Cross2D")
      numGlobalElements = nx * ny;
    else if (matrixType == "Elasticity2D")
      numGlobalElements = 2 * nx * ny;
    else if (matrixType == "Laplace3D" || matrixType == "Brick3D")
      numGlobalElements = nx * ny * nz;
    else if (matrixType == "Elasticity3D")
      numGlobalElements = 3 * nx * ny * nz;
    else {
      std::string msg = matrixType + " is unsupported (in unit testing)";
      throw(MueLu::Exceptions::RuntimeError(msg));
    }

    RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, 0, comm);
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixType, map, matrixList);
    RCP<Matrix> Op = Pr->BuildMatrix();

    return Op;
  }  // BuildMatrix()

  // Create a tridiagonal matrix (stencil = [b,a,c]) with the specified number of rows
  // dofMap: row map of matrix
  static RCP<Matrix> BuildTridiag(RCP<const Map> dofMap, Scalar a, Scalar b, Scalar c, Xpetra::UnderlyingLib lib = Xpetra::NotSpecified) {  // global_size_t

    if (lib == Xpetra::NotSpecified)
      lib = TestHelpers::Parameters::getLib();

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Teuchos::RCP<Matrix> mtx = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(dofMap, 3);

    LocalOrdinal NumMyElements                               = dofMap->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = dofMap->getLocalElementList();
    GlobalOrdinal indexBase                                  = dofMap->getIndexBase();

    GlobalOrdinal NumEntries;
    LocalOrdinal nnz = 3;
    std::vector<Scalar> Values(nnz);
    std::vector<GlobalOrdinal> Indices(nnz);

    // access information from strided block map
    std::vector<size_t> strInfo = std::vector<size_t>();
    size_t blockSize            = 1;
    LocalOrdinal blockId        = -1;
    GlobalOrdinal offset        = 0;

    Teuchos::RCP<const StridedMap> strdofMap = Teuchos::rcp_dynamic_cast<const StridedMap>(dofMap);
    if (strdofMap != Teuchos::null) {
      strInfo   = strdofMap->getStridingData();
      blockSize = strdofMap->getFixedBlockSize();
      blockId   = strdofMap->getStridedBlockId();
      offset    = strdofMap->getOffset();
      TEUCHOS_TEST_FOR_EXCEPTION(blockId > -1 && strInfo[blockId] == 1, MueLu::Exceptions::RuntimeError,
                                 "MueLu::TestHelpers::BuildTridiag: strInfo block size must be > 1.");
      // todo: write one more special case for a row being first and last bock row
    } else {
      // no striding information. emulate block matrix
      blockId = 0;
      strInfo.push_back(blockSize);  // default block size = 1
    }

    GlobalOrdinal rrmax = (dofMap->getMaxAllGlobalIndex() - offset - indexBase) / blockSize;

    // loop over all rows
    for (LocalOrdinal i = 0; i < NumMyElements; i++) {
      GlobalOrdinal rr = (MyGlobalElements[i] - offset - indexBase) / blockSize + indexBase;  // node index

      // distinguish 5 different cases

      GlobalOrdinal blockOffset = 0;
      for (LocalOrdinal k = 0; k < blockId; k++)
        blockOffset += Teuchos::as<GlobalOrdinal>(strInfo[k]);

      if (MyGlobalElements[i] == blockOffset + offset + indexBase) {
        // very first row
        Indices[0] = MyGlobalElements[i];
        Values[0]  = a;
        Indices[1] = MyGlobalElements[i] + 1;
        Values[1]  = c;
        NumEntries = 2;
      } else if (MyGlobalElements[i] == rrmax * Teuchos::as<GO>(blockSize) + blockOffset + Teuchos::as<GO>(strInfo[blockId]) - 1 + offset + indexBase) {
        // very last row
        Indices[0] = MyGlobalElements[i] - 1;
        Values[0]  = b;
        Indices[1] = MyGlobalElements[i];
        Values[1]  = a;
        NumEntries = 2;
      } else if (MyGlobalElements[i] == rr * Teuchos::as<GO>(blockSize) + blockOffset + Teuchos::as<GO>(strInfo[blockId]) - 1 + offset + indexBase) {
        // last row in current node block
        Indices[0] = MyGlobalElements[i] - 1;
        Values[0]  = b;
        Indices[1] = MyGlobalElements[i];
        Values[1]  = a;
        Indices[2] = (rr + 1) * blockSize + blockOffset + offset + indexBase;
        Values[2]  = c;
        NumEntries = 3;
      } else if (MyGlobalElements[i] == rr * Teuchos::as<GO>(blockSize) + blockOffset + offset + indexBase) {
        // first row in current node block
        Indices[0] = (rr - 1) * blockSize + blockOffset + strInfo[blockId] - 1 + offset + indexBase;
        Values[0]  = b;
        Indices[1] = MyGlobalElements[i];
        Values[1]  = a;
        Indices[2] = MyGlobalElements[i] + 1;
        Values[2]  = c;
        NumEntries = 3;
      } else {
        // usual row entries in block rows
        Indices[0] = MyGlobalElements[i] - 1;
        Values[0]  = b;
        Indices[1] = MyGlobalElements[i];
        Values[1]  = a;
        Indices[2] = MyGlobalElements[i] + 1;
        Values[2]  = c;
        NumEntries = 3;
      }

      // debug output
      /*std::cout << "--------> Proc " << comm->getRank() << " lrow " << i << " grow " << MyGlobalElements[i] << " node " << rr << " values = " ;
      for (size_t k = 0; k < NumEntries; k++) {
        std::cout << " " << Indices[k] << "->" << Values[k] << "   ";
      }
      std::cout << std::endl;*/

      // put the off-diagonal entries
      // Xpetra wants ArrayViews (sigh)
      Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
      Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
      mtx->insertGlobalValues(MyGlobalElements[i], iv, av);
    }

    mtx->fillComplete();

    return mtx;
  }  // BuildTridiag()

  // Create a 1D Poisson matrix with the specified number of rows
  // nx: global number of rows
  static RCP<Matrix> Build1DPoisson(GO nx, Xpetra::UnderlyingLib lib = Xpetra::NotSpecified) {  // global_size_t
    Teuchos::ParameterList matrixList;
    matrixList.set("nx", nx);
    matrixList.set("matrixType", "Laplace1D");
    RCP<Matrix> A = BuildMatrix(matrixList, lib);
    return A;
  }  // Build1DPoisson()

  // Create a 2D Poisson matrix with the specified number of rows
  // nx: global number of rows
  // ny: global number of rows
  static RCP<Matrix> Build2DPoisson(GO nx, GO ny = -1, Xpetra::UnderlyingLib lib = Xpetra::NotSpecified) {  // global_size_t
    Teuchos::ParameterList matrixList;
    if (ny == -1) ny = nx;
    matrixList.set("nx", nx);
    matrixList.set("ny", ny);
    matrixList.set("matrixType", "Laplace2D");
    RCP<Matrix> A = BuildMatrix(matrixList, lib);
    return A;
  }  // Build2DPoisson()

  static RCP<RealValuedMultiVector>
  BuildGeoCoordinates(const int numDimensions, const Array<GO> gNodesPerDir,
                      Array<LO>& lNodesPerDir, Array<GO>& meshData,
                      const std::string meshLayout = "Global Lexicographic") {
    // Get MPI infos
    Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();
    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    LO numRanks                         = comm->getSize();
    LO myRank                           = comm->getRank();

    meshData.resize(10 * numRanks);

    ////////////////////////////////////
    //                                //
    //   Step 1: Compute map layout   //
    //                                //
    ////////////////////////////////////
    if (numRanks == 1) {
      if (meshLayout == "Local Lexicographic") {
        meshData[0] = 0;  // Local Proc
        meshData[1] = 0;  // Local Proc?
        meshData[2] = 0;  // Parent Proc
      }
      for (int dim = 0; dim < 3; ++dim) {
        lNodesPerDir[dim] = gNodesPerDir[dim];
        if (meshLayout == "Local Lexicographic") {
          meshData[3 + 2 * dim]     = 0;                      // Lowest  index in direction dim
          meshData[3 + 2 * dim + 1] = gNodesPerDir[dim] - 1;  // Highest index in direction dim
        }
      }
      if (meshLayout == "Local Lexicographic") {
        meshData[9] = 0;
      }
    } else if (numRanks == 4) {
      if (meshLayout == "Local Lexicographic") {
        // First fill the rank related info in meshData
        for (int rank = 0; rank < numRanks; ++rank) {
          meshData[10 * rank + 0] = rank;  // Local Proc
          meshData[10 * rank + 1] = rank;  // Local Proc?
          meshData[10 * rank + 2] = 0;     // Parent Proc
        }
      }

      if (numDimensions == 1) {
        LO numNodesPerRank = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / numRanks);
        if (myRank == numRanks - 1) {
          lNodesPerDir[0] = gNodesPerDir[0] - myRank * numNodesPerRank;
        } else {
          lNodesPerDir[0] = numNodesPerRank;
        }
        lNodesPerDir[1] = gNodesPerDir[1];

        if (meshLayout == "Local Lexicographic") {
          for (int rank = 0; rank < numRanks; ++rank) {
            meshData[10 * rank + 3] = rank * numNodesPerRank;
            meshData[10 * rank + 4] =
                (rank == numRanks - 1) ? (gNodesPerDir[0] - 1) : ((rank + 1) * numNodesPerRank - 1);
            meshData[10 * rank + 5] = 0;
            meshData[10 * rank + 6] = 0;
          }
        }
      } else {
        if (myRank == 0 || myRank == 2) {
          lNodesPerDir[0] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
        } else {
          lNodesPerDir[0] = std::floor(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
        }
        if (myRank == 0 || myRank == 1) {
          lNodesPerDir[1] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
        } else {
          lNodesPerDir[1] = std::floor(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
        }
        if (meshLayout == "Local Lexicographic") {
          for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
              int rank = 2 * j + i;
              if (i == 0) {
                meshData[10 * rank + 3] = 0;
                meshData[10 * rank + 4] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2) - 1;
              } else if (i == 1) {
                meshData[10 * rank + 3] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
                meshData[10 * rank + 4] = gNodesPerDir[0] - 1;
              }
              if (j == 0) {
                meshData[10 * rank + 5] = 0;
                meshData[10 * rank + 6] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2) - 1;
              } else if (j == 1) {
                meshData[10 * rank + 5] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
                meshData[10 * rank + 6] = gNodesPerDir[1] - 1;
              }
            }
          }
        }
      }
      lNodesPerDir[2] = gNodesPerDir[2];
      if (meshLayout == "Local Lexicographic") {
        for (int rank = 0; rank < numRanks; ++rank) {
          meshData[10 * rank + 7] = 0;
          meshData[10 * rank + 8] = gNodesPerDir[2] - 1;
        }
        meshData[9] = 0;
        for (int rank = 1; rank < numRanks; ++rank) {
          meshData[10 * rank + 9] = meshData[10 * (rank - 1) + 9] + (meshData[10 * (rank - 1) + 4] - meshData[10 * (rank - 1) + 3] + 1) * (meshData[10 * (rank - 1) + 6] - meshData[10 * (rank - 1) + 5] + 1) * (meshData[10 * (rank - 1) + 8] - meshData[10 * (rank - 1) + 7] + 1);
        }
      }
    }

    GO gNumNodes = 1;
    LO lNumNodes = 1;
    for (int dim = 0; dim < 3; ++dim) {
      gNumNodes = gNumNodes * gNodesPerDir[dim];
      lNumNodes = lNumNodes * lNodesPerDir[dim];
    }

    GO myGIDOffset = 0;
    Array<GO> myLowestTuple(3);
    if (numDimensions == 1) {
      myGIDOffset      = myRank * std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / numRanks);
      myLowestTuple[0] = myGIDOffset;
    } else {
      if (myRank == 1) {
        myLowestTuple[0] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
        if (meshLayout == "Global Lexicographic") {
          myGIDOffset = myLowestTuple[0];
        } else if (meshLayout == "Local Lexicographic") {
          myGIDOffset = myLowestTuple[0] * std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2) * gNodesPerDir[2];
        }
      } else if (myRank == 2) {
        myLowestTuple[1] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
        if (meshLayout == "Global Lexicographic") {
          myGIDOffset = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2) * gNodesPerDir[0];
        } else if (meshLayout == "Local Lexicographic") {
          myGIDOffset = myLowestTuple[1] * gNodesPerDir[0] * gNodesPerDir[2];
        }
      } else if (myRank == 3) {
        myLowestTuple[0] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
        myLowestTuple[1] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
        if (meshLayout == "Global Lexicographic") {
          myGIDOffset = myLowestTuple[1] * gNodesPerDir[0] + myLowestTuple[0];
        } else if (meshLayout == "Local Lexicographic") {
          myGIDOffset = (myLowestTuple[0] * (gNodesPerDir[1] - myLowestTuple[1]) + myLowestTuple[1] * gNodesPerDir[0]) * gNodesPerDir[2];
        }
      }
    }

    ////////////////////////////////////
    //                                //
    //    Step 2: Compute map GIDs    //
    //                                //
    ////////////////////////////////////
    Array<GO> myGIDs(lNumNodes);
    if (meshLayout == "Global Lexicographic") {
      for (LO k = 0; k < lNodesPerDir[2]; ++k) {
        for (LO j = 0; j < lNodesPerDir[1]; ++j) {
          for (LO i = 0; i < lNodesPerDir[0]; ++i) {
            myGIDs[k * lNodesPerDir[1] * lNodesPerDir[0] + j * lNodesPerDir[0] + i] =
                myGIDOffset + k * gNodesPerDir[1] * gNodesPerDir[0] + j * gNodesPerDir[0] + i;
          }
        }
      }
    } else if (meshLayout == "Local Lexicographic") {
      for (LO nodeIdx = 0; nodeIdx < lNumNodes; ++nodeIdx) {
        myGIDs[nodeIdx] = myGIDOffset + nodeIdx;
      }
    }

    RCP<const Map> coordMap = MapFactory::Build(lib, gNumNodes, myGIDs(), 0, comm);

    ///////////////////////////////////////
    //                                   //
    //    Step 3: Compute coordinates    //
    //                                   //
    ///////////////////////////////////////
    RCP<RealValuedMultiVector> Coordinates = RealValuedMultiVectorFactory::Build(coordMap,
                                                                                 numDimensions);
    Array<ArrayRCP<real_type> > myCoords(numDimensions);
    for (int dim = 0; dim < numDimensions; ++dim) {
      myCoords[dim] = Coordinates->getDataNonConst(dim);
    }

    LO nodeIdx = 0;
    Array<LO> ijk(3);
    for (ijk[2] = 0; ijk[2] < lNodesPerDir[2]; ++ijk[2]) {
      for (ijk[1] = 0; ijk[1] < lNodesPerDir[1]; ++ijk[1]) {
        for (ijk[0] = 0; ijk[0] < lNodesPerDir[0]; ++ijk[0]) {
          nodeIdx = ijk[2] * lNodesPerDir[1] * lNodesPerDir[0] + ijk[1] * lNodesPerDir[0] + ijk[0];
          for (int dim = 0; dim < numDimensions; ++dim) {
            if (gNodesPerDir[dim] == 1) {
              myCoords[dim][nodeIdx] = 0.0;
            } else {
              myCoords[dim][nodeIdx] =
                  Teuchos::as<real_type>(myLowestTuple[dim] + ijk[dim]) / (gNodesPerDir[dim] - 1);
            }
          }
        }
      }
    }

    return Coordinates;
  }  // BuildGeoCoordinates

  // Xpetra version of CreateMap
  static RCP<Map> BuildMap(Xpetra::UnderlyingLib lib, const std::set<GlobalOrdinal>& gids, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
    Teuchos::Array<GlobalOrdinal> mapvec;
    mapvec.reserve(gids.size());
    mapvec.assign(gids.begin(), gids.end());
    GlobalOrdinal count = Teuchos::as<GlobalOrdinal>(mapvec.size());
    GlobalOrdinal gcount;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, count, Teuchos::outArg(gcount));

    Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
        MapFactory::Build(lib, gcount, mapvec(), 0, comm);
    mapvec.clear();
    return map;
  }

  // Xpetra version of SplitMap
  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > SplitMap(Xpetra::UnderlyingLib lib, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& Amap, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& Agiven) {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Amap.getComm();

    GlobalOrdinal count = 0;
    Teuchos::Array<GlobalOrdinal> myaugids(Amap.getLocalNumElements());
    for (size_t i = 0; i < Amap.getLocalNumElements(); ++i) {
      const GlobalOrdinal gid = Amap.getGlobalElement(i);
      if (Agiven.isNodeGlobalElement(gid)) continue;
      myaugids[Teuchos::as<GlobalOrdinal>(count)] = gid;
      ++count;
    }
    myaugids.resize(count);
    GlobalOrdinal gcount;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &count, &gcount);
    return MapFactory::Build(lib, gcount, myaugids(), 0, comm);
  }

  static Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockDiagonalExampleMatrix(Xpetra::UnderlyingLib lib, int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
    GlobalOrdinal nOverallDOFGidsPerProc = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, noBlocks - 2)) * 10;

    GlobalOrdinal procOffset = comm->getRank() * nOverallDOFGidsPerProc;

    std::set<GlobalOrdinal> myDOFGids;
    for (GlobalOrdinal i = 0; i < nOverallDOFGidsPerProc; i++)
      myDOFGids.insert(i + procOffset);

    Teuchos::RCP<Map> fullmap = TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(lib, myDOFGids, comm);

    std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);
    GlobalOrdinal nPartGIDs            = nOverallDOFGidsPerProc;
    Teuchos::RCP<Map> remainingpartmap = fullmap;
    for (int it = 0; it < noBlocks; it++) {
      if (it == noBlocks - 1) {
        maps[0] = remainingpartmap;
        break;
      }
      // collect first half of GIDs
      nPartGIDs = nPartGIDs / 2;
      std::set<GlobalOrdinal> myHalfGIDs;
      for (GlobalOrdinal j = 0; j < nPartGIDs; j++)
        myHalfGIDs.insert(j + procOffset);

      Teuchos::RCP<Map> halfmap = TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(lib, myHalfGIDs, comm);

      Teuchos::RCP<Map> secondmap = TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SplitMap(lib, *remainingpartmap, *halfmap);
      remainingpartmap            = halfmap;

      maps[noBlocks - 1 - it] = secondmap;
    }

    // create diagonal blocks
    std::vector<Teuchos::RCP<CrsMatrix> > blocks(noBlocks, Teuchos::null);
    for (int it = 0; it < noBlocks; it++) {
      blocks[it] = CrsMatrixFactory::Build(maps[it], 1);

      LocalOrdinal NumMyElements                               = maps[it]->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = maps[it]->getLocalElementList();

      for (LocalOrdinal i = 0; i < NumMyElements; i++)
        blocks[it]->insertGlobalValues(MyGlobalElements[i],
                                       Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                       Teuchos::tuple<Scalar>(it + 1));
      blocks[it]->fillComplete();
    }

    // create map extractor
    Teuchos::RCP<const MapExtractor> rgMapExtractor = Teuchos::rcp(new MapExtractor(fullmap, maps, false));
    Teuchos::RCP<const MapExtractor> doMapExtractor = Teuchos::rcp(new MapExtractor(fullmap, maps, false));

    // build blocked operator
    Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 1));

    for (int it = 0; it < noBlocks; it++) {
      Teuchos::RCP<CrsMatrixWrap> csrwrap =
          Teuchos::rcp(new CrsMatrixWrap(blocks[it]));
      bop->setMatrix(Teuchos::as<size_t>(it), Teuchos::as<size_t>(it), csrwrap);
    }
    bop->fillComplete();
    return bop;
  }

  static Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlockDiagonalExampleMatrixThyra(Xpetra::UnderlyingLib lib, int noBlocks, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
    std::vector<Teuchos::RCP<const Map> > maps(noBlocks, Teuchos::null);

    maps[0] = MapFactory::Build(lib, comm->getSize() * 5, 5, 0, comm);
    for (int it = 1; it < noBlocks; it++) {
      GlobalOrdinal localDofs = Teuchos::as<GlobalOrdinal>(Teuchos::ScalarTraits<GlobalOrdinal>::pow(2, it - 1) * 5);
      maps[it]                = MapFactory::Build(lib, comm->getSize() * localDofs, localDofs, 0, comm);
    }

    // create diagonal blocks
    std::vector<Teuchos::RCP<CrsMatrix> > blocks(noBlocks, Teuchos::null);
    for (int it = 0; it < noBlocks; it++) {
      blocks[it] = CrsMatrixFactory::Build(maps[it], 1);

      LocalOrdinal NumMyElements                               = maps[it]->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = maps[it]->getLocalElementList();

      for (LocalOrdinal i = 0; i < NumMyElements; i++)
        blocks[it]->insertGlobalValues(MyGlobalElements[i],
                                       Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                       Teuchos::tuple<Scalar>(it + 1));
      blocks[it]->fillComplete();
    }

    // create map extractor
    // To generate the Thyra style map extractor we do not need a full map but only the
    // information about the Map details (i.e. lib and indexBase). We can extract this
    // information from maps[0]
    Teuchos::RCP<const MapExtractor> rgMapExtractor =
        Teuchos::rcp(new MapExtractor(maps[0], maps, true));
    Teuchos::RCP<const MapExtractor> doMapExtractor =
        Teuchos::rcp(new MapExtractor(maps[0], maps, true));
    // build blocked operator
    Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 1));

    for (int it = 0; it < noBlocks; it++) {
      Teuchos::RCP<CrsMatrixWrap> csrwrap = Teuchos::rcp(new CrsMatrixWrap(blocks[it]));
      bop->setMatrix(Teuchos::as<size_t>(it), Teuchos::as<size_t>(it), csrwrap);
    }
    bop->fillComplete();
    return bop;
  }

  static Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > CreateBlocked3x3MatrixThyra(const Teuchos::Comm<int>& comm, Xpetra::UnderlyingLib lib) {
    std::vector<RCP<const Map> > maps = std::vector<RCP<const Map> >(3, Teuchos::null);
    maps[0]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
    maps[1]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
    maps[2]                           = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMap(100);
    RCP<Matrix> A00                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], 4.0, -1.0, -1.0, lib);
    RCP<Matrix> A01                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[0], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A10                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A11                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], 4.0, -1.0, -1.0, lib);
    RCP<Matrix> A12                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[1], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A21                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], -1.0, 0.0, 0.0, lib);
    RCP<Matrix> A22                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildTridiag(maps[2], 4.0, -1.0, -1.0, lib);

    // create map extractor
    // To generate the Thyra style map extractor we do not need a full map but only the
    // information about the Map details (i.e. lib and indexBase). We can extract this
    // information from maps[0]
    Teuchos::RCP<const MapExtractor> rgMapExtractor =
        Teuchos::rcp(new MapExtractor(maps[0], maps, true));
    Teuchos::RCP<const MapExtractor> doMapExtractor =
        Teuchos::rcp(new MapExtractor(maps[0], maps, true));
    // build blocked operator
    Teuchos::RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExtractor, doMapExtractor, 5));
    bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(0), A00);
    bop->setMatrix(Teuchos::as<size_t>(0), Teuchos::as<size_t>(1), A01);
    bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(0), A10);
    bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(1), A11);
    bop->setMatrix(Teuchos::as<size_t>(1), Teuchos::as<size_t>(2), A12);
    bop->setMatrix(Teuchos::as<size_t>(2), Teuchos::as<size_t>(1), A21);
    bop->setMatrix(Teuchos::as<size_t>(2), Teuchos::as<size_t>(2), A22);
    bop->fillComplete();
    return bop;
  }

  // Needed to initialize correctly a level used for testing SingleLevel factory Build() methods.
  // This method initializes LevelID and linked list of level
  static void createSingleLevelHierarchy(Level& currentLevel) {
    RCP<FactoryManager> factoryHandler = rcp(new FactoryManager());
    factoryHandler->SetKokkosRefactor(false);
    currentLevel.SetFactoryManager(factoryHandler);

    currentLevel.SetLevelID(0);
  }

  // Needed to initialize correctly levels used for testing TwoLevel factory Build() methods.
  // This method initializes LevelID and linked list of level
  static void createTwoLevelHierarchy(Level& fineLevel, Level& coarseLevel) {
    RCP<FactoryManager> factoryHandler = rcp(new FactoryManager());
    factoryHandler->SetKokkosRefactor(false);
    fineLevel.SetFactoryManager(factoryHandler);
    coarseLevel.SetFactoryManager(factoryHandler);

    coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));

    fineLevel.SetLevelID(0);
    coarseLevel.SetLevelID(1);
  }

  static RCP<SmootherPrototype> createSmootherPrototype(const std::string& type = "Gauss-Seidel", LO sweeps = 1) {
    std::string ifpackType = "RELAXATION";
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: type", type);
    ifpackList.set("relaxation: sweeps", (LO)sweeps);
    ifpackList.set("relaxation: damping factor", (SC)1.0);
    return Teuchos::rcp(new TrilinosSmoother(ifpackType, ifpackList));
  }

  // Create a matrix as specified by parameter list options
  static RCP<Matrix> BuildBlockMatrixAsPoint(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) {
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    GO GO_INVALID                       = Teuchos::OrdinalTraits<GO>::invalid();
    RCP<Matrix> Op;

    if (lib == Xpetra::NotSpecified)
      lib = TestHelpers::Parameters::getLib();

    // Make the base graph
    RCP<Matrix> old_matrix        = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMatrix(matrixList, lib);
    RCP<const CrsGraph> old_graph = old_matrix->getCrsGraph();
    RCP<const Map> old_rowmap     = old_graph->getRowMap();
    RCP<const Map> old_colmap     = old_graph->getColMap();
    int blocksize                 = 3;

    // Block Map
    LO orig_num_rows = (LO)old_graph->getRowMap()->getLocalNumElements();
    Teuchos::Array<GlobalOrdinal> owned_rows(blocksize * orig_num_rows);
    for (LO i = 0; i < orig_num_rows; i++) {
      GO old_gid = old_rowmap->getGlobalElement(i);
      for (int j = 0; j < blocksize; j++) {
        owned_rows[i * blocksize + j] = old_gid * blocksize + j;
      }
    }
    RCP<Map> new_map = Xpetra::MapFactory<LO, GO, NO>::Build(lib, GO_INVALID, owned_rows(), 0, comm);
    if (new_map.is_null()) throw std::runtime_error("BuildBlockMatrixAsPoint: Map constructor failed");

    // Block Graph / Matrix
    RCP<CrsMatrix> new_matrix = Xpetra::CrsMatrixFactory<SC, LO, GO, NO>::Build(new_map, blocksize * old_graph->getLocalMaxNumRowEntries());
    if (new_matrix.is_null()) throw std::runtime_error("BuildBlockMatrixAsPoint: Matrix constructor failed");
    for (LO i = 0; i < orig_num_rows; i++) {
      Teuchos::ArrayView<const LO> old_indices;
      Teuchos::ArrayView<const SC> old_values;
      Teuchos::Array<GO> new_indices(1);
      Teuchos::Array<SC> new_values(1);
      old_matrix->getLocalRowView(i, old_indices, old_values);
      for (int ii = 0; ii < blocksize; ii++) {
        GO GRID = new_map->getGlobalElement(i * blocksize + ii);
        for (LO j = 0; j < (LO)old_indices.size(); j++) {
          for (int jj = 0; jj < blocksize; jj++) {
            new_indices[0] = old_colmap->getGlobalElement(old_indices[j]) * blocksize + jj;
            new_values[0]  = old_values[j] * (SC)((ii == jj && i == old_indices[j]) ? blocksize * blocksize : 1);
            new_matrix->insertGlobalValues(GRID, new_indices(), new_values);
          }
        }
      }
    }
    new_matrix->fillComplete();
    Op = rcp(new CrsMatrixWrap(new_matrix));
    if (new_map.is_null()) throw std::runtime_error("BuildBlockMatrixAsPoint: CrsMatrixWrap constructor failed");
    Op->SetFixedBlockSize(blocksize);

    return Op;
  }  // BuildBlockMatrixAsPoint()

};  // class TestFactory

// Helper class which has some Tpetra specific code inside
// We put this into an extra helper class as we need partial specializations and
// do not want to introduce partial specializations for the full TestFactory class
//
// The BuildBlockMatrix is only available with Tpetra. However, if both Epetra
// and Tpetra are enabled it may be that Tpetra is not instantiated on either
// GO=int/long long and/or Node=Serial/OpenMP. We need partial specializations
// with an empty BuildBlockMatrix routine for all instantiations Teptra is not
// enabled for, but are existing in Xpetra due to Epetra enabled.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraTestFactory {
#include "MueLu_UseShortNames.hpp"
 public:
  // Create a matrix as specified by parameter list options
  static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) {
    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
    RCP<Matrix> Op;

    if (lib == Xpetra::NotSpecified)
      lib = TestHelpers::Parameters::getLib();

    // This only works for Tpetra
    if (lib != Xpetra::UseTpetra) return Op;

    // Thanks for the code, Travis!

    // Make the graph
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > FirstMatrix = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMatrix(matrixList, lib);
    RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > FGraph      = FirstMatrix->getCrsGraph();

    int blocksize                                                                = 3;
    RCP<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> > TGraph = rcp_dynamic_cast<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> >(FGraph);
    RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > TTGraph      = TGraph->getTpetra_CrsGraph();

    RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bcrsmatrix = rcp(new Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*TTGraph, blocksize));

    const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& meshRowMap = *bcrsmatrix->getRowMap();
    const Scalar zero                                                = Teuchos::ScalarTraits<Scalar>::zero();
    const Scalar one                                                 = Teuchos::ScalarTraits<Scalar>::one();
    const Scalar two                                                 = one + one;
    const Scalar three                                               = two + one;

    Teuchos::Array<Scalar> basematrix(blocksize * blocksize, zero);
    basematrix[0] = two;
    basematrix[2] = three;
    basematrix[3] = three;
    basematrix[4] = two;
    basematrix[7] = three;
    basematrix[8] = two;
    Teuchos::Array<Scalar> offmatrix(blocksize * blocksize, zero);
    offmatrix[0] = offmatrix[4] = offmatrix[8] = -1;

    Teuchos::Array<LocalOrdinal> lclColInds(1);
    for (LocalOrdinal lclRowInd = meshRowMap.getMinLocalIndex(); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
      lclColInds[0] = lclRowInd;
      bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &basematrix[0], 1);

      // Off diagonals
      if (lclRowInd > meshRowMap.getMinLocalIndex()) {
        lclColInds[0] = lclRowInd - 1;
        bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &offmatrix[0], 1);
      }
      if (lclRowInd < meshRowMap.getMaxLocalIndex()) {
        lclColInds[0] = lclRowInd + 1;
        bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &offmatrix[0], 1);
      }
    }

    RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > temp = rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bcrsmatrix));
    Op                                                                      = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(temp));
    return Op;
  }  // BuildBlockMatrix()

 private:
  TpetraTestFactory() {}  // static class

};  // class TpetraTestFactory

// TAW: 3/14/2016: If both Epetra and Tpetra are enabled we need partial specializations
//                 on GO=int/long long as well as NO=EpetraNode to disable BuildBlockMatrix
#ifdef HAVE_MUELU_EPETRA
// partial specializations (GO=int not enabled with Tpetra)
#if !defined(HAVE_TPETRA_INST_INT_INT)
template <class Scalar, class LocalOrdinal, class Node>
class TpetraTestFactory<Scalar, LocalOrdinal, int, Node> {
  typedef int GlobalOrdinal;
#include "MueLu_UseShortNames.hpp"
 public:
  static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }
  static RCP<Matrix> BuildBlockMatrixAsPoint(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }

 private:
  TpetraTestFactory() {}  // static class
};                        // class TpetraTestFactory
#endif

// partial specializations (GO=long long not enabled with Tpetra)
#if !defined(HAVE_TPETRA_INST_INT_LONG_LONG)
template <class Scalar, class LocalOrdinal, class Node>
class TpetraTestFactory<Scalar, LocalOrdinal, long long, Node> {
  typedef long long GlobalOrdinal;
#include "MueLu_UseShortNames.hpp"
 public:
  static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }
  static RCP<Matrix> BuildBlockMatrixAsPoint(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }

 private:
  TpetraTestFactory() {}  // static class
};                        // class TpetraTestFactory
#endif

// partial specializations (NO=EpetraNode not enabled with Tpetra)
#if ((defined(EPETRA_HAVE_OMP) && !(defined(HAVE_TPETRA_INST_OPENMP))) || \
     (!defined(EPETRA_HAVE_OMP) && !(defined(HAVE_TPETRA_INST_SERIAL))))

template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
class TpetraTestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> {
  typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"
 public:
  static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }
  static RCP<Matrix> BuildBlockMatrixAsPoint(Teuchos::ParameterList& matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }

 private:
  TpetraTestFactory() {}  // static class
};                        // class TpetraTestFactory
#endif
#endif  // endif HAVE_MUELU_EPETRA

//! Return the list of files in the directory. Only files that are matching '*filter*' are returned.
ArrayRCP<std::string> GetFileList(const std::string& dirPath, const std::string& filter);

}  // namespace TestHelpers

}  // namespace MueLuTests

// Macro to skip a test when UnderlyingLib==Epetra or Tpetra
#define MUELU_TEST_ONLY_FOR(UnderlyingLib)                                                                                        \
  if (TestHelpers::Parameters::getLib() != UnderlyingLib) {                                                                       \
    out << "Skipping test for " << ((TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) ? "Epetra" : "Tpetra") << std::endl; \
    return;                                                                                                                       \
  }

// Macro to skip a test when Epetra is used with Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal) \
  if (!(TestHelpers::Parameters::getLib() == Xpetra::UseEpetra && (Teuchos::OrdinalTraits<LocalOrdinal>::name() != string("int") || Teuchos::OrdinalTraits<GlobalOrdinal>::name() != string("int"))))

// Macro to skip a test when Epetra is used with Scalar != double or Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal)                                        \
  if (!(TestHelpers::Parameters::getLib() == Xpetra::UseEpetra && Teuchos::ScalarTraits<Scalar>::name() != string("double"))) \
  MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal)

//

// TODO: add directly to Teuchos ?
//#include "../xpetra/test/Xpetra_UnitTestHelpers.hpp" // declaration of TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL

//

//! Namespace for MueLu test classes
namespace MueLuTests {

using namespace TestHelpers;
}

#endif  // ifndef MUELU_TEST_HELPERS_H

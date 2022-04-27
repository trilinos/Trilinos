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
#ifndef MUELU_TEST_HELPERS_KOKKOS_H
#define MUELU_TEST_HELPERS_KOKKOS_H

#include <string>
#ifndef _MSC_VER
#include <dirent.h>
#endif

// Teuchos
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

// Xpetra
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsGraph.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Level.hpp"

// Galeri
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraMatrixTypes.hpp>

#include "MueLu_NoFactory.hpp"



// Conditional Tpetra stuff
#ifdef HAVE_MUELU_TPETRA
#include <TpetraCore_config.h>
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Xpetra_TpetraRowMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#endif

#include <MueLu_TestHelpers_Common_kokkos.hpp>

namespace MueLuTests {
  using Teuchos::arcp;
  using Teuchos::arcpFromArrayView;
  using Teuchos::arcp_reinterpret_cast;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;
  using Teuchos::rcp_implicit_cast;

# include <MueLu_TestHelpers_Common_kokkos.hpp>

  namespace TestHelpers_kokkos {

    using Xpetra::global_size_t;

    class Parameters {

    private:
      Parameters() {} // static class

    public:

      static Xpetra::Parameters xpetraParameters;

      inline static RCP<const Teuchos::Comm<int> > getDefaultComm() {
        return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
      }

      inline static Xpetra::UnderlyingLib getLib() {
        return TestHelpers_kokkos::Parameters::xpetraParameters.GetLib();
      }
    };

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    class TestFactory {
#include "MueLu_UseShortNames.hpp"

    private:
      TestFactory() {} // static class

    public:

      typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
      typedef Xpetra::MultiVector<real_type,LO,GO,NO> RealValuedMultiVector;
      typedef Xpetra::MultiVectorFactory<real_type,LO,GO,NO> RealValuedMultiVectorFactory;

      // Create a map containing a specified number of local elements per process.
      static const RCP<const Map> BuildMap(LO numElementsPerProc) {
        RCP<const Teuchos::Comm<int> > comm   = TestHelpers_kokkos::Parameters::getDefaultComm();
        const global_size_t           INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();

        return MapFactory::Build(TestHelpers_kokkos::Parameters::getLib(), INVALID, numElementsPerProc, 0, comm);

      } // BuildMap()

      // Create a matrix as specified by parameter list options
      static RCP<Matrix> BuildMatrix(ParameterList& matrixList, Xpetra::UnderlyingLib lib) {
        RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();

        if (lib == Xpetra::NotSpecified)
          lib = TestHelpers_kokkos::Parameters::getLib();

        GO nx, ny, nz;
        nx = ny = nz = 5;
        nx = matrixList.get("nx", nx);
        ny = matrixList.get("ny", ny);
        nz = matrixList.get("nz", nz);

        std::string matrixType = matrixList.get("matrixType", "Laplace1D");
        GO numGlobalElements; //global_size_t
        if (matrixType == "Laplace1D")
          numGlobalElements = nx;
        else if (matrixType == "Laplace2D" || matrixType == "Star2D")
          numGlobalElements = nx*ny;
        else if (matrixType == "Laplace3D")
          numGlobalElements = nx*ny*nz;
        else {
          std::string msg = matrixType + " is unsupported (in unit testing)";
          throw(MueLu::Exceptions::RuntimeError(msg));
        }

        RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, 0, comm);
        RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixType, map, matrixList);
        RCP<Matrix> Op = Pr->BuildMatrix();

        return Op;
      }

      // Create a 1D Poisson matrix with the specified number of rows
      // nx: global number of rows
      static RCP<Matrix> Build1DPoisson(GO nx, Xpetra::UnderlyingLib lib = Xpetra::NotSpecified) { //global_size_t
        ParameterList matrixList;
        matrixList.set("nx",            nx);
        matrixList.set("matrixType",    "Laplace1D");
        return BuildMatrix(matrixList,lib);
      }

      // Create a 2D Poisson matrix with the specified number of rows
      // nx: global number of rows
      // ny: global number of rows
      static RCP<Matrix> Build2DPoisson(GO nx, GO ny = -1, Xpetra::UnderlyingLib lib = Xpetra::NotSpecified) { //global_size_t
        if (ny == -1) ny = nx;

        ParameterList matrixList;
        matrixList.set("nx",            nx);
        matrixList.set("ny",            ny);
        matrixList.set("matrixType",    "Laplace2D");

        return BuildMatrix(matrixList,lib);
      }

      // Create a tridiagonal matrix (stencil = [b,a,c]) with the specified number of rows
      static RCP<Matrix> BuildTridiag(RCP<const Map> rowMap, SC a, SC b, SC c, Xpetra::UnderlyingLib lib = Xpetra::NotSpecified) {
        if (lib == Xpetra::NotSpecified)
          lib = TestHelpers_kokkos::Parameters::getLib();

        RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

        RCP<Matrix> mtx = MatrixFactory::Build(rowMap, 3);

        ArrayView<const GO> MyGlobalElements = rowMap->getLocalElementList();
        GO indexBase = rowMap->getIndexBase();

        GO NumEntries;
        LO nnz = 3;
        std::vector<SC> Values(nnz);
        std::vector<GO> Indices(nnz);

        // access information from strided block map
        std::vector<size_t> strInfo = std::vector<size_t>();
        size_t blockSize = 1;
        LO blockId = -1;
        GO offset = 0;

        RCP<const StridedMap> strrowMap = rcp_dynamic_cast<const StridedMap>(rowMap);
        if (strrowMap != Teuchos::null) {
          strInfo   = strrowMap->getStridingData();
          blockSize = strrowMap->getFixedBlockSize();
          blockId   = strrowMap->getStridedBlockId();
          offset    = strrowMap->getOffset();
          TEUCHOS_TEST_FOR_EXCEPTION(blockId > -1 && strInfo[blockId]==1, MueLu::Exceptions::RuntimeError,
                                     "MueLu::TestHelpers_kokkos::BuildTridiag: strInfo block size must be > 1.");
          // todo: write one more special case for a row being first and last bock row
        } else {
          // no striding information. emulate block matrix
          blockId = 0;
          strInfo.push_back(blockSize); // default block size = 1
        }

        GO rrmax = (rowMap->getMaxAllGlobalIndex()-offset-indexBase) / blockSize;

        // loop over all rows
        LO numMyElements = rowMap->getLocalNumElements();
        for (LO i = 0; i < numMyElements; i++) {
          GO rr = (MyGlobalElements[i]-offset-indexBase) / blockSize + indexBase;  // node index

          // distinguish 5 different cases

          GO blockOffset = 0;
          for (LO k = 0; k < blockId; k++)
            blockOffset += Teuchos::as<GO>(strInfo[k]);

          if (MyGlobalElements[i] == blockOffset + offset + indexBase) {
            // very first row
            Indices[0] = MyGlobalElements[i];
            Values [0] = a;
            Indices[1] = MyGlobalElements[i] + 1;
            Values [1] = c;
            NumEntries = 2;

          } else if (MyGlobalElements[i] == rrmax * Teuchos::as<GO>(blockSize) + blockOffset + Teuchos::as<GO>(strInfo[blockId]) - 1 + offset + indexBase) {
            // very last row
            Indices[0] = MyGlobalElements[i] - 1;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            NumEntries = 2;

          } else if (MyGlobalElements[i] == rr * Teuchos::as<GO>(blockSize) + blockOffset + Teuchos::as<GO>(strInfo[blockId]) - 1 + offset + indexBase) {
            // last row in current node block
            Indices[0] = MyGlobalElements[i] - 1;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            Indices[2] = (rr+1)*blockSize + blockOffset + offset + indexBase;
            Values [2] = c;
            NumEntries = 3;

          } else if (MyGlobalElements[i] == rr * Teuchos::as<GO>(blockSize) + blockOffset + offset + indexBase) {
            // first row in current node block
            Indices[0] = (rr-1)*blockSize + blockOffset + strInfo[blockId] - 1 + offset + indexBase;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            Indices[2] = MyGlobalElements[i] + 1;
            Values [2] = c;
            NumEntries = 3;

          } else {
            // usual row entries in block rows
            Indices[0] = MyGlobalElements[i] - 1;
            Values [0] = b;
            Indices[1] = MyGlobalElements[i];
            Values [1] = a;
            Indices[2] = MyGlobalElements[i] + 1;
            Values [2] = c;
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
          ArrayView<SC> av(&Values [0], NumEntries);
          ArrayView<GO> iv(&Indices[0], NumEntries);
          mtx->insertGlobalValues(MyGlobalElements[i], iv, av);
        }

        mtx->fillComplete();

        return mtx;
      }

      static RCP<RealValuedMultiVector>
      BuildGeoCoordinates(const int numDimensions, const Array<GO> gNodesPerDir,
                          Array<LO>& lNodesPerDir, Array<GO>& meshData,
                          const std::string meshLayout = "Global Lexicographic") {

        // Get MPI infos
        Xpetra::UnderlyingLib lib = TestHelpers_kokkos::Parameters::getLib();
        RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();
        LO numRanks = comm->getSize();
        LO myRank   = comm->getRank();

        meshData.resize(10*numRanks);

        ////////////////////////////////////
        //                                //
        //   Step 1: Compute map layout   //
        //                                //
        ////////////////////////////////////
        if(numRanks == 1) {
          if(meshLayout == "Local Lexicographic") {
            meshData[0] = 0; // Local Proc
            meshData[1] = 0; // Local Proc?
            meshData[2] = 0; // Parent Proc
          }
          for(int dim = 0; dim < 3; ++dim) {
            lNodesPerDir[dim] = gNodesPerDir[dim];
            if(meshLayout == "Local Lexicographic") {
              meshData[3 + 2*dim]    = 0;                        // Lowest  index in direction dim
              meshData[3 + 2*dim +1] = gNodesPerDir[dim] - 1;    // Highest index in direction dim
            }
          }
          if(meshLayout == "Local Lexicographic") {
            meshData[9] = 0;
          }
        } else if(numRanks == 4) {
          if(meshLayout == "Local Lexicographic") {
            // First fill the rank related info in meshData
            for(int rank = 0; rank < numRanks; ++rank) {
              meshData[10*rank + 0] = rank; // Local Proc
              meshData[10*rank + 1] = rank; // Local Proc?
              meshData[10*rank + 2] = 0;        // Parent Proc
            }
          }

          if(numDimensions == 1) {
            LO numNodesPerRank = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / numRanks);
            if(myRank == numRanks - 1) {
              lNodesPerDir[0] = gNodesPerDir[0] - myRank*numNodesPerRank;
            } else {
              lNodesPerDir[0] = numNodesPerRank;
            }
            lNodesPerDir[1] = gNodesPerDir[1];

            if(meshLayout == "Local Lexicographic") {
              for(int rank = 0; rank < numRanks; ++rank) {
                meshData[10*rank + 3] = rank*numNodesPerRank;
                meshData[10*rank + 4] =
                  (rank == numRanks - 1) ? (gNodesPerDir[0] - 1) : ((rank + 1)*numNodesPerRank - 1);
                meshData[10*rank + 5] = 0;
                meshData[10*rank + 6] = 0;
              }
            }
          } else {
            if(myRank == 0 || myRank == 2) {
              lNodesPerDir[0] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
            } else {
              lNodesPerDir[0] = std::floor(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
            }
            if(myRank == 0 || myRank == 1) {
              lNodesPerDir[1] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
            } else {
              lNodesPerDir[1] = std::floor(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
            }
            if(meshLayout == "Local Lexicographic") {
              for(int j = 0; j < 2; ++j) {
                for(int i = 0; i < 2; ++i) {
                  int rank = 2*j + i;
                  if(i == 0) {
                    meshData[10*rank + 3] = 0;
                    meshData[10*rank + 4] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2) - 1;
                  } else if(i == 1) {
                    meshData[10*rank + 3] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
                    meshData[10*rank + 4] = gNodesPerDir[0] - 1;
                  }
                  if(j == 0) {
                    meshData[10*rank + 5] = 0;
                    meshData[10*rank + 6] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2) - 1;
                  } else if(j == 1) {
                    meshData[10*rank + 5] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
                    meshData[10*rank + 6] = gNodesPerDir[1] - 1;
                  }
                }
              }
            }
          }
          lNodesPerDir[2] = gNodesPerDir[2];
          if(meshLayout == "Local Lexicographic") {
            for(int rank = 0; rank < numRanks; ++rank) {
              meshData[10*rank + 7] = 0;
              meshData[10*rank + 8] = gNodesPerDir[2] - 1;
            }
            meshData[9] = 0;
            for(int rank = 1; rank < numRanks; ++rank) {
              meshData[10*rank + 9] = meshData[10*(rank - 1) + 9]
                + (meshData[10*(rank - 1) + 4] - meshData[10*(rank - 1) + 3] + 1)
                *(meshData[10*(rank - 1) + 6] - meshData[10*(rank - 1) + 5] + 1)
                *(meshData[10*(rank - 1) + 8] - meshData[10*(rank - 1) + 7] + 1);
            }
          }
        }

        GO gNumNodes = 1;
        LO lNumNodes = 1;
        for(int dim = 0; dim < 3; ++dim) {
          gNumNodes = gNumNodes*gNodesPerDir[dim];
          lNumNodes = lNumNodes*lNodesPerDir[dim];
        }

        GO myGIDOffset = 0;
        Array<GO> myLowestTuple(3);
        if(numDimensions == 1) {
          myGIDOffset = myRank*std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / numRanks);
          myLowestTuple[0] = myGIDOffset;
        } else {
          if(myRank == 1) {
            myLowestTuple[0] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
            if(meshLayout == "Global Lexicographic") {
              myGIDOffset = myLowestTuple[0];
            } else if(meshLayout == "Local Lexicographic") {
              myGIDOffset = myLowestTuple[0]*std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2)
                *gNodesPerDir[2];
            }
          } else if(myRank == 2) {
            myLowestTuple[1] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
            if(meshLayout == "Global Lexicographic") {
              myGIDOffset = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2)*gNodesPerDir[0];
            } else if(meshLayout == "Local Lexicographic") {
              myGIDOffset = myLowestTuple[1]*gNodesPerDir[0]*gNodesPerDir[2];
            }
          } else if(myRank == 3) {
            myLowestTuple[0] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[0]) / 2);
            myLowestTuple[1] = std::ceil(Teuchos::as<real_type>(gNodesPerDir[1]) / 2);
            if(meshLayout == "Global Lexicographic") {
              myGIDOffset = myLowestTuple[1]*gNodesPerDir[0] + myLowestTuple[0];
            } else if(meshLayout == "Local Lexicographic") {
              myGIDOffset = (myLowestTuple[0]*(gNodesPerDir[1] - myLowestTuple[1])
                             + myLowestTuple[1]*gNodesPerDir[0])*gNodesPerDir[2];
            }
          }
        }

        ////////////////////////////////////
        //                                //
        //    Step 2: Compute map GIDs    //
        //                                //
        ////////////////////////////////////
        Array<GO> myGIDs(lNumNodes);
        if(meshLayout == "Global Lexicographic") {
          for(LO k = 0; k < lNodesPerDir[2]; ++k) {
            for(LO j = 0; j < lNodesPerDir[1]; ++j) {
              for(LO i = 0; i < lNodesPerDir[0]; ++i) {
                myGIDs[k*lNodesPerDir[1]*lNodesPerDir[0] + j*lNodesPerDir[0] + i] =
                  myGIDOffset + k*gNodesPerDir[1]*gNodesPerDir[0] + j*gNodesPerDir[0] + i;
              }
            }
          }
        } else if(meshLayout == "Local Lexicographic") {
          for(LO nodeIdx = 0; nodeIdx < lNumNodes; ++nodeIdx) {
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
        for(int dim = 0; dim < numDimensions; ++dim) {
          myCoords[dim] = Coordinates->getDataNonConst(dim);
        }

        LO nodeIdx = 0;
        Array<LO> ijk(3);
        for(ijk[2] = 0; ijk[2] < lNodesPerDir[2]; ++ijk[2]) {
          for(ijk[1] = 0; ijk[1] < lNodesPerDir[1]; ++ijk[1]) {
            for(ijk[0] = 0; ijk[0] < lNodesPerDir[0]; ++ijk[0]) {
              nodeIdx = ijk[2]*lNodesPerDir[1]*lNodesPerDir[0] + ijk[1]*lNodesPerDir[0] + ijk[0];
              for(int dim = 0; dim < numDimensions; ++dim) {
                 if(gNodesPerDir[dim] == 1) {
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
      }

      // Needed to initialize correctly a level used for testing SingleLevel factory Build() methods.
      // This method initializes LevelID and linked list of level
      static void createSingleLevelHierarchy(Level& currentLevel) {
        RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
        currentLevel.SetFactoryManager(factoryHandler);

        currentLevel.SetLevelID(0);
#ifdef HAVE_MUELU_TIMER_SYNCHRONIZATION
        currentLevel.SetComm(TestHelpers_kokkos::Parameters::getDefaultComm());
#endif
      }

      // Needed to initialize correctly levels used for testing TwoLevel factory Build() methods.
      // This method initializes LevelID and linked list of level
      static void createTwoLevelHierarchy(Level& fineLevel, Level& coarseLevel) {
        RCP<MueLu::FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
        fineLevel  .SetFactoryManager(factoryHandler);
        coarseLevel.SetFactoryManager(factoryHandler);

        coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));

        fineLevel  .SetLevelID(0);
        coarseLevel.SetLevelID(1);
#ifdef HAVE_MUELU_TIMER_SYNCHRONIZATION
        fineLevel  .SetComm(TestHelpers_kokkos::Parameters::getDefaultComm());
        coarseLevel.SetComm(TestHelpers_kokkos::Parameters::getDefaultComm());
#endif
      }

#if 0
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
      static RCP<SmootherPrototype> createSmootherPrototype(const std::string& type="Gauss-Seidel", LO sweeps=1) {
        Teuchos::ParameterList  ifpackList;
        ifpackList.set("relaxation: type", type);
        ifpackList.set("relaxation: sweeps", (LO) sweeps);
        ifpackList.set("relaxation: damping factor", (SC) 1.0);
        return rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
      }
#endif
#endif
    }; // class TestFactory

    // Helper class which has some Tpetra specific code inside
    // We put this into an extra helper class as we need partial specializations and
    // do not want to introduce partial specializations for the full TestFactory class
    //
    // The BuildBlockMatrix is only available with Teptra. However, if both Epetra
    // and Tpetra are enabled it may be that Tpetra is not instantiated on either
    // GO=int/long long and/or Node=Serial/OpenMP. We need partial specializations
    // with an empty BuildBlockMatrix routine for all instantiations Teptra is not
    // enabled for, but are existing in Xpetra due to Epetra enabled.
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    class TpetraTestFactory {
#include "MueLu_UseShortNames.hpp"
    public:

      // Create a matrix as specified by parameter list options
      static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList &matrixList, Xpetra::UnderlyingLib lib) {
        RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();
        RCP<Matrix> Op;

         if (lib == Xpetra::NotSpecified)
           lib = TestHelpers_kokkos::Parameters::getLib();

         // This only works for Tpetra
         if (lib!=Xpetra::UseTpetra) return Op;

#if defined(HAVE_MUELU_TPETRA)
         // Thanks for the code, Travis!

         // Make the graph
         RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > FirstMatrix = TestHelpers_kokkos::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildMatrix(matrixList,lib);
         RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > Graph = FirstMatrix->getCrsGraph();

         int blocksize = 3;
         RCP<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> > TGraph = rcp_dynamic_cast<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> >(Graph);
         RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > TTGraph = TGraph->getTpetra_CrsGraph();

         RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bcrsmatrix = rcp(new Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (*TTGraph, blocksize));

         const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& meshRowMap = *bcrsmatrix->getRowMap();
         const Scalar zero   = Teuchos::ScalarTraits<Scalar>::zero();
         const Scalar one   = Teuchos::ScalarTraits<Scalar>::one();
         const Scalar two   = one+one;
         const Scalar three = two+one;

         Teuchos::Array<Scalar> basematrix(blocksize*blocksize, zero);
         basematrix[0] = two;
         basematrix[2] = three;
         basematrix[3] = three;
         basematrix[4] = two;
         basematrix[7] = three;
         basematrix[8] = two;
         Teuchos::Array<LocalOrdinal> lclColInds(1);
         for (LocalOrdinal lclRowInd = meshRowMap.getMinLocalIndex (); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
           lclColInds[0] = lclRowInd;
           bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &basematrix[0], 1);
         }

         RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > temp = rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bcrsmatrix));
         Op = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(temp));
#endif
         return Op;
      } // BuildBlockMatrix()

    private:
      TpetraTestFactory() {} // static class

    }; // class TpetraTestFactory

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
      static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList &matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }
    private:
      TpetraTestFactory() {} // static class
    }; // class TpetraTestFactory
#endif

    // partial specializations (GO=long long not enabled with Tpetra)
#if !defined(HAVE_TPETRA_INST_INT_LONG_LONG)
    template <class Scalar, class LocalOrdinal, class Node>
    class TpetraTestFactory<Scalar, LocalOrdinal, long long, Node> {
      typedef long long GlobalOrdinal;
#include "MueLu_UseShortNames.hpp"
    public:
      static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList &matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }
    private:
      TpetraTestFactory() {} // static class
    }; // class TpetraTestFactory
#endif

    // partial specializations (NO=EpetraNode not enabled with Tpetra)
#if ((defined(EPETRA_HAVE_OMP) && !(defined(HAVE_TPETRA_INST_OPENMP))) || \
    (!defined(EPETRA_HAVE_OMP) && !(defined(HAVE_TPETRA_INST_SERIAL))))

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
    class TpetraTestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> {
      typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"
    public:
      static RCP<Matrix> BuildBlockMatrix(Teuchos::ParameterList &matrixList, Xpetra::UnderlyingLib lib) { return Teuchos::null; }
    private:
      TpetraTestFactory() {} // static class
    }; // class TpetraTestFactory
#endif
#endif // endif HAVE_MUELU_EPETRA

    //! Return the list of files in the directory. Only files that are matching '*filter*' are returned.
    ArrayRCP<std::string> GetFileList(const std::string& dirPath, const std::string & filter);

  } // namespace TestHelpers_kokkos

} // namespace MueLu


// Macro to skip a test when UnderlyingLib==Epetra or Tpetra
#define MUELU_TEST_ONLY_FOR(UnderlyingLib) \
  if (TestHelpers_kokkos::Parameters::getLib() != UnderlyingLib) { \
    out << "Skipping test for " << ((TestHelpers_kokkos::Parameters::getLib()==Xpetra::UseEpetra) ? "Epetra" : "Tpetra") << std::endl; \
    return; \
  }

// Macro to skip a test when Epetra is used with Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal) \
  if (!(TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra && (Teuchos::OrdinalTraits<LocalOrdinal>::name() != string("int") || Teuchos::OrdinalTraits<GlobalOrdinal>::name() != string("int"))))

// Macro to skip a test when Epetra is used with Scalar != double or Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) \
  if (!(TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra && Teuchos::ScalarTraits<Scalar>::name() != string("double"))) \
    MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal)


//! Namespace for MueLu test classes
namespace MueLuTests {

  using namespace TestHelpers_kokkos;
}

#endif // ifndef MUELU_TEST_HELPERS_KOKKOS_H

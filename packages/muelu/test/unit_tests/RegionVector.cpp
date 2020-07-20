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
#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Galeri_XpetraParameters.hpp>

// Region MG headers
#include "SetupRegionUtilities.hpp"
#include "SetupRegionVector_def.hpp"
#include "SetupRegionMatrix_def.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RegionVector, RegionCompositeVector, Scalar, LocalOrdinal, GlobalOrdinal, Node)
{
#   include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  // Get MPI parameter
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  LO numRanks = comm->getSize();
  LO myRank = comm->getRank();

  GO nx = 5, ny = 5, nz = 1;
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");
  Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
  std::string   matrixType = galeriParameters.GetMatrixType();

  // Build maps for the problem
  const LO numDofsPerNode = 1;
  RCP<Map> nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(TestHelpers::Parameters::getLib(),
                                                             "Cartesian2D", comm, galeriList);
  RCP<Map> dofMap  = Xpetra::MapFactory<LO,GO,Node>::Build(nodeMap, numDofsPerNode);

  // Build the Xpetra problem
  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
    Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), dofMap, galeriList);

  // Generate the operator
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(numDofsPerNode);

  // Set global geometric data
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  Array<GO> procsPerDim(3);
  gNodesPerDir[0] = galeriList.get<GO>("nx");
  gNodesPerDir[1] = galeriList.get<GO>("ny");
  gNodesPerDir[2] = 1;
  lNodesPerDir[0] = galeriList.get<LO>("lnx");
  lNodesPerDir[1] = galeriList.get<LO>("lny");
  lNodesPerDir[2] = 1;
  procsPerDim[0] = galeriList.get<GO>("mx");
  procsPerDim[1] = galeriList.get<GO>("my");
  procsPerDim[2] = 1;

  Array<int> boundaryConditions;
  int maxRegPerGID = 0;
  int numInterfaces = 0;
  LO numLocalRegionNodes = 0;
  Array<GO>  sendGIDs;
  Array<int> sendPIDs;
  Array<LO>  rNodesPerDim(3);
  Array<LO>  compositeToRegionLIDs(nodeMap->getNodeNumElements()*numDofsPerNode);
  Array<GO>  quasiRegionGIDs;
  Array<GO>  quasiRegionCoordGIDs;
  Array<GO>  interfaceGIDs;
  Array<LO>  interfaceLIDsData;
  createRegionData(2, false, numDofsPerNode,
                   gNodesPerDir(), lNodesPerDir(), procsPerDim(),
                   nodeMap, dofMap,
                   maxRegPerGID, numLocalRegionNodes, boundaryConditions,
                   sendGIDs, sendPIDs, numInterfaces, rNodesPerDim,
                   quasiRegionGIDs, quasiRegionCoordGIDs, compositeToRegionLIDs,
                   interfaceGIDs, interfaceLIDsData);

  const int maxRegPerProc = 1;
  std::vector<RCP<Map> > rowMapPerGrp(maxRegPerProc),        colMapPerGrp(maxRegPerProc);
  rowMapPerGrp[0] = Xpetra::MapFactory<LO,GO,Node>::Build(A->getRowMap()->lib(),
                                                          Teuchos::OrdinalTraits<GO>::invalid(),
                                                          quasiRegionGIDs(),
                                                          A->getRowMap()->getIndexBase(),
                                                          A->getRowMap()->getComm());
  colMapPerGrp[0] = rowMapPerGrp[0];

  std::vector<RCP<Map> > revisedRowMapPerGrp(maxRegPerProc), revisedColMapPerGrp(maxRegPerProc);
  revisedRowMapPerGrp[0] = Xpetra::MapFactory<LO,GO,Node>::Build(A->getRowMap()->lib(),
                                                                 Teuchos::OrdinalTraits<GO>::invalid(),
                                                                 quasiRegionGIDs.size()*numDofsPerNode,
                                                                 A->getRowMap()->getIndexBase(),
                                                                 A->getRowMap()->getComm());
  revisedColMapPerGrp[0] = revisedRowMapPerGrp[0];

  std::vector<RCP<Import> > rowImportPerGrp(maxRegPerProc);
  std::vector<RCP<Import> > colImportPerGrp(maxRegPerProc);
  rowImportPerGrp[0] = ImportFactory::Build(dofMap, rowMapPerGrp[0]);
  colImportPerGrp[0] = ImportFactory::Build(dofMap, colMapPerGrp[0]);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  { // test compositeToRegional
    // Create a composite vector
    RCP<Vector> compVec = VectorFactory::Build(dofMap);
    TEUCHOS_ASSERT(!compVec.is_null());
    compVec->putScalar(0.0);
    const size_t localLength = compVec->getLocalLength();
    for (size_t k = 0; k < localLength; ++k) {
      compVec->replaceLocalValue(k, static_cast<Scalar>(k) + static_cast<magnitude_type>(myRank)/10);
    }

    // Create a region vector
    Array<RCP<Vector> > regVec(maxRegPerProc);
    Array<Teuchos::RCP<Vector> > quasiRegVec(maxRegPerProc);

    compositeToRegional(compVec, quasiRegVec, regVec, revisedRowMapPerGrp, rowImportPerGrp);

    if(numRanks == 1){
      TEST_EQUALITY(regVec[0]->getLocalLength(),  25);
      TEST_EQUALITY(regVec[0]->getGlobalLength(), 25);
      TEST_EQUALITY(quasiRegVec[0]->getLocalLength(),  25);
      TEST_EQUALITY(quasiRegVec[0]->getGlobalLength(), 25);
      ArrayRCP<SC> refValues;
      Teuchos::ArrayRCP<const SC> myValues;
      refValues.deepCopy(ArrayView<const SC>({0.0, 1, 2, 3, 4,
                                               5, 6, 7, 8, 9,
                                              10, 11, 12, 13, 14,
                                              15, 16, 17, 18, 19,
                                              20, 21, 22, 23, 24}));
      for (int j = 0; j < maxRegPerProc; j++) {
        myValues = regVec[j]->getData(0);
        for(size_t idx = 0; idx < regVec[j]->getLocalLength(); ++idx) {
          TEST_FLOATING_EQUALITY(myValues[idx], refValues[idx], 100*TMT::eps());
        }
        myValues = quasiRegVec[j]->getData(0);
        for(size_t idx = 0; idx < quasiRegVec[j]->getLocalLength(); ++idx) {
          TEST_FLOATING_EQUALITY(myValues[idx], refValues[idx], 100*TMT::eps());
        }

      }

    } else if(numRanks == 4) {
      // All ranks will have the same number of rows/cols/entries
      TEST_EQUALITY(regVec[0]->getLocalLength(),  9);
      TEST_EQUALITY(regVec[0]->getGlobalLength(), 36);
      TEST_EQUALITY(quasiRegVec[0]->getLocalLength(),  9);
      TEST_EQUALITY(quasiRegVec[0]->getGlobalLength(), 36);

      ArrayRCP<SC> refValues;
      Teuchos::ArrayRCP<const SC> myValues;
      if(myRank == 0) {
        refValues.deepCopy(ArrayView<const SC>({0.0, 1, 2, 3, 4, 5, 6, 7, 8}));

      } else if(myRank == 1) {
        refValues.deepCopy(ArrayView<const SC>({2.0, 0.1, 1.1, 5.0, 2.1, 3.1, 8.0, 4.1, 5.1}));

      } else if(myRank == 2) {
        refValues.deepCopy(ArrayView<const SC>({6, 7, 8, 0.2, 1.2, 2.2, 3.2, 4.2, 5.2}));

      } else if(myRank == 3) {
        refValues.deepCopy(ArrayView<const SC>({8, 4.1, 5.1, 2.2, 0.3, 1.3, 5.2, 2.3, 3.3}));
      }
      // Loop over region matrix data and compare it to ref data
      for (int j = 0; j < maxRegPerProc; j++){
        myValues = regVec[j]->getData(0);
        for(size_t idx = 0; idx < regVec[j]->getLocalLength(); ++idx) {
          TEST_FLOATING_EQUALITY(myValues[idx], refValues[idx], 100*TMT::eps());
        }
        myValues = quasiRegVec[j]->getData(0);
        for(size_t idx = 0; idx < quasiRegVec[j]->getLocalLength(); ++idx) {
          TEST_FLOATING_EQUALITY(myValues[idx], refValues[idx], 100*TMT::eps());
        }
      }
    }

  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  {
    // Create a region vector
    Array<RCP<Vector> > regVec(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      regVec[j] = VectorFactory::Build(revisedRowMapPerGrp[j]);
      regVec[j]->putScalar( myRank / 10.0 );
      for (size_t k = 0; k < regVec[j]->getLocalLength(); ++k){
        regVec[j]->sumIntoLocalValue(k, k);
      }
    }

    // Create a composite vector
    RCP<Vector> compVec = VectorFactory::Build(dofMap, true);

    regionalToComposite(regVec, compVec, rowImportPerGrp);

    if(numRanks == 1) {
      TEST_EQUALITY(compVec->getLocalLength(),  25);
      TEST_EQUALITY(compVec->getGlobalLength(), 25);

      ArrayRCP<SC> refValues;
      Teuchos::ArrayRCP<const SC> myValues;
      refValues.deepCopy(ArrayView<const SC>({0.0, 1.0, 2.0, 3.0, 4.0,
                                              5.0, 6.0, 7.0, 8.0, 9.0,
                                              10.0, 11.0, 12.0, 13.0, 14.0,
                                              15.0, 16.0, 17.0, 18.0, 19.0,
                                              20.0, 21.0, 22.0, 23.0, 24.0}));
      myValues = compVec->getData(0);
      for(size_t idx = 0; idx < compVec->getLocalLength(); ++idx) {
        TEST_FLOATING_EQUALITY(myValues[idx], refValues[idx], 100*TMT::eps());
      }

    } else if(numRanks == 4) {
      // All ranks will have the same global length of the vector, but local length differs.
      // Due to region numbering and region interface duplications,
      // lower rank IDs have more local entries then higher rank IDs.
      if( myRank == 0){
        TEST_EQUALITY(compVec->getLocalLength(),  9);
        TEST_EQUALITY(compVec->getGlobalLength(), 25);
      } else if (myRank == 1) {
        TEST_EQUALITY(compVec->getLocalLength(),  6);
        TEST_EQUALITY(compVec->getGlobalLength(), 25);
      } else if (myRank == 2) {
        TEST_EQUALITY(compVec->getLocalLength(),  6);
        TEST_EQUALITY(compVec->getGlobalLength(), 25);
      } else if (myRank == 3) {
        TEST_EQUALITY(compVec->getLocalLength(),  4);
        TEST_EQUALITY(compVec->getGlobalLength(), 25);
      }

      ArrayRCP<SC> refValues;
      Teuchos::ArrayRCP<const SC> myValues;
      if(myRank == 0) {
        refValues.deepCopy(ArrayView<const SC>({0.0, 1.0, 2.1, 3.0, 4.0, 8.1, 6.2, 8.2, 16.6}));

      } else if(myRank == 1) {
        refValues.deepCopy(ArrayView<const SC>({1.1, 2.1, 4.1, 5.1, 8.4, 10.4}));

      } else if(myRank == 2) {
        refValues.deepCopy(ArrayView<const SC>({3.2, 4.2, 8.5, 6.2, 7.2, 14.5}));

      } else if(myRank == 3) {
        refValues.deepCopy(ArrayView<const SC>({4.3, 5.3, 7.3, 8.3}));
      }
      // Loop over region vector data and compare it to reference data
      myValues = compVec->getData(0);
      for(size_t idx = 0; idx < compVec->getLocalLength(); ++idx) {
        TEST_FLOATING_EQUALITY(myValues[idx], refValues[idx], 100*TMT::eps());
      }
    }
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  {
    // initialize region vector with all ones.
    Array<RCP<Vector> > interfaceScaling(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      interfaceScaling[j] = VectorFactory::Build(revisedRowMapPerGrp[j]);
      interfaceScaling[j]->putScalar(1.0);
    }

    // transform to composite layout while adding interface values via the Export() combine mode
    RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(dofMap, true);
    regionalToComposite(interfaceScaling, compInterfaceScalingSum, rowImportPerGrp);

    /* transform composite layout back to regional layout. Now, GIDs associated
     * with region interface should carry a scaling factor (!= 1).
     */
    Array<RCP<Vector> > regVec(maxRegPerProc);
    Array<RCP<Vector> > quasiRegInterfaceScaling(maxRegPerProc);
    compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
                        interfaceScaling, revisedRowMapPerGrp, rowImportPerGrp);

    if(numRanks == 1) {
      TEST_EQUALITY(interfaceScaling[0]->getLocalLength(),    25);
      TEST_EQUALITY(interfaceScaling[0]->getGlobalLength(),   25);

      // No scaling on one rank, so all values are 1.0
      Teuchos::ArrayRCP<const SC> myScaling;
      for (int j = 0; j < maxRegPerProc; j++){
        myScaling = interfaceScaling[j]->getData(0);
        for(size_t idx = 0; idx < interfaceScaling[j]->getLocalLength(); ++idx) {
          TEST_FLOATING_EQUALITY(myScaling[idx], TST::one(), 100*TMT::eps());
        }
      }

    } else if(numRanks == 4) {
      // All ranks will have the same number of rows/cols/entries
      TEST_EQUALITY(interfaceScaling[0]->getLocalLength(),    9);
      TEST_EQUALITY(interfaceScaling[0]->getGlobalLength(),    36);

      ArrayRCP<SC> refValues;
      Teuchos::ArrayRCP<const SC> myScaling;
      if(myRank == 0) {
        refValues.deepCopy(ArrayView<const SC>({1, 1, 2, 1, 1, 2, 2, 2, 4}));

      } else if(myRank == 1) {
        refValues.deepCopy(ArrayView<const SC>({2, 1, 1, 2, 1, 1, 4, 2, 2}));

      } else if(myRank == 2) {
        refValues.deepCopy(ArrayView<const SC>({2, 2, 4, 1, 1, 2, 1, 1, 2}));

      } else if(myRank == 3) {
        refValues.deepCopy(ArrayView<const SC>({4, 2, 2, 2, 1, 1, 2, 1, 1}));
      }
      // Loop over region vector data and compare it to reference data
      for (int j = 0; j < maxRegPerProc; j++){
        myScaling = interfaceScaling[j]->getData(0);
        for(size_t idx = 0; idx < interfaceScaling[j]->getLocalLength(); ++idx) {
          TEST_FLOATING_EQUALITY(myScaling[idx], refValues[idx], 100*TMT::eps());
        }
      }
    }

  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


} // RegionCompositeVector

#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RegionVector,RegionCompositeVector,Scalar,LO,GO,Node)

#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests

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

#include <Teuchos_UnitTestHarness.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_LineDetectionFactory.hpp"
#include "MueLu_SemiCoarsenPFactory_kokkos.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SemiCoarsenPFactory_kokkos, TestSemiCoarsenP, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  typedef TestHelpers_kokkos::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Global Lexicographic";
  const std::string coupling   = "uncoupled";
  const LO numDimensions       = 3;
  const LO numSize             = 5;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(numDimensions);
  Array<GO> gNodesPerDir(numDimensions);
  for (int dim = 0; dim < numDimensions; ++dim)
    gNodesPerDir[dim] = numSize;

  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coord_type;
  typedef Xpetra::MultiVector<coord_type, LO, GO, NO> CoordMV;
  RCP<CoordMV> fineCoords =
      TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                           lNodesPerDir, meshData,
                                                                           meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("ny", gNodesPerDir[1]);
  matrixList.set("nz", gNodesPerDir[2]);
  matrixList.set("matrixType", "Laplace3D");
  matrixList.set("left boundary", "Neumann");
  matrixList.set("right boundary", "Neumann");
  matrixList.set("front boundary", "Neumann");
  matrixList.set("back boundary", "Neumann");
  matrixList.set("bottom boundary", "Dirichlet");
  matrixList.set("top boundary", "Dirichlet");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", fineCoords->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();

  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);

  fineLevel.Set("A", A);
  fineLevel.Set("Nullspace", nullSpace);
  fineLevel.Set("Coordinates", fineCoords);
  fineLevel.Set("CoarseNumZLayers", numSize);

  RCP<LineDetectionFactory> LineDetectionFact = rcp(new LineDetectionFactory());
  LineDetectionFact->SetParameter("linedetection: orientation",
                                  Teuchos::ParameterEntry(std::string("coordinates")));
  LineDetectionFact->SetParameter("linedetection: num layers",
                                  Teuchos::ParameterEntry(numSize));

  RCP<SemiCoarsenPFactory_kokkos> SemiCoarsenPFact = rcp(new SemiCoarsenPFactory_kokkos());
  SemiCoarsenPFact->SetParameter("semicoarsen: coarsen rate", Teuchos::ParameterEntry(2));
  SemiCoarsenPFact->SetFactory("LineDetection_VertLineIds", LineDetectionFact);
  SemiCoarsenPFact->SetFactory("LineDetection_Layers", LineDetectionFact);
  SemiCoarsenPFact->SetFactory("CoarseNumZLayers", LineDetectionFact);
  coarseLevel.Request("P", SemiCoarsenPFact.get());
  coarseLevel.Request("Coordinates", SemiCoarsenPFact.get());
  SemiCoarsenPFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, SemiCoarsenPFact.get());
  RCP<CoordMV> coarseCoords = coarseLevel.Get<RCP<CoordMV>>("Coordinates", SemiCoarsenPFact.get());

  coarseLevel.Release("P", SemiCoarsenPFact.get());
  coarseLevel.Release("Coordinates", SemiCoarsenPFact.get());

  // Prolongate coarse coords and compute difference
  using STS                       = Teuchos::ScalarTraits<SC>;
  const auto one                  = STS::one();
  RCP<MultiVector> coarseCoordsSC = Utilities::RealValuedToScalarMultiVector(coarseCoords);
  RCP<MultiVector> fineCoordsDiff = Utilities::RealValuedToScalarMultiVector(fineCoords);
  P->apply(*coarseCoordsSC, *fineCoordsDiff, Teuchos::NO_TRANS, one, -one);

  // check prolongation of coarse coordinates
  // in this special case, the third layer will have the correct coordinates
  using impl_SC                 = typename Kokkos::ArithTraits<SC>::val_type;
  using impl_ATS                = Kokkos::ArithTraits<impl_SC>;
  const auto fineCoordsDiffView = fineCoordsDiff->getDeviceLocalView(Xpetra::Access::ReadOnly);
  const int numNodes            = fineCoordsDiff->getLocalLength();
  const int numVectors          = fineCoordsDiff->getNumVectors();
  const auto fineMap            = fineCoordsDiff->getMap()->getLocalMap();
  const GO gStartLayer3         = 2 * numSize * numSize;
  const GO gEndLayer3           = 3 * numSize * numSize;
  using range_policy            = Kokkos::RangePolicy<LO, typename SemiCoarsenPFactory_kokkos::execution_space>;
  LO numBadCoords               = 0;
  Kokkos::parallel_reduce(
      "Checking zeros on third layer",
      range_policy(0, numNodes), KOKKOS_LAMBDA(const LO node, LO& lNumBadCoords) {
        const GO gnode = fineMap.getGlobalElement(node);
        if (gnode >= gStartLayer3 && gnode < gEndLayer3)
          for (int k = 0; k < numVectors; ++k)
            if (impl_ATS::magnitude(fineCoordsDiffView(node, k)) > 100 * impl_ATS::eps())
              lNumBadCoords += 1;
      },
      numBadCoords);
  LO gNumBadCoords;
  const auto comm = Teuchos::DefaultComm<int>::getComm();
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &numBadCoords, &gNumBadCoords);
  TEST_EQUALITY(gNumBadCoords, 0);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SemiCoarsenPFactory_kokkos, TestSemiCoarsenP, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

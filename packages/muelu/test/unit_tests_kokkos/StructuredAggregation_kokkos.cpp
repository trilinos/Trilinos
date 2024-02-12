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
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_IndexManager_kokkos.hpp"
#include "MueLu_StructuredAggregationFactory_kokkos.hpp"
#include "MueLu_NullspaceFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation_kokkos, CreateCrsGraphConstant, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef typename Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  typedef Teuchos::ScalarTraits<Scalar> TST;
  typedef TestHelpers_kokkos::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;

  out << "version: " << MueLu::Version() << std::endl;

  Level currentLevel;
  test_factory::createSingleLevelHierarchy(currentLevel);
  currentLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  const std::string coupling   = "uncoupled";
  LO numDimensions             = 3;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 6;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                           lNodesPerDir, meshData,
                                                                           meshLayout);

  // Since we are doing uncoupled aggregation we fill meshData with -1
  for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
    meshData[idx] = -1;
  }

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  currentLevel.Request("A");
  currentLevel.Set("A", A);
  currentLevel.Set("numDimensions", numDimensions);
  currentLevel.Set("Coordinates", Coordinates);
  currentLevel.Set("lNodesPerDim", lNodesPerDir);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  currentLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory_kokkos> StructuredAggFact = rcp(new StructuredAggregationFactory_kokkos());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));
  StructuredAggFact->SetParameter("aggregation: output type",
                                  Teuchos::ParameterEntry(std::string("CrsGraph")));

  currentLevel.Request("prolongatorGraph", StructuredAggFact.get());  // request prolongatorGraph
  currentLevel.Request(*StructuredAggFact);
  StructuredAggFact->Build(currentLevel);
}  // CreateCrsGraphConstant

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation_kokkos, CreateCrsGraphLinear, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef typename Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  typedef Teuchos::ScalarTraits<Scalar> TST;
  typedef TestHelpers_kokkos::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;

  out << "version: " << MueLu::Version() << std::endl;

  Level currentLevel;
  test_factory::createSingleLevelHierarchy(currentLevel);
  currentLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  const std::string coupling   = "uncoupled";
  LO numDimensions             = 3;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 6;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                           lNodesPerDir, meshData,
                                                                           meshLayout);

  // Since we are doing uncoupled aggregation we fill meshData with -1
  for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
    meshData[idx] = -1;
  }

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  currentLevel.Request("A");
  currentLevel.Set("A", A);
  currentLevel.Set("numDimensions", numDimensions);
  currentLevel.Set("Coordinates", Coordinates);
  currentLevel.Set("lNodesPerDim", lNodesPerDir);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  currentLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory_kokkos> StructuredAggFact = rcp(new StructuredAggregationFactory_kokkos());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(1));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));
  StructuredAggFact->SetParameter("aggregation: output type",
                                  Teuchos::ParameterEntry(std::string("CrsGraph")));

  currentLevel.Request("prolongatorGraph", StructuredAggFact.get());  // request prolongatorGraph
  currentLevel.Request(*StructuredAggFact);
  StructuredAggFact->Build(currentLevel);

  RCP<CrsGraph> prolongatorGraph;
  currentLevel.Get("prolongatorGraph", prolongatorGraph, StructuredAggFact.get());

}  // CreateCrsGraphLinear

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation_kokkos, UncoupledTentative3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }

  using real_type             = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using RealValuedMultiVector = typename Xpetra::MultiVector<real_type, LO, GO, NO>;
  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using test_factory          = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  const std::string coupling   = "uncoupled";
  LO numDimensions             = 3;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 6;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                           lNodesPerDir, meshData,
                                                                           meshLayout);

  // Since we are doing uncoupled aggregation we fill meshData with -1
  for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
    meshData[idx] = -1;
  }

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory_kokkos> StructuredAggFact = rcp(new StructuredAggregationFactory_kokkos());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory_kokkos> TentativePFact = rcp(new TentativePFactory_kokkos());
  TentativePFact->SetFactory("Aggregates", StructuredAggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());
  coarseLevel.Request(*TentativePFact);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  RCP<MultiVector> coarseNullSpace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  coarseLevel.Release("P", TentativePFact.get());  // release Ptent
  coarseLevel.Release("Nullspace", TentativePFact.get());

  // check normalization and orthogonality of prolongator columns
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, Node> > PtentTPtent = Xpetra::MatrixMatrix<SC, LO, GO, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<SC, LO, GO, Node> > diagVec     = Xpetra::VectorFactory<SC, LO, GO, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 1000 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 1000 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // UncoupledTentative3D

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation_kokkos, CreateCrsGraphConstant, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation_kokkos, CreateCrsGraphLinear, SC, LO, GO, NO)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation_kokkos, UncoupledTentative3D, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

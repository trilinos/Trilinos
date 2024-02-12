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

#include "MueLu_StructuredAggregationFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_LocalLexicographicIndexManager.hpp"
#include "MueLu_GlobalLexicographicIndexManager.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, GlobalLexiTentative1D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Global Lexicographic";
  LO numDimensions             = 1;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 20;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace1D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // GlobalLexiTentative1D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, GlobalLexiTentative2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Global Lexicographic";
  LO numDimensions             = 2;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      gNodesPerDir[dim] = 12;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // GlobalLexiTentative2D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, GlobalLexiTentative3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Global Lexicographic";
  LO numDimensions             = 3;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      gNodesPerDir[dim] = 6;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // GlobalLexiTentative3D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, LocalLexiTentative1D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  LO numDimensions             = 1;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 20;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace1D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // LocalLexiTentative1D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, LocalLexiTentative2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  LO numDimensions             = 2;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 12;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // LocalLexiTentative2D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, LocalLexiTentative3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
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
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace1D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // LocalLexiTentative3D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, UncoupledLocalLexiTentative1D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  const std::string coupling   = "uncoupled";
  LO numDimensions             = 1;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 20;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
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
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace1D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: mode",
                                  Teuchos::ParameterEntry(coupling));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // UncoupledLocalLexiTentative1D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, UncoupledLocalLexiTentative2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  const std::string coupling   = "uncoupled";
  LO numDimensions             = 2;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 12;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
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
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: mode",
                                  Teuchos::ParameterEntry(coupling));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // UncoupledLocalLexiTentative2D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, UncoupledLocalLexiTentative3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

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
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
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
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: mode",
                                  Teuchos::ParameterEntry(coupling));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // UncoupledLocalLexiTentative3D

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, UncoupledMultilevelScalar, Scalar,
                                  LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using magnitude_type        = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using real_type             = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

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
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  // Since we are doing uncoupled aggregation we fill meshData with -1
  for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
    meshData[idx] = -1;
  }

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("ny", gNodesPerDir[1]);
  matrixList.set("nz", gNodesPerDir[2]);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), 1);
  nullSpace->putScalar((Scalar)1.0);
  Teuchos::Array<magnitude_type> norms(1);
  nullSpace->norm1(norms);
  if (Coordinates->getMap()->getComm()->getRank() == 0) {
    out << "||NS|| = " << norms[0] << std::endl;
  }

  // fill hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  // create the factory manager and the factories
  FactoryManager M;

  RCP<Factory> AmalgFact     = rcp(new AmalgamationFactory());
  RCP<Factory> CDropfact     = rcp(new CoalesceDropFactory());
  RCP<Factory> Aggfact       = rcp(new StructuredAggregationFactory());
  RCP<Factory> coarseMapFact = rcp(new CoarseMapFactory());
  RCP<Factory> Pfact         = rcp(new TentativePFactory());
  RCP<Factory> Rfact         = rcp(new TransPFactory());
  RCP<Factory> Tfact         = rcp(new CoordinatesTransferFactory());
  RCP<RAPFactory> Acfact     = rcp(new RAPFactory());
  RCP<Factory> NSfact        = rcp(new NullspaceFactory());

  // Set paramters needed by the factories
  CDropfact->SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  CDropfact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.0));

  Aggfact->SetParameter("aggregation: mode", Teuchos::ParameterEntry(std::string("uncoupled")));

  Tfact->SetParameter("structured aggregation", Teuchos::ParameterEntry(true));

  // Set interfactory dependencies
  CDropfact->SetFactory("UnAmalgamationInfo", AmalgFact);
  Aggfact->SetFactory("Graph", CDropfact);
  Aggfact->SetFactory("DofsPerNode", CDropfact);
  coarseMapFact->SetFactory("Aggregates", Aggfact);
  Pfact->SetFactory("UnAmalgamationInfo", AmalgFact);
  Pfact->SetFactory("Aggregates", Aggfact);
  Pfact->SetFactory("CoarseMap", coarseMapFact);
  NSfact->SetFactory("Nullspace", Pfact);

  Acfact->AddTransferFactory(Tfact);

  // Set default factories in the manager
  M.SetFactory("P", Pfact);
  M.SetFactory("R", Rfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Nullspace", NSfact);
  M.SetFactory("Graph", CDropfact);
  M.SetFactory("lNodesPerDim", Tfact);
  M.SetFactory("numDimensions", Tfact);
  M.SetFactory("lCoarseNodesPerDim", Aggfact);

  // setup smoothers
  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LocalOrdinal)1);
  smootherParamList.set("relaxation: damping factor", (Scalar)1.0);
  RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
  RCP<SmootherFactory> SmooFact    = rcp(new SmootherFactory(smooProto));
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

  // Set smoothers and coarse solver in the manager
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarseSolveFact);

  // Populate data on the finest level of the hierarchy
  RCP<Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A", A);                  // set fine level matrix
  Finest->Set("Nullspace", nullSpace);  // set null space information for finest level
  // Finest->Set("Coordinates", Coordinates);    // set fine level coordinates
  Finest->Set("numDimensions", numDimensions);  // set GeneralGeometricPFactory specific info
  Finest->Set("lNodesPerDim", lNodesPerDir);    // set GeneralGeometricPFactory specific info

  // Setup the hierarchy
  LO maxLevels = 5;
  H->SetMaxCoarseSize(10);
  H->Setup(M, 0, maxLevels);

  // // Extract the prolongator operator
  // RCP<Level> lvl1 = H->GetLevel(1);
  // RCP<Xpetra::Matrix<SC,LO,GO,Node> > P = lvl1->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
  // RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // // Construct vectors to check that a linear vector remains linear after projection
  // RCP<Xpetra::MultiVector<SC,LO,GO,NO> > vector0
  //   = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(PCrs->getRangeMap(),1);
  // RCP<Xpetra::MultiVector<SC,LO,GO,NO> > vector1
  //   = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(PCrs->getDomainMap(),1);
  // ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
  // if(comm->getRank() == 0) {
  //   for(LO i = 0; i < 3; ++i) {
  //     for(LO j = 0; j < 3; ++j) {
  //       coarse_data[3*i + j] = 2.0*i + j;
  //     }
  //   }
  // }

  // PCrs->apply(*vector1, *vector0, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
  //             Teuchos::ScalarTraits<SC>::zero());

  // ArrayRCP<const SC> fine_data   = vector0->getData(0);
  // bool is_constant = true, is_injected = true;
  // Array<LO> fine_inds(9);
  // Array<LO> coarse_inds(9);
  // if(comm->getRank() == 0) {
  //   fine_inds[0] = 0;
  //   fine_inds[1] = 2;
  //   fine_inds[2] = 3;
  //   fine_inds[3] = 8;
  //   fine_inds[4] = 10;
  //   fine_inds[5] = 11;
  //   fine_inds[6] = 12;
  //   fine_inds[7] = 14;
  //   fine_inds[8] = 15;
  //   for(LO i = 0; i < 9; ++i) {
  //     if(std::abs(fine_data[fine_inds[i]] - coarse_data[i]) > 1.0e-10) { is_injected = false; }
  //   }
  //   if( (std::abs(fine_data[1] - coarse_data[0]) > 1.0e-10) ||
  //       (std::abs(fine_data[4] - coarse_data[0]) > 1.0e-10) ||
  //       (std::abs(fine_data[5] - coarse_data[0]) > 1.0e-10)) { is_constant = false; }
  // }

  TEST_EQUALITY(A != Teuchos::null, true);
  TEST_EQUALITY(H != Teuchos::null, true);
  TEST_EQUALITY(nullSpace != Teuchos::null, true);
  // TEST_EQUALITY(is_injected, true);
  // TEST_EQUALITY(is_constant, true);

}  // UncoupledMultilevelScalar

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, ProlongatorGraphUncoupled, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, Node>;

  out << "version: " << MueLu::Version() << std::endl;

  for (int interpolationOrder = 0; interpolationOrder < 2; ++interpolationOrder) {
    out << "Tesing 7x7 grid with piece-wise "
        << (interpolationOrder == 0 ? "constant" : "linear") << " interpolation" << std::endl;
    Level fineLevel;
    test_factory::createSingleLevelHierarchy(fineLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test

    // Set global geometric data
    const std::string meshLayout = "Local Lexicographic";
    const std::string coupling   = "uncoupled";
    LO numDimensions             = 2;
    Array<GO> meshData;
    Array<LO> lNodesPerDir(3);
    Array<GO> gNodesPerDir = {{7, 7, 1}};

    RCP<RealValuedMultiVector> Coordinates =
        TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                      lNodesPerDir, meshData,
                                                                      meshLayout);

    // Since we are doing uncoupled aggregation we fill meshData with -1
    for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
      meshData[idx] = -1;
    }

    Teuchos::ParameterList matrixList;
    matrixList.set("nx", gNodesPerDir[0]);
    matrixList.set("ny", gNodesPerDir[1]);
    matrixList.set("matrixType", "Laplace2D");
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
        BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", Coordinates->getMap(),
                                                                  matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();
    fineLevel.Request("A");
    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates", Coordinates);
    fineLevel.Set("numDimensions", numDimensions);
    fineLevel.Set("gNodesPerDim", gNodesPerDir);
    fineLevel.Set("lNodesPerDim", lNodesPerDir);
    fineLevel.Set("aggregation: mesh data", meshData);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
    RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
    StructuredAggFact->SetFactory("Graph", dropFact);
    StructuredAggFact->SetFactory("DofsPerNode", dropFact);
    StructuredAggFact->SetParameter("aggregation: mesh layout",
                                    Teuchos::ParameterEntry(meshLayout));
    StructuredAggFact->SetParameter("aggregation: mode",
                                    Teuchos::ParameterEntry(coupling));
    StructuredAggFact->SetParameter("aggregation: output type",
                                    Teuchos::ParameterEntry(std::string("Graph")));
    StructuredAggFact->SetParameter("aggregation: coarsening order",
                                    Teuchos::ParameterEntry(interpolationOrder));
    StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                    Teuchos::ParameterEntry(std::string("{3}")));

    fineLevel.Request("prolongatorGraph", StructuredAggFact.get());
    fineLevel.Request(*StructuredAggFact);
    StructuredAggFact->Build(fineLevel);

    RCP<CrsGraph> prolongatorGraph;
    fineLevel.Get("prolongatorGraph", prolongatorGraph, StructuredAggFact.get());

    TEST_EQUALITY(prolongatorGraph != Teuchos::null, true);

    int numErrors      = 0;
    const int numRanks = Coordinates->getMap()->getComm()->getSize();
    TEST_EQUALITY(prolongatorGraph->getGlobalNumRows() == 49, true);
    if (interpolationOrder == 0) {
      TEST_EQUALITY(prolongatorGraph->getGlobalNumEntries() == 49, true);
    } else {
      TEST_EQUALITY(prolongatorGraph->getGlobalNumEntries() == (numRanks == 1 ? 169 : 148), true);
    }
    if (numRanks == 1) {
      TEST_EQUALITY(prolongatorGraph->getGlobalNumCols() == 9, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumRows() == 49, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumCols() == 9, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumEntries() == (interpolationOrder == 1 ? 169 : 49), true);

      ArrayView<const LO> rowIndices;
      const LO indicesConstant[49] = {0, 0, 1, 1, 1, 2, 2,
                                      0, 0, 1, 1, 1, 2, 2,
                                      3, 3, 4, 4, 4, 5, 5,
                                      3, 3, 4, 4, 4, 5, 5,
                                      3, 3, 4, 4, 4, 5, 5,
                                      6, 6, 7, 7, 7, 8, 8,
                                      6, 6, 7, 7, 7, 8, 8};
      const LO indicesLinear[169]  = {0, 0, 1, 3, 4, 0, 1, 3, 4, 1, 1, 2, 4, 5, 1, 2, 4, 5, 2,
                                      0, 1, 3, 4, 0, 1, 3, 4, 0, 1, 3, 4, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5,
                                      0, 1, 3, 4, 0, 1, 3, 4, 0, 1, 3, 4, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5,
                                      3, 3, 4, 6, 7, 3, 4, 6, 7, 4, 4, 5, 7, 8, 4, 5, 7, 8, 5,
                                      3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 4, 5, 7, 8, 4, 5, 7, 8, 4, 5, 7, 8, 4, 5, 7, 8,
                                      3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 4, 5, 7, 8, 4, 5, 7, 8, 4, 5, 7, 8, 4, 5, 7, 8,
                                      6, 3, 4, 6, 7, 3, 4, 6, 7, 7, 4, 5, 7, 8, 4, 5, 7, 8, 8};
      const LO rowOffset[50]       = {0, 1, 5, 9, 10, 14, 18, 19,
                                      23, 27, 31, 35, 39, 43, 47,
                                      51, 55, 59, 63, 67, 71, 75,
                                      76, 80, 84, 85, 89, 93, 94,
                                      98, 102, 106, 110, 114, 118, 122,
                                      126, 130, 134, 138, 142, 146, 150,
                                      151, 155, 159, 160, 164, 168, 169};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);

        if (interpolationOrder == 0) {
          if (rowIndices[0] != indicesConstant[rowIdx]) {
            ++numErrors;
          }
        } else {
          for (size_t entryIdx = 0; entryIdx < prolongatorGraph->getNumEntriesInLocalRow(rowIdx); ++entryIdx) {
            if (rowIndices[entryIdx] != indicesLinear[rowOffset[rowIdx] + entryIdx]) {
              ++numErrors;
            }
          }
        }
      }
    } else if (numRanks == 4) {
      const int myRank             = Coordinates->getMap()->getComm()->getRank();
      Array<size_t> nodeNumRows    = {{16, 12, 12, 9}};
      Array<int> rankOffsets       = {{0, 16, 28, 40}};
      Array<size_t> nodeNumEntries = {{52, 36, 36, 24}};

      TEST_EQUALITY(prolongatorGraph->getGlobalNumCols() == 16, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumRows() == nodeNumRows[myRank], true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumCols() == 4, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumEntries() ==
                        ((interpolationOrder == 0) ? nodeNumRows[myRank] : nodeNumEntries[myRank]),
                    true);

      ArrayView<const LO> rowIndices;
      const LO indicesConstant[49] = {0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3, 2, 2, 3, 3,
                                      0, 0, 1, 0, 0, 1, 2, 2, 3, 2, 2, 3,
                                      0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3,
                                      0, 0, 1, 0, 0, 1, 2, 2, 3};
      const LO indicesLinear[148]  = {0, 0, 1, 2, 3, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 0, 1, 2, 3, 3,
                                      0, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 3,
                                      0, 0, 1, 2, 3, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 0, 1, 2, 3, 3,
                                      0, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 3};
      const LO rowOffset[50]       = {0, 1, 5, 9, 10, 14, 18, 22, 26, 30, 34, 38, 42, 43, 47, 51, 52,
                                      53, 57, 58, 62, 66, 70, 74, 78, 82, 83, 87, 88,
                                      89, 93, 97, 98, 102, 106, 110, 114, 115, 119, 123, 124,
                                      125, 129, 130, 134, 138, 142, 143, 147, 148};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);

        if (interpolationOrder == 0) {
          if (rowIndices[0] != indicesConstant[rankOffsets[myRank] + rowIdx]) {
            ++numErrors;
          }
        } else {
          for (size_t entryIdx = 0; entryIdx < prolongatorGraph->getNumEntriesInLocalRow(rowIdx); ++entryIdx) {
            if (rowIndices[entryIdx] != indicesLinear[rowOffset[rankOffsets[myRank] + rowIdx] + entryIdx]) {
              ++numErrors;
            }
          }
        }
      }
    }
    TEST_EQUALITY(numErrors == 0, true);
  }  // Loop over interpolationOrder

  for (int interpolationOrder = 0; interpolationOrder < 2; ++interpolationOrder) {
    out << "Tesing 7x6 grid with piece-wise "
        << (interpolationOrder == 0 ? "constant" : "linear") << " interpolation" << std::endl;
    Level fineLevel;
    test_factory::createSingleLevelHierarchy(fineLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test

    // Set global geometric data
    const std::string meshLayout = "Local Lexicographic";
    const std::string coupling   = "uncoupled";
    LO numDimensions             = 2;
    Array<GO> meshData;
    Array<LO> lNodesPerDir(3);
    Array<GO> gNodesPerDir = {{7, 6, 1}};

    RCP<RealValuedMultiVector> Coordinates =
        TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                      lNodesPerDir, meshData,
                                                                      meshLayout);

    // Since we are doing uncoupled aggregation we fill meshData with -1
    for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
      meshData[idx] = -1;
    }

    Teuchos::ParameterList matrixList;
    matrixList.set("nx", gNodesPerDir[0]);
    matrixList.set("ny", gNodesPerDir[1]);
    matrixList.set("matrixType", "Laplace2D");
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
        BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", Coordinates->getMap(),
                                                                  matrixList);
    RCP<Matrix> A = Pr->BuildMatrix();
    fineLevel.Request("A");
    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates", Coordinates);
    fineLevel.Set("numDimensions", numDimensions);
    fineLevel.Set("gNodesPerDim", gNodesPerDir);
    fineLevel.Set("lNodesPerDim", lNodesPerDir);
    fineLevel.Set("aggregation: mesh data", meshData);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
    RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
    StructuredAggFact->SetFactory("Graph", dropFact);
    StructuredAggFact->SetFactory("DofsPerNode", dropFact);
    StructuredAggFact->SetParameter("aggregation: mesh layout",
                                    Teuchos::ParameterEntry(meshLayout));
    StructuredAggFact->SetParameter("aggregation: mode",
                                    Teuchos::ParameterEntry(coupling));
    StructuredAggFact->SetParameter("aggregation: output type",
                                    Teuchos::ParameterEntry(std::string("Graph")));
    StructuredAggFact->SetParameter("aggregation: coarsening order",
                                    Teuchos::ParameterEntry(interpolationOrder));
    StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                    Teuchos::ParameterEntry(std::string("{3}")));

    fineLevel.Request("prolongatorGraph", StructuredAggFact.get());
    fineLevel.Request(*StructuredAggFact);
    StructuredAggFact->Build(fineLevel);

    RCP<CrsGraph> prolongatorGraph;
    fineLevel.Get("prolongatorGraph", prolongatorGraph, StructuredAggFact.get());

    TEST_EQUALITY(prolongatorGraph != Teuchos::null, true);

    int numErrors      = 0;
    const int numRanks = Coordinates->getMap()->getComm()->getSize();
    TEST_EQUALITY(prolongatorGraph->getGlobalNumRows() == 42, true);
    if (interpolationOrder == 0) {
      TEST_EQUALITY(prolongatorGraph->getGlobalNumEntries() == 42, true);
    } else {
      TEST_EQUALITY(prolongatorGraph->getGlobalNumEntries() == (numRanks == 1 ? 141 : 120), true);
    }
    if (numRanks == 1) {
      TEST_EQUALITY(prolongatorGraph->getGlobalNumCols() == 9, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumRows() == 42, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumCols() == 9, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumEntries() == (interpolationOrder == 1 ? 141 : 42), true);

      ArrayView<const LO> rowIndices;
      const LO indicesConstant[42] = {0, 0, 1, 1, 1, 2, 2,
                                      0, 0, 1, 1, 1, 2, 2,
                                      3, 3, 4, 4, 4, 5, 5,
                                      3, 3, 4, 4, 4, 5, 5,
                                      3, 3, 4, 4, 4, 5, 5,
                                      6, 6, 7, 7, 7, 8, 8};
      const LO indicesLinear[141]  = {0, 0, 1, 3, 4, 0, 1, 3, 4, 1, 1, 2, 4, 5, 1, 2, 4, 5, 2,
                                      0, 1, 3, 4, 0, 1, 3, 4, 0, 1, 3, 4, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5,
                                      0, 1, 3, 4, 0, 1, 3, 4, 0, 1, 3, 4, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5, 1, 2, 4, 5,
                                      3, 3, 4, 6, 7, 3, 4, 6, 7, 4, 4, 5, 7, 8, 4, 5, 7, 8, 5,
                                      3, 4, 6, 7, 3, 4, 6, 7, 3, 4, 6, 7, 4, 5, 7, 8, 4, 5, 7, 8, 4, 5, 7, 8, 4, 5, 7, 8,
                                      6, 3, 4, 6, 7, 3, 4, 6, 7, 7, 4, 5, 7, 8, 4, 5, 7, 8, 8};
      const LO rowOffset[43]       = {0, 1, 5, 9, 10, 14, 18, 19,
                                      23, 27, 31, 35, 39, 43, 47,
                                      51, 55, 59, 63, 67, 71, 75,
                                      76, 80, 84, 85, 89, 93, 94,
                                      98, 102, 106, 110, 114, 118, 122,
                                      123, 127, 131, 132, 136, 140, 141};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);

        if (interpolationOrder == 0) {
          if (rowIndices[0] != indicesConstant[rowIdx]) {
            ++numErrors;
          }
        } else {
          for (size_t entryIdx = 0; entryIdx < prolongatorGraph->getNumEntriesInLocalRow(rowIdx); ++entryIdx) {
            if (rowIndices[entryIdx] != indicesLinear[rowOffset[rowIdx] + entryIdx]) {
              ++numErrors;
            }
          }
        }
      }
    } else if (numRanks == 4) {
      const int myRank             = Coordinates->getMap()->getComm()->getRank();
      Array<size_t> nodeNumRows    = {{12, 9, 12, 9}};
      Array<int> rankOffsets       = {{0, 12, 21, 33}};
      Array<size_t> nodeNumEntries = {{36, 24, 36, 24}};

      TEST_EQUALITY(prolongatorGraph->getGlobalNumCols() == 16, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumRows() == nodeNumRows[myRank], true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumCols() == 4, true);
      TEST_EQUALITY(prolongatorGraph->getLocalNumEntries() ==
                        ((interpolationOrder == 0) ? nodeNumRows[myRank] : nodeNumEntries[myRank]),
                    true);

      ArrayView<const LO> rowIndices;
      const LO indicesConstant[42] = {0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3,
                                      0, 0, 1, 0, 0, 1, 2, 2, 3,
                                      0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3,
                                      0, 0, 1, 0, 0, 1, 2, 2, 3};
      const LO indicesLinear[120]  = {0, 0, 1, 2, 3, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 0, 1, 2, 3, 3,
                                      0, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 3,
                                      0, 0, 1, 2, 3, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 0, 1, 2, 3, 3,
                                      0, 0, 1, 2, 3, 1,
                                      0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                                      2, 0, 1, 2, 3, 3};
      const LO rowOffset[43]       = {0, 1, 5, 9, 10, 14, 18, 22, 26, 27, 31, 35, 36,
                                      37, 41, 42, 46, 50, 54, 55, 59, 60,
                                      61, 65, 69, 70, 74, 78, 82, 86, 87, 91, 95, 96,
                                      97, 101, 102, 106, 110, 114, 115, 119, 120};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);

        if (interpolationOrder == 0) {
          if (rowIndices[0] != indicesConstant[rankOffsets[myRank] + rowIdx]) {
            ++numErrors;
          }
        } else {
          for (size_t entryIdx = 0; entryIdx < prolongatorGraph->getNumEntriesInLocalRow(rowIdx); ++entryIdx) {
            if (rowIndices[entryIdx] != indicesLinear[rowOffset[rankOffsets[myRank] + rowIdx] + entryIdx]) {
              ++numErrors;
            }
          }
        }
      }
    }
    TEST_EQUALITY(numErrors == 0, true);
  }  // Loop over interpolationOrder

}  // ProlongatorGraphUncoupled

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, UncoupledAggSingleCoarseNode, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  const std::string coupling   = "uncoupled";
  LO numDimensions             = 2;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 5;
    } else {
      gNodesPerDir[dim] = 1;
    }
  }

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  // Since we are doing uncoupled aggregation we fill meshData with -1
  for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
    meshData[idx] = -1;
  }

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("matrixType", "Laplace2D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  // only one NS vector -> exercises manual orthogonalization
  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: mode",
                                  Teuchos::ParameterEntry(coupling));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(0));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{2}")));
  StructuredAggFact->SetParameter("aggregation: single coarse point",
                                  Teuchos::ParameterEntry(true));

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", StructuredAggFact);
  RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
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
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Ptent, true, *Ptent, false, out);
  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > diagVec     = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (TST::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), Teuchos::as<magnitude_type>(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

}  // UncoupledAggSingleCoarseNode

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredAggregation, UncoupledGraphSingleCoarseNode, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, Node>;

  out << "version: " << MueLu::Version() << std::endl;
  int interpolationOrder = 0;
  out << "Tesing 7x6 grid with piece-wise "
      << (interpolationOrder == 0 ? "constant" : "linear") << " interpolation" << std::endl;
  Level fineLevel;
  test_factory::createSingleLevelHierarchy(fineLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test

  // Set global geometric data
  const std::string meshLayout = "Local Lexicographic";
  const std::string coupling   = "uncoupled";
  LO numDimensions             = 2;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(3);
  Array<GO> gNodesPerDir = {{7, 7, 1}};

  RCP<RealValuedMultiVector> Coordinates =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  // Since we are doing uncoupled aggregation we fill meshData with -1
  for (size_t idx = 0; idx < (size_t)meshData.size(); ++idx) {
    meshData[idx] = -1;
  }

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("ny", gNodesPerDir[1]);
  matrixList.set("matrixType", "Laplace2D");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace2D", Coordinates->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();
  fineLevel.Request("A");
  fineLevel.Set("A", A);
  fineLevel.Set("Coordinates", Coordinates);
  fineLevel.Set("numDimensions", numDimensions);
  fineLevel.Set("gNodesPerDim", gNodesPerDir);
  fineLevel.Set("lNodesPerDim", lNodesPerDir);
  fineLevel.Set("aggregation: mesh data", meshData);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  RCP<StructuredAggregationFactory> StructuredAggFact = rcp(new StructuredAggregationFactory());
  StructuredAggFact->SetFactory("Graph", dropFact);
  StructuredAggFact->SetFactory("DofsPerNode", dropFact);
  StructuredAggFact->SetParameter("aggregation: mesh layout",
                                  Teuchos::ParameterEntry(meshLayout));
  StructuredAggFact->SetParameter("aggregation: mode",
                                  Teuchos::ParameterEntry(coupling));
  StructuredAggFact->SetParameter("aggregation: output type",
                                  Teuchos::ParameterEntry(std::string("Graph")));
  StructuredAggFact->SetParameter("aggregation: coarsening order",
                                  Teuchos::ParameterEntry(interpolationOrder));
  StructuredAggFact->SetParameter("aggregation: coarsening rate",
                                  Teuchos::ParameterEntry(std::string("{3}")));
  StructuredAggFact->SetParameter("aggregation: single coarse point",
                                  Teuchos::ParameterEntry(true));

  fineLevel.Request("prolongatorGraph", StructuredAggFact.get());
  fineLevel.Request(*StructuredAggFact);
  StructuredAggFact->Build(fineLevel);

  RCP<CrsGraph> prolongatorGraph;
  fineLevel.Get("prolongatorGraph", prolongatorGraph, StructuredAggFact.get());

  TEST_EQUALITY(prolongatorGraph != Teuchos::null, true);

  int numErrors      = 0;
  const int numRanks = Coordinates->getMap()->getComm()->getSize();
  TEST_EQUALITY(prolongatorGraph->getGlobalNumRows() == 49, true);
  TEST_EQUALITY(prolongatorGraph->getGlobalNumEntries() == 49, true);
  if (numRanks == 1) {
    TEST_EQUALITY(prolongatorGraph->getGlobalNumCols() == 9, true);
    TEST_EQUALITY(prolongatorGraph->getLocalNumRows() == 49, true);
    TEST_EQUALITY(prolongatorGraph->getLocalNumCols() == 9, true);
    TEST_EQUALITY(prolongatorGraph->getLocalNumEntries() == 49, true);

    ArrayView<const LO> rowIndices;
    const LO indicesConstant[49] = {0, 0, 1, 1, 1, 2, 2,
                                    0, 0, 1, 1, 1, 2, 2,
                                    3, 3, 4, 4, 4, 5, 5,
                                    3, 3, 4, 4, 4, 5, 5,
                                    3, 3, 4, 4, 4, 5, 5,
                                    6, 6, 7, 7, 7, 8, 8,
                                    6, 6, 7, 7, 7, 8, 8};
    for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
      prolongatorGraph->getLocalRowView(rowIdx, rowIndices);
      if (rowIndices[0] != indicesConstant[rowIdx]) {
        ++numErrors;
      }
    }
  } else if (numRanks == 4) {
    const int myRank                = Coordinates->getMap()->getComm()->getRank();
    const Array<size_t> nodeNumRows = {{16, 12, 12, 9}};
    const Array<size_t> nodeNumCols = {{4, 2, 2, 1}};

    TEST_EQUALITY(prolongatorGraph->getGlobalNumCols() == 9, true);
    TEST_EQUALITY(prolongatorGraph->getLocalNumRows() == nodeNumRows[myRank], true);
    TEST_EQUALITY(prolongatorGraph->getLocalNumCols() == nodeNumCols[myRank], true);
    TEST_EQUALITY(prolongatorGraph->getLocalNumEntries() == nodeNumRows[myRank], true);

    ArrayView<const LO> rowIndices;
    if (myRank == 0) {
      const LO indicesConstant[16] = {0, 0, 1, 1,
                                      0, 0, 1, 1,
                                      2, 2, 3, 3,
                                      2, 2, 3, 3};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);
        if (rowIndices[0] != indicesConstant[rowIdx]) {
          ++numErrors;
        }
      }
    } else if (myRank == 1) {
      const LO indicesConstant[12] = {0, 0, 0,
                                      0, 0, 0,
                                      1, 1, 1,
                                      1, 1, 1};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);
        if (rowIndices[0] != indicesConstant[rowIdx]) {
          ++numErrors;
        }
      }
    } else if (myRank == 2) {
      const LO indicesConstant[12] = {0, 0, 1, 1,
                                      0, 0, 1, 1,
                                      0, 0, 1, 1};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);
        if (rowIndices[0] != indicesConstant[rowIdx]) {
          ++numErrors;
        }
      }
    } else if (myRank == 3) {
      const LO indicesConstant[9] = {0, 0, 0,
                                     0, 0, 0,
                                     0, 0, 0};
      for (size_t rowIdx = 0; rowIdx < prolongatorGraph->getLocalNumRows(); ++rowIdx) {
        prolongatorGraph->getLocalRowView(rowIdx, rowIndices);
        if (rowIndices[0] != indicesConstant[rowIdx]) {
          ++numErrors;
        }
      }
    }
  }
  TEST_EQUALITY(numErrors == 0, true);

}  // UncoupledGraphSingleCoarseNode

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, GlobalLexiTentative1D, Scalar, LO, GO, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, GlobalLexiTentative2D, Scalar, LO, GO, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, GlobalLexiTentative3D, Scalar, LO, GO, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, LocalLexiTentative1D, Scalar, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, LocalLexiTentative2D, Scalar, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, LocalLexiTentative3D, Scalar, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, UncoupledLocalLexiTentative1D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, UncoupledLocalLexiTentative2D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, UncoupledLocalLexiTentative3D, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, UncoupledMultilevelScalar, Scalar, LO, GO, Node)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, ProlongatorGraphUncoupled, Scalar, LO, GO, Node)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, UncoupledAggSingleCoarseNode, Scalar, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredAggregation, UncoupledGraphSingleCoarseNode, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

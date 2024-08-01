// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory_kokkos, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }
  out << "version: " << MueLu::Version() << std::endl;

  RCP<TentativePFactory_kokkos> tentPFact = rcp(new TentativePFactory_kokkos);
  TEST_EQUALITY(tentPFact != Teuchos::null, true);
}

#if 0
  TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentative_LapackQR)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    out << "Test QR with user-supplied nullspace" << std::endl;

    using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
    using TMT            = Teuchos::ScalarTraits<magnitude_type>;

    Level fineLevel, coarseLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);

    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(/*199*/29);
    A->SetFixedBlockSize(1);
    fineLevel.Set("A",A);

    LO NSdim = 2;
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    nullSpace->randomize();
    fineLevel.Set("Nullspace",nullSpace);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
    UncoupledAggFact->SetFactory("Graph", dropFact);

    UncoupledAggFact->SetMinNodesPerAggregate(3);
    UncoupledAggFact->SetMaxNeighAlreadySelected(0);
    UncoupledAggFact->SetOrdering("natural");

    RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
    coarseMapFact->SetFactory("Aggregates", UncoupledAggFact);

    RCP<TentativePFactory> TentativePFact = rcp(new TentativePFactory());
    TentativePFact->SetFactory("Aggregates", UncoupledAggFact);
    TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
    TentativePFact->SetFactory("CoarseMap", coarseMapFact);

    coarseLevel.Request("P",TentativePFact.get());         // request Ptent
    coarseLevel.Request("Nullspace",TentativePFact.get()); // request coarse nullspace
    coarseLevel.Request(*TentativePFact);
    TentativePFact->Build(fineLevel,coarseLevel);

    RCP<Matrix> Ptent;
    coarseLevel.Get("P",Ptent,TentativePFact.get());

    RCP<MultiVector> coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

    //check interpolation
    RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(),NSdim);
    Ptent->apply(*coarseNullspace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

    RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(),NSdim);
    diff->putScalar(0.0);

    coarseLevel.Release("P",TentativePFact.get()); // release Ptent
    coarseLevel.Release("Nullspace",TentativePFact.get());   // release coarse nullspace

    //diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
    diff->update(1.0,*nullSpace,-1.0,*PtN,0.0);

    Teuchos::Array<magnitude_type> norms(NSdim);
    diff->norm2(norms);
    for (LO i=0; i<NSdim; ++i) {
      out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
      TEST_COMPARE_CONST(norms[i], <, 100*TMT::eps());
    }

    Teuchos::ArrayRCP<const double> col1 = coarseNullspace->getData(0);
    Teuchos::ArrayRCP<const double> col2 = coarseNullspace->getData(1);
    TEST_EQUALITY(col1.size() == col2.size(), true);

  } // MakeTentative  Lapack QR
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory_kokkos, MakeTentative, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO);

  using STS            = Teuchos::ScalarTraits<Scalar>;
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test QR with user-supplied nullspace" << std::endl;

  Level fineLevel, coarseLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(199);
  fineLevel.Request("A");
  fineLevel.Set("A", A);

  // Only one NS vector -> exercises manual orthogonalization
  LO NSdim                   = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());

  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
  ParameterList aggParams;
  aggParams.set("aggregation: ordering", "natural");
  aggParams.set("aggregation: min agg size", 3);
  aggParams.set("aggregation: max selected neighbors", 0);
  aggFact->SetParameterList(aggParams);
  aggFact->SetFactory("DofsPerNode", dropFact);
  aggFact->SetFactory("Graph", dropFact);

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", aggFact);

  RCP<TentativePFactory_kokkos> TentativePFact = rcp(new TentativePFactory_kokkos());
  TentativePFact->SetFactory("Aggregates", aggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());
  coarseLevel.Request(*TentativePFact);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  RCP<MultiVector> coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  // Check interpolation by computing ||fineNS - P*coarseNS||
  RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);
  Ptent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  coarseLevel.Release("P", TentativePFact.get());  // release Ptent
  coarseLevel.Release("Nullspace", TentativePFact.get());

  // diff = fineNS - (P*coarseNS)
  diff->update(1.0, *nullSpace, -1.0, *PtN, 0.0);

  Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(NSdim);
  diff->norm2(norms);
  for (LO i = 0; i < NSdim; ++i) {
    out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
    TEST_COMPARE_CONST(norms[i], <, 100 * TMT::eps());
  }

  // Check normalization and orthogonality of prolongator columns by
  // making sure that P^T*P = I
  RCP<Matrix> PtentTPtent = MatrixMatrix::Multiply(*Ptent, true, *Ptent, false, out);
  RCP<Vector> diagVec     = VectorFactory::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (STS::name().find("complex") == std::string::npos) {  // skip check for Scalar=complex
    TEST_EQUALITY(diagVec->norm1(), diagVec->getGlobalLength());
    TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  }
  if (STS::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(STS::magnitude(diagVec->meanValue()), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());
}  // MakeTentative

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory_kokkos, MakeTentativeVectorBasedUsingDefaultNullspace, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO);

  using STS            = Teuchos::ScalarTraits<Scalar>;
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test QR when nullspace isn't supplied by user" << std::endl;

  Level fineLevel, coarseLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  auto A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(200);

  A->SetFixedBlockSize(2);

  fineLevel.Set("A", A);

  auto amalgFact = rcp(new AmalgamationFactory());

  auto dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  auto aggFact = rcp(new UncoupledAggregationFactory());
  ParameterList aggParams;
  aggParams.set("aggregation: ordering", "natural");
  aggParams.set("aggregation: min agg size", 3);
  aggParams.set("aggregation: max selected neighbors", 0);
  aggFact->SetParameterList(aggParams);
  aggFact->SetFactory("DofsPerNode", dropFact);
  aggFact->SetFactory("Graph", dropFact);

  auto coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", aggFact);

  auto TentativePFact = rcp(new TentativePFactory_kokkos());
  TentativePFact->SetFactory("Aggregates", aggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());          // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());  // request coarse nullspace
  coarseLevel.Request(*TentativePFact);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  auto coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  size_t NSdim = coarseNullspace->getNumVectors();
  TEST_EQUALITY(NSdim, 2);

  // coarseNullspace->describe(out, Teuchos::VERB_EXTREME);

  // Check interpolation by computing ||fineNS - P*coarseNS||
  auto PtN = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);
  Ptent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  auto diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  coarseLevel.Release("P", TentativePFact.get());  // release Ptent
  coarseLevel.Release("Nullspace", TentativePFact.get());

  auto nspFact = Teuchos::rcp(new NullspaceFactory());
  fineLevel.Request("Nullspace", nspFact.get());

  nspFact->Build(fineLevel);

  auto fineNullspace = fineLevel.Get<RCP<MultiVector> >("Nullspace", nspFact.get());

  TEST_EQUALITY(fineNullspace->getNumVectors(), 2);

  // diff = fineNS - (P*coarseNS)
  diff->update(1.0, *fineNullspace, -1.0, *PtN, 0.0);

  Array<magnitude_type> norms(NSdim);
  diff->norm2(norms);
  for (decltype(NSdim) i = 0; i < NSdim; ++i)
    TEST_COMPARE_CONST(norms[i], <, 100 * TMT::eps());

  // Check normalization and orthogonality of prolongator columns by
  // making sure that P^T*P = I
  RCP<Matrix> PtentTPtent = MatrixMatrix::Multiply(*Ptent, true, *Ptent, false, out);
  RCP<Vector> diagVec     = VectorFactory::Build(PtentTPtent->getRowMap());
  PtentTPtent->getLocalDiagCopy(*diagVec);
  if (STS::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(diagVec->norm1(), STS::magnitude(diagVec->getGlobalLength()), 100 * TMT::eps());
  TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100 * TMT::eps());
  if (STS::name().find("complex") == std::string::npos)  // skip check for Scalar=complex
    TEST_FLOATING_EQUALITY(STS::magnitude(diagVec->meanValue()), TMT::one(), 100 * TMT::eps());
  TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());
}  // MakeTentativeVectorBasedUsingDefaultNullspace

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TentativePFactory_kokkos, MakeTentativeUsingDefaultNullspace, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO);

  using STS            = Teuchos::ScalarTraits<Scalar>;
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test QR when nullspace isn't supplied by user" << std::endl;

  Level fineLevel, coarseLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  auto A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(199);

  fineLevel.Set("A", A);

  auto amalgFact = rcp(new AmalgamationFactory());

  auto dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  auto aggFact = rcp(new UncoupledAggregationFactory());
  ParameterList aggParams;
  aggParams.set("aggregation: ordering", "natural");
  aggParams.set("aggregation: min agg size", 3);
  aggParams.set("aggregation: max selected neighbors", 0);
  aggFact->SetParameterList(aggParams);
  aggFact->SetFactory("DofsPerNode", dropFact);
  aggFact->SetFactory("Graph", dropFact);

  auto coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", aggFact);

  auto TentativePFact = rcp(new TentativePFactory_kokkos());
  TentativePFact->SetFactory("Aggregates", aggFact);
  TentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  TentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", TentativePFact.get());          // request Ptent
  coarseLevel.Request("Nullspace", TentativePFact.get());  // request coarse nullspace
  coarseLevel.Request(*TentativePFact);
  TentativePFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> Ptent;
  coarseLevel.Get("P", Ptent, TentativePFact.get());

  auto coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());

  coarseLevel.Release("P", TentativePFact.get());          // release Ptent
  coarseLevel.Release("Nullspace", TentativePFact.get());  // release coarse nullspace

  // grab default fine level nullspace (vector of all ones)
  auto nullSpace = MultiVectorFactory::Build(A->getRowMap(), 1);
  nullSpace->putScalar(1.0);

  // check interpolation
  LO NSdim = 1;
  auto PtN = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  Ptent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  auto diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  // diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
  diff->update(1.0, *nullSpace, -1.0, *PtN, 0.0);

  Teuchos::Array<magnitude_type> norms(NSdim);
  diff->norm2(norms);
  for (LO i = 0; i < NSdim; ++i)
    TEST_COMPARE_CONST(norms[i], <, 100 * TMT::eps());

}  // MakeTentativeUsingDefaultNullspace

#if 0
  TEUCHOS_UNIT_TEST(TentativePFactory, NonStandardMaps)
  {

    //#warning Unit test PgPFactory NonStandardMaps disabled
    //  return;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    // generate problem
    LO maxLevels = 3;
    //LO its=10;
    GO nEle = 63;
    GO nIndexBase = 10;
    const RCP<const Map> map = MapFactory::Build(lib, nEle, nIndexBase, comm);
    RCP<Matrix> mtx = Galeri::Xpetra::MatrixTraits<Map,CrsMatrixWrap>::Build(map, 3);

    LocalOrdinal NumMyElements = map->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
    GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
    assert(NumGlobalElements == nEle);

    GlobalOrdinal NumEntries;
    LocalOrdinal nnz=2;
    std::vector<Scalar> Values(nnz);
    std::vector<GlobalOrdinal> Indices(nnz);

    Scalar a = 2.0;
    Scalar b = -1.0;
    Scalar c = -1.0;

    for (LocalOrdinal i = 0; i < NumMyElements; ++i)
    {
      if (MyGlobalElements[i] == nIndexBase)
      {
        // off-diagonal for first row
        Indices[0] = nIndexBase;
        NumEntries = 1;
        Values[0] = c;
      }
      else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1)
      {
        // off-diagonal for last row
        Indices[0] = nIndexBase + NumGlobalElements - 2;
        NumEntries = 1;
        Values[0] = b;
      }
      else
      {
        // off-diagonal for internal row
        Indices[0] = MyGlobalElements[i] - 1;
        Values[1] = b;
        Indices[1] = MyGlobalElements[i] + 1;
        Values[0] = c;
        NumEntries = 2;
      }

      // put the off-diagonal entries
      // Xpetra wants ArrayViews (sigh)
      Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
      Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
      mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

      // Put in the diagonal entry
      mtx->insertGlobalValues(MyGlobalElements[i],
          Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
          Teuchos::tuple<Scalar>(a) );

    } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)


    mtx->fillComplete(map,map);

    std::cout << map->getIndexBase() << std::endl;

    RCP<Matrix> Op = Teuchos::rcp_dynamic_cast<Matrix>(mtx);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (SC) 1.0);

    // fill hierarchy
    RCP<Hierarchy> H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    RCP<Level> Finest = H->GetLevel(); // first associate level with hierarchy (for defaultFactoryHandler!)

    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",Op);                      // set fine level matrix
    Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

    // define transfer operators
    RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
    UncoupledAggFact->SetMinNodesPerAggregate(3);
    UncoupledAggFact->SetMaxNeighAlreadySelected(0);
    UncoupledAggFact->SetOrdering("natural");

    RCP<TentativePFactory> Pfact = rcp(new TentativePFactory());
    RCP<Factory>          Rfact = rcp( new TransPFactory() );
    RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
    H->SetMaxCoarseSize(1);

    // setup smoothers
    Teuchos::ParameterList smootherParamList;
    smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smootherParamList.set("relaxation: sweeps", (LO) 1);
    smootherParamList.set("relaxation: damping factor", (SC) 1.0);
    RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    FactoryManager M;
    M.SetFactory("P", Pfact);
    M.SetFactory("R", Rfact);
    M.SetFactory("A", Acfact);
    M.SetFactory("Ptent", Pfact);
    M.SetFactory("Aggregates", UncoupledAggFact);
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", SmooFact);

    H->Setup(M, 0, maxLevels);

    RCP<Level> coarseLevel = H->GetLevel(1);
    TEST_EQUALITY(coarseLevel->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("A",MueLu::NoFactory::get()), true);

    TEST_EQUALITY(coarseLevel->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsRequested("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("P",Pfact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel->GetKeepFlag("R",Rfact.get()), 0);
    RCP<Level> coarseLevel2 = H->GetLevel(2);
    TEST_EQUALITY(coarseLevel2->IsRequested("A",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("P",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("R",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("PreSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("PostSmoother",MueLu::NoFactory::get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("A",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("P",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",MueLu::NoFactory::get()), false); // coarse level only has presmoother = coarse solver
    TEST_EQUALITY(coarseLevel2->IsAvailable("R",MueLu::NoFactory::get()), true);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("A",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother",MueLu::NoFactory::get()), MueLu::Final);
    // TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",MueLu::NoFactory::get()), MueLu::Final); // coarse level only has presmoother = coarse solver
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("R",MueLu::NoFactory::get()), MueLu::Final);
    TEST_EQUALITY(coarseLevel2->IsRequested("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsRequested("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("P",Pfact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PreSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("PostSmoother",SmooFact.get()), false);
    TEST_EQUALITY(coarseLevel2->IsAvailable("R",Rfact.get()), false);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("P",Pfact.get()), 0);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("PreSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("PostSmoother",SmooFact.get()), 0);
    TEST_EQUALITY(coarseLevel2->GetKeepFlag("R",Rfact.get()), 0);

  }

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_IFPACK) && defined(HAVE_MUELU_IFPACK2)
  TEUCHOS_UNIT_TEST(TentativePFactory, EpetraVsTpetra)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "Test QR when nullspace isn't supplied by user" << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    using TST            = Teuchos::ScalarTraits<SC>;
    using magnitude_type = typename TST::magnitudeType;
    using TMT            = Teuchos::ScalarTraits<magnitudeType>;

    Teuchos::Array<magnitude_type> results(2);

    // run test only on 1 proc
    if(comm->getSize() == 1)
    {
      Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

      // run Epetra and Tpetra test
      for (int run = 0; run < 2; run++)
      {
        if (run == 0) lib = Xpetra::UseEpetra;
        else lib = Xpetra::UseTpetra;

        // generate problem
        LO maxLevels = 3;
        LO its=10;
        LO nEle = 63;
        const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
        Teuchos::ParameterList matrixParameters;
        matrixParameters.set("nx", Teuchos::as<GO>(nEle));
        RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap,MultiVector>("Laplace1D", map, matrixParameters);
        RCP<Matrix> Op = Pr->BuildMatrix();

        // build nullspace
        RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
        nullSpace->putScalar( (SC) 1.0);
        Teuchos::Array<magnitude_type> norms(1);
        nullSpace->norm1(norms);
        if (comm->getRank() == 0)
          out << "||NS|| = " << norms[0] << std::endl;

        // fill hierarchy
        RCP<Hierarchy> H = rcp( new Hierarchy() );
        H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
        RCP<Level> Finest = H->GetLevel(); // first associate level with hierarchy (for defaultFactoryHandler!)

        Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
        Finest->Set("A",Op);                      // set fine level matrix
        Finest->Set("Nullspace",nullSpace);       // set null space information for finest level

        // define transfer operators
        RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
        UncoupledAggFact->SetMinNodesPerAggregate(3);
        UncoupledAggFact->SetMaxNeighAlreadySelected(0);
        UncoupledAggFact->SetOrdering("natural");

        RCP<TentativePFactory> Pfact = rcp(new TentativePFactory());
        RCP<Factory>          Rfact = rcp( new TransPFactory() );
        RCP<RAPFactory>        Acfact = rcp( new RAPFactory() );
        H->SetMaxCoarseSize(1);

        // setup smoothers
        Teuchos::ParameterList smootherParamList;
        smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
        smootherParamList.set("relaxation: sweeps", (LO) 1);
        smootherParamList.set("relaxation: damping factor", (SC) 1.0);
        RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother("RELAXATION", smootherParamList) );
        RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
        Acfact->setVerbLevel(Teuchos::VERB_HIGH);

        RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

        FactoryManager M;
        M.SetFactory("P", Pfact);
        M.SetFactory("R", Rfact);
        M.SetFactory("A", Acfact);
        M.SetFactory("Ptent", Pfact);
        M.SetFactory("Aggregates", UncoupledAggFact);
        M.SetFactory("Smoother", SmooFact);
        M.SetFactory("CoarseSolver", coarseSolveFact);

        H->Setup(M, 0, maxLevels);

        // test some basic multgrid data
        RCP<Level> coarseLevel = H->GetLevel(1);
        RCP<Matrix> P1 = coarseLevel->Get< RCP<Matrix> >("P");
        RCP<Matrix> R1 = coarseLevel->Get< RCP<Matrix> >("R");
        TEST_EQUALITY(P1->getGlobalNumRows(), 63);
        TEST_EQUALITY(P1->getGlobalNumCols(), 21);
        TEST_EQUALITY(R1->getGlobalNumRows(), 21);
        TEST_EQUALITY(R1->getGlobalNumCols(), 63);
        RCP<Level> coarseLevel2 = H->GetLevel(2);
        RCP<Matrix> P2 = coarseLevel2->Get< RCP<Matrix> >("P");
        RCP<Matrix> R2 = coarseLevel2->Get< RCP<Matrix> >("R");
        TEST_EQUALITY(P2->getGlobalNumRows(), 21);
        TEST_EQUALITY(P2->getGlobalNumCols(), 7);
        TEST_EQUALITY(R2->getGlobalNumRows(), 7);
        TEST_EQUALITY(R2->getGlobalNumCols(), 21);

        Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > PtentTPtent = Xpetra::MatrixMatrix<Scalar,LO,GO,Node>::Multiply(*P1,true,*P1,false,out);
        Teuchos::RCP<Xpetra::Vector<Scalar,LO,GO,Node> > diagVec = Xpetra::VectorFactory<Scalar,LO,GO,Node>::Build(PtentTPtent->getRowMap());
        PtentTPtent->getLocalDiagCopy(*diagVec);
        TEST_FLOATING_EQUALITY(diagVec->norm1(), TST::magnitude(diagVec->getGlobalLength()), 100*TMT:eps());
        TEST_FLOATING_EQUALITY(diagVec->normInf(), TMT::one(), 100*TMT:eps());
        TEST_FLOATING_EQUALITY(diagVec->meanValue(), TMT::one(), 100*TMT:eps());
        TEST_EQUALITY(PtentTPtent->getGlobalNumEntries(), diagVec->getGlobalLength());

        // Define RHS
        RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

        X->putScalar(1.0);
        X->norm2(norms);
        if (comm->getRank() == 0)
          out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

        Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

        // Use AMG directly as an iterative method
        {
          X->putScalar( (SC) 0.0);

          H->Iterate(*RHS,*X,its);

          X->norm2(norms);
          if (comm->getRank() == 0)
            out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
          results[run] = norms[0];
        }
      }

      TEST_FLOATING_EQUALITY(results[0], results[1], 100*TMT::eps()); // check results of EPETRA vs TPETRA
    } // comm->getSize == 1

  } // TentativePFactory_EpetraVsTpetra
#endif
#endif

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory_kokkos, Constructor, SC, LO, GO, NO)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory_kokkos, MakeTentative, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory_kokkos, MakeTentativeUsingDefaultNullspace, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TentativePFactory_kokkos, MakeTentativeVectorBasedUsingDefaultNullspace, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

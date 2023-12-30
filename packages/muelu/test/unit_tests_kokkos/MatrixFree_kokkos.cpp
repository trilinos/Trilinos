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

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_MatrixFreeTentativePFactory.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"
#include "MueLu_NullspaceFactory_kokkos.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MatrixFreeTentativePFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  if (TestHelpers_kokkos::Parameters::getLib() == Xpetra::UseEpetra) {
    out << "skipping test for linAlgebra==UseEpetra" << std::endl;
    return;
  }
  out << "version: " << MueLu::Version() << std::endl;

  RCP<MatrixFreeTentativePFactory> MFtentPFact = rcp(new MatrixFreeTentativePFactory);
  TEST_EQUALITY(MFtentPFact != Teuchos::null, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MatrixFreeTentativePFactory, MakeTentative, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  out << "Test with user-supplied nullspace" << std::endl;

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
  nullSpace->putScalar(1.0);
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());

  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  ParameterList aggParams;
  aggParams.set("aggregation: ordering", "natural");
  aggParams.set("aggregation: deterministic", true);
  aggParams.set("aggregation: min agg size", 3);
  aggParams.set("aggregation: max selected neighbors", 0);
  aggFact->SetParameterList(aggParams);
  aggFact->SetFactory("DofsPerNode", dropFact);
  aggFact->SetFactory("Graph", dropFact);

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", aggFact);

  RCP<MatrixFreeTentativePFactory> MFTentativePFact = rcp(new MatrixFreeTentativePFactory());
  MFTentativePFact->SetFactory("Aggregates", aggFact);
  MFTentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  MFTentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", MFTentativePFact.get());  // request Ptent
  coarseLevel.Request("Nullspace", MFTentativePFact.get());
  coarseLevel.Request(*MFTentativePFact);
  MFTentativePFact->Build(fineLevel, coarseLevel);

  RCP<Operator> Ptent;
  coarseLevel.Get("P", Ptent, MFTentativePFact.get());

  RCP<MultiVector> coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", MFTentativePFact.get());

  // Check interpolation by computing ||fineNS - P*coarseNS||
  // PtN = P*coarseNS
  RCP<MultiVector> PtN = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);
  Ptent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

  RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  coarseLevel.Release("P", MFTentativePFact.get());  // release Ptent
  coarseLevel.Release("Nullspace", MFTentativePFact.get());

  // diff = fineNS - PtN = fineNS - P*coarseNS
  diff->update(1.0, *nullSpace, -1.0, *PtN, 0.0);

  Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(NSdim);
  diff->norm2(norms);
  for (LO i = 0; i < NSdim; ++i) {
    out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
    TEST_COMPARE_CONST(norms[i], <, 100 * TMT::eps());
  }

}  // MakeTentative

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MatrixFreeTentativePFactory, MakeTentativeVectorBasedUsingDefaultNullspace, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  out << "Test when nullspace isn't supplied by user" << std::endl;

  Level fineLevel, coarseLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  auto A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(200);

  A->SetFixedBlockSize(2);

  fineLevel.Set("A", A);

  auto amalgFact = rcp(new AmalgamationFactory());

  auto dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  auto aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  ParameterList aggParams;
  aggParams.set("aggregation: ordering", "natural");
  aggParams.set("aggregation: deterministic", true);
  aggParams.set("aggregation: min agg size", 3);
  aggParams.set("aggregation: max selected neighbors", 0);
  aggFact->SetParameterList(aggParams);
  aggFact->SetFactory("DofsPerNode", dropFact);
  aggFact->SetFactory("Graph", dropFact);

  auto coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", aggFact);

  auto MFTentativePFact = rcp(new MatrixFreeTentativePFactory());
  MFTentativePFact->SetFactory("Aggregates", aggFact);
  MFTentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  MFTentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", MFTentativePFact.get());          // request Ptent
  coarseLevel.Request("Nullspace", MFTentativePFact.get());  // request coarse nullspace
  coarseLevel.Request(*MFTentativePFact);
  MFTentativePFact->Build(fineLevel, coarseLevel);

  RCP<Operator> Ptent;
  coarseLevel.Get("P", Ptent, MFTentativePFact.get());

  auto coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", MFTentativePFact.get());

  size_t NSdim = coarseNullspace->getNumVectors();
  TEST_EQUALITY(NSdim, 2);

  // coarseNullspace->describe(out, Teuchos::VERB_EXTREME);

  // Check interpolation by computing ||fineNS - P*coarseNS||
  auto PtN = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);
  // Ptent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0); // TODO: don't support strides :)

  auto diff = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  diff->putScalar(0.0);

  coarseLevel.Release("P", MFTentativePFact.get());  // release Ptent
  coarseLevel.Release("Nullspace", MFTentativePFact.get());

  auto nspFact = Teuchos::rcp(new NullspaceFactory_kokkos());
  fineLevel.Request("Nullspace", nspFact.get());

  nspFact->Build(fineLevel);

  auto fineNullspace = fineLevel.Get<RCP<MultiVector> >("Nullspace", nspFact.get());

  TEST_EQUALITY(fineNullspace->getNumVectors(), 2);

  // diff = fineNS - (P*coarseNS)
  diff->update(1.0, *fineNullspace, -1.0, *PtN, 0.0);

  Array<magnitude_type> norms(NSdim);
  diff->norm2(norms);
  // for (decltype(NSdim) i = 0; i < NSdim; ++i)
  //   TEST_COMPARE_CONST(norms[i], <, 100*TMT::eps());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MatrixFreeTentativePFactory, MakeTentativeUsingDefaultNullspace, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  out << "Test when nullspace isn't supplied by user" << std::endl;

  Level fineLevel, coarseLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

  auto A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(199);

  fineLevel.Set("A", A);

  auto amalgFact = rcp(new AmalgamationFactory());

  auto dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  auto aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  ParameterList aggParams;
  aggParams.set("aggregation: ordering", "natural");
  aggParams.set("aggregation: deterministic", true);
  aggParams.set("aggregation: min agg size", 3);
  aggParams.set("aggregation: max selected neighbors", 0);
  aggFact->SetParameterList(aggParams);
  aggFact->SetFactory("DofsPerNode", dropFact);
  aggFact->SetFactory("Graph", dropFact);

  auto coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", aggFact);

  auto MFTentativePFact = rcp(new MatrixFreeTentativePFactory());
  MFTentativePFact->SetFactory("Aggregates", aggFact);
  MFTentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
  MFTentativePFact->SetFactory("CoarseMap", coarseMapFact);

  coarseLevel.Request("P", MFTentativePFact.get());          // request Ptent
  coarseLevel.Request("Nullspace", MFTentativePFact.get());  // request coarse nullspace
  coarseLevel.Request(*MFTentativePFact);
  MFTentativePFact->Build(fineLevel, coarseLevel);

  RCP<Operator> Ptent;
  coarseLevel.Get("P", Ptent, MFTentativePFact.get());

  auto coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", MFTentativePFact.get());

  coarseLevel.Release("P", MFTentativePFact.get());          // release Ptent
  coarseLevel.Release("Nullspace", MFTentativePFact.get());  // release coarse nullspace

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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MatrixFreeTentativePFactory, MatrixVsMatrixFree, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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
  out << "Test matvec speed in matrix and matrix-free versions and compare results" << std::endl;

  unsigned int numMatVecs = 100;

  RCP<Matrix> A      = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(199);
  RCP<MultiVector> X = MultiVectorFactory::Build(A->getDomainMap(), 1);
  RCP<MultiVector> Xcoarse;
  RCP<MultiVector> MFXcoarse;
  RCP<MultiVector> diff;
  X->randomize();

  // Matrix version
  {
    Level fineLevel, coarseLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);
    fineLevel.Request("A");
    fineLevel.Set("A", A);

    // Only one NS vector -> exercises manual orthogonalization
    LO NSdim                   = 1;
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
    nullSpace->putScalar(1.0);
    fineLevel.Set("Nullspace", nullSpace);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());

    RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
    ParameterList aggParams;
    aggParams.set("aggregation: ordering", "natural");
    aggParams.set("aggregation: deterministic", true);
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
    Xcoarse = MultiVectorFactory::Build(Ptent->getDomainMap(), 1);
    Ptent->apply(*X, *Xcoarse, Teuchos::TRANS);

    RCP<MultiVector> coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", TentativePFact.get());
    RCP<MultiVector> PtN             = MultiVectorFactory::Build(Ptent->getRangeMap(), NSdim);

    // warm-up apply operations
    for (unsigned int i = 0; i < 5; ++i)
      Ptent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

    {
      RCP<Teuchos::TimeMonitor> tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MatVec original")));
      for (unsigned int i = 0; i < numMatVecs; ++i)
        Ptent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);
    }
  }

  // Matrix-free version
  {
    Level fineLevel, coarseLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);

    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);
    fineLevel.Request("A");
    fineLevel.Set("A", A);

    // Only one NS vector -> exercises manual orthogonalization
    LO NSdim                   = 1;
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
    nullSpace->putScalar(1.0);
    fineLevel.Set("Nullspace", nullSpace);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());

    RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
    ParameterList aggParams;
    aggParams.set("aggregation: ordering", "natural");
    aggParams.set("aggregation: deterministic", true);
    aggParams.set("aggregation: min agg size", 3);
    aggParams.set("aggregation: max selected neighbors", 0);
    aggFact->SetParameterList(aggParams);
    aggFact->SetFactory("DofsPerNode", dropFact);
    aggFact->SetFactory("Graph", dropFact);

    RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
    coarseMapFact->SetFactory("Aggregates", aggFact);

    RCP<MatrixFreeTentativePFactory> MFTentativePFact = rcp(new MatrixFreeTentativePFactory());
    MFTentativePFact->SetFactory("Aggregates", aggFact);
    MFTentativePFact->SetFactory("UnAmalgamationInfo", amalgFact);
    MFTentativePFact->SetFactory("CoarseMap", coarseMapFact);

    coarseLevel.Request("P", MFTentativePFact.get());  // request Ptent
    coarseLevel.Request("Nullspace", MFTentativePFact.get());
    coarseLevel.Request(*MFTentativePFact);
    MFTentativePFact->Build(fineLevel, coarseLevel);

    RCP<Operator> MFPtent;
    coarseLevel.Get("P", MFPtent, MFTentativePFact.get());
    MFXcoarse = MultiVectorFactory::Build(MFPtent->getDomainMap(), 1);
    MFPtent->apply(*X, *MFXcoarse, Teuchos::TRANS);

    RCP<MultiVector> coarseNullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", MFTentativePFact.get());
    RCP<MultiVector> PtN             = MultiVectorFactory::Build(MFPtent->getRangeMap(), NSdim);
    // warm-up apply operations
    for (unsigned int i = 0; i < 5; ++i)
      MFPtent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);

    {
      RCP<Teuchos::TimeMonitor> tm = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("MatVec matrix-free")));
      for (unsigned int i = 0; i < numMatVecs; ++i)
        MFPtent->apply(*coarseNullspace, *PtN, Teuchos::NO_TRANS, 1.0, 0.0);
    }

    // need scope
    diff = MultiVectorFactory::Build(MFPtent->getDomainMap(), 1);
    diff->putScalar(0.0);
  }

  Teuchos::TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
  Teuchos::TimeMonitor::zeroOutTimers();

  // Compare Xcoarse and MFXcoarse (diff = Xcoarse - MFXcoarse + 0*diff)
  diff->update(1.0, *Xcoarse, -1.0, *MFXcoarse, 0.0);
  Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
  diff->norm2(norms);
  out << "||diff||_2 = " << norms[0] << std::endl;
  TEST_COMPARE_CONST(norms[0], <, 100 * TMT::eps());

}  // MatrixVsMatrixFree

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MatrixFreeTentativePFactory, Constructor, SC, LO, GO, NO)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MatrixFreeTentativePFactory, MakeTentative, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MatrixFreeTentativePFactory, MakeTentativeUsingDefaultNullspace, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MatrixFreeTentativePFactory, MatrixVsMatrixFree, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

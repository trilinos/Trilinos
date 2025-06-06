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

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_ConstraintFactory.hpp>
#include <MueLu_Constraint.hpp>
#include <MueLu_PatternFactory.hpp>
#include <MueLu_EminPFactory.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, NullspaceConstraint, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

  std::vector<std::string> matrixTypes = {
      "Laplace1D",
      "Laplace2D",
      "Laplace3D",
      "Brick3D",
      "Elasticity2D",
      "Elasticity3D",
  };
  for (auto &matrixType : matrixTypes) {
    out << "\n\n==================================================\nTesting " << matrixType << "\n\n";

    Level fineLevel, coarseLevel;
    test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);

    GlobalOrdinal nx, ny = 1, nz = 1;
    if (matrixType == "Laplace1D") {
      nx = 200;
    } else if ((matrixType == "Laplace2D") || (matrixType == "Brick2D") || (matrixType == "Elasticity2D")) {
      nx = 20;
      ny = 20;
    } else if ((matrixType == "Laplace3D") || (matrixType == "Brick3D") || (matrixType == "Elasticity3D")) {
      nx = 20;
      ny = 20;
      nz = 20;
    }

    Teuchos::ParameterList galeriList;
    galeriList.set("matrixType", matrixType);
    galeriList.set("nx", nx);
    galeriList.set("ny", ny);
    galeriList.set("nz", nz);

    auto [A, coordinates, nullSpace, DofsPerNode] = test_factory::BuildMatrixCoordsNullspace(galeriList);

    fineLevel.Request("A");
    fineLevel.Set("A", A);
    fineLevel.Set("Coordinates", coordinates);
    fineLevel.Set("Nullspace", nullSpace);
    fineLevel.Set("DofsPerNode", DofsPerNode);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
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

    RCP<PatternFactory> patternFact = rcp(new PatternFactory());
    patternFact->SetFactory("P", TentativePFact);

    RCP<ConstraintFactory> constraintFact = rcp(new ConstraintFactory());
    constraintFact->SetFactory("CoarseNullspace", TentativePFact);
    constraintFact->SetFactory("Ppattern", patternFact);

    RCP<EminPFactory> eminFact = rcp(new EminPFactory());
    eminFact->SetFactory("P", TentativePFact);
    eminFact->SetFactory("Constraint", constraintFact);

    coarseLevel.Request("P", TentativePFact.get());  // request Ptent
    coarseLevel.Request("P", eminFact.get());        // request P
    coarseLevel.Request("Nullspace", TentativePFact.get());
    coarseLevel.Request("Constraint", constraintFact.get());
    coarseLevel.Request(*eminFact);
    TentativePFact->Build(fineLevel, coarseLevel);

    RCP<Matrix> Ptent, P;
    coarseLevel.Get("P", Ptent, TentativePFact.get());
    coarseLevel.Get("P", P, eminFact.get());

    RCP<Constraint> constraint;
    coarseLevel.Get("Constraint", constraint, constraintFact.get());

    using Magnitude = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
    const auto eps  = Teuchos::ScalarTraits<Magnitude>::eps();

    // Test that Ptent satisfies the constraint.
    TEST_COMPARE(constraint->ResidualNorm(Ptent), <, 400 * eps);

    // Test that P satisfies the constraint.
    TEST_COMPARE(constraint->ResidualNorm(P), <, 20000 * eps);

    // Test that P has lower energy norm than Ptent.
    auto energyNormPtent = EminPFactory::ComputeProlongatorEnergyNorm(A, Ptent, out);
    auto energyNormP     = EminPFactory::ComputeProlongatorEnergyNorm(A, P, out);
    TEST_COMPARE(energyNormP, <, energyNormPtent);
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

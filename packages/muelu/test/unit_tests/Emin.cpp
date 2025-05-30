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
#include <Xpetra_IO.hpp>

#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_ConstraintFactory.hpp>
#include <MueLu_Constraint.hpp>
#include <MueLu_DenseConstraint.hpp>
#include <MueLu_PatternFactory.hpp>
#include <MueLu_EdgeProlongatorPatternFactory.hpp>
#include <MueLu_EminPFactory.hpp>
#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include "MueLu_NoFactory.hpp"
#include "MueLu_ReitzingerPFactory_decl.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_VerbosityLevel.hpp"

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

    Level fineLevel;
    Level coarseLevel;
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

    // Test that both Ptent satisfies the constraint after converting it to a vector and back to a matrix.
    auto vecP = MultiVectorFactory::Build(constraint->getDomainMap(), 1);
    constraint->AssignMatrixEntriesToVector(*Ptent, *vecP);
    auto Ptent2 = constraint->GetMatrixWithEntriesFromVector(*vecP);
    TEST_COMPARE(constraint->ResidualNorm(Ptent2), <, 400 * eps);

    // Teuchos::rcp_const_cast<CrsGraph>(Ptent->getCrsGraph())->computeGlobalConstants();
    // Teuchos::rcp_const_cast<CrsGraph>(Ptent2->getCrsGraph())->computeGlobalConstants();
    // Ptent->describe(out, Teuchos::VERB_EXTREME);
    // Ptent2->describe(out, Teuchos::VERB_EXTREME);

    // Test that P satisfies the constraint.
    TEST_COMPARE(constraint->ResidualNorm(P), <, 20000 * eps);

    // Test that P has lower energy norm than Ptent.
    auto energyNormPtent = EminPFactory::ComputeProlongatorEnergyNorm(A, Ptent, out);
    auto energyNormP     = EminPFactory::ComputeProlongatorEnergyNorm(A, P, out);
    TEST_COMPARE(energyNormP, <, energyNormPtent);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
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

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  int numProcs                       = comm->getSize();
  if (numProcs != 1) {
    std::cout << "\nThis test must be run on 1 processor!\n"
              << std::endl;
    return;
  }

  std::string scalarName = Teuchos::ScalarTraits<Scalar>::name();
  out << "scalar type = " << scalarName << std::endl;
  if (scalarName.find("complex") != std::string::npos) {
    out << "Skipping Test for SC=complex" << std::endl;
    return;
  }

  LocalOrdinal indexBase = 0;

  GlobalOrdinal global_num_fine_edges   = 56;
  GlobalOrdinal global_num_coarse_edges = 5;
  GlobalOrdinal global_num_fine_nodes   = 25;
  GlobalOrdinal global_num_coarse_nodes = 4;

  LocalOrdinal local_num_fine_edges   = 56;
  LocalOrdinal local_num_coarse_edges = 5;
  LocalOrdinal local_num_fine_nodes   = 25;
  LocalOrdinal local_num_coarse_nodes = 4;

  auto lib = TestHelpers::Parameters::getLib();

  auto fine_edge_map    = MapFactory::Build(lib, global_num_fine_edges, local_num_fine_edges, indexBase, comm);
  auto coarse_edge_map  = MapFactory::Build(lib, global_num_coarse_edges, local_num_coarse_edges, indexBase, comm);
  auto fine_nodal_map   = MapFactory::Build(lib, global_num_fine_nodes, local_num_fine_nodes, indexBase, comm);
  auto coarse_nodal_map = MapFactory::Build(lib, global_num_coarse_nodes, local_num_coarse_nodes, indexBase, comm);

  auto A = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("emin_matrices/A.mm", fine_edge_map, fine_edge_map, fine_edge_map, fine_edge_map);
  auto D = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("emin_matrices/D0h.mm", fine_edge_map, fine_nodal_map, fine_nodal_map, fine_edge_map);

  // Auxiliary nodal hierarchy
  RCP<Matrix> NodeAggMatrix, NodeAggMatrixCoarse, Pnodal, Ptentnodal;
  {
    auto A_D0     = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *D, false, out, true, true);
    NodeAggMatrix = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*D, true, *A_D0, false, out, true, true);

    RCP<Hierarchy> H = rcp(new Hierarchy("Nodal Hierarchy"));
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<Level> fineLevelNodal = H->GetLevel();
    fineLevelNodal->Set("A", NodeAggMatrix);
    FactoryManager M;
    M.SetKokkosRefactor(false);
    RCP<TentativePFactory> Ptentfact = rcp(new TentativePFactory());
    Ptentfact->SetParameter("sa: keep tentative prolongator", Teuchos::ParameterEntry(true));
    Ptentfact->SetParameter("tentative: calculate qr", Teuchos::ParameterEntry(false));
    Ptentfact->SetParameter("tentative: constant column sums", Teuchos::ParameterEntry(false));
    M.SetFactory("Ptent", Ptentfact);
    H->SetMaxCoarseSize(1);
    H->Setup(M, 0, 2);

    RCP<Level> coarseLevelNodal = H->GetLevel(1);
    Ptentnodal                  = coarseLevelNodal->Get<RCP<Matrix>>("Ptent");
    Pnodal                      = coarseLevelNodal->Get<RCP<Matrix>>("P");
    NodeAggMatrixCoarse         = coarseLevelNodal->Get<RCP<Matrix>>("A");
  }

  // Edge hierarchy
  RCP<Matrix> P0, P;
  RCP<Constraint> constraint;
  {
    Level fineLevel, coarseLevel;
    test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
    coarseLevel.SetFactoryManager(Teuchos::null);

    fineLevel.Set("A", A);
    fineLevel.Set("D0", D);

    // Ptentnodal is used to construct coarse D0
    coarseLevel.Set("Pnodal", Ptentnodal);

    RCP<ReitzingerPFactory> reitzingerFact = rcp(new ReitzingerPFactory());

    // PnodalEmin is used to construct pattern for P
    coarseLevel.Set("PnodalEmin", Pnodal);

    fineLevel.Set("NodeAggMatrix", NodeAggMatrix);
    coarseLevel.Set("NodeAggMatrix", NodeAggMatrixCoarse);

    RCP<EdgeProlongatorPatternFactory> patternFact = rcp(new EdgeProlongatorPatternFactory());
    patternFact->SetFactory("CoarseD0", reitzingerFact);

    RCP<ConstraintFactory> constraintFact = rcp(new ConstraintFactory());
    constraintFact->SetFactory("Ppattern", patternFact);
    constraintFact->SetParameter("emin: constraint type", Teuchos::ParameterEntry(std::string("maxwell")));
    constraintFact->SetFactory("CoarseD0", reitzingerFact);

    RCP<EminPFactory> eminFact = rcp(new EminPFactory());
    eminFact->SetFactory("Constraint", constraintFact);
    eminFact->SetFactory("P", constraintFact);

    // coarseLevel.Request("D0", reitzingerFact.get());

    coarseLevel.Request("Constraint", constraintFact.get());
    coarseLevel.Request("P", constraintFact.get());
    coarseLevel.Request("P", eminFact.get());
    coarseLevel.Request(*eminFact);

    coarseLevel.Get("Constraint", constraint, constraintFact.get());

    // This is the initial guess used for emin.
    coarseLevel.Get("P", P0, constraintFact.get());

    // This is the result after running the minimization.
    coarseLevel.Get("P", P, eminFact.get());
  }

  const auto eps = Teuchos::ScalarTraits<magnitude_type>::eps();

  // Test that P0 satisfies the constraint.
  TEST_COMPARE(constraint->ResidualNorm(P0), <, 4000 * eps);

  // Test that P satisfies the constraint.
  TEST_COMPARE(constraint->ResidualNorm(P), <, 4000 * eps);

  // Test that P has lower energy norm than P0.
  auto energyNormP0 = EminPFactory::ComputeProlongatorEnergyNorm(A, P0, out);
  auto energyNormP  = EminPFactory::ComputeProlongatorEnergyNorm(A, P, out);
  TEST_COMPARE(energyNormP, <, energyNormP0);
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

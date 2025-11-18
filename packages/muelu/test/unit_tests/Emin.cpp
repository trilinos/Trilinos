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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void testNullspaceConstraint(const std::string &matrixType, Teuchos::FancyOStream &out, bool &success) {
#include "MueLu_UseShortNames.hpp"

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real                  = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real, LO, GO, NO>;
  using test_factory          = TestHelpers::TestFactory<SC, LO, GO, NO>;

  out << "version: " << MueLu::Version() << std::endl;

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
  } else if ((matrixType == "Laplace3D") || (matrixType == "Brick3D")) {
    nx = 20;
    ny = 20;
    nz = 20;
  } else if (matrixType == "Elasticity3D") {
    nx = 10;
    ny = 10;
    nz = 10;
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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, NullspaceConstraint_Laplace1D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testNullspaceConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>("Laplace1D", out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, NullspaceConstraint_Laplace2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testNullspaceConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>("Laplace2D", out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, NullspaceConstraint_Laplace3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testNullspaceConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>("Laplace3D", out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, NullspaceConstraint_Brick3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testNullspaceConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>("Brick3D", out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, NullspaceConstraint_Elasticity2D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testNullspaceConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>("Elasticity2D", out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, NullspaceConstraint_Elasticity3D, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testNullspaceConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>("Elasticity3D", out, success);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void testMaxwellConstraint(const std::string &inputDir,
                           GlobalOrdinal global_num_fine_edges,
                           GlobalOrdinal global_num_coarse_edges,
                           GlobalOrdinal global_num_fine_nodes,
                           GlobalOrdinal global_num_coarse_nodes,
                           const bool readWriteForTesting,
                           Teuchos::FancyOStream &out, bool &success) {
#include "MueLu_UseShortNames.hpp"

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

  // #define GrindEmin

  LocalOrdinal indexBase = 0;

  auto lib = TestHelpers::Parameters::getLib();

  auto fine_edge_map    = MapFactory::Build(lib, global_num_fine_edges, indexBase, comm);
  auto coarse_edge_map  = MapFactory::Build(lib, global_num_coarse_edges, indexBase, comm);
  auto fine_nodal_map   = MapFactory::Build(lib, global_num_fine_nodes, indexBase, comm);
  auto coarse_nodal_map = MapFactory::Build(lib, global_num_coarse_nodes, indexBase, comm);

  auto A = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(inputDir + "A.dat", fine_edge_map, fine_edge_map, fine_edge_map, fine_edge_map);
  auto D = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(inputDir + "D0h.dat", fine_edge_map, fine_nodal_map, fine_nodal_map, fine_edge_map);

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
  bool useExternalP0 = false;
  if (readWriteForTesting) {
    // Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("Ptent.code", *Ptentnodal);
    // Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("Pnodal.code", *Pnodal);
    // Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("NodeAggMatrixCoarse.code", *NodeAggMatrixCoarse);
    Ptentnodal = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(inputDir + "Ptent.dat", fine_nodal_map, coarse_nodal_map, coarse_nodal_map, fine_nodal_map);
    Pnodal     = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(inputDir + "Pn.dat", fine_nodal_map, coarse_nodal_map, coarse_nodal_map, fine_nodal_map);

    if (useExternalP0)
      P0 = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(inputDir + "Pe.dat", fine_edge_map, coarse_edge_map, coarse_edge_map, fine_edge_map);
  }

  {
    Level fineLevel, coarseLevel;
    if (P0 != Teuchos::null) coarseLevel.Set("P0", P0);
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
#ifdef GrindEmin
    eminFact->SetParameter("emin: num iterations", Teuchos::ParameterEntry(110));
#endif

    if (readWriteForTesting)
      coarseLevel.Request("D0", reitzingerFact.get());

    coarseLevel.Request("Constraint", constraintFact.get());
    coarseLevel.Request("P", constraintFact.get());
    coarseLevel.Request("P", eminFact.get());
    coarseLevel.Request(*eminFact);

    coarseLevel.Get("Constraint", constraint, constraintFact.get());

    // The initial guess used for emin starts up being call P

    if (coarseLevel.IsAvailable("P0")) {
      coarseLevel.Get("P0", P0);
    } else {
      coarseLevel.Get("P", P0, constraintFact.get());
    }

    // This is the result after running the minimization.
    coarseLevel.Get("P", P, eminFact.get());
    if (readWriteForTesting) {
      RCP<Matrix> D0H;
      D0H = coarseLevel.Get<RCP<Matrix>>("D0");
      Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("D0H.code", *D0H);
      Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("Pe.code", *P);
    }
  }

  const auto eps = Teuchos::ScalarTraits<magnitude_type>::eps();

  // Test that P0 satisfies the constraint.
  TEST_COMPARE(constraint->ResidualNorm(P0), <, 400000 * eps);

  // Test that P satisfies the constraint.
  TEST_COMPARE(constraint->ResidualNorm(P), <, 40000000000000 * eps);

  // Test that P has lower energy norm than P0.
  auto energyNormP0 = EminPFactory::ComputeProlongatorEnergyNorm(A, P0, out);
  auto energyNormP  = EminPFactory::ComputeProlongatorEnergyNorm(A, P, out);
  TEST_COMPARE(energyNormP, <, energyNormP0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_1, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/",
                                                                   /*global_num_fine_edges=*/56,
                                                                   /*global_num_coarse_edges=*/5,
                                                                   /*global_num_fine_nodes=*/25,
                                                                   /*global_num_coarse_nodes=*/4,
                                                                   /*readWriteForTesting=*/false,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_Tris, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/tris/",
                                                                   /*global_num_fine_edges=*/12416,
                                                                   /*global_num_coarse_edges=*/1365,
                                                                   /*global_num_fine_nodes=*/4225,
                                                                   /*global_num_coarse_nodes=*/484,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_TrisWithDir, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/tris/withDir/",
                                                                   /*global_num_fine_edges=*/12352,
                                                                   /*global_num_coarse_edges=*/1387,
                                                                   /*global_num_fine_nodes=*/4160,
                                                                   /*global_num_coarse_nodes=*/484,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_Quads, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/quads/",
                                                                   /*global_num_fine_edges=*/1512,
                                                                   /*global_num_coarse_edges=*/180,
                                                                   /*global_num_fine_nodes=*/784,
                                                                   /*global_num_coarse_nodes=*/100,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_QuadsWithDir, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/quads/withDir/",
                                                                   /*global_num_fine_edges=*/1485,
                                                                   /*global_num_coarse_edges=*/171,
                                                                   /*global_num_fine_nodes=*/756,
                                                                   /*global_num_coarse_nodes=*/90,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_Tets, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/tets/",
                                                                   /*global_num_fine_edges=*/5859,
                                                                   /*global_num_coarse_edges=*/279,
                                                                   /*global_num_fine_nodes=*/1000,
                                                                   /*global_num_coarse_nodes=*/64,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_TetsWithDir, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/tets/withDir/",
                                                                   /*global_num_fine_edges=*/5850,
                                                                   /*global_num_coarse_edges=*/342,
                                                                   /*global_num_fine_nodes=*/990,
                                                                   /*global_num_coarse_nodes=*/62,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_Hexes, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/hexes/",
                                                                   /*global_num_fine_edges=*/2700,
                                                                   /*global_num_coarse_edges=*/144,
                                                                   /*global_num_fine_nodes=*/1000,
                                                                   /*global_num_coarse_nodes=*/64,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(EminPFactory, MaxwellConstraint_HexesWithDir, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  testMaxwellConstraint<Scalar, LocalOrdinal, GlobalOrdinal, Node>(/*inputDir=*/"emin_matrices/hexes/withDir/",
                                                                   /*global_num_fine_edges=*/2691,
                                                                   /*global_num_coarse_edges=*/155,
                                                                   /*global_num_fine_nodes=*/990,
                                                                   /*global_num_coarse_nodes=*/60,
                                                                   /*readWriteForTesting=*/true,
                                                                   out, success);
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint_Laplace1D, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint_Laplace2D, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint_Laplace3D, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint_Brick3D, SC, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint_Elasticity2D, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, NullspaceConstraint_Elasticity3D, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_1, SC, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_Tris, SC, LO, GO, Node)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_TrisWithDir, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_Quads, SC, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_QuadsWithDir, SC, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_Tets, SC, LO, GO, Node)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_TetsWithDir, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_Hexes, SC, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(EminPFactory, MaxwellConstraint_HexesWithDir, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

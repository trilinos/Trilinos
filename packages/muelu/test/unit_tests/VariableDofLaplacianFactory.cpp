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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_VariableDofLaplacianFactory.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(VariableDofLaplacianFactory, VarLaplConstructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  if (!TYPE_EQUAL(GO, int)) {
    out << "Skipping test for GO != int" << std::endl;
    return;
  }
  if (!TYPE_EQUAL(SC, double)) {
    out << "Skipping test for SC != double; no support for complex SC in Galeri::Xpetra branch" << std::endl;
    return;
  }
  out << "version: " << MueLu::Version() << std::endl;

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  GlobalOrdinal nEle       = 64;
  const RCP<const Map> map = MapFactory::Build(lib, nEle, 0, comm);
  Teuchos::ParameterList matrixParameters;
  matrixParameters.set("nx", 8);
  matrixParameters.set("ny", 8);

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Laplace2D", map, matrixParameters);
  RCP<Matrix> A = Pr->BuildMatrix();

  // TODO coords should use SC = double
  RCP<RealValuedMultiVector> coords = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("2D", A->getRowMap(), matrixParameters);

  // build hierarchy
  RCP<Level> l = rcp(new Level());
  l->SetLevelID(0);
  l->SetComm(comm);
  l->Set("A", A);
  l->Set("Coordinates", coords);

  Teuchos::ArrayRCP<LocalOrdinal> dofPresent(A->getRowMap()->getLocalNumElements(), 1);
  l->Set<Teuchos::ArrayRCP<LocalOrdinal> >("DofPresent", dofPresent);

  VariableDofLaplacianFactory lapFact;

  l->Request("A", &lapFact);

  RCP<Matrix> lapA = l->Get<RCP<Matrix> >("A", &lapFact);

  // lapA->describe(out, Teuchos::VERB_EXTREME);

  TEST_EQUALITY(lapA->getRowMap()->isSameAs(*A->getRowMap()), true);

  Teuchos::RCP<Vector> oneVec = VectorFactory::Build(A->getRowMap());
  Teuchos::RCP<Vector> res    = VectorFactory::Build(A->getRowMap());
  oneVec->putScalar(Teuchos::ScalarTraits<Scalar>::one());
  res->putScalar(Teuchos::as<Scalar>(27) * Teuchos::ScalarTraits<Scalar>::one());
  Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(lapA)->apply(*oneVec, *res);
  TEST_COMPARE(res->normInf(), <, 100 * TMT::eps());
}  // VarLaplConstructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(VariableDofLaplacianFactory, VarLaplConstructor2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  if (!TYPE_EQUAL(GO, int)) {
    out << "Skipping test for GO != int" << std::endl;
    return;
  }
  out << "version: " << MueLu::Version() << std::endl;

  using TST                   = Teuchos::ScalarTraits<SC>;
  using magnitude_type        = typename TST::magnitudeType;
  using TMT                   = Teuchos::ScalarTraits<magnitude_type>;
  using real_type             = typename TST::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // TAW 04/21: test is crashing on 4 procs with Epetra due to an unknown reason in the Epetra_BlockMap constructor (MPI communication)
  if (comm->getSize() > 2 && lib == Xpetra::UseEpetra) {
    out << "Skipping test for more than 2 procs when using Epetra" << std::endl;
    return;
  }

  GlobalOrdinal nx = 6, ny = 6;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Map> nodeMap = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(lib, "Cartesian2D", comm, galeriList);

  // build coordinates before expanding map (nodal coordinates, not dof-based)
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);
  RCP<const Map> dofMap                  = MapFactory::Build(nodeMap, 2);  // expand map for 2 DOFs per node

  galeriList.set("right boundary", "Neumann");
  galeriList.set("bottom boundary", "Neumann");
  galeriList.set("top boundary", "Neumann");
  galeriList.set("front boundary", "Neumann");
  galeriList.set("back boundary", "Neumann");
  galeriList.set("keepBCs", false);

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", dofMap, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(2);

  //->describe(out, Teuchos::VERB_EXTREME);

  TEST_EQUALITY(dofMap->getLocalNumElements(), 2 * nodeMap->getLocalNumElements());

  // build hierarchy
  RCP<Level> l = rcp(new Level());
  l->SetLevelID(0);
  l->SetComm(comm);
  l->Set("A", A);
  l->Set("Coordinates", coordinates);

  Teuchos::ArrayRCP<LocalOrdinal> dofPresent(A->getRowMap()->getLocalNumElements(), 1);
  l->Set<Teuchos::ArrayRCP<LocalOrdinal> >("DofPresent", dofPresent);

  // A->getColMap()->describe(out,Teuchos::VERB_EXTREME);

  VariableDofLaplacianFactory lapFact;
  lapFact.SetParameter("maxDofPerNode", Teuchos::ParameterEntry(2));
  l->Request("A", &lapFact);

  RCP<Matrix> lapA = l->Get<RCP<Matrix> >("A", &lapFact);

  // lapA->describe(out, Teuchos::VERB_EXTREME);

  Teuchos::RCP<Vector> oneVec = VectorFactory::Build(lapA->getRowMap());
  Teuchos::RCP<Vector> res    = VectorFactory::Build(lapA->getRowMap());
  oneVec->putScalar(Teuchos::ScalarTraits<Scalar>::one());
  res->putScalar(Teuchos::as<Scalar>(27) * Teuchos::ScalarTraits<Scalar>::one());
  Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(lapA)->apply(*oneVec, *res);
  TEST_COMPARE(res->normInf(), <, 100 * TMT::eps());

  Teuchos::ArrayRCP<LocalOrdinal> dofPresent2(3 * lapA->getRowMap()->getLocalNumElements(), 1);
  for (decltype(dofPresent2.size()) i = 2; i < dofPresent2.size(); i = i + 3) {
    dofPresent2[i] = 0;
  }
  l->Set<Teuchos::ArrayRCP<LocalOrdinal> >("DofPresent", dofPresent2);

  // A->getColMap()->describe(out,Teuchos::VERB_EXTREME);

  VariableDofLaplacianFactory lapFact2;
  lapFact2.SetParameter("maxDofPerNode", Teuchos::ParameterEntry(3));
  l->Request("A", &lapFact2);

  RCP<Matrix> lapA2 = l->Get<RCP<Matrix> >("A", &lapFact2);

  Teuchos::RCP<Vector> oneVec2 = VectorFactory::Build(lapA2->getRowMap());
  Teuchos::RCP<Vector> res2    = VectorFactory::Build(lapA2->getRowMap());
  oneVec2->putScalar(Teuchos::ScalarTraits<Scalar>::one());
  res2->putScalar(Teuchos::as<Scalar>(27) * Teuchos::ScalarTraits<Scalar>::one());
  Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(lapA2)->apply(*oneVec2, *res2);
  TEST_COMPARE(res2->normInf(), <, 100 * TMT::eps());
  TEST_EQUALITY(res->getMap()->isSameAs(*(res2->getMap())), true);

  // lapA2->describe(out, Teuchos::VERB_EXTREME);
}  // VarLaplConstructor2

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(VariableDofLaplacianFactory, VarLaplPtent, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  if (!TYPE_EQUAL(GO, int)) {
    out << "Skipping test for GO != int" << std::endl;
    return;
  }
  out << "version: " << MueLu::Version() << std::endl;

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // TAW 04/21: test is crashing on 4 procs with Epetra due to an unknown reason in the Epetra_BlockMap constructor (MPI communication)
  if (comm->getSize() > 2 && lib == Xpetra::UseEpetra) {
    out << "Skipping test for more than 2 procs when using Epetra" << std::endl;
    return;
  }

  GlobalOrdinal nx = 6, ny = 6;

  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> mv_type_double;
  typedef Xpetra::MultiVectorFactory<double, LocalOrdinal, GlobalOrdinal, Node> MVFactory_double;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Map> nodeMap = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(lib, "Cartesian2D", comm, galeriList);

  // build coordinates before expanding map (nodal coordinates, not dof-based)
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double, LocalOrdinal, GlobalOrdinal, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);
  RCP<const Map> dofMap                  = MapFactory::Build(nodeMap, 2);  // expand map for 2 DOFs per node

  galeriList.set("right boundary", "Neumann");
  galeriList.set("bottom boundary", "Neumann");
  galeriList.set("top boundary", "Neumann");
  galeriList.set("front boundary", "Neumann");
  galeriList.set("back boundary", "Neumann");
  galeriList.set("keepBCs", false);

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", dofMap, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(2);

  //->describe(out, Teuchos::VERB_EXTREME);

  TEST_EQUALITY(dofMap->getLocalNumElements(), 2 * nodeMap->getLocalNumElements());

  Teuchos::ArrayRCP<LocalOrdinal> dofPresent(A->getRowMap()->getLocalNumElements(), 1);

  // build hierarchy
  typedef Teuchos::ScalarTraits<Scalar> TST;
  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;

  // generate laplacian matrix using level l
  Level l;
  l.Set("A", A);
  l.Set("Coordinates", coordinates);
  l.Set<Teuchos::ArrayRCP<LocalOrdinal> >("DofPresent", dofPresent);

  VariableDofLaplacianFactory lapFact;
  lapFact.SetParameter("maxDofPerNode", Teuchos::ParameterEntry(2));
  l.Request("A", &lapFact);

  RCP<Matrix> lapA = l.Get<RCP<Matrix> >("A", &lapFact);
  lapA->describe(out, Teuchos::VERB_EXTREME);

  // create a new two-level hierarchy and set lapA as input "A"
  // This way we avoid dealing with factory managers
  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);

  fineLevel.Set("Coordinates", coordinates);
  fineLevel.Set("A", lapA);

  TentativePFactory PFact;

  coarseLevel.Request("P", &PFact);
  RCP<Matrix> pMat = coarseLevel.Get<RCP<Matrix> >("P", &PFact);

  TEST_EQUALITY(pMat->getRowMap()->isSameAs(*(lapA->getRowMap())), true);
  TEST_EQUALITY(pMat->getLocalNumEntries(), pMat->getRowMap()->getLocalNumElements());
}  // VarLaplPtent

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(VariableDofLaplacianFactory, VarLaplConstructor, SC, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(VariableDofLaplacianFactory, VarLaplConstructor2, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(VariableDofLaplacianFactory, VarLaplPtent, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests

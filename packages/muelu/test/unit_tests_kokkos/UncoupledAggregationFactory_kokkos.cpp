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
// // Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
//#include <Galeri_XpetraMaps.hpp>

#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory_kokkos.hpp>
#include <MueLu_config.hpp>
#include <MueLu_TestHelpers_kokkos.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_UncoupledAggregationFactory_kokkos.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory_kokkos, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>

  out << "version: " << MueLu::Version() << std::endl;
  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  TEST_EQUALITY(aggFact != Teuchos::null, true);
  TEST_THROW(aggFact->SetOrdering("unknown_ordering"), Teuchos::Exceptions::InvalidParameterValue);
  aggFact->SetOrdering("natural");
  TEST_EQUALITY(aggFact->GetOrdering() == "natural", true);
  aggFact->SetOrdering("graph");
  TEST_EQUALITY(aggFact->GetOrdering() == "graph", true);
  aggFact->SetOrdering("random");
  TEST_EQUALITY(aggFact->GetOrdering() == "random", true);

  aggFact->SetMaxNeighAlreadySelected(100);
  TEST_EQUALITY(aggFact->GetMaxNeighAlreadySelected() == 100, true);
  aggFact->SetMinNodesPerAggregate(10);
  TEST_EQUALITY(aggFact->GetMinNodesPerAggregate() == 10, true);

  auto validParamList = aggFact->GetValidParameterList();
  TEST_EQUALITY(validParamList != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory_kokkos, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();

  // Make a Matrix with multiple degrees of freedom per node
  GlobalOrdinal nx = 20, ny = 20;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);

  const GO nxx       = 200;
  using test_factory = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;
  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                    // manage parameters of xpetra

  RCP<Matrix> A      = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nxx);
  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node
  A->SetFixedBlockSize(2);

  Level aLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);
  aLevel.Request("A");
  aLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);

  aLevel.Request(*aggFact);

  aggFact->Build(aLevel);

}  // Build

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory_kokkos, Build_MIS2_Coarsen, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();

  // Make a Matrix with multiple degrees of freedom per node
  GlobalOrdinal nx = 20, ny = 20;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);

  const GO nxx       = 200;
  using test_factory = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;
  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                    // manage parameters of xpetra

  RCP<Matrix> A      = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nxx);
  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node
  A->SetFixedBlockSize(2);

  Level aLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);
  aLevel.Request("A");
  aLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: coloring algorithm",
                        Teuchos::ParameterEntry(std::string("mis2 coarsening")));

  aLevel.Request(*aggFact);

  aggFact->Build(aLevel);

}  // Build_MIS2_Coarsen

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory_kokkos, Build_MIS2_Aggregate, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();

  // Make a Matrix with multiple degrees of freedom per node
  GlobalOrdinal nx = 20, ny = 20;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);

  const GO nxx       = 200;
  using test_factory = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;
  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                    // manage parameters of xpetra

  RCP<Matrix> A      = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nxx);
  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node
  A->SetFixedBlockSize(2);

  Level aLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);
  aLevel.Request("A");
  aLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: coloring algorithm",
                        Teuchos::ParameterEntry(std::string("mis2 aggregation")));

  aLevel.Request(*aggFact);

  aggFact->Build(aLevel);

}  // Build_MIS2_Aggregate

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory_kokkos, Constructor, Scalar, LO, GO, Node)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory_kokkos, Build, Scalar, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory_kokkos, Build_MIS2_Coarsen, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory_kokkos, Build_MIS2_Aggregate, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

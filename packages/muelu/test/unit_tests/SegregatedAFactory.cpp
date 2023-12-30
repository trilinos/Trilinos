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
#include <Teuchos_FancyOStream.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_IO.hpp>

#include <MueLu_SegregatedAFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SegregatedAFactory, Basic, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  using TST            = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename TST::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  // test with simple example matrix and segregate first 4 rows from rest
  //     x x x x                x x x x
  //     x x x x                x x x x
  //     x x x x x x            x x x x
  //     x x x x x x      ==>   x x x x
  //         x x x x x                 x x x
  //         x x x x x                 x x x
  //             x x x                 x x x
  {
    RCP<Map> rowMap     = MapFactory::Build(lib, 7, 0, comm);
    RCP<CrsGraph> graph = CrsGraphFactory::Build(rowMap, 6);

    graph->insertGlobalIndices(0, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3));
    graph->insertGlobalIndices(1, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3));
    graph->insertGlobalIndices(2, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(3, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(4, Teuchos::tuple<GlobalOrdinal>(2, 3, 4, 5, 6));
    graph->insertGlobalIndices(5, Teuchos::tuple<GlobalOrdinal>(2, 3, 4, 5, 6));
    graph->insertGlobalIndices(6, Teuchos::tuple<GlobalOrdinal>(4, 5, 6));
    graph->fillComplete();

    RCP<Matrix> A = MatrixFactory::Build(graph.getConst());
    A->setAllToScalar(TST::one());
    A->fillComplete();

    TEST_EQUALITY(A.is_null(), false);
    TEST_EQUALITY(A->getGlobalNumEntries(), 33);
    TEST_FLOATING_EQUALITY(A->getFrobeniusNorm(), 5.744562646538029, 2e1 * TMT::eps());

    RCP<Map> map = MapFactory::Build(lib, 4, 0, comm);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    const std::string mapName = "Segregate Map";
    level.Set(mapName, map);

    TEST_ASSERT(level.IsAvailable("A", MueLu::NoFactory::get()));
    TEST_ASSERT(level.IsAvailable(mapName, MueLu::NoFactory::get()));

    RCP<SegregatedAFactory> segregateFact = rcp(new SegregatedAFactory());
    segregateFact->SetFactory("A", MueLu::NoFactory::getRCP());
    segregateFact->SetParameter("map: name", Teuchos::ParameterEntry(mapName));

    // request segregated operator
    level.Request("A", segregateFact.get());
    segregateFact->Build(level);

    RCP<Matrix> Aout = level.Get<RCP<Matrix> >("A", segregateFact.get());

    // Output (segregated Operator)
    // Node ID   Global Row  Num Entries(Index,Value)
    // 0         0           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         1           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         2           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         3           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         4           3 (4, 1)  (5, 1)  (6, 1)
    // 0         5           3 (4, 1)  (5, 1)  (6, 1)
    // 0         6           3 (4, 1)  (5, 1)  (6, 1)

    TEST_EQUALITY(Aout.is_null(), false);
    TEST_EQUALITY(Aout->getGlobalNumEntries(), 25);
    TEST_FLOATING_EQUALITY(Aout->getFrobeniusNorm(), 5.0, 2e1 * TMT::eps());
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SegregatedAFactory, Basic, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests

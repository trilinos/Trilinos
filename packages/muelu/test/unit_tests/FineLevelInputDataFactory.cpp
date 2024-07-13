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

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_FineLevelInputDataFactory.hpp"

namespace MueLuTests {

/////////////////////////
// helper function

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
GenerateProblemMatrix(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map,
                      Scalar a = 2.0, Scalar b = -1.0, Scalar c = -1.0) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::RCP<CrsMatrixWrap> mtx = Galeri::Xpetra::MatrixTraits<Map, CrsMatrixWrap>::Build(map, 3);

  LocalOrdinal NumMyElements                               = map->getLocalNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();
  GlobalOrdinal NumGlobalElements                          = map->getGlobalNumElements();
  GlobalOrdinal nIndexBase                                 = map->getIndexBase();

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz = 2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  for (LocalOrdinal i = 0; i < NumMyElements; ++i) {
    if (MyGlobalElements[i] == nIndexBase) {
      // off-diagonal for first row
      Indices[0] = nIndexBase;
      NumEntries = 1;
      Values[0]  = c;
    } else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1) {
      // off-diagonal for last row
      Indices[0] = nIndexBase + NumGlobalElements - 2;
      NumEntries = 1;
      Values[0]  = b;
    } else {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1]  = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0]  = c;
      NumEntries = 2;
    }

    // put the off-diagonal entries
    // Xpetra wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0], NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0], NumEntries);
    mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    mtx->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(a));

  }  // for (LocalOrdinal i = 0; i < NumMyElements; ++i)

  mtx->fillComplete(map, map);

  return mtx;
}

// dummy helper class (placeholder)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class DummyFactory : public MueLu::SingleLevelFactoryBase {
 public:
  DummyFactory() {}
  virtual ~DummyFactory() {}
  RCP<const Teuchos::ParameterList> GetValidParameterList() const { return Teuchos::null; };
  void DeclareInput(MueLu::Level& currentLevel) const {};
  void Build(MueLu::Level& currentLevel) const {};
};  // class DummyFactory

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FineLevelInputDataFactory, InputData, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  RCP<const Map> map;
  GO numElements = 100;

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  map = MapFactory::Build(lib, numElements, 0, comm);

  RCP<Matrix> A = GenerateProblemMatrix<Scalar, LO, GO, Node>(map, 2, -1, -1);

  // build hierarchy
  RCP<Level> levelOne = rcp(new Level());
  levelOne->SetLevelID(0);

  FineLevelInputDataFactory inputData;

#ifdef HAVE_MUELU_DEBUG
  inputData.DisableMultipleCallCheck();
#endif

  // Test 1: different fine level variable name and coarse level variable name, default fine level factory (=NoFactory)
  levelOne->Set("A-output", A);
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("A-output")));

  levelOne->Request("A-output", &inputData);
  RCP<Matrix> AA = levelOne->Get<RCP<Matrix> >("A-output", &inputData);

  TEST_EQUALITY(AA.get(), A.get());

  levelOne->Release("A-output", &inputData);

  // Test 2: set fine level factory
  MueLuTests::DummyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> dummyFact;
  levelOne->Set("A-output2", A);
  inputData.SetFactory("Fine level factory", MueLu::NoFactory::getRCP());
  inputData.SetFactory("Coarse level factory", Teuchos::rcpFromRef(dummyFact));
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("A-output2")));

  levelOne->Request("A-output2", &inputData);
  AA = levelOne->Get<RCP<Matrix> >("A-output2", &inputData);
  TEST_EQUALITY(AA.get(), A.get());
  levelOne->Release("A-output2", &inputData);

  // Test 3: same as test 2 on coarse level
  levelOne->SetLevelID(1);
  levelOne->Request("A-output4", &dummyFact);
  levelOne->Set("A-output4", A, &dummyFact);

  inputData.SetFactory("Coarse level factory", Teuchos::rcpFromRef(dummyFact));
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("A-output4")));

  levelOne->Request("A-output4", &inputData);

  AA = levelOne->Get<RCP<Matrix> >("A-output4", &inputData);

  TEST_EQUALITY(AA.get(), A.get());

  // Test 4: same as test 2 on coarse level
  levelOne->Request("A-output5", &dummyFact);
  levelOne->Set("A-output5", A, &dummyFact);

  inputData.SetFactory("Coarse level factory", Teuchos::rcpFromRef(dummyFact));
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("A-output5")));
  levelOne->Request("A-output5", &inputData);
  AA = levelOne->Get<RCP<Matrix> >("A-output5", &inputData);
  TEST_EQUALITY(AA.get(), A.get());
}  // InputData

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FineLevelInputDataFactory, InputDataMap, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  RCP<Map> map = MapFactory::Build(lib, 100, 0, comm);

  // build hierarchy
  RCP<Level> levelOne = rcp(new Level());
  levelOne->SetLevelID(0);

  FineLevelInputDataFactory inputData;

#ifdef HAVE_MUELU_DEBUG
  inputData.DisableMultipleCallCheck();
#endif

  // Test 1: different fine level variable name and coarse level variable name, default fine level factory (=NoFactory)
  levelOne->Set("m", map);
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("m")));

  levelOne->Request("m", &inputData);

  TEST_THROW(levelOne->Get<RCP<Map> >("m", &inputData), Teuchos::bad_any_cast);

  inputData.SetParameter("Variable type", Teuchos::ParameterEntry(std::string("Map")));
  RCP<Map> mm = levelOne->Get<RCP<Map> >("m", &inputData);
  TEST_EQUALITY(mm.get(), map.get());

  levelOne->Release("m", &inputData);

  // Test 2: set fine level factory
  MueLuTests::DummyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> dummyFact;
  levelOne->Set("m", map);
  inputData.SetFactory("Fine level factory", MueLu::NoFactory::getRCP());
  inputData.SetFactory("Coarse level factory", Teuchos::rcpFromRef(dummyFact));
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("m")));

  levelOne->Request("m", &inputData);
  mm = levelOne->Get<RCP<Map> >("m", &inputData);
  TEST_EQUALITY(mm.get(), map.get());
  levelOne->Release("m", &inputData);

  // Test 3: same as test 2 on coarse level
  levelOne->SetLevelID(1);
  levelOne->Request("m", &dummyFact);
  levelOne->Set("m", map, &dummyFact);

  inputData.SetFactory("Coarse level factory", Teuchos::rcpFromRef(dummyFact));
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("m")));

  levelOne->Request("m", &inputData);

  mm = levelOne->Get<RCP<Map> >("m", &inputData);

  TEST_EQUALITY(mm.get(), map.get());

  // Test 4: same as test 2 on coarse level
  levelOne->Request("m", &dummyFact);
  levelOne->Set("m", map, &dummyFact);

  inputData.SetFactory("Coarse level factory", Teuchos::rcpFromRef(dummyFact));
  inputData.SetParameter("Variable", Teuchos::ParameterEntry(std::string("m")));
  levelOne->Request("m", &inputData);
  mm = levelOne->Get<RCP<Map> >("m", &inputData);
  TEST_EQUALITY(mm.get(), map.get());
}  // InputDataMap

// helper class for testing functionality of FineLevelInputDataFactory
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class FineLevelInputDataFactoryTester {
#include "MueLu_UseShortNames.hpp"
 public:
  void test_testfunc(const FineLevelInputDataFactory& fac) {
    // std::cout << "FineLevelInputDataFactoryTester" << std::endl;
    fac.test();
  }
};

// unit test just for demonstration purposes
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(FineLevelInputDataFactory, TestFunc, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  FineLevelInputDataFactory fac;

  FineLevelInputDataFactoryTester<Scalar, LocalOrdinal, GlobalOrdinal, Node> tester;

  tester.test_testfunc(fac);

  // TEST_EQUALITY(AA.get(), A.get());
}  // TestFunc

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FineLevelInputDataFactory, InputData, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FineLevelInputDataFactory, InputDataMap, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(FineLevelInputDataFactory, TestFunc, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests

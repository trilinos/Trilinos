// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Galeri_XpetraUtils.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_MLParameterListInterpreter.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MLParameterListInterpreter, SetParameterList, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  if (!TYPE_EQUAL(SC, double)) {
    out << "Skipping for SC != double" << std::endl;
    return;
  }
  out << "version: " << MueLu::Version() << std::endl;

  // TODO: this test can be done at compilation time
#if !defined(HAVE_MUELU_EPETRA) or !defined(HAVE_MUELU_IFPACK) or !defined(HAVE_MUELU_AMESOS)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseEpetra, "Epetra, Ifpack, Amesos");
#endif

#if !defined(HAVE_MUELU_IFPACK2) or !defined(HAVE_MUELU_AMESOS2)
  MUELU_TESTING_DO_NOT_TEST(Xpetra::UseTpetra, "Tpetra, Ifpack2, Amesos2");
#endif

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(99);
  Teuchos::ParameterList galeriParameters;
  galeriParameters.set<GO>("nx", 99);
  RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", A->getRowMap(), galeriParameters);

  ArrayRCP<std::string> fileList = TestHelpers::GetFileList(std::string("ParameterList/MLParameterListInterpreter/"), std::string(".xml"));

  for (int i = 0; i < fileList.size(); i++) {
    out << "Processing file: " << fileList[i] << std::endl;
    Teuchos::ParameterList myList;
    myList.set("xml parameter file", "ParameterList/MLParameterListInterpreter/" + fileList[i]);

    Teuchos::ArrayRCP<typename MultiVector::scalar_type> xcoord = coordinates->getDataNonConst(0);
    myList.set("x-coordinates", xcoord.get());

    MLParameterListInterpreter mueluFactory(myList, A->getRowMap()->getComm());

    RCP<Hierarchy> H = mueluFactory.CreateHierarchy();
    H->GetLevel(0)->template Set<RCP<Matrix> >("A", A);

    mueluFactory.SetupHierarchy(*H);

    // TODO: check no unused parameters
    // TODO: check results of Iterate()
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MLParameterListInterpreter, SetParameterList, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

// TODO: some tests of the parameter list parser can be done without building the Hierarchy.

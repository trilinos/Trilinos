// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <MueLu_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#ifdef HAVE_MUELU_ML
#include <ml_epetra_utils.h>
#endif

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_ParameterListUtils.hpp"

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

TEUCHOS_UNIT_TEST(MueLu_CreateSublists, SetParameterList) {
  std::string dir("ParameterList/CreateSublists/");

  ArrayRCP<std::string> fileList = TestHelpers::GetFileList(dir, std::string(".xml"));

  for (int i = 0; i < fileList.size(); i++) {
    out << "Processing file: " << fileList[i] << std::endl;

    Teuchos::RCP<const Teuchos::ParameterList> inputList = Teuchos::getParametersFromXmlFile(dir + fileList[i]);
    Teuchos::ParameterList outputList;

    MueLu::CreateSublists(*inputList, outputList);

    // Test against reference output (replace '.xml' by '.output' to get the filename)
    Teuchos::RCP<Teuchos::ParameterList> refOutputList = Teuchos::getParametersFromXmlFile(dir + fileList[i].substr(0, fileList[i].find_last_of(".")) + ".output");
    TEST_EQUALITY(outputList, *refOutputList);
  }
}

#ifdef HAVE_MUELU_ML
TEUCHOS_UNIT_TEST(ML_CreateSublists, SetParameterList) {
  std::string dir("ParameterList/CreateSublists/");

  ArrayRCP<std::string> fileList = TestHelpers::GetFileList(dir, std::string(".xml"));

  for (int i = 0; i < fileList.size(); i++) {
    out << "Processing file: " << fileList[i] << std::endl;

    Teuchos::RCP<const Teuchos::ParameterList> inputList = Teuchos::getParametersFromXmlFile(dir + fileList[i]);
    Teuchos::ParameterList outputList;

    ML_CreateSublists(*inputList, outputList);

    // Test against reference output (replace '.xml' by '.output' to get the filename)
    Teuchos::RCP<Teuchos::ParameterList> refOutputList = Teuchos::getParametersFromXmlFile(dir + fileList[i].substr(0, fileList[i].find_last_of(".")) + ".output");
    TEST_EQUALITY(outputList, *refOutputList);
  }
}
#endif

}  // namespace MueLuTests

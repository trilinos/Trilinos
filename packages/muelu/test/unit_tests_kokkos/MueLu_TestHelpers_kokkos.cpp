// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestRepository.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"

namespace MueLuTests {

// static members initialization of the class TestHelpers::Parameters
Xpetra::Parameters TestHelpers_kokkos::Parameters::xpetraParameters = Xpetra::Parameters(Teuchos::UnitTestRepository::getCLP());

}  // namespace MueLuTests

namespace MueLuTests {
namespace TestHelpers_kokkos {

ArrayRCP<std::string> GetFileList(const std::string& dirPath, const std::string& filter) {
  RCP<std::vector<std::string> > files = rcp(new std::vector<std::string>());

  DIR* dir = opendir(dirPath.c_str());
  TEUCHOS_TEST_FOR_EXCEPTION(dir == NULL, MueLu::Exceptions::RuntimeError, "GetFileList(\"" + dirPath + "\") : " + strerror(errno));

  struct dirent* dirEntry;
  while ((dirEntry = readdir(dir)) != NULL) {
    std::string dirEntryS(dirEntry->d_name);

    size_t pos = dirEntryS.rfind(filter);
    if (pos != std::string::npos && pos + filter.size() == dirEntryS.size())
      files->push_back(std::string(dirEntryS));
  }

  closedir(dir);

  return arcp(files);
}

}  // namespace TestHelpers_kokkos
}  // namespace MueLuTests

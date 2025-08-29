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
#include <filesystem>

namespace MueLuTests {

// static members initialization of the class TestHelpers::Parameters
Xpetra::Parameters TestHelpers_kokkos::Parameters::xpetraParameters = Xpetra::Parameters(Teuchos::UnitTestRepository::getCLP());

}  // namespace MueLuTests

namespace MueLuTests::TestHelpers_kokkos {

// Use <filesystem> to list files in the directory matching the filter
ArrayRCP<std::string> GetFileList(const std::string& dirPath, const std::string& filter) {
  RCP<std::vector<std::string> > files = rcp(new std::vector<std::string>());

  // Use filesystem to iterate over the directory
  try {
    for (const auto& entry : std::filesystem::directory_iterator(dirPath)) {
      if (entry.is_regular_file()) {
        std::string fileName = entry.path().filename().string();
        // Apply the filter
        size_t pos = fileName.rfind(filter);
        if (pos != std::string::npos && pos + filter.size() == fileName.size()) {
          files->push_back(fileName);
        }
      }
    }
  } catch (const std::filesystem::filesystem_error& e) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
                               "GetFileList(\"" + dirPath + "\") : " + e.what());
  }

  return arcp(files);
}

}  // namespace MueLuTests::TestHelpers_kokkos

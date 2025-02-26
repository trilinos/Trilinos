// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestRepository.hpp>
#include <filesystem>
#include "MueLu_TestHelpers.hpp"

namespace MueLuTests {

// static members initialization of the class TestHelpers::Parameters
Xpetra::Parameters TestHelpers::Parameters::xpetraParameters = Xpetra::Parameters(Teuchos::UnitTestRepository::getCLP());

}  // namespace MueLuTests

namespace MueLuTests::TestHelpers {

ArrayRCP<std::string> GetFileList(const std::string& dirPath, const std::string& filter) {
  namespace fs                        = std::filesystem;
  RCP<std::vector<std::string>> files = rcp(new std::vector<std::string>());

  try {
    if (!fs::exists(dirPath) || !fs::is_directory(dirPath)) {
      throw MueLu::Exceptions::RuntimeError("GetFileList(\"" + dirPath + "\") : Directory does not exist or is not a directory.");
    }

    for (const auto& entry : fs::directory_iterator(dirPath)) {
      const auto& path = entry.path();
      if (path.has_filename()) {
        std::string filename = path.filename().string();
        size_t pos           = filename.rfind(filter);
        if (pos != std::string::npos && pos + filter.size() == filename.size()) {
          files->push_back(filename);
        }
      }
    }
  } catch (const fs::filesystem_error& e) {
    throw MueLu::Exceptions::RuntimeError("GetFileList(\"" + dirPath + "\") : " + e.what());
  }

  return arcp(files);
}

}  // namespace MueLuTests::TestHelpers

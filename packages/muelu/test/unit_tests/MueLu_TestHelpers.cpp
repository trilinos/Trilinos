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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestRepository.hpp>

#include "MueLu_TestHelpers.hpp"

namespace MueLuTests {

  // static members initialization of the class TestHelpers::Parameters
  Xpetra::Parameters TestHelpers::Parameters::xpetraParameters = Xpetra::Parameters(Teuchos::UnitTestRepository::getCLP());

}

namespace MueLuTests {
  namespace TestHelpers {

  ArrayRCP<std::string> GetFileList(const std::string & dirPath, const std::string & filter) {

    RCP<std::vector<std::string> > files = rcp(new std::vector<std::string>());

    DIR *dir;
    struct dirent *dirEntry;

    dir = opendir(dirPath.c_str());
    TEUCHOS_TEST_FOR_EXCEPTION(dir == NULL, MueLu::Exceptions::RuntimeError, "GetFileList(\"" + dirPath + "\") : " + strerror(errno));

    while ((dirEntry = readdir(dir)) != NULL) {
      std::string dirEntryS(dirEntry->d_name);

      if (dirEntryS.rfind(filter) != string::npos)
        files->push_back(std::string(dirEntryS));
    }

    closedir(dir);

    return arcp(files);
  }

  }
}

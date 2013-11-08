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

#include <cstdlib>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <MueLu.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_TestHelpers.hpp>

#include <MueLu_EasyParameterListInterpreter.hpp>

// These files must be included last
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // Set seed
  std::srand(12345);

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int numProc = comm->getSize();
  int myRank  = comm->getRank();

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  Teuchos::CommandLineProcessor clp(false);
  ::Xpetra::Parameters xpetraParameters(clp);

  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  // =========================================================================
  // Problem construction
  // =========================================================================
  RCP<Matrix> A = MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::Build1DPoisson(9999, lib);

  Teuchos::ArrayRCP<std::string> fileList;
  if (numProc == 1) {
    // Run all xml configs in serial/single mpi mode
    fileList = MueLuTests::TestHelpers::GetFileList(std::string("EasyParameterListInterpreter/"), std::string(".xml"));
  } else {
    // In addition, rerun some files in parallel mode
    fileList = MueLuTests::TestHelpers::GetFileList(std::string("EasyParameterListInterpreter/"), std::string("_np" + Teuchos::toString(numProc) + ".xml"));
  }

  bool failed = false;
  for (int i = 0; i < fileList.size(); i++) {
    std::string xmlFile = "EasyParameterListInterpreter/" + fileList[i];
    std::string baseFile = xmlFile.substr(0, xmlFile.find_last_of('.')) + (lib == Xpetra::UseEpetra ? "_epetra" : "_tpetra");

    // Redirect output
    std::filebuf buffer;
    buffer.open((baseFile + ".out").c_str(), std::ios::out);
    std::streambuf* oldbuffer = std::cout.rdbuf(&buffer);

    // NOTE: we cannot use EasyParameterListInterpreter(xmlFile, comm), because we want to update the ParameterList
    // first to include verbosity
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFile, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);
    paramList.set("verbosity", "test");

    try {
      EasyParameterListInterpreter mueluFactory(paramList);

      RCP<Hierarchy> H = mueluFactory.CreateHierarchy();

      H->GetLevel(0)->Set<RCP<Matrix> >("A", A);
      // H->GetLevel(0)->Set("Nullspace",   nullspace);
      // H->GetLevel(0)->Set("Coordinates", coordinates);

      mueluFactory.SetupHierarchy(*H);

    } catch (Teuchos::ExceptionBase& e) {
      if (!myRank) {
        std::string msg = e.what();
        msg = msg.substr(msg.find_last_of('\n')+1);
        std::cout << "Caught exception: " << msg << std::endl;
      }
    }

    // Redirect output back
    std::cout.rdbuf(oldbuffer);
    buffer.close();

    {
      // Tpetra produces different eigenvalues in Chebyshev due to using
      // std::rand() for generating random vectors, which may be initialized
      // using different seed, and may have different algorithm from one
      // gcc version to another, or to anogther compiler (like clang)
      // This leads to us always failing this test.
      // NOTE: Epetra, on the other hand, rolls out its out random number
      // generator, which always produces same results
      std::string sed_cmd = "sed -i 's/lambdaMax\\ =\\ [0-9]*\\.[0-9]*/lambdaMax\\ =\\ <ignored>/' ";
      system((sed_cmd + baseFile + ".res").c_str());
      system((sed_cmd + baseFile + ".out").c_str());

      sed_cmd = "sed -i 's/lambdaMin\\ =\\ [0-9]*\\.[0-9]*/lambdaMin\\ =\\ <ignored>/' ";
      system((sed_cmd + baseFile + ".res").c_str());
      system((sed_cmd + baseFile + ".out").c_str());
    }

    // Run comparison
    std::string cmd = "diff -w " + baseFile + ".res " + baseFile + ".out";
    int ret = system(cmd.c_str());
    if (ret)
      failed = true;
    if (!myRank)
      std::cout << xmlFile << ": " << (ret ? "failed" : "passed") << std::endl;
  }

  if (comm->getRank() == 0)
    std::cout << std::endl << "End Result: TEST " << (failed ? "FAILED" : "PASSED") << std::endl;

  return 0;
}

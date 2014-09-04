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

#include <MueLu_ParameterListInterpreter.hpp>

// These files must be included last
#include <MueLu_UseDefaultTypes.hpp>

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

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
  ParameterList matrixParameters;
  matrixParameters.set("nx",         9999);
  matrixParameters.set("matrixType", "Laplace1D");
  RCP<Matrix>      A           = MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(matrixParameters.get<int>("nx"), lib);
  RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", A->getRowMap(), matrixParameters);

  const int numLists = 2;

  std::string dirList[numLists];
  dirList[0] = "EasyParameterListInterpreter/";
  dirList[1] = "FactoryParameterListInterpreter/";

  Teuchos::ArrayRCP<std::string> fileLists[numLists];
  if (numProc == 1) {
    // Run all xml configs in serial/single mpi mode
    fileLists[0] = MueLuTests::TestHelpers::GetFileList(dirList[0], std::string(".xml"));
    fileLists[1] = MueLuTests::TestHelpers::GetFileList(dirList[1], std::string(".xml"));
  } else {
    // In addition, rerun some files in parallel mode
    fileLists[0] = MueLuTests::TestHelpers::GetFileList(dirList[0], std::string("_np" + Teuchos::toString(numProc) + ".xml"));
    fileLists[1] = MueLuTests::TestHelpers::GetFileList(dirList[1], std::string("_np" + Teuchos::toString(numProc) + ".xml"));
  }

  bool failed = false;
  for (int k = 0; k < numLists; k++) {
    Teuchos::ArrayRCP<std::string> fileList = fileLists[k];

    for (int i = 0; i < fileList.size(); i++) {
      // Set seed
      std::srand(12345);

      // Reset (potentially) cached value of the estimate
      A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());

      std::string xmlFile  = dirList[k] + fileList[i];
      std::string baseFile = xmlFile.substr(0, xmlFile.find_last_of('.'));
      std::size_t found = baseFile.find("_np");
      if (numProc == 1 && found != std::string::npos) {
#ifdef HAVE_MPI
        baseFile = baseFile.substr(0, found);
#else
        std::cout << "Skipping \"" << xmlFile << "\" as MPI is not enabled" << std::endl;
        continue;
#endif
      }
      baseFile = baseFile + (lib == Xpetra::UseEpetra ? "_epetra" : "_tpetra");

      std::filebuf buffer;
      std::streambuf* oldbuffer = NULL;
      if (myRank == 0) {
        // Redirect output
        buffer.open((baseFile + ".out").c_str(), std::ios::out);
        oldbuffer = std::cout.rdbuf(&buffer);
      }

      // NOTE: we cannot use ParameterListInterpreter(xmlFile, comm), because we want to update the ParameterList
      // first to include "test" verbosity
      Teuchos::ParameterList paramList;
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFile, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);
      paramList.set("verbosity", "test");

      try {
        ParameterListInterpreter mueluFactory(paramList);

        RCP<Hierarchy> H = mueluFactory.CreateHierarchy();

        H->GetLevel(0)->Set<RCP<Matrix> >("A", A);
        // H->GetLevel(0)->Set("Nullspace",   nullspace);
        H->GetLevel(0)->Set("Coordinates", coordinates);

        mueluFactory.SetupHierarchy(*H);

      } catch (Teuchos::ExceptionBase& e) {
        std::string msg = e.what();
        msg = msg.substr(msg.find_last_of('\n')+1);

        if (myRank == 0) {
          std::cout << "Caught exception: " << msg << std::endl;

          // Redirect output back
          std::cout.rdbuf(oldbuffer);
          buffer.close();
        }

        if (msg == "Zoltan interface is not available" ||
            msg == "Zoltan2 interface is not available") {

          if (myRank == 0)
            std::cout << xmlFile << ": skipped (missing library)" << std::endl;

          continue;
        }
      }

      if (myRank == 0) {
        // Redirect output back
        std::cout.rdbuf(oldbuffer);
        buffer.close();

        // Tpetra produces different eigenvalues in Chebyshev due to using
        // std::rand() for generating random vectors, which may be initialized
        // using different seed, and may have different algorithm from one
        // gcc version to another, or to anogther compiler (like clang)
        // This leads to us always failing this test.
        // NOTE1 : Epetra, on the other hand, rolls out its out random number
        // generator, which always produces same results
        // NOTE2 : sed behaviour differs between Mac and Linux
        // You can run "sed -i 's//' " in Linux, but you always have to specify
        // "sed -i "<smth,could be empty>" 's//'" in Mac. Both, however, take '-i<extension>'

        // Ignore the value of "lambdaMax"
        std::string sed_cmd = "sed -iorig 's/lambdaMax:\\ [0-9]*\\.[0-9]*/lambdaMax\\ =\\ <ignored>/' ";
        system((sed_cmd + baseFile + ".res").c_str());
        system((sed_cmd + baseFile + ".out").c_str());

        std::string sed_pref = "sed -i ";
#ifdef __APPLE__
        sed_pref = sed_pref +  "\"\" ";
#endif

        // Ignore the value of "lambdaMin"
        sed_cmd = sed_pref + "'s/lambdaMin:\\ [0-9]*\\.[0-9]*/lambdaMin\\ =\\ <ignored>/' ";
        system((sed_cmd + baseFile + ".res").c_str());
        system((sed_cmd + baseFile + ".out").c_str());

        // Ignore the value of "chebyshev: max eigenvalue"
        // NOTE: we skip lines with default value ([default])
        sed_cmd = sed_pref + "'/[default]/! s/chebyshev:\\ max\\ eigenvalue\\ =\\ [0-9]*\\.[0-9]*/chebyshev:\\ max\\ eigenvalue\\ =\\ <ignored>/' ";
        system((sed_cmd + baseFile + ".res").c_str());
        system((sed_cmd + baseFile + ".out").c_str());

        // Stript template args for some classes
        std::vector<std::string> classes;
        classes.push_back("Xpetra::Matrix");
        classes.push_back("MueLu::Constraint");
        for (size_t q = 0; q < classes.size(); q++) {
          sed_cmd = sed_pref + "'s/" + classes[q] + "<.*>/" + classes[q] + "<ignored> >/' ";
          system((sed_cmd + baseFile + ".res").c_str());
          system((sed_cmd + baseFile + ".out").c_str());
        }

#ifdef __APPLE__
        // Some Macs print outs ptrs as 0x0 instead of 0, fix that
        sed_cmd = sed_pref + "'/RCP/ s/=0x0/=0/g' ";
        system((sed_cmd + baseFile + ".out").c_str());
#endif

        // Run comparison
        std::string cmd = "diff -w -I\"^\\s*$\" " + baseFile + ".res " + baseFile + ".out";
        int ret = system(cmd.c_str());
        if (ret)
          failed = true;

        std::cout << xmlFile << ": " << (ret ? "failed" : "passed") << std::endl;
      }
    }
  }

  if (myRank == 0)
    std::cout << std::endl << "End Result: TEST " << (failed ? "FAILED" : "PASSED") << std::endl;

  return 0;
}

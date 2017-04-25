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
#include <fstream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <MueLu.hpp>

#if defined (HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h> // needed for check whether KLU2 is available
#endif

#include <MueLu_Exceptions.hpp>
#include <MueLu_TestHelpers.hpp>

#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_ML2MueLuParameterTranslator.hpp>

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif

void run_sed(const std::string& pattern, const std::string& baseFile);

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
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

  bool runHeavyTests = false;
  clp.setOption("heavytests", "noheavytests",  &runHeavyTests, "whether to exercise tests that take a long time to run");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  // =========================================================================
  // Problem construction
  // =========================================================================
  Teuchos::ParameterList matrixParameters;
  matrixParameters.set("nx",         Teuchos::as<GO>(9999));
  matrixParameters.set("matrixType", "Laplace1D");
  RCP<Matrix>      A           = MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(matrixParameters.get<GO>("nx"), lib);
  RCP<MultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", A->getRowMap(), matrixParameters);

  std::string outDir = "Output/";

  std::vector<std::string> dirList;
  if (runHeavyTests) {
    dirList.push_back("EasyParameterListInterpreter-heavy/");
    dirList.push_back("FactoryParameterListInterpreter-heavy/");
  } else {
    dirList.push_back("EasyParameterListInterpreter/");
    dirList.push_back("FactoryParameterListInterpreter/");
  }
#if defined(HAVE_MPI) && defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_AMESOS2_KLU2)
  // The ML interpreter have internal ifdef, which means that the resulting
  // output would depend on configuration (reguarl interpreter does not have
  // that). Therefore, we need to stabilize the configuration here.
  // In addition, we run ML parameter list tests only if KLU is available
  dirList.push_back("MLParameterListInterpreter/");
  dirList.push_back("MLParameterListInterpreter2/");
#endif
  int numLists = dirList.size();

  bool failed = false;
  Teuchos::Time timer("Interpreter timer");
  //double lastTime = timer.wallTime();
  for (int k = 0; k < numLists; k++) {
    Teuchos::ArrayRCP<std::string> fileList = MueLuTests::TestHelpers::GetFileList(dirList[k],
      (numProc == 1 ? std::string(".xml") : std::string("_np" + Teuchos::toString(numProc) + ".xml")));

    for (int i = 0; i < fileList.size(); i++) {
      // Set seed
      std::srand(12345);

      // Reset (potentially) cached value of the estimate
      A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());

      std::string xmlFile  = dirList[k] + fileList[i];
      std::string outFile  = outDir     + fileList[i];
      std::string baseFile = outFile.substr(0, outFile.find_last_of('.'));
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
      std::string goldFile = baseFile + ".gold";
      std::ifstream f(goldFile.c_str());
      if (!f.good()) {
        if (myRank == 0)
          std::cout << "Warning: comparison file " << goldFile << " not found.  Skipping test" << std::endl;
        continue;
      }

      std::filebuf    buffer;
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
      if      (dirList[k] == "EasyParameterListInterpreter/" || dirList[k] == "EasyParameterListInterpreter-heavy/")
        paramList.set("verbosity", "test");
      else if (dirList[k] == "FactoryParameterListInterpreter/" || dirList[k] == "FactoryParameterListInterpreter-heavy/")
        paramList.sublist("Hierarchy").set("verbosity", "Test");
      else if (dirList[k] == "MLParameterListInterpreter/")
        paramList.set("ML output",     42);
      else if (dirList[k] == "MLParameterListInterpreter2/")
        paramList.set("ML output",     10);

      try {
        timer.start();
        Teuchos::RCP<HierarchyManager> mueluFactory;

        // create parameter list interpreter
        // here we have to distinguish between the general MueLu parameter list interpreter
        // and the ML parameter list interpreter. Note that the ML paramter interpreter also
        // works with Tpetra matrices.
        if (dirList[k] == "EasyParameterListInterpreter/"         ||
            dirList[k] == "EasyParameterListInterpreter-heavy/"   ||
            dirList[k] == "FactoryParameterListInterpreter/"      ||
            dirList[k] == "FactoryParameterListInterpreter-heavy/") {
          mueluFactory = Teuchos::rcp(new ParameterListInterpreter(paramList));

        } else if (dirList[k] == "MLParameterListInterpreter/") {
          mueluFactory = Teuchos::rcp(new MLParameterListInterpreter(paramList));

        } else if (dirList[k] == "MLParameterListInterpreter2/") {
          //std::cout << "ML ParameterList: " << std::endl;
          //std::cout << paramList << std::endl;
          RCP<Teuchos::ParameterList> mueluParamList = Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(paramList,"SA"));
          //std::cout << "MueLu ParameterList: " << std::endl;
          //std::cout << *mueluParamList << std::endl;
          mueluFactory = Teuchos::rcp(new ParameterListInterpreter(*mueluParamList));
        }

        RCP<Hierarchy> H = mueluFactory->CreateHierarchy();

        H->GetLevel(0)->template Set<RCP<Matrix> >("A", A);

        if (dirList[k] == "MLParameterListInterpreter/") {
          // MLParameterInterpreter needs the nullspace information if rebalancing is active!
          // add default constant null space vector
          RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getRowMap(), 1);
          nullspace->putScalar(1.0);
          H->GetLevel(0)->Set("Nullspace", nullspace);
        }

        H->GetLevel(0)->Set("Coordinates", coordinates);

        mueluFactory->SetupHierarchy(*H);

        if (strncmp(fileList[i].c_str(), "reuse", 5) == 0) {
          // Build the Hierarchy the second time
          // Should be faster if we actually do the reuse
          A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
          mueluFactory->SetupHierarchy(*H);
        }

        timer.stop();
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
            msg == "Zoltan2 interface is not available" ||
            msg == "MueLu::FactoryFactory:BuildFactory(): Cannot create a Zoltan2Interface object: Zoltan2 is disabled: HAVE_MUELU_ZOLTAN2 && HAVE_MPI == false.") {

          if (myRank == 0)
            std::cout << xmlFile << ": skipped (missing library)" << std::endl;

          continue;
        }
      }

      std::string cmd;
      if (myRank == 0) {
        // Redirect output back
        std::cout.rdbuf(oldbuffer);
        buffer.close();

        // Create a copy of outputs
        cmd = "cp -f ";
        system((cmd + baseFile + ".gold " + baseFile + ".gold_filtered").c_str());
        system((cmd + baseFile + ".out " + baseFile + ".out_filtered").c_str());

        // Tpetra produces different eigenvalues in Chebyshev due to using
        // std::rand() for generating random vectors, which may be initialized
        // using different seed, and may have different algorithm from one
        // gcc version to another, or to anogther compiler (like clang)
        // This leads to us always failing this test.
        // NOTE1 : Epetra, on the other hand, rolls out its out random number
        // generator, which always produces same results

        // Ignore the value of "lambdaMax"
        run_sed("'s/lambdaMax: [0-9]*.[0-9]*/lambdaMax = <ignored>/'", baseFile);

        // Ignore the value of "lambdaMin"
        run_sed("'s/lambdaMin: [0-9]*.[0-9]*/lambdaMin = <ignored>/'", baseFile);

        // Ignore the value of "chebyshev: max eigenvalue"
        // NOTE: we skip lines with default value ([default])
        run_sed("'/[default]/! s/chebyshev: max eigenvalue = [0-9]*.[0-9]*/chebyshev: max eigenvalue = <ignored>/'", baseFile);

        // Ignore the exact type of direct solver (it is selected semi-automatically
        // depending on how Trilinos was configured
        run_sed("'s/Amesos\\([2]*\\)Smoother{type = .*}/Amesos\\1Smoother{type = <ignored>}/'", baseFile);
        run_sed("'s/SuperLU solver interface, direct solve/<Direct> solver interface/'", baseFile);
        run_sed("'s/KLU2 solver interface/<Direct> solver interface/'", baseFile);
        run_sed("'s/Basker solver interface/<Direct> solver interface/'", baseFile);

        // Strip template args for some classes
        std::vector<std::string> classes;
        classes.push_back("Xpetra::Matrix");
        classes.push_back("MueLu::Constraint");
        classes.push_back("MueLu::SmootherPrototype");
        for (size_t q = 0; q < classes.size(); q++)
          run_sed("'s/" + classes[q] + "<.*>/" + classes[q] + "<ignored> >/'", baseFile);

        // Strip ParamterList pointers
        run_sed("'s/Teuchos::RCP<Teuchos::ParameterList const>{.*}/Teuchos::RCP<Teuchos::ParameterList const>{<ignored>} >/'", baseFile);

#ifdef __APPLE__
        // Some Macs print outs ptrs as 0x0 instead of 0, fix that
        run_sed("'/RCP/ s/=0x0/=0/g'", baseFile);
#endif

        // Run comparison (ignoring whitespaces)
        cmd = "diff -u -w -I\"^\\s*$\" " + baseFile + ".gold_filtered " + baseFile + ".out_filtered";
        int ret = system(cmd.c_str());
        if (ret)
          failed = true;

        //std::ios_base::fmtflags ff(std::cout.flags());
        //std::cout.precision(2);
        //std::cout << xmlFile << " (" << std::setiosflags(std::ios::fixed)
        //          << timer.wallTime() - lastTime << " sec.) : " << (ret ? "failed" : "passed") << std::endl;
        //lastTime = timer.wallTime();
        //std::cout.flags(ff); // reset flags to whatever they were prior to printing time
        std::cout << xmlFile << " : " << (ret ? "failed" : "passed") << std::endl;
      }
    }
  }

  if (myRank == 0)
    std::cout << std::endl << "End Result: TEST " << (failed ? "FAILED" : "PASSED") << std::endl;

  return (failed ? EXIT_FAILURE : EXIT_SUCCESS);
}

int main(int argc, char* argv[]) {
  bool success = false;
  bool verbose = true;

  try {
    const bool throwExceptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions);
    Xpetra::Parameters xpetraParameters(clp);

    clp.recogniseAllOptions(false);
    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:               return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double,int,int,Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
      return main_<double,int,long,Node>(clp, lib, argc, argv);
#else
#  if defined(HAVE_MUELU_INST_DOUBLE_INT_INT)           && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
      return main_<double,int,int,Node> (clp, lib, argc, argv);
#  elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)     && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
      return main_<double,int,long,Node>(clp, lib, argc, argv);
#  elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
      return main_<double,int,long long,Node>(clp, lib, argc, argv);
#  else
      throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation for Tpetra");
#  endif
#endif

#else
      throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

void run_sed(const std::string& pattern, const std::string& baseFile) {
  // sed behaviour differs between Mac and Linux
  // You can run "sed -i 's//' " in Linux, but you always have to specify
  // "sed -i "<smth,could be empty>" 's//'" in Mac. Both, however, take '-i<extension>'
  std::string sed_pref = "sed -i ";
#ifdef __APPLE__
  sed_pref = sed_pref +  "\"\" ";
#endif

  system((sed_pref + pattern + " " + baseFile + ".gold_filtered").c_str());
  system((sed_pref + pattern + " " + baseFile + ".out_filtered").c_str());
}

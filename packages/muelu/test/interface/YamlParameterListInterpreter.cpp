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
#include <MueLu_YamlParser_def.hpp>

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif

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
  int matched = 0;
  int total = 0;
  Teuchos::Time timer("Interpreter timer");
  //double lastTime = timer.wallTime();
  for (int k = 0; k < numLists; k++) {
    Teuchos::ArrayRCP<std::string> fileList = MueLuTests::TestHelpers::GetFileList(dirList[k],
      (numProc == 1 ? std::string(".xml") : std::string("_np" + Teuchos::toString(numProc) + ".xml")));

    for (int i = 0; i < fileList.size(); i++) {
      // Set seed
      std::srand(12345);
      std::string xmlFile  = dirList[k] + fileList[i];
      std::size_t found = xmlFile.find("_np");
      if (numProc == 1 && found != std::string::npos) {
#ifndef HAVE_MPI
        std::cout << "Skipping \"" << xmlFile << "\" as MPI is not enabled" << std::endl;
        continue;
#endif
      }
      RCP<Teuchos::ParameterList> xmlList = Teuchos::getParametersFromXmlFile(xmlFile);
      RCP<Teuchos::ParameterList> yamlList = rcp(new ParameterList);
      if(comm->getRank() == 0)
      {
        MueLu::convertXmlToYaml(xmlFile);
      }
      std::string yamlFile = xmlFile.substr(0, xmlFile.length() - 4) + ".yaml";
      MueLu::updateParametersFromYamlFileAndBroadcast(yamlFile, yamlList.ptr(), *comm);
      if(comm->getRank() == 0)
      {
        //remove(yamlFile.c_str());
      }
      if(!MueLu::haveSameValuesUnordered(*xmlList, *yamlList))
      {
        std::cout << yamlFile << " : failed\n";
        std::cout << "XML version: \n\n" << *xmlList << "\n\n";
        std::cout << "YAML version: \n\n" << *yamlList << "\n\n";
        failed = true;
      }
      else
      {
        std::cout << yamlFile << " : passed\n";
        matched++;
      }
      total++;
    }
  }
  std::cout << matched << " out of " << total << " matched after conversion to YAML." << std::endl;
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

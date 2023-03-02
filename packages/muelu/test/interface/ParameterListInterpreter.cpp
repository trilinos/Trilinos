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
#include <algorithm>
#include <cstdio>

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
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib& lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  using real_type             = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================

  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int numProc = comm->getSize();
  int myRank  = comm->getRank();

  // =========================================================================
  // Parameters initialization
  // =========================================================================

  bool runHeavyTests = false;
  std::string xmlForceFile = "";
  bool useKokkos = false;
  if(lib == Xpetra::UseTpetra) {
# ifdef HAVE_MUELU_SERIAL
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosSerialWrapperNode).name())
      useKokkos = false;
# endif
# ifdef HAVE_MUELU_OPENMP
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosOpenMPWrapperNode).name())
      useKokkos = true;
# endif
# ifdef HAVE_MUELU_CUDA
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosCudaWrapperNode).name())
      useKokkos = true;
# endif
# ifdef HAVE_MUELU_HIP
    if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosHIPWrapperNode).name())
      useKokkos = true;
# endif
  }
  bool compareWithGold = true;
#ifdef KOKKOS_ENABLE_CUDA
  if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosCudaWrapperNode).name())
    // Behavior of some algorithms on Cuda is non-deterministic, so we won't check the output.
    compareWithGold = false;
#endif
#ifdef KOKKOS_ENABLE_HIP
  if (typeid(Node).name() == typeid(Kokkos::Compat::KokkosHIPWrapperNode).name())
    // Behavior of some algorithms on HIP is non-deterministic, so we won't check the output.
    compareWithGold = false;
#endif
  clp.setOption("useKokkosRefactor", "noKokkosRefactor", &useKokkos, "use kokkos refactor");
  clp.setOption("heavytests", "noheavytests",  &runHeavyTests, "whether to exercise tests that take a long time to run");
  clp.setOption("xml", &xmlForceFile, "xml input file (useful for debugging)");
  clp.setOption("compareWithGold", "skipCompareWithGold", &compareWithGold, "compare runs against gold files");
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
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<real_type,LO,GO,Map,RealValuedMultiVector>("1D", A->getRowMap(), matrixParameters);
  RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getRowMap(), 1);
  nullspace->putScalar(1.0);

  std::string prefix;
  if (useKokkos) {
    prefix = "kokkos/";
  } else {
    prefix = "default/";
  }
  std::string outDir = prefix+"Output/";

  std::vector<std::string> dirList;
  if (runHeavyTests) {
    dirList.push_back(prefix+"EasyParameterListInterpreter-heavy/");
    if (!useKokkos) {
      // commented since extended xml interface does not support kokkos factories
      dirList.push_back(prefix+"FactoryParameterListInterpreter-heavy/");
    }
  } else {
    dirList.push_back(prefix+"EasyParameterListInterpreter/");
    if (!useKokkos) {
      // commented since extended xml interface does not support kokkos factories
      dirList.push_back(prefix+"FactoryParameterListInterpreter/");
    }
  }
#if defined(HAVE_MPI) && defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_AMESOS2_KLU2)
  // The ML interpreter have internal ifdef, which means that the resulting
  // output would depend on configuration (reguarl interpreter does not have
  // that). Therefore, we need to stabilize the configuration here.
  // In addition, we run ML parameter list tests only if KLU is available
  dirList.push_back(prefix+"MLParameterListInterpreter/");
  dirList.push_back(prefix+"MLParameterListInterpreter2/");
#endif
  int numLists = dirList.size();

  bool failed = false;
  bool jumpOut = false;
  Teuchos::Time timer("Interpreter timer");
  //double lastTime = timer.wallTime();
  for (int k = 0; k < numLists; k++) {
    Teuchos::ArrayRCP<std::string> fileList = MueLuTests::TestHelpers::GetFileList(dirList[k],
      (numProc == 1 ? std::string(".xml") : std::string("_np" + Teuchos::toString(numProc) + ".xml")));

    std::sort(fileList.begin(),fileList.end());

    for (int i = 0; i < fileList.size(); i++) {
      // Set seed
      std::srand(12345);

      // Reset (potentially) cached value of the estimate
      A->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());

      std::string xmlFile;
      std::string outFile;
      std::string baseFile;
      std::size_t found;
      if (xmlForceFile==""){
        xmlFile  = dirList[k] + fileList[i];
        outFile  = outDir     + fileList[i];
        baseFile = outFile.substr(0, outFile.find_last_of('.'));
        found = baseFile.find("_np");
      } else {
        xmlFile = xmlForceFile;
        std::string dir = prefix + xmlForceFile.substr(0, xmlForceFile.find_last_of('/')+1);
        dirList[k] = dir;
        outFile = outDir + xmlForceFile.substr(xmlForceFile.find_last_of('/')+1, xmlForceFile.size());
        baseFile = outFile.substr(0, outFile.find_last_of('.'));
        found = baseFile.find("_np");
        jumpOut = true;
        xmlFile = prefix + xmlFile;
        std::cout << "Test dir: " << dirList[k] << std::endl;
      }

      if (numProc == 1 && found != std::string::npos) {
#ifdef HAVE_MPI
        baseFile = baseFile.substr(0, found);
#else
        std::cout << "Skipping \"" << xmlFile << "\" as MPI is not enabled" << std::endl;
        continue;
#endif
      }
      if (myRank == 0)
        std::cout << "Testing: "<< xmlFile << std::endl;

      baseFile = baseFile + (lib == Xpetra::UseEpetra ? "_epetra" : "_tpetra");
      std::string goldFile = baseFile + ".gold";
      std::ifstream f(goldFile.c_str());
      if (!f.good()) {
        if (myRank == 0)
          std::cout << "Warning: comparison file " << goldFile << " not found.  Skipping test" << std::endl;
        continue;
      }

      std::stringbuf buffer;
      std::streambuf* oldbuffer = NULL;
      //   // Redirect output
      oldbuffer = std::cout.rdbuf(&buffer);

      // NOTE: we cannot use ParameterListInterpreter(xmlFile, comm), because we want to update the ParameterList
      // first to include "test" verbosity
      Teuchos::ParameterList paramList;
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFile, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);
      if      (dirList[k] == prefix+"EasyParameterListInterpreter/" || dirList[k] == prefix+"EasyParameterListInterpreter-heavy/")
        paramList.set("verbosity", "interfacetest");
      else if (dirList[k] == prefix+"FactoryParameterListInterpreter/" || dirList[k] == prefix+"FactoryParameterListInterpreter-heavy/")
        paramList.sublist("Hierarchy").set("verbosity", "InterfaceTest");
      else if (dirList[k] == prefix+"MLParameterListInterpreter/" || dirList[k] == prefix+"MLParameterListInterpreter2/")
        paramList.set("ML output",     666);

      try {
        timer.start();
        Teuchos::RCP<HierarchyManager> mueluFactory;

        // create parameter list interpreter
        // here we have to distinguish between the general MueLu parameter list interpreter
        // and the ML parameter list interpreter. Note that the ML paramter interpreter also
        // works with Tpetra matrices.
        bool coordsSet = false;
        if (dirList[k] == prefix+"EasyParameterListInterpreter/"         ||
            dirList[k] == prefix+"EasyParameterListInterpreter-heavy/") {
          paramList.set("use kokkos refactor", useKokkos);
          mueluFactory = Teuchos::rcp(new ParameterListInterpreter(paramList,comm));
        } else if (dirList[k] == prefix+"FactoryParameterListInterpreter/"      ||
                   dirList[k] == prefix+"FactoryParameterListInterpreter-heavy/") {
          paramList.sublist("Hierarchy").set("use kokkos refactor", useKokkos);
          mueluFactory = Teuchos::rcp(new ParameterListInterpreter(paramList,comm));

        } else if (dirList[k] == prefix+"MLParameterListInterpreter/") {
          paramList.set("use kokkos refactor", useKokkos);

          paramList.set("x-coordinates",coordinates->getDataNonConst(0).get());
          if (coordinates->getNumVectors() > 1)
            paramList.set("y-coordinates",coordinates->getDataNonConst(1).get());
          if (coordinates->getNumVectors() > 2)
            paramList.set("z-coordinates",coordinates->getDataNonConst(2).get());
          coordsSet = true;

          // MLParameterInterpreter needs the nullspace information if rebalancing is active!
          // add default constant null space vector
          paramList.set("null space: type", "pre-computed");
          paramList.set("null space: dimension", Teuchos::as<int>(nullspace->getNumVectors()));
          paramList.set("null space: vectors", nullspace->getDataNonConst(0).get());

          mueluFactory = Teuchos::rcp(new MLParameterListInterpreter(paramList));

        } else if (dirList[k] == prefix+"MLParameterListInterpreter2/") {
          RCP<Teuchos::ParameterList> mueluParamList = Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(paramList,"SA"));
          mueluParamList->set("use kokkos refactor", useKokkos);
          mueluFactory = Teuchos::rcp(new ParameterListInterpreter(*mueluParamList));
        } else
          TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Not a test directory");

        RCP<Hierarchy> H = mueluFactory->CreateHierarchy();

        H->GetLevel(0)->template Set<RCP<Matrix> >("A", A);

        if (!coordsSet)
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

        if (myRank == 0)
          std::cout << "Caught exception: " << msg << std::endl;

        if (msg == "Zoltan interface is not available" ||
            msg == "Zoltan2 interface is not available" ||
            msg == "MueLu::FactoryFactory:BuildFactory(): Cannot create a Zoltan2Interface object: Zoltan2 is disabled: HAVE_MUELU_ZOLTAN2 && HAVE_MPI == false.") {

          if (myRank == 0)
            std::cout << xmlFile << ": skipped (missing library)" << std::endl;

          continue;
        }
      }

      // Redirect output back
      std::cout.rdbuf(oldbuffer);
#ifdef HAVE_MPI
      std::string logStr = buffer.str();
      if (myRank == 0)
        remove((baseFile + ".out").c_str());
      RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm)->getRawMpiComm();
      MPI_File logfile;
      comm->barrier();
      MPI_File_open((*mpiComm)(), (baseFile + ".out").c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &logfile);
      MPI_File_set_atomicity(logfile, true);
      const char* msg = logStr.c_str();
      int err = MPI_File_write_ordered(logfile, msg, logStr.size(), MPI_CHAR, MPI_STATUS_IGNORE);
      TEUCHOS_ASSERT(err == MPI_SUCCESS);
      MPI_File_close(&logfile);
#else
      std::ofstream outStream;
      outStream.open((baseFile + ".out").c_str(), std::ofstream::out);
      outStream << buffer.str();
      outStream.close();
#endif

      std::string cmd;
      if (myRank == 0) {
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

        // make sure complex tests pass
        run_sed("'s/relaxation: damping factor = (1,0)/relaxation: damping factor = 1/'", baseFile);
        run_sed("'s/damping factor: (1,0)/damping factor: 1/'", baseFile);
        run_sed("'s/relaxation: min diagonal value = (0,0)/relaxation: min diagonal value = 0/'", baseFile);

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
        run_sed("'s/Teuchos::RCP<Teuchos::ParameterList const>{.*}/Teuchos::RCP<Teuchos::ParameterList const>{<ignored>}/'", baseFile);

#ifdef __APPLE__
        // Some Macs print outs ptrs as 0x0 instead of 0, fix that
        run_sed("'/RCP/ s/=0x0/=0/g'", baseFile);
#endif

        // Run comparison (ignoring whitespaces)
        cmd = "diff -u -w -I\"^\\s*$\" " + baseFile + ".gold_filtered " + baseFile + ".out_filtered";
        int ret = 0;
        if (compareWithGold)
          ret = system(cmd.c_str());
        else
          std::cout << "Skipping comparison with gold file\n";
        if (ret) {
          failed = true;
        }

        //std::ios_base::fmtflags ff(std::cout.flags());
        //std::cout.precision(2);
        //std::cout << xmlFile << " (" << std::setiosflags(std::ios::fixed)
        //          << timer.wallTime() - lastTime << " sec.) : " << (ret ? "failed" : "passed") << std::endl;
        //lastTime = timer.wallTime();
        //std::cout.flags(ff); // reset flags to whatever they were prior to printing time
        std::cout << xmlFile << " : " << (ret ? "failed" : "passed") << std::endl;
      }
      if (jumpOut)
        break;
    }
    if (jumpOut)
      break;
  }

  if (myRank == 0)
    std::cout << std::endl << "End Result: TEST " << (failed ? "FAILED" : "PASSED") << std::endl;

  return (failed ? EXIT_FAILURE : EXIT_SUCCESS);
}


//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
#ifdef KOKKOS_ENABLE_OPENMP
  omp_set_num_threads(1);
#endif
  return Automatic_Test_ETI(argc,argv);
}


void run_sed(const std::string& pattern, const std::string& baseFile) {
  // sed behaviour differs between Mac and Linux
  // You can run "sed -i 's//' " in Linux, but you always have to specify
  // "sed -i "<smth,could be empty>" 's//'" in Mac. Both, however, take '-i<extension>'
  std::string sed_pref = "sed -i ";
#ifdef __APPLE__
  sed_pref = sed_pref +  "\"\" ";
#endif
  int ret;
  ret = system((sed_pref + pattern + " " + baseFile + ".gold_filtered").c_str());
  TEUCHOS_ASSERT_EQUALITY(ret,0);
  ret = system((sed_pref + pattern + " " + baseFile + ".out_filtered").c_str());
  TEUCHOS_ASSERT_EQUALITY(ret,0);
}

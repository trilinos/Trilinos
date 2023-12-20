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
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_CrsMatrixFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

#include <MueLu_CreateTpetraPreconditioner.hpp>
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_CreateEpetraPreconditioner.hpp>
#include <EpetraExt_MMHelpers.h>
#include <EpetraExt_RowMatrixOut.h>
#endif

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > build_map_for_transfer(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > oldMap) {
#include <MueLu_UseShortNames.hpp>
  // 10 unknowns per proc except N-1: Assume Ids run from 0 to MaxGlobalIndex
  int Nproc = oldMap->getComm()->getSize();
  int MyPID = oldMap->getComm()->getRank();
  int N     = (int)oldMap->getGlobalNumElements();

  Xpetra::UnderlyingLib lib = oldMap->lib();

  if (Nproc == 1) return oldMap;

  int start = MyPID == Nproc - 1 ? 0 : std::max(N - (MyPID + 1) * 10, 0);
  int stop  = std::max(N - MyPID * 10, 0);

  Teuchos::Array<GlobalOrdinal> elems;
  if (stop - start > 0) {
    elems.resize(stop - start);
    for (int i = 0; i < stop - start; i++)
      elems[i] = start + i;
  }

#if 0
  printf("[%d] elems(%d) = ",MyPID,(int)elems.size());
  for(int i=0; i<(int)elems.size(); i++)
    printf("%d ",elems[i]);
  printf("\n");
  fflush(stdout);
#endif

  return Xpetra::MapFactory<LO, GO, Node>::Build(lib, oldMap->getGlobalNumElements(), elems(), oldMap->getIndexBase(), oldMap->getComm());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > build_map_for_transfer_repartition(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > oldMap) {
#include <MueLu_UseShortNames.hpp>
  // Turn P procs into ~P/9 procs
  int Nproc = oldMap->getComm()->getSize();
  int MyPID = oldMap->getComm()->getRank();
  int N     = (int)oldMap->getGlobalNumElements();

  Xpetra::UnderlyingLib lib = oldMap->lib();

  if (Nproc == 1) return oldMap;
  int ideal_new_unknowns_per_proc = N / (Nproc / 9.0);

  int start = ideal_new_unknowns_per_proc * MyPID < N ? ideal_new_unknowns_per_proc * MyPID : 0;
  int stop  = ideal_new_unknowns_per_proc * MyPID < N ? std::min(N, ideal_new_unknowns_per_proc * (MyPID + 1)) : 0;

  Teuchos::Array<GlobalOrdinal> elems;
  if (stop - start > 0) {
    elems.resize(stop - start);
    for (int i = 0; i < stop - start; i++)
      elems[i] = start + i;
  }

#if 0
  printf("[%d] elems(%d) = ",MyPID,(int)elems.size());
  for(int i=0; i<(int)elems.size(); i++)
    printf("%d ",elems[i]);
  printf("\n");
  fflush(stdout);
#endif

  int i_am_active = (elems.size() > 0);
  int num_active  = 0;
  Teuchos::reduce(&i_am_active, &num_active, 1, Teuchos::REDUCE_SUM, 0, *oldMap->getComm());
  if (MyPID == 0) printf("Repartitioning to %d/%d processors\n", num_active, Nproc);

  return Xpetra::MapFactory<LO, GO, Node>::Build(lib, oldMap->getGlobalNumElements(), elems(), oldMap->getIndexBase(), oldMap->getComm());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > toXpetraCrs(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &A) {
#include <MueLu_UseShortNames.hpp>
  RCP<const CrsMatrixWrap> crswrapA = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
  TEUCHOS_TEST_FOR_EXCEPTION(crswrapA == Teuchos::null, MueLu::Exceptions::BadCast,
                             "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  RCP<CrsMatrix> crsA = crswrapA->getCrsMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(crsA == Teuchos::null, std::runtime_error,
                             "Xpetra::CrsMatrixWrap doesn't have a CrsMatrix");
  return crsA;
}

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TestTransferAndFillComplete(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &A, Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &importer) {
#include <MueLu_UseShortNames.hpp>
  //  Xpetra::UnderlyingLib lib = A->getRowMap()->lib();
  RCP<TimeMonitor> tm;

  // It makes me sad to do this
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > constA = Teuchos::rcp<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(&*A, false);

  // Only makes sense in parallel
  if (A->getRowMap()->getComm()->getSize() == 1) return;

  // ==================
  // Optimized Transfer
  // ==================
  {
    tm                                              = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("OptimizedTransfer")));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > B = Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(toXpetraCrs(A), *importer);
  }

  A->getRowMap()->getComm()->barrier();

  // ==================
  // Naive Transfer
  // ==================
  {
    tm                                                                   = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("NaiveTransfer")));
    RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > B = Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(importer->getTargetMap(), 0);
    B->doImport(*toXpetraCrs(A), *importer, Xpetra::ADD);
    B->fillComplete();
  }
  A->getRowMap()->getComm()->barrier();
}

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 200, ny = 200, nz = 10;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "import.xml";
  clp.setOption("xml", &xmlFileName, "read parameters from a file");
  bool printTimings = true;
  clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  std::string timingsFormat = "table-fixed";
  clp.setOption("time-format", &timingsFormat, "timings format (table-fixed | table-scientific | yaml)");
  int numImports = 100;
  clp.setOption("numImport", &numImports, "#times to test");
  std::string mapmode = "small";
  clp.setOption("mapMode", &mapmode, "map to use ('small' or 'repartition')");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  ParameterList paramList;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  bool isDriver = paramList.isSublist("Run1");
  if (isDriver) {
    // update galeriParameters with the values from the XML file
    ParameterList &realParams = galeriParameters.GetParameterList();

    for (ParameterList::ConstIterator it = realParams.begin(); it != realParams.end(); it++) {
      const std::string &name = realParams.name(it);
      if (paramList.isParameter(name))
        realParams.setEntry(name, paramList.getEntry(name));
    }
  }

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  std::ostringstream galeriStream;
  comm->barrier();
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Matrix Build")));

  RCP<Xpetra::Matrix<Scalar, LO, GO, Node> > A;
  RCP<const Map> map;
  RCP<RealValuedMultiVector> coordinates;
  typedef typename RealValuedMultiVector::scalar_type Real;
  RCP<MultiVector> nullspace;
  galeriStream << "========================================================\n"
               << xpetraParameters << galeriParameters;

  // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
  //                                 d1  d2  d3
  //                                 d4  d5  d6
  //                                 d7  d8  d9
  //                                 d10 d11 d12
  // A perfect distribution is only possible when the #processors is a perfect square.
  // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
  // size. For example, np=14 will give a 7-by-2 distribution.
  // If you don't want Galeri to do this, specify mx or my on the galeriList.
  std::string matrixType = galeriParameters.GetMatrixType();

  // Create map and coordinates
  // In the future, we hope to be able to first create a Galeri problem, and then request map and coordinates from it
  // At the moment, however, things are fragile as we hope that the Problem uses same map and coordinates inside
  if (matrixType == "Laplace1D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Real, LO, GO, Map, RealValuedMultiVector>("1D", map, galeriList);

  } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
             matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Real, LO, GO, Map, RealValuedMultiVector>("2D", map, galeriList);

  } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
    map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<Real, LO, GO, Map, RealValuedMultiVector>("3D", map, galeriList);
  }

  // Expand map to do multiple DOF per node for block problems
  if (matrixType == "Elasticity2D")
    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 2);
  if (matrixType == "Elasticity3D")
    map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

  galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
               << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
               << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
               << "========================================================" << std::endl;

  if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
    // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
    galeriList.set("right boundary", "Neumann");
    galeriList.set("bottom boundary", "Neumann");
    galeriList.set("top boundary", "Neumann");
    galeriList.set("front boundary", "Neumann");
    galeriList.set("back boundary", "Neumann");
  }

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
  A = Pr->BuildMatrix();

  comm->barrier();
  tm = Teuchos::null;

  galeriStream << "Galeri complete.\n========================================================" << std::endl;

  int numReruns = 1;
  if (paramList.isParameter("number of reruns"))
    numReruns = paramList.get<int>("number of reruns");

  const bool mustAlreadyExist = true;
  for (int rerunCount = 1; rerunCount <= numReruns; rerunCount++) {
    ParameterList mueluList, runList;
    bool stop = false;
    if (isDriver) {
      runList   = paramList.sublist("Run1", mustAlreadyExist);
      mueluList = runList.sublist("MueLu", mustAlreadyExist);
    } else {
      mueluList = paramList;
      stop      = true;
    }

    int runCount = 1;
    do {
      int savedOut    = -1;
      FILE *openedOut = NULL;
      if (isDriver) {
        if (runList.isParameter("filename")) {
          // Redirect all output into a filename We have to redirect all output,
          // including printf's, therefore we cannot simply replace C++ cout
          // buffers, and have to use heavy machinary (dup2)
          std::string filename = runList.get<std::string>("filename");
          if (numReruns > 1)
            filename += "_run" + MueLu::toString(rerunCount);
          filename += (lib == Xpetra::UseEpetra ? ".epetra" : ".tpetra");

          savedOut  = dup(STDOUT_FILENO);
          openedOut = fopen(filename.c_str(), "w");
          dup2(fileno(openedOut), STDOUT_FILENO);
        }
      }

      // Instead of checking each time for rank, create a rank 0 stream
      RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      Teuchos::FancyOStream &out       = *fancy;
      out.setOutputToRootOnly(0);
      out << galeriStream.str();

      // Build the target map for the importing
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > importMap;
      if (mapmode == "small")
        importMap = build_map_for_transfer<Scalar, LO, GO, Node>(A->getRowMap());
      else if (mapmode == "repartition")
        importMap = build_map_for_transfer_repartition<Scalar, LO, GO, Node>(A->getRowMap());
      else
        throw std::runtime_error("Invalid map mode");

      Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > importer = Xpetra::ImportFactory<LO, GO, Node>::Build(A->getRowMap(), importMap);

      for (int i = 0; i < numImports; i++) {
        // =========================================================================
        // Optimized transfer & fill complete loop
        // =========================================================================
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 4 - TransferAndFillComplete")));
        TestTransferAndFillComplete<Scalar, LO, GO, Node>(A, importer);
        comm->barrier();
        tm = Teuchos::null;
      }

      // Cleanup
      globalTimeMonitor = Teuchos::null;

      // =========================================================================
      // Timing stuff
      // =========================================================================
      if (printTimings) {
        RCP<ParameterList> reportParams = rcp(new ParameterList);
        if (timingsFormat == "yaml") {
          reportParams->set("Report format", "YAML");  // "Table" or "YAML"
          reportParams->set("YAML style", "compact");  // "spacious" or "compact"
        }
        reportParams->set("How to merge timer sets", "Union");
        reportParams->set("alwaysWriteLocal", false);
        reportParams->set("writeGlobalStats", true);
        reportParams->set("writeZeroTimers", false);
        // FIXME: no "ignoreZeroTimers"

        const std::string filter = "";

        std::ios_base::fmtflags ff(out.flags());
        if (timingsFormat == "table-fixed")
          out << std::fixed;
        else
          out << std::scientific;
        TimeMonitor::report(comm.ptr(), out, filter, reportParams);
        out << std::setiosflags(ff);
      }

      TimeMonitor::clearCounters();

      if (isDriver) {
        if (openedOut != NULL) {
          TEUCHOS_ASSERT(savedOut >= 0);
          dup2(savedOut, STDOUT_FILENO);
          fclose(openedOut);
          openedOut = NULL;
        }
        try {
          runList   = paramList.sublist("Run" + MueLu::toString(++runCount), mustAlreadyExist);
          mueluList = runList.sublist("MueLu", mustAlreadyExist);
        } catch (Teuchos::Exceptions::InvalidParameterName &) {
          stop = true;
        }
      }

    } while (!stop);
  }

  return EXIT_SUCCESS;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}

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
#include <unistd.h>
#include <iostream>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp

#include <MueLu_Utilities.hpp>

#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 100, ny = 100, nz = 100;
    Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

    std::string xmlFileName = "viztest.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file");
    bool printTimings = true;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    int writeMatricesOPT = -2;
    clp.setOption("write", &writeMatricesOPT, "write matrices to file (-1 means all; i>=0 means level i)");
    std::string dsolveType = "cg", solveType;
    clp.setOption("solver", &dsolveType, "solve type: (none | cg | gmres | standalone)");
    double dtol = 1e-12, tol;
    clp.setOption("tol", &dtol, "solver convergence tolerance");

    std::string mapFile;
    clp.setOption("map", &mapFile, "map data file");
    std::string colMapFile;
    clp.setOption("colmap", &colMapFile, "colmap data file");
    std::string domainMapFile;
    clp.setOption("domainmap", &domainMapFile, "domainmap data file");
    std::string rangeMapFile;
    clp.setOption("rangemap", &rangeMapFile, "rangemap data file");
    std::string matrixFile;
    clp.setOption("matrix", &matrixFile, "matrix data file");
    std::string coordFile;
    clp.setOption("coords", &coordFile, "coordinates data file");
    std::string nullFile;
    clp.setOption("nullspace", &nullFile, "nullspace data file");
    int numRebuilds = 0;
    clp.setOption("rebuild", &numRebuilds, "#times to rebuild hierarchy");
    int maxIts = 200;
    clp.setOption("its", &maxIts, "maximum number of solver iterations");
    bool scaleResidualHist = true;
    clp.setOption("scale", "noscale", &scaleResidualHist, "scaled Krylov residual history");

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
      ParameterList& realParams = galeriParameters.GetParameterList();

      for (ParameterList::ConstIterator it = realParams.begin(); it != realParams.end(); it++) {
        const std::string& name = realParams.name(it);
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
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
    RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

    RCP<Matrix> A;
    RCP<const Map> map;
    RCP<MultiVector> coordinates;
    RCP<MultiVector> nullspace;
    int ndims = 2;
    if (matrixFile.empty()) {
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
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", map, galeriList);

      } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
                 matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
        map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("2D", map, galeriList);

      } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
        map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("3D", map, galeriList);
        ndims       = 3;
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

      if (matrixType == "Elasticity2D" ||
          matrixType == "Elasticity3D") {
        nullspace = Pr->BuildNullspace();
        A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
      }

    } else {
      if (!mapFile.empty())
        map = Xpetra::IO<SC, LO, GO, Node>::ReadMap(mapFile, lib, comm);
      comm->barrier();

      const bool binaryFormat = false;

      if (!binaryFormat && !map.is_null()) {
        RCP<const Map> colMap    = (!colMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(colMapFile, lib, comm) : Teuchos::null);
        RCP<const Map> domainMap = (!domainMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(domainMapFile, lib, comm) : Teuchos::null);
        RCP<const Map> rangeMap  = (!rangeMapFile.empty() ? Xpetra::IO<SC, LO, GO, Node>::ReadMap(rangeMapFile, lib, comm) : Teuchos::null);
        A                        = Xpetra::IO<SC, LO, GO, Node>::Read(matrixFile, map, colMap, domainMap, rangeMap);

      } else {
        A = Xpetra::IO<SC, LO, GO, Node>::Read(matrixFile, lib, comm, binaryFormat);

        if (!map.is_null()) {
          RCP<Matrix> newMatrix = MatrixFactory::Build(map, 1);
          RCP<Import> importer  = ImportFactory::Build(A->getRowMap(), map);
          newMatrix->doImport(*A, *importer, Xpetra::INSERT);
          newMatrix->fillComplete();

          A.swap(newMatrix);
        }
      }
      map = A->getMap();

      comm->barrier();

      if (!coordFile.empty()) {
        // NOTE: currently we only allow reading scalar matrices, thus coordinate
        // map is same as matrix map
        coordinates = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(coordFile, map);
      }

      if (!nullFile.empty())
        nullspace = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(nullFile, map);
    }

    comm->barrier();
    tm = Teuchos::null;

    galeriStream << "Galeri complete.\n========================================================" << std::endl;

    if (ndims == 2)
      paramList.set<std::string>("aggregation: output filename", "MPI-Viz-Output-2D-Level%LEVELID-Proc%PROCID");
    else if (ndims == 3)
      paramList.set<std::string>("aggregation: output filename", "MPI-Viz-Output-3D-Level%LEVELID-Proc%PROCID");
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

      mueluList.set("use kokkos refactor", false);

      if (nullspace.is_null()) {
        int blkSize = 1;
        if (mueluList.isSublist("Matrix")) {
          // Factory style parameter list
          const Teuchos::ParameterList& operatorList = paramList.sublist("Matrix");
          if (operatorList.isParameter("PDE equations"))
            blkSize = operatorList.get<int>("PDE equations");

        } else if (paramList.isParameter("number of equations")) {
          // Easy style parameter list
          blkSize = paramList.get<int>("number of equations");
        }

        nullspace = MultiVectorFactory::Build(map, blkSize);
        for (int i = 0; i < blkSize; i++) {
          RCP<const Map> domainMap = A->getDomainMap();
          GO indexBase             = domainMap->getIndexBase();

          ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
          for (int j = 0; j < nsData.size(); j++) {
            GO GID = domainMap->getGlobalElement(j) - indexBase;

            if ((GID - i) % blkSize == 0)
              nsData[j] = Teuchos::ScalarTraits<SC>::one();
          }
        }
      }

      int runCount = 1;
      do {
        solveType = dsolveType;
        tol       = dtol;

        int savedOut    = -1;
        FILE* openedOut = NULL;
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
          if (runList.isParameter("solver")) solveType = runList.get<std::string>("solver");
          if (runList.isParameter("tol")) tol = runList.get<double>("tol");
        }

        // Instead of checking each time for rank, create a rank 0 stream
        RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        Teuchos::FancyOStream& out       = *fancy;
        out.setOutputToRootOnly(0);

        out << galeriStream.str();

        // =========================================================================
        // Preconditioner construction
        // =========================================================================
        comm->barrier();
        tm.reset();
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));

        RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(mueluList));

        comm->barrier();
        tm.reset();
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

        RCP<Hierarchy> H;
        for (int i = 0; i <= numRebuilds; i++) {
          A->SetMaxEigenvalueEstimate(-one);

          H = mueLuFactory->CreateHierarchy();
          H->GetLevel(0)->Set("A", A);
          H->GetLevel(0)->Set("Nullspace", nullspace);
          if (!coordinates.is_null())
            H->GetLevel(0)->Set("Coordinates", coordinates);
          mueLuFactory->SetupHierarchy(*H);
        }

        comm->barrier();
        tm = Teuchos::null;

        // =========================================================================
        // System solution (Ax = b)
        // =========================================================================
        comm->barrier();
        tm.reset();
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));

        RCP<Vector> X = VectorFactory::Build(map);
        RCP<Vector> B = VectorFactory::Build(map);

        {
          // we set seed for reproducibility
          Utilities::SetRandomSeed(*comm);
          X->randomize();
          A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

          Teuchos::Array<typename STS::magnitudeType> norms(1);
          B->norm2(norms);
          B->scale(one / norms[0]);
          X->putScalar(zero);
        }
        tm = Teuchos::null;

        if (writeMatricesOPT > -2) {
          tm.reset();
          tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3.5 - Matrix output")));
          H->Write(writeMatricesOPT, writeMatricesOPT);
          tm = Teuchos::null;
        }

        comm->barrier();
        if (solveType == "none") {
          // Do not perform a solve

        } else if (solveType == "standalone") {
          tm.reset();
          tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Fixed Point Solve")));

          H->IsPreconditioner(false);
          H->Iterate(*B, *X, maxIts);

        } else if (solveType == "cg" || solveType == "gmres") {
#ifdef HAVE_MUELU_BELOS
          tm.reset();
          tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 5 - Belos Solve")));

          // Operator and Multivector type that will be used with Belos
          typedef MultiVector MV;
          typedef Belos::OperatorT<MV> OP;

          H->IsPreconditioner(true);

          // Define Operator and Preconditioner
          Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));  // Turns a Xpetra::Matrix object into a Belos operator
          Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));   // Turns a MueLu::Hierarchy object into a Belos operator

          // Construct a Belos LinearProblem object
          RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
          belosProblem->setRightPrec(belosPrec);

          bool set = belosProblem->setProblem();
          if (set == false) {
            out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
            return EXIT_FAILURE;
          }

          // Belos parameter list
          Teuchos::ParameterList belosList;
          belosList.set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
          belosList.set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
          belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
          belosList.set("Output Frequency", 1);
          belosList.set("Output Style", Belos::Brief);
          if (!scaleResidualHist)
            belosList.set("Implicit Residual Scaling", "None");

          // Create an iterative solver manager
          RCP<Belos::SolverManager<SC, MV, OP> > solver;
          if (solveType == "cg") {
            solver = rcp(new Belos::PseudoBlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
          } else if (solveType == "gmres") {
            solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
          }

          // Perform solve
          Belos::ReturnType ret = Belos::Unconverged;
          ret                   = solver->solve();

          // Get the number of iterations for this solve.
          out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

          // Check convergence
          if (ret != Belos::Converged)
            out << std::endl
                << "ERROR:  Belos did not converge! " << std::endl;
          else
            out << std::endl
                << "SUCCESS:  Belos converged!" << std::endl;
#endif  // ifdef HAVE_MUELU_BELOS
        } else {
          throw MueLu::Exceptions::RuntimeError("Unknown solver type: \"" + solveType + "\"");
        }
        comm->barrier();
        tm                = Teuchos::null;
        globalTimeMonitor = Teuchos::null;

        if (printTimings) {
          const bool alwaysWriteLocal = false;
          const bool writeGlobalStats = true;
          const bool writeZeroTimers  = false;
          const bool ignoreZeroTimers = true;
          const std::string filter    = "";
          TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                                 writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
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
          } catch (Teuchos::Exceptions::InvalidParameterName& e) {
            stop = true;
          }
        }

      } while (!stop);
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

int main(int argc, char* argv[]) {
  bool success = false;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try {
    const bool throwExceptions     = false;
    const bool recogniseAllOptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
    Xpetra::Parameters xpetraParameters(clp);

    std::string node = "";
    clp.setOption("node", &node, "node type (serial | openmp | cuda | hip)");

    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double, int, int, Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
      if (node == "") {
        typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_ENABLE_SERIAL
        typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_ENABLE_OPENMP
        typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_ENABLE_CUDA
        typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("CUDA node type is disabled");
#endif
      } else if (node == "hip") {
#ifdef KOKKOS_ENABLE_HIP
        typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_HIP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("HIP node type is disabled");
#endif
      } else if (node == "sycl") {
#ifdef KOKKOS_ENABLE_SYCL
        typedef Tpetra::KokkosCompat::KokkosSYCLWrapperNode Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#else
#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double, int, int, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double, int, long, Node>(clp, lib, argc, argv);
#elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SYCL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double, int, long long, Node>(clp, lib, argc, argv);
#else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#endif
#endif
#else
        throw MueLu::Exceptions::RuntimeError("SYCL node type is disabled");
#endif
      } else {
        throw MueLu::Exceptions::RuntimeError("Unrecognized node type");
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

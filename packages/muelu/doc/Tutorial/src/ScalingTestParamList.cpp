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

#include <Tpetra_KokkosCompat_DefaultNode.hpp>  // For Epetra only runs this points to FakeKokkos in Xpetra

#include "Xpetra_ConfigDefs.hpp"
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>

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

// Define default data types
typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

int main(int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  Teuchos::CommandLineProcessor clp(false);

  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "scalingTest.xml";
  clp.setOption("xml", &xmlFileName, "read parameters from a file [default = 'scalingTest.xml']");
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
  bool scaleResidualHistory = true;
  clp.setOption("scale", "noscale", &scaleResidualHistory, "scaled Krylov residual history");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

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
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D")
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 2);
    if (matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO, GO, Node>::Build(map, 3);

    galeriStream << "Processor subdomains in x direction: " << galeriList.get<int>("mx") << std::endl
                 << "Processor subdomains in y direction: " << galeriList.get<int>("my") << std::endl
                 << "Processor subdomains in z direction: " << galeriList.get<int>("mz") << std::endl
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
      map = Utils2::ReadMap(mapFile, xpetraParameters.GetLib(), comm);
    comm->barrier();

    if (lib == Xpetra::UseEpetra) {
      A = Utils::Read(matrixFile, map);

    } else {
      // Tpetra matrix reader is still broken, so instead we read in
      // a matrix in a binary format and then redistribute it
      const bool binaryFormat = true;
      A                       = Utils::Read(matrixFile, lib, comm, binaryFormat);

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

    if (!coordFile.empty())
      coordinates = Utils2::ReadMultiVector(coordFile, map);

    if (!nullFile.empty())
      nullspace = Utils2::ReadMultiVector(nullFile, map);
  }

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
      A->SetMaxEigenvalueEstimate(-one);

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
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));
      //============================================ SPLIT
      RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(mueluList));
      //============================================ SPLIT
      comm->barrier();
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));
      //============================================ SPLIT
      RCP<Hierarchy> H;
      //============================================ SPLIT
      for (int i = 0; i <= numRebuilds; i++) {
        A->SetMaxEigenvalueEstimate(-one);
        //============================================ SPLIT
        H = mueLuFactory->CreateHierarchy();
        H->GetLevel(0)->Set("A", A);
        H->GetLevel(0)->Set("Nullspace", nullspace);
        if (!coordinates.is_null())
          H->GetLevel(0)->Set("Coordinates", coordinates);
        mueLuFactory->SetupHierarchy(*H);
        //============================================ SPLIT
      }

      comm->barrier();
      tm = Teuchos::null;

      // =========================================================================
      // System solution (Ax = b)
      // =========================================================================
      comm->barrier();
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));

      RCP<Vector> X = VectorFactory::Build(map);
      RCP<Vector> B = VectorFactory::Build(map);

      {
        // we set seed for reproducibility
        Utils::SetRandomSeed(*comm);
        X->randomize();
        A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

        Teuchos::Array<STS::magnitudeType> norms(1);
        B->norm2(norms);
        B->scale(one / norms[0]);
        X->putScalar(zero);
      }
      tm = Teuchos::null;

      if (writeMatricesOPT > -2) {
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3.5 - Matrix output")));
        H->Write(writeMatricesOPT, writeMatricesOPT);
        tm = Teuchos::null;
      }

      comm->barrier();
      if (solveType == "none") {
        // Do not perform a solve

      } else if (solveType == "standalone") {
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Fixed Point Solve")));

        H->IsPreconditioner(false);
        H->Iterate(*B, *X, maxIts);

      } else if (solveType == "cg" || solveType == "gmres") {
#ifdef HAVE_MUELU_BELOS
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
        if (!scaleResidualHistory)
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
        try {
          ret = solver->solve();

          // Get the number of iterations for this solve.
          out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

        } catch (const std::exception& ex) {
          out << std::endl
              << "ERROR:  Belos threw an error!  The exception message is:" << std::endl;
          std::cout << ex.what() << std::endl;
        }

        catch (...) {
          out << std::endl
              << "ERROR:  Belos threw an unknown error! " << std::endl;
        }

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

      if (printTimings)
        TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);

      TimeMonitor::clearCounters();

      if (isDriver) {
        if (openedOut != NULL) {
          dup2(savedOut, STDOUT_FILENO);
          fclose(openedOut);
          openedOut = NULL;
        }
        try {
          runList   = paramList.sublist("Run" + MueLu::toString(++runCount), mustAlreadyExist);
          mueluList = runList.sublist("MueLu", mustAlreadyExist);
        } catch (std::exception) {
          stop = true;
        }
      }

    } while (stop == false);
  }

  return 0;
}  // main

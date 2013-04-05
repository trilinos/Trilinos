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
#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp

#include <MueLu_Utilities.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false); // Note:

  GO nx,ny,nz;
  nx=100;
  ny=100;
  nz=100;
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName = "scalingTest.xml"; clp.setOption("xml",                   &xmlFileName,  "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");
  int  amgAsPrecond = 1;                       clp.setOption("precond",               &amgAsPrecond, "apply multigrid as preconditioner");
  int   amgAsSolver = 0;                       clp.setOption("fixPoint",              &amgAsSolver,  "apply multigrid as solver");
  bool printTimings = true;                    clp.setOption("timings", "notimings",  &printTimings, "print timings to screen");
  int  writeMatricesOPT = -2;                  clp.setOption("write",                  &writeMatricesOPT, "write matrices to file (-1 means all; i>=0 means level i)");

  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  if (comm->getRank() == 0) {
    std::cout << "========================================================" << std::endl
              << xpetraParameters << matrixParameters;
  }

  //
  // Construct the problem
  //

  RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
  RCP<TimeMonitor> tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

  RCP<const Map> map;
  RCP<MultiVector> coordinates;

  // Retrieve matrix parameters (they may have been changed on the command line), and pass them to Galeri.
  // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
  //                                 d1  d2  d3
  //                                 d4  d5  d6
  //                                 d7  d8  d9
  //                                 d10 d11 d12
  // A perfect distribution is only possible when the #processors is a perfect square.
  // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
  // size. For example, np=14 will give a 7-by-2 distribution.
  // If you don't want Galeri to do this, specify mx or my on the galeriList.
  Teuchos::ParameterList pl = matrixParameters.GetParameterList();
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", pl.get("nx",nx));
  galeriList.set("ny", pl.get("ny",ny));
  galeriList.set("nz", pl.get("nz",nz));
  //galeriList.set("mx", comm->getSize());
  //galeriList.set("my", 1);

  if (matrixParameters.GetMatrixType() == "Laplace1D") {
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",map,matrixParameters.GetParameterList());
  }
  else if (matrixParameters.GetMatrixType() == "Laplace2D" || matrixParameters.GetMatrixType() == "Star2D" || matrixParameters.GetMatrixType() == "Elasticity2D") {
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D",map,matrixParameters.GetParameterList());
  }
  else if (matrixParameters.GetMatrixType() == "Laplace3D" || matrixParameters.GetMatrixType() == "Elasticity3D") {
    map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D",map,matrixParameters.GetParameterList());
  }
  // Expand map to do multiple DOF per node
  if (matrixParameters.GetMatrixType() == "Elasticity2D")
    map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 2);
  if (matrixParameters.GetMatrixType() == "Elasticity3D")
    map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 3);

  if (comm->getRank() == 0) {
    GO mx = galeriList.get("mx", -1);
    GO my = galeriList.get("my", -1);
    GO mz = galeriList.get("mz", -1);
    std::cout << "Processor subdomains in x direction: " << mx << std::endl
              << "Processor subdomains in y direction: " << my << std::endl
              << "Processor subdomains in z direction: " << mz << std::endl
              << "========================================================" << std::endl;
  }

  Teuchos::ParameterList matrixParams = matrixParameters.GetParameterList();
  matrixParams.set("mx", galeriList.get("mx", -1));
  matrixParams.set("my", galeriList.get("my", -1));
  matrixParams.set("mz", galeriList.get("mz", -1));
  if (matrixParameters.GetMatrixType() == "Elasticity2D" || matrixParameters.GetMatrixType() == "Elasticity3D") {
    matrixParams.set("right boundary" , "Neumann");
    matrixParams.set("bottom boundary", "Neumann");
    matrixParams.set("top boundary"   , "Neumann");
    matrixParams.set("front boundary" , "Neumann");
    matrixParams.set("back boundary"  , "Neumann");
  }

  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParams);
  RCP<Matrix> A = Pr->BuildMatrix();

  RCP<MultiVector> nullspace = MultiVectorFactory::Build(map,1);
  if (matrixParameters.GetMatrixType() == "Elasticity2D" ||
      matrixParameters.GetMatrixType() == "Elasticity3D") {
    nullspace = Pr->BuildNullspace();
    A->SetFixedBlockSize((matrixParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
  } else {
    nullspace->putScalar( (SC) 1.0);
  }

  comm->barrier();

  tm = Teuchos::null;

  //
  // Construct a multigrid preconditioner
  //

  // Multigrid Hierarchy
  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));
  ParameterListInterpreter mueLuFactory(xmlFileName,*comm);

  comm->barrier();

  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

  RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();

  H->SetDefaultVerbLevel(MueLu::Extreme);


  H->GetLevel(0)->Set("A", A);

  H->GetLevel(0)->Set("Nullspace", nullspace);
  H->GetLevel(0)->Set("Coordinates", coordinates);

  mueLuFactory.SetupHierarchy(*H);

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(-1);
  H->print(*fos);

  tm = Teuchos::null;
  //
  // Solve Ax = b
  //

  tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));
  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);

  {
    X->setSeed(846930886);
    X->randomize();
    A->apply(*X,*B,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Teuchos::Array<ST::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(1.0/norms[0]);
    X->putScalar( (SC) 0.0);
  }
  tm = Teuchos::null;

  if (amgAsSolver) {
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Fixed Point Solve")));

    H->IsPreconditioner(false);
    H->Iterate(*B,25,*X);


    tm = Teuchos::null;
  }

  if (writeMatricesOPT > -2) {
    H->Write(writeMatricesOPT,writeMatricesOPT);
  }

#ifdef HAVE_MUELU_BELOS
  if (amgAsPrecond) {
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 5 - Belos Solve")));
    // Operator and Multivector type that will be used with Belos
    typedef MultiVector          MV;
    typedef Belos::OperatorT<MV> OP;
    H->IsPreconditioner(true);

    // Define Operator and Preconditioner
    Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(A)); // Turns a Xpetra::Matrix object into a Belos operator
    Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO, LMO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setLeftPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      if (comm->getRank() == 0)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    int maxIts = 100;
    double tol = 1e-12; // FIXME: use command line option
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
    //belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency",1);
    belosList.set("Output Style",Belos::Brief);

    // Create an iterative solver manager
    RCP< Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

    // Perform solve
    Belos::ReturnType ret=Belos::Unconverged;
    try {
      ret = solver->solve();

      // Get the number of iterations for this solve.
      if (comm->getRank() == 0)
        std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

    } //try

    catch(...) {
      if (comm->getRank() == 0)
        std::cout << std::endl << "ERROR:  Belos threw an error! " << std::endl;
    }

    // Check convergence
    if (ret != Belos::Converged) {
      if (comm->getRank() == 0) std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    } else {
      if (comm->getRank() == 0) std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
    }
    tm = Teuchos::null;
  } //if (amgAsPrecond)
#endif //ifdef HAVE_MUELU_BELOS

  globalTimeMonitor = Teuchos::null;

  if (printTimings)
    TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union);

  MueLu::MutuallyExclusiveTime<MueLu::BaseClass>::PrintParentChildPairs();


} //main

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
#include <unistd.h>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_PermutedTransferFactory.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosXpetraAdapter.hpp" // this header defines Belos::XpetraOp()
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuOp()
#endif

//
typedef double Scalar;
typedef int    LocalOrdinal;
//FIXME we need a HAVE_MUELU_LONG_LONG_INT option
// #ifdef HAVE_TEUCHOS_LONG_LONG_INT
// typedef long long int GlobalOrdinal;
// #else
typedef int GlobalOrdinal;
//#endif
//
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
//
#include "MueLu_UseShortNames.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  //using Galeri::Xpetra::CreateCartesianCoordinates;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  //out->setOutputToRootOnly(-1);
  //out->precision(12);

//FIXME we need a HAVE_MUELU_LONG_LONG_INT option
//#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
//#endif

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  // custom parameters
  LO maxLevels = 10; //was 3 in simple
  LO its=10;
  std::string smooType="sgs";
  int pauseForDebugger=0;
  int amgAsSolver=1;
  int amgAsPrecond=1;
  int useExplicitR=1;
  int sweeps=2;
  int maxCoarseSize=50;  //FIXME clp doesn't like long long int
  Scalar SADampingFactor=4./3;
  double tol = 1e-7;
  std::string aggOrdering = "natural";
  int minPerAgg=2; //was 3 in simple
  int maxNbrAlreadySelected=0;
  int writeMatrix=0;
  int printTimings=0;
  LO minRowsPerProc=2000;
  double nonzeroImbalance=1.2;

  clp.setOption("aggOrdering",&aggOrdering,"aggregation ordering strategy (natural,random,graph)");
    clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");
  clp.setOption("dump",&writeMatrix,"write matrix to file");
  clp.setOption("explicitR",&useExplicitR,"restriction will be explicitly stored as transpose of prolongator");
    clp.setOption("fixPoint",&amgAsSolver,"apply multigrid as solver");
  clp.setOption("its",&its,"number of multigrid cycles");
  clp.setOption("maxCoarseSize",&maxCoarseSize,"maximum #dofs in coarse operator");
  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  clp.setOption("maxNbrSel",&maxNbrAlreadySelected,"maximum # of nbrs allowed to be in other aggregates");
  clp.setOption("minPerAgg",&minPerAgg,"minimum #DOFs per aggregate");
  clp.setOption("minRowsPerProc",&minRowsPerProc,"min #rows allowable per proc before repartitioning occurs");
  clp.setOption("nnzImbalance",&nonzeroImbalance,"max allowable nonzero imbalance before repartitioning occurs");
    clp.setOption("precond",&amgAsPrecond,"apply multigrid as preconditioner");
    clp.setOption("saDamping",&SADampingFactor,"prolongator damping factor");
  clp.setOption("smooType",&smooType,"smoother type ('l1-sgs', 'sgs 'or 'cheby')");
    clp.setOption("sweeps",&sweeps,"sweeps to be used in SGS (or Chebyshev degree)");
  clp.setOption("timings",&printTimings,"print timings to screen");
    clp.setOption("tol",&tol,"stopping tolerance for Krylov method");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));

  if (pauseForDebugger) {
    Utils::PauseForDebugger();
  }

  matrixParameters.check();
  xpetraParameters.check();
  // TODO: check custom parameters
  std::transform(smooType.begin(), smooType.end(), smooType.begin(), ::tolower);
  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  if (comm->getRank() == 0) {
    std::cout << xpetraParameters << matrixParameters;
    // TODO: print custom parameters // Or use paramList::print()!
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  RCP<const Map> map;
  RCP<Matrix> A;

  RCP<MultiVector> Coordinates;
  {
    TimeMonitor tm(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build"));

    map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);
    A = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Matrix vs. CrsMatrixWrap

    if (matrixParameters.GetMatrixType() == "Laplace1D") {
      Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",map,matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace2D") {
      Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D",map,matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace3D") {
      Coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D",map,matrixParameters.GetParameterList());
    }
  }
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  // dump matrix to file
  if (writeMatrix) {
    std::string fileName = "Amat.mm";
    Utils::Write(fileName,*A);
  }

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);

  nullSpace->norm1(norms);
  if (comm->getRank() == 0)
    std::cout << "||NS|| = " << norms[0] << std::endl;

  Teuchos::ParameterList status;
  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H;
  {
    TimeMonitor tm(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup"));

    H = rcp( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

    RCP<MueLu::Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",           A);
    Finest->Set("Nullspace",   nullSpace);
    Finest->Set("Coordinates", Coordinates); //FIXME: XCoordinates, YCoordinates, ..

    FactoryManager M;
    M.SetFactory("Graph", rcp(new CoalesceDropFactory()));        /* do not use the permuted nullspace (otherwise, circular dependencies) */

    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    *out << "========================= Aggregate option summary  =========================" << std::endl;
    *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
    *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
    UCAggFact->SetMinNodesPerAggregate(minPerAgg);  //TODO should increase if run anything other than 1D
    UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
    std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
    if (aggOrdering == "natural") {
      *out << "aggregate ordering :                    NATURAL" << std::endl;
      UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
    } else if (aggOrdering == "random") {
      *out << "aggregate ordering :                    RANDOM" << std::endl;
      UCAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
    } else if (aggOrdering == "graph") {
      *out << "aggregate ordering :                    GRAPH" << std::endl;
      UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
    } else {
      std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
      throw(MueLu::Exceptions::RuntimeError(msg));
    }
    UCAggFact->SetPhase3AggCreation(0.5);
    M.SetFactory("Aggregates", UCAggFact);
    *out << "=============================================================================" << std::endl;

    //  M.SetFactory("Ptent", rcp(new TentativePFactory()));

    RCP<SaPFactory> SaPfact = rcp(new SaPFactory() );
    SaPfact->SetDampingFactor(SADampingFactor);
    M.SetFactory("P", SaPfact);
    RCP<FactoryBase2> rfact = rcp(new TransPFactory());
    rfact->SetFactory("P", M.GetFactory("P"));
    M.SetFactory("R", rfact);

    RCP<RAPFactory> Acfact = rcp(new RAPFactory());
    Acfact->SetFactory("P", M.GetFactory("P"));
    Acfact->SetFactory("R", M.GetFactory("R"));

    Acfact->setVerbLevel(Teuchos::VERB_HIGH);
    M.SetFactory("A", Acfact);

    RCP<RAPFactory>       AcfactFinal;

    RCP<PermutedTransferFactory> permPFactory, permRFactory;
    RCP<MultiVectorTransferFactory> mvTransFact;
    if (useExplicitR) {

//FIXME we need a HAVE_MUELU_LONG_LONG_INT option
// #if defined(HAVE_TEUCHOS_LONG_LONG_INT)
//       // for long long
//       AcfactFinal = Acfact;
// #else
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
      //Matrix used to transfer coordinates to coarse grid
      RCP<const FactoryBase> Rtentfact = M.GetFactory("R"); //for projecting coordinates
      //Factory that will invoke the coordinate transfer. This factory associates data and operator.
      mvTransFact = rcp(new MultiVectorTransferFactory("Coordinates","R"));
      mvTransFact->SetFactory("R", Rtentfact);
      //Register it with the coarse operator factory
      Acfact->AddTransferFactory(mvTransFact);
      //Set up repartitioning
      RCP<ZoltanInterface>      zoltan = rcp(new ZoltanInterface());
      zoltan->SetFactory("A", Acfact);
      zoltan->SetFactory("Coordinates", mvTransFact);
      RCP<RepartitionFactory> RepartitionFact = rcp(new RepartitionFactory(minRowsPerProc,nonzeroImbalance));
      RepartitionFact->SetFactory("Partition", zoltan);
      RepartitionFact->SetFactory("A", Acfact);
      permPFactory = rcp( new PermutedTransferFactory(MueLu::INTERPOLATION));
      permPFactory->SetFactory("A", Acfact);
      permPFactory->SetFactory("P", SaPfact);
      permPFactory->SetFactory("Importer", RepartitionFact);

      permRFactory = rcp( new PermutedTransferFactory(MueLu::RESTRICTION));
      permRFactory->SetFactory("A", Acfact);
      permPFactory->SetFactory("P", M.GetFactory("Ptent"));
      permRFactory->SetFactory("R", M.GetFactory("R"));
      permRFactory->SetFactory("Importer", RepartitionFact);
      permRFactory->SetFactory("Nullspace", mvTransFact);

      AcfactFinal = rcp(new RAPFactory());
      AcfactFinal->SetFactory("P", permPFactory);
      AcfactFinal->SetFactory("R", permRFactory);
      M.SetFactory("A", AcfactFinal);
#else
      AcfactFinal = Acfact;
#endif
      //#endif // TEUCHOS_LONG_LONG_INT
    } else {

        H->SetImplicitTranspose(true);
        Acfact->SetImplicitTranspose(true);
      AcfactFinal = Acfact;
        if (comm->getRank() == 0) std::cout << "\n\n* ***** USING IMPLICIT RESTRICTION OPERATOR ***** *\n" << std::endl;
    } //if (useExplicitR)

    H->SetMaxCoarseSize((GO) maxCoarseSize);

    RCP<SmootherPrototype> smooProto;
    std::string ifpackType;
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (LO) sweeps);
    ifpackList.set("relaxation: damping factor", (SC) 1.0);
    if (smooType == "sgs") {
      ifpackType = "RELAXATION";
      ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    }
    else if (smooType == "l1-sgs") {
      ifpackType = "RELAXATION";
      ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
      ifpackList.set("relaxation: use l1",true);
    } else if (smooType == "cheby") {
      ifpackType = "CHEBYSHEV";
      ifpackList.set("chebyshev: degree", (LO) sweeps);

      if (matrixParameters.GetMatrixType() == "Laplace1D") {
	ifpackList.set("chebyshev: ratio eigenvalue", (SC) 3);
      }
      else if (matrixParameters.GetMatrixType() == "Laplace2D") {
	ifpackList.set("chebyshev: ratio eigenvalue", (SC) 7);
      }
      else if (matrixParameters.GetMatrixType() == "Laplace3D") {
	ifpackList.set("chebyshev: ratio eigenvalue", (SC) 20);
      }
      // ifpackList.set("chebyshev: max eigenvalue", (double) -1.0);
      // ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
    }

    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );

    /* test direct solve */
    /*    if (maxLevels == 1) {
          Teuchos::ParameterList amesosList;
          amesosList.set("PrintTiming", true);
          smooProto = rcp( new DirectSolver("", amesosList) );
          }
    */

    RCP<SmootherFactory> SmooFact;
    SmooFact = rcp( new SmootherFactory(smooProto) );
    AcfactFinal->setVerbLevel(Teuchos::VERB_HIGH);

    /*
    for (int i=0; i<comm->getSize(); ++i) {
      if (comm->getRank() == i) {
        std::cout << "pid " << i << ": address(Acfact)      = " << Acfact.get() << std::endl;
        std::cout << "pid " << i << ": address(AcfactFinal) = " << AcfactFinal.get() << std::endl;
      }
      comm->barrier();
    }
    */

    if (maxLevels > 1) {
      M.SetFactory("A",AcfactFinal);
      M.SetFactory("Smoother",SmooFact);
      M.SetFactory("Aggregates", UCAggFact);

      if (useExplicitR) {
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
        M.SetFactory("P",permPFactory);
        M.SetFactory("R",permRFactory);
        M.SetFactory("Nullspace", permRFactory);
#endif
      }

      int startLevel=0;
      //      std::cout << startLevel << " " << maxLevels << std::endl;
      H->Setup(M,startLevel,maxLevels);
    }//maxLevels>1


    if (maxLevels == 1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Mue!");
      // H->SetCoarsestSolver(*SmooFact, MueLu::BOTH);
    }

  } // end of Setup TimeMonitor

  //   *out  << "======================\n Multigrid statistics \n======================" << std::endl;
  //   status.print(*out,Teuchos::ParameterList::PrintOptions().indent(2));

  // Define B
  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> B = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  A->apply(*X,*B,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  B->norm2(norms);
  B->scale(1.0/norms[0]);

#define AMG_SOLVER
#ifdef AMG_SOLVER
  // Use AMG directly as an iterative method
  if (amgAsSolver) {
    //*out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    {
      X->putScalar( (SC) 0.0);

      TimeMonitor tm(*TimeMonitor::getNewTimer("ScalingTest: 3 - Fixed Point Solve"));

      H->IsPreconditioner(false);
      H->Iterate(*B,its,*X);

      //X->norm2(norms);
      //*out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
    }
  } // if (fixedPt)
#endif //ifdef AMG_SOLVER

  // Use AMG as a preconditioner in Belos
#ifdef HAVE_MUELU_BELOS
  if (amgAsPrecond) {
    RCP<TimeMonitor> tm;
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 5 - Belos Solve")));
    // Operator and Multivector type that will be used with Belos
    typedef MultiVector          MV;
    typedef Belos::OperatorT<MV> OP;
    H->IsPreconditioner(true);

    // Define Operator and Preconditioner
    Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(A)); // Turns a Xpetra::Operator object into a Belos operator
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
      {
        TimeMonitor tm2(*TimeMonitor::getNewTimer("ScalingTest: 5bis - Belos Internal Solve"));
        ret = solver->solve();
      } // end of TimeMonitor

      // Get the number of iterations for this solve.
      if (comm->getRank() == 0)
        std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

      // Compute actual residuals.
      int numrhs = 1;
      std::vector<double> actual_resids( numrhs ); //TODO: double?
      std::vector<double> rhs_norm( numrhs );
      RCP<MultiVector> resid = MultiVectorFactory::Build(map, numrhs);

      typedef Belos::OperatorTraits<SC,MV,OP>  OPT;
      typedef Belos::MultiVecTraits<SC,MV>     MVT;

      OPT::Apply( *belosOp, *X, *resid );
      MVT::MvAddMv( -1.0, *resid, 1.0, *B, *resid );
      MVT::MvNorm( *resid, actual_resids );
      MVT::MvNorm( *B, rhs_norm );
      *out<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        *out<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        //if (actRes > tol) { badRes = true; }
      }

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
#endif // HAVE_MUELU_BELOS

  // Timer final summaries
  globalTimeMonitor = Teuchos::null; // stop this timer before summary

  if (printTimings)
    TimeMonitor::summarize();

  return EXIT_SUCCESS;
}

// TODO: add warning if:
// DEBUG_MODE, LONG_LONG or KLU

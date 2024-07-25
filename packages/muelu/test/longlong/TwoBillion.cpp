// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <unistd.h>
#include <iostream>

// Teuchos
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_TestingHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include "Xpetra_ConfigDefs.hpp"
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

// MueLu
#include <MueLu_ConfigDefs.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosXpetraAdapter.hpp"  // this header defines Belos::XpetraOp()
#include "BelosMueLuAdapter.hpp"   // this header defines Belos::MueLuOp()
#endif

typedef double Scalar;
typedef int LocalOrdinal;
typedef long long GlobalOrdinal;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

using Xpetra::global_size_t;

int main(int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);

    Teuchos::CommandLineProcessor clp(false);
    std::string xmlFileName = "TwoBillion.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file");
    int num_per_proc = 1000;
    clp.setOption("dpc", &num_per_proc, "DOFs per core");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    int NumProcs = comm->getSize();
    int MyPID    = comm->getRank();

    if (!MyPID) printf("TwoBillion: Running Test\n");

    const long long FIRST_GID = 3000000000L;
    //  const long long FIRST_GID  = 0L;

    const long long IndexBase = 0L;
    //  const long long IndexBase  = 3000000000L;

    global_size_t NumGlobalElements = NumProcs * num_per_proc;

    // Create Map w/ GIDs starting at > 2 billion
    RCP<const Map> map;
    RCP<CrsMatrix> Acrs;
    Teuchos::Array<GlobalOrdinal> mygids(num_per_proc);

    for (int i = 0; i < num_per_proc; i++)
      mygids[i] = FIRST_GID + MyPID * num_per_proc + i;

    for (int i = 0; i < NumProcs; ++i) {
      if (i == MyPID)
        std::cout << "pid " << i << " : 1st GID = " << mygids[0] << std::endl;
    }

    // for(int i=0;i<num_per_proc;i++)
    //   printf("[%d] mygids[%d] = %lld\n",MyPID,i,mygids[i]);

    map = MapFactory::Build(Xpetra::UseTpetra, Teuchos::OrdinalTraits<global_size_t>::invalid(), mygids(), IndexBase, comm);

    //  RCP<Teuchos::FancyOStream> fox = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    //  fox->setOutputToRootOnly(-1);
    //  map->describe(*fox,Teuchos::VERB_EXTREME);

    // Create 1D Laplacian w/ GIDs starting at > 2 billion
    Teuchos::Array<Scalar> myvals(3);
    Teuchos::Array<GlobalOrdinal> mycols(3);
    Teuchos::ArrayView<Scalar> ValView;
    Teuchos::ArrayView<GlobalOrdinal> ColView;

    Acrs = CrsMatrixFactory::Build(map, 3);
    for (int i = 0; i < num_per_proc; i++) {
      if (mygids[i] == FIRST_GID) {
        mycols[0] = mygids[i];
        myvals[0] = 2;
        mycols[1] = mygids[i] + 1;
        myvals[1] = -1;
        ValView   = myvals.view(0, 2);
        ColView   = mycols.view(0, 2);
        //      printf("[%d %lld] cols %lld %lld\n",MyPID,mygids[i],mycols[0],mycols[1]);
      } else if (mygids[i] == FIRST_GID + (long long)NumGlobalElements - 1) {
        mycols[0] = mygids[i] - 1;
        myvals[0] = -1;
        mycols[1] = mygids[i];
        myvals[1] = 2;
        ValView   = myvals.view(0, 2);
        ColView   = mycols.view(0, 2);
        //      printf("[%d %lld] cols %lld %lld\n",MyPID,mygids[i],mycols[0],mycols[1]);
      } else {
        mycols[0] = mygids[i] - 1;
        myvals[0] = -1;
        mycols[1] = mygids[i];
        myvals[1] = -2;
        mycols[2] = mygids[i] + 1;
        myvals[1] = -1;
        ValView   = myvals();
        ColView   = mycols();
        //      printf("[%d %lld] cols %lld %lld %lld\n",MyPID,mygids[i],mycols[0],mycols[1],mycols[2]);
      }
      Acrs->insertGlobalValues(mygids[i], ColView, ValView);
    }
    Acrs->fillComplete();

    RCP<Matrix> A = rcp(new CrsMatrixWrap(Acrs));

    RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
    nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());
    RCP<MultiVector> coordinates        = MultiVectorFactory::Build(map, 1);
    Teuchos::ArrayRCP<Scalar> coordVals = coordinates->getDataNonConst(0);
    double h                            = 1.0 / NumGlobalElements;
    for (LocalOrdinal i = 0; i < num_per_proc; ++i)
      coordVals[i] = MyPID * num_per_proc * h + i * h;
    coordVals = Teuchos::null;

    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), *comm);

    RCP<HierarchyManager> mueLuFactory;
    mueLuFactory = rcp(new ParameterListInterpreter(xmlFileName, *comm));

    RCP<Hierarchy> H;
    H = mueLuFactory->CreateHierarchy();
    H->GetLevel(0)->Set("A", A);
    H->GetLevel(0)->Set("Nullspace", nullspace);
    H->GetLevel(0)->Set("Coordinates", coordinates);
    mueLuFactory->SetupHierarchy(*H);

    //
    //
    // SOLVE
    //
    //

    // Define X, B
    RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
    RCP<MultiVector> B = MultiVectorFactory::Build(map, 1);
    Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);

    X->setSeed(846930886);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);
    B->norm2(norms);
    B->scale(1.0 / norms[0]);

    //
    // Use AMG as a preconditioner in Belos
    //
#ifdef HAVE_MUELU_BELOS
    // Operator and Multivector type that will be used with Belos
    typedef MultiVector MV;
    typedef Belos::OperatorT<MV> OP;
    H->IsPreconditioner(true);

    // Define Operator and Preconditioner
    Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));  // Turns a Xpetra::Operator object into a Belos operator
    Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));   // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setLeftPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      if (comm->getRank() == 0)
        std::cout << std::endl
                  << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    int maxIts    = 100;
    double optTol = 1e-8;
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations", maxIts);     // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", optTol);  // Relative convergence tolerance requested
    // belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency", 1);
    belosList.set("Output Style", Belos::Brief);

    // Create an iterative solver manager
    RCP<Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

    // Perform solve
    Belos::ReturnType ret = Belos::Unconverged;
    try {
      ret = solver->solve();

      // Get the number of iterations for this solve.
      if (comm->getRank() == 0)
        std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

      if (solver->getNumIters() > 6) {
        if (comm->getRank() == 0) std::cout << std::endl
                                            << "ERROR:  Belos did not converge! " << std::endl;
        return (EXIT_FAILURE);
      }

      // Compute actual residuals.
      int numrhs = 1;
      std::vector<double> actual_resids(numrhs);  // TODO: double?
      std::vector<double> rhs_norm(numrhs);
      RCP<MultiVector> resid = MultiVectorFactory::Build(map, numrhs);

      typedef Belos::OperatorTraits<SC, MV, OP> OPT;
      typedef Belos::MultiVecTraits<SC, MV> MVT;

      OPT::Apply(*belosOp, *X, *resid);
      MVT::MvAddMv(-1.0, *resid, 1.0, *B, *resid);
      MVT::MvNorm(*resid, actual_resids);
      MVT::MvNorm(*B, rhs_norm);
      *out << "---------- Actual Residuals (normalized) ----------" << std::endl
           << std::endl;
      for (int i = 0; i < numrhs; i++) {
        double actRes = actual_resids[i] / rhs_norm[i];
        *out << "Problem " << i << " : \t" << actRes << std::endl;
        // if (actRes > tol) { badRes = true; }
      }

    }  // try

    catch (...) {
      if (comm->getRank() == 0)
        std::cout << std::endl
                  << "ERROR:  Belos threw an error! " << std::endl;
    }

    success = (ret == Belos::Converged);
    // Check convergence
    if (success) {
      if (comm->getRank() == 0) std::cout << std::endl
                                          << "SUCCESS:  Belos converged!" << std::endl;
    } else {
      if (comm->getRank() == 0) std::cout << std::endl
                                          << "ERROR:  Belos did not converge! " << std::endl;
    }
#endif  // HAVE_MUELU_BELOS
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

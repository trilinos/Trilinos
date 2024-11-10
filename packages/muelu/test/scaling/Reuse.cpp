// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_IO.hpp>

#include <Teuchos_Time.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp

#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#endif
//- -- --------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  typedef Teuchos::ScalarTraits<SC> STS;

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int mypid                           = comm->getRank();

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(cout));
    Teuchos::FancyOStream &out       = *fancy;
    out.setOutputToRootOnly(0);

    SC zero = STS::zero(), one = STS::one();

    //
    // Parameters
    //

    Xpetra::Parameters xpetraParameters(clp);  // manage parameters of Xpetra

    std::string xmlFileName = "reuse.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'reuse.xml'");
    std::string matrixPrefix = "jac";
    clp.setOption("matrix", &matrixPrefix, "prefix for matrix file names.  Default = 'jac'");
    bool binary = false;
    clp.setOption("binary", "nobinary", &binary, "matrix files are binary. Default = false");
    std::string rhsPrefix = "rhs";
    clp.setOption("rhs", &rhsPrefix, "prefix for rhs file names. Default = 'rhs'");
    bool printTimings = true;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    int first_matrix = 0;
    clp.setOption("firstMatrix", &first_matrix, "first matrix in the sequence to use");
    int last_matrix = 1;
    clp.setOption("lastMatrix", &last_matrix, "last matrix in the sequence to use");
    bool inMemory = false;
    clp.setOption("inmemory", "noinmemory", &inMemory, "load all matrices in memory");
    bool doReuse = true;
    clp.setOption("reuse", "noreuse", &doReuse, "if you want to try reuse");
    const int numPDEs = 2;

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Global Time"))), tm;
    RCP<Time> timer;

    // Operator and Multivector type that will be used with Belos
    typedef MultiVector MV;
    typedef Belos::OperatorT<MV> OP;

    // Stats tracking
    int numSteps = last_matrix - first_matrix + 1;
    Array<Array<int> > iteration_counts(numSteps);
    Array<Array<double> > iteration_times(numSteps);
    Array<Array<double> > setup_times(numSteps);
    Array<RCP<Matrix> > matrices(numSteps);
    Array<RCP<MultiVector> > rhss(numSteps);

    for (int i = 0; i < numSteps; i++) {
      iteration_counts[i].resize(numSteps);
      iteration_times[i].resize(numSteps);
      setup_times[i].resize(numSteps);
    }

    ParameterListInterpreter mueLuFactory(xmlFileName, *comm);

    for (int i = first_matrix; i <= last_matrix; i++) {
      out << "==================================================================" << std::endl;

      char matrixFileName[80];
      char rhsFileName[80];
      char timerName[80];

      sprintf(matrixFileName, "%s%d.mm", matrixPrefix.c_str(), i);

      // Load the matrix
      RCP<Matrix> Aprecond = matrices[i];
      if (Aprecond.is_null()) {
        out << "[" << i << "] Loading matrix \"" << matrixFileName << "\"... ";
        Aprecond = Xpetra::IO<SC, LO, GO, Node>::Read(std::string(matrixFileName), xpetraParameters.GetLib(), comm, binary);
        out << "done" << std::endl;

        Aprecond->SetFixedBlockSize(numPDEs);

        if (inMemory)
          matrices[i] = Aprecond;

      } else {
        // Reset estimate
        Aprecond->SetMaxEigenvalueEstimate(-one);
      }

      // Build the nullspace
      RCP<MultiVector> nullspace = MultiVectorFactory::Build(Aprecond->getRowMap(), numPDEs);
      nullspace->putScalar(zero);

      Teuchos::ArrayRCP<SC> data0, data1;
      data0 = nullspace->getDataNonConst(0);
      data1 = nullspace->getDataNonConst(1);
      for (size_t k = 0; k < Aprecond->getRowMap()->getLocalNumElements(); k += 2)
        data0[k + 0] = data1[k + 1] = one;

      // Build the preconditioner
      out << "[" << i << "] Building preconditioner \"" << matrixFileName << "\" ..." << std::endl;
      sprintf(timerName, "Reuse: Preconditioner Setup i=%d", i);
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(timerName)));

      RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
      H->SetDefaultVerbLevel(MueLu::Extreme);
      H->GetLevel(0)->Set("A", Aprecond);
      H->GetLevel(0)->Set("Nullspace", nullspace);
      H->IsPreconditioner(true);

      sprintf(timerName, "Reuse: Setup i=%d j=%d", i, i);
      timer = TimeMonitor::getNewTimer(timerName);
      timer->start();
      mueLuFactory.SetupHierarchy(*H);
      setup_times[i - first_matrix][i - first_matrix] = timer->stop();
      timer                                           = Teuchos::null;

      Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator
      tm                         = Teuchos::null;

      // Loop over all future matrices
      int j = i;
      for (; j <= last_matrix; j++) {
        out << "------------------------------------------------------------------" << std::endl;

        sprintf(matrixFileName, "%s%d.mm", matrixPrefix.c_str(), j);
        sprintf(rhsFileName, "%s%d.mm", rhsPrefix.c_str(), j);

        RCP<Matrix> Amatvec = matrices[j];
        if (j != i) {
          if (Amatvec.is_null()) {
            // Load the matrix
            out << "[" << j << "]<-[" << i << "] Loading matrix \"" << matrixFileName << "\"... ";
            Amatvec = Xpetra::IO<SC, LO, GO, Node>::Read(std::string(matrixFileName), xpetraParameters.GetLib(), comm, binary);
            out << "done" << std::endl;

            Amatvec->SetFixedBlockSize(numPDEs);

            if (inMemory)
              matrices[j] = Amatvec;

          } else {
            // Reset estimate
            Amatvec->SetMaxEigenvalueEstimate(-one);
          }

          // Preconditioner update
          sprintf(timerName, "Reuse: Preconditioner Update i=%d", i);
          tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(timerName)));
          // No-op at present

          sprintf(timerName, "Reuse: Setup i=%d j=%d", i, j);
          timer = TimeMonitor::getNewTimer(timerName);
          timer->start();
          out << "[" << j << "]<-[" << i << "] Updating preconditioner \"" << matrixFileName << "\" ..." << std::endl;

          H->GetLevel(0)->Set("A", Amatvec);
          mueLuFactory.SetupHierarchy(*H);

          setup_times[i - first_matrix][j - first_matrix] = timer->stop();
          timer                                           = Teuchos::null;

          tm = Teuchos::null;

        } else {
          Amatvec = Aprecond;
        }

        // Load the RHS
        RCP<MultiVector> rhs = rhss[j];
        if (rhs.is_null()) {
          out << "[" << j << "] Loading rhs " << rhsFileName << "\"... ";
          rhs = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(std::string(rhsFileName), Amatvec->getRowMap());
          out << "done" << std::endl;

          if (inMemory)
            rhss[j] = rhs;
        }

        // Create an LHS
        RCP<Vector> X = VectorFactory::Build(Amatvec->getRowMap());
        X->putScalar(0.0);

        // Define Operator and Preconditioner
        Teuchos::RCP<OP> belosOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(Amatvec));  // Turns a Xpetra::Matrix object into a Belos operator

        // Construct a Belos LinearProblem object
        RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, rhs));
        belosProblem->setLeftPrec(belosPrec);
        belosProblem->setProblem();

        // Belos parameter list
        int maxIts = 100;
        double tol = 1e-4;
        Teuchos::ParameterList belosList;
        belosList.set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
        belosList.set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
        // belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
        belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
        belosList.set("Output Frequency", 1);
        belosList.set("Output Style", Belos::Brief);

        // Create an iterative solver manager
        RCP<Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

        // Perform solve
        sprintf(timerName, "Reuse: Solve i=%d j=%d", i, j);
        timer = TimeMonitor::getNewTimer(timerName);
        timer->start();
        Belos::ReturnType ret = Belos::Unconverged;
        ret                   = solver->solve();

        double my_time = timer->stop();
        timer          = Teuchos::null;

        // Get the number of iterations for this solve.

        // Check convergence
        if (ret != Belos::Converged) {
          out << "ERROR  :  Belos did not converge! " << std::endl;
          break;
        } else {
          out << "SUCCESS:  Belos converged!" << std::endl;
          out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
          iteration_counts[i - first_matrix][j - first_matrix] = solver->getNumIters();
          iteration_times[i - first_matrix][j - first_matrix]  = my_time;
        }

        if (doReuse == 0) {
          j++;
          break;
        }
      }  // end j
      for (; j <= last_matrix; j++)
        setup_times[i - first_matrix][j - first_matrix] = Teuchos::ScalarTraits<double>::nan();
    }  // end i

    globalTimeMonitor = Teuchos::null;

    if (printTimings)
      // TimeMonitor::summarize(comm.ptr(), out, false, true, false, Teuchos::Union);
      TimeMonitor::summarize(comm.ptr(), out);

    if (!mypid) {
      printf("************************* Iteration Counts ***********************\n");
      for (int i = 0; i < numSteps; i++) {
        for (int j = 0; j < i; j++) printf("        ");
        for (int j = i; j < numSteps; j++) {
          if (STS::isnaninf(setup_times[i][j])) {
            if (i == j)
              printf("       -");
            break;
          }
          printf(" %7d", iteration_counts[i][j]);
        }
        printf("\n");
      }
      // For convenince, print data without reuse in a single line
      for (int i = 0; i < numSteps; i++)
        printf(" %7d", iteration_counts[i][i]);
      printf("\n");

      printf("************************* Iteration Times ***********************\n");
      for (int i = 0; i < numSteps; i++) {
        for (int j = 0; j < i; j++) printf("        ");
        for (int j = i; j < numSteps; j++) {
          if (STS::isnaninf(setup_times[i][j]))
            break;
          printf(" %7.2f", iteration_times[i][j]);
        }
        printf("\n");
      }
      // For convenince, print data without reuse in a single line
      for (int i = 0; i < numSteps; i++)
        printf(" %7.2f", iteration_times[i][i]);
      printf("\n");

      printf("************************* Setup Times ***********************\n");
      for (int i = 0; i < numSteps; i++) {
        for (int j = 0; j < i; j++) printf("        ");
        for (int j = i; j < numSteps; j++) {
          if (STS::isnaninf(setup_times[i][j]))
            break;
          printf(" %7.2f", setup_times[i][j]);
        }
        printf("\n");
      }
      // For convenince, print data without reuse in a single line
      for (int i = 0; i < numSteps; i++)
        printf(" %7.2f", setup_times[i][i]);
      printf("\n");

      printf("************************* Total Times ***********************\n");
      for (int i = 0; i < numSteps; i++) {
        for (int j = 0; j < i; j++) printf("        ");
        for (int j = i; j < numSteps; j++) {
          if (STS::isnaninf(setup_times[i][j]))
            break;
          printf(" %7.2f", setup_times[i][j] + iteration_times[i][j]);
        }
        printf("\n");
      }
      // For convenince, print data without reuse in a single line
      for (int i = 0; i < numSteps; i++)
        printf(" %7.2f", setup_times[i][i] + iteration_times[i][i]);
      printf("\n");
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}

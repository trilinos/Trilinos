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
#include <Xpetra_DefaultPlatform.hpp>

#include <Teuchos_Time.hpp>

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
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int mypid = comm->getRank();

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(cout));
  Teuchos::FancyOStream& cout = *fancy;
  cout.setOutputToRootOnly(0);

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false);
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName  = "reuse.xml"; clp.setOption("xml",                   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'reuse.xml'");
  std::string matrixPrefix = "jac";       clp.setOption("matrix",               &matrixPrefix, "prefix for matrix file names.  Default = 'jac'");
  bool        binary       = false;       clp.setOption("binary", "nobinary",         &binary, "matrix files are binary. Default = false");
  std::string rhsPrefix    = "rhs";       clp.setOption("rhs" ,                    &rhsPrefix, "prefix for rhs file names. Default = 'rhs'");
  bool        printTimings = true;        clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  int         first_matrix = 0;           clp.setOption("firstMatrix",          &first_matrix, "first matrix in the sequence to use");
  int         last_matrix  = 1;           clp.setOption("lastMatrix",            &last_matrix, "last matrix in the sequence to use");
  string      do_reuse_str = "none";      clp.setOption("doReuse",              &do_reuse_str, "if you want to try reuse");
  const int   numPDEs      = 2;

  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  int do_reuse = 0;
  if      (!strcmp(do_reuse_str.c_str(), "none"))   do_reuse = 0;
  else if (!strcmp(do_reuse_str.c_str(), "simple")) do_reuse = 1;
  else if (!strcmp(do_reuse_str.c_str(), "fast"))   do_reuse = 2;
  else return EXIT_FAILURE;


  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time"))), tm;
  RCP<Time> timer;

  // Operator and Multivector type that will be used with Belos
  typedef MultiVector          MV;
  typedef Belos::OperatorT<MV> OP;

  // Stats tracking
  int ArraySize = last_matrix - first_matrix + 1;
  Array<Array<int>    > iteration_counts(ArraySize);
  Array<Array<double> > iteration_times(ArraySize);
  Array<Array<double> > setup_times(ArraySize);

  for (int i = 0; i < ArraySize; i++) {
    iteration_counts[i].resize(ArraySize);
    iteration_times[i] .resize(ArraySize);
    setup_times[i]     .resize(ArraySize);
  }

  ParameterListInterpreter mueLuFactory(xmlFileName, *comm);

  SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();
  for (int i = first_matrix; i <= last_matrix; i++) {
    cout << "==================================================================" << std::endl;

    char matrixFileName[80];
    char rhsFileName[80];
    char timerName[80];

    sprintf(matrixFileName,"%s%d.mm", matrixPrefix.c_str(), i);

    // Load the matrix
    RCP<Matrix> Aprecond;
    cout << "[" << i << "] Loading matrix \"" << matrixFileName << "\"... ";
    try {
      Aprecond = Utils::Read(string(matrixFileName), xpetraParameters.GetLib(), comm, binary);
      cout << "done" << std::endl;
    } catch (...) {
      cout << "failed" << std::endl;
      return 1;
    }
    Aprecond->SetFixedBlockSize(numPDEs);

    // Build the nullspace
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(Aprecond->getRowMap(), numPDEs);
    nullspace->putScalar(zero);

    Teuchos::ArrayRCP<SC> data0, data1;
    data0 = nullspace->getDataNonConst(0);
    data1 = nullspace->getDataNonConst(1);
    for (size_t k = 0; k < Aprecond->getRowMap()->getNodeNumElements(); k += 2)
      data0[k+0] = data1[k+1] = one;

    // Build the preconditioner
    cout << "[" << i << "] Building preconditioner \"" << matrixFileName << "\" ..." << std::endl;
    sprintf(timerName, "Reuse: Preconditioner Setup i=%d", i);
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(timerName)));

    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->SetDefaultVerbLevel(MueLu::Extreme);
    H->GetLevel(0)->Set("A",         Aprecond);
    H->GetLevel(0)->Set("Nullspace", nullspace);
    H->IsPreconditioner(true);

    sprintf(timerName, "Reuse: Setup i=%d j=%d", i, i);
    timer = TimeMonitor::getNewTimer(timerName);
    timer->start();
    mueLuFactory.SetupHierarchy(*H);
    setup_times[i-first_matrix][i-first_matrix] = timer->stop();
    timer = Teuchos::null;

    Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO, LMO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator
    tm = Teuchos::null;

    // Loop over all future matrices
    int j = i;
    for (; j <= last_matrix; j++) {
      cout << "------------------------------------------------------------------" << std::endl;

      sprintf(matrixFileName, "%s%d.mm", matrixPrefix.c_str(), j);
      sprintf(rhsFileName,    "%s%d.mm", rhsPrefix.c_str(),    j);

      RCP<Matrix> Amatvec;
      if (j != i) {
        // Load the matrix
        cout << "[" << j << "]<-[" << i << "] Loading matrix \"" << matrixFileName << "\"... ";
        try {
          Amatvec = Utils::Read(string(matrixFileName), xpetraParameters.GetLib(), comm, binary);
          cout << "done" << std::endl;
        } catch (...) {
          cout << "failed" << std::endl;
          return 1;
        }

      } else {
        Amatvec = Aprecond;
      }
      Amatvec->SetFixedBlockSize(numPDEs);

      // Preconditioner update
      sprintf(timerName, "Reuse: Preconditioner Update i=%d", i);
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(timerName)));
      // No-op at present

      sprintf(timerName, "Reuse: Setup i=%d j=%d", i, j);
      timer = TimeMonitor::getNewTimer(timerName);
      timer->start();
      if (do_reuse == 0 && j != i) {
        // No reuse: Do a full recompute
        cout << "[" << j << "]<-[" << i << "] Building preconditioner \"" << matrixFileName << "\" ..." << std::endl;

        H = mueLuFactory.CreateHierarchy();

        H->SetDefaultVerbLevel(MueLu::Extreme);
        H->GetLevel(0)->Set("A",         Amatvec);
        H->GetLevel(0)->Set("Nullspace", nullspace);
        H->IsPreconditioner(true);

        mueLuFactory.SetupHierarchy(*H);

      } else if (do_reuse == 2 && j != i) {
        // "Fast" reuse
        // NTS: This isn't quite a real recompute yet.
        cout << "[" << j << "]<-[" << i << "] Updating preconditioner \"" << matrixFileName << "\" ..." << std::endl;

        H->GetLevel(0)->Set("A", Amatvec);
        mueLuFactory.SetupHierarchy(*H);
      }
      setup_times[i-first_matrix][j-first_matrix] = timer->stop();
      timer = Teuchos::null;

      tm = Teuchos::null;

      // Load the RHS
      cout << "[" << j << "] Loading rhs " << rhsFileName << std::endl;
      RCP<MultiVector> rhs = Utils::Read(string(rhsFileName), Amatvec->getRowMap());

      // Create an LHS
      RCP<Vector> X = VectorFactory::Build(Amatvec->getRowMap());
      X->putScalar(0.0);

      // Define Operator and Preconditioner
      Teuchos::RCP<OP> belosOp = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(Amatvec)); // Turns a Xpetra::Matrix object into a Belos operator

      // Construct a Belos LinearProblem object
      RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, rhs));
      belosProblem->setLeftPrec(belosPrec);
      belosProblem->setProblem();

      // Belos parameter list
      int    maxIts = 100;
      double tol    = 1e-12;
      Teuchos::ParameterList belosList;
      belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
      belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
      // belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
      belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList.set("Output Frequency",      1);
      belosList.set("Output Style",          Belos::Brief);

      // Create an iterative solver manager
      RCP< Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

      // Perform solve
      sprintf(timerName, "Reuse: Solve i=%d j=%d", i, j);
      timer = TimeMonitor::getNewTimer(timerName);
      timer->start();
      Belos::ReturnType ret = Belos::Unconverged;
      try {
        ret = solver->solve();

      } catch(...) {
        cout << std::endl << "ERROR:  Belos threw an error! " << std::endl;
      }
      double my_time = timer->stop();
      timer = Teuchos::null;

      // Get the number of iterations for this solve.

      // Check convergence
      if (ret != Belos::Converged) {
        cout << "ERROR  :  Belos did not converge! " << std::endl;
        break;
      } else {
        cout << "SUCCESS:  Belos converged!" << std::endl;
        cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
        iteration_counts[i-first_matrix][j-first_matrix] = solver->getNumIters();
        iteration_times [i-first_matrix][j-first_matrix] = my_time;
      }

    } //end j
    for (; j <= last_matrix; j++)
      setup_times[i-first_matrix][j-first_matrix] = Teuchos::ScalarTraits<SC>::nan();
  } // end i

  globalTimeMonitor = Teuchos::null;

  if (printTimings)
    // TimeMonitor::summarize(comm.ptr(), cout, false, true, false, Teuchos::Union);
    TimeMonitor::summarize(comm.ptr(), cout);

  if (!mypid) {
    printf("************************* Iteration Counts ***********************\n");
    for (int i = 0; i < ArraySize; i++) {
      for (int j = 0; j < ArraySize; j++)
        printf("%3d ", iteration_counts[i][j]);
      printf(";\n");
    }

    printf("************************* Iteration Times ***********************\n");
    for (int i = 0; i < ArraySize; i++) {
      for (int j = 0; j < ArraySize; j++)
        printf("%10.2f ", iteration_times[i][j]);
      printf(";\n");
    }

    printf("************************* Setup Times ***********************\n");
    for (int i = 0; i < ArraySize; i++) {
      for (int j = 0; j < ArraySize; j++)
        printf("%10.2f ", setup_times[i][j]);
      printf(";\n");
    }
  }

  //MueLu::MutuallyExclusiveTime<MueLu::BaseClass>::PrintParentChildPairs();

  return EXIT_SUCCESS;
} //main

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
#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_IO.hpp>

// Teuchos
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
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp

#include <MueLu_Utilities.hpp>

#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>  // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>   // => This header defines Belos::MueLuOp
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream &fancyout  = *fancy;
    fancyout.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 100, ny = 100, nz = 100;
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

    std::string xmlFileName = "scalingTest.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'scalingTest.xml'");
    int amgAsPrecond = 1;
    clp.setOption("precond", &amgAsPrecond, "apply multigrid as preconditioner");
    int amgAsSolver = 0;
    clp.setOption("fixPoint", &amgAsSolver, "apply multigrid as solver");
    bool printTimings = true;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    int writeMatricesOPT = -2;
    clp.setOption("write", &writeMatricesOPT, "write matrices to file (-1 means all; i>=0 means level i)");
    double tol = 1e-12;
    clp.setOption("tol", &tol, "solver convergence tolerance");
    std::string krylovMethod = "cg";
    clp.setOption("krylov", &krylovMethod, "outer Krylov method");
    int maxIts = 100;
    clp.setOption("maxits", &maxIts, "maximum number of Krylov iterations");
    int output = 1;
    clp.setOption("output", &output, "how often to print Krylov residual history");
    std::string matrixFileName = "A.mm";
    clp.setOption("matrixfile", &matrixFileName, "matrix market file containing matrix");
    std::string rhsFileName = "";
    clp.setOption("rhsfile", &rhsFileName, "matrix market file containing right-hand side");
    int nPDE = 1;
    clp.setOption("numpdes", &nPDE, "number of PDE equations");
    std::string convType = "r0";
    clp.setOption("convtype", &convType, "convergence type (r0 or none)");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    fancyout << "========================================================\n"
             << xpetraParameters << matrixParameters;

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time"))), tm;

    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

    RCP<Matrix> A              = Xpetra::IO<SC, LO, GO, Node>::Read(std::string(matrixFileName), lib, comm);
    RCP<const Map> map         = A->getRowMap();
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getDomainMap(), nPDE);
    // RCP<MultiVector> fakeCoordinates = MultiVectorFactory::Build(A->getDomainMap(),1);
    A->SetFixedBlockSize(nPDE);
    std::cout << "#pdes = " << A->GetFixedBlockSize() << std::endl;
    if (nPDE == 1)
      nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());
    else {
      for (int i = 0; i < nPDE; ++i) {
        Teuchos::ArrayRCP<SC> nsData = nullspace->getDataNonConst(i);
        for (int j = 0; j < nsData.size(); ++j) {
          GO gel = A->getDomainMap()->getGlobalElement(j) - A->getDomainMap()->getIndexBase();
          if ((gel - i) % nPDE == 0)
            nsData[j] = Teuchos::ScalarTraits<SC>::one();
        }
      }
    }

    comm->barrier();
    tm = Teuchos::null;

    fancyout << "Galeri complete.\n========================================================" << std::endl;

    // =========================================================================
    // Preconditioner construction
    // =========================================================================
    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));
    ParameterListInterpreter mueLuFactory(xmlFileName, *comm);

    comm->barrier();
    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();

    // By default, we use Extreme. However, typically the xml file contains verbosity parameter
    // which is used instead
    H->SetDefaultVerbLevel(MueLu::Extreme);

    H->GetLevel(0)->Set("A", A);
    H->GetLevel(0)->Set("Nullspace", nullspace);
    // H->GetLevel(0)->Set("Coordinates", fakeCoordinates);

    mueLuFactory.SetupHierarchy(*H);

    comm->barrier();
    tm = Teuchos::null;

    // Print out the hierarchy stats. We should not need this line, but for some reason the
    // print out in the hierarchy construction does not work.
    H->print(fancyout);

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================
    comm->barrier();
    tm.reset();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));

    RCP<Vector> X      = VectorFactory::Build(map, 1);
    RCP<MultiVector> B = VectorFactory::Build(map, 1);

    if (rhsFileName != "")
      B = Xpetra::IO<SC, LO, GO, Node>::ReadMultiVector(std::string(rhsFileName), A->getRowMap());
    else {
      // we set seed for reproducibility
      X->setSeed(846930886);
      bool useSameRandomGen = false;
      X->randomize(useSameRandomGen);
      A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

      Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(1);
      B->norm2(norms);
      // B->scale(1.0/norms[0]);
    }
    X->putScalar(zero);
    tm = Teuchos::null;

    if (writeMatricesOPT > -2)
      H->Write(writeMatricesOPT, writeMatricesOPT);

    comm->barrier();
    if (amgAsSolver) {
      tm.reset();
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Fixed Point Solve")));

      H->IsPreconditioner(false);
      Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(1);
      norms = Utilities::ResidualNorm(*A, *X, *B);
      std::cout << "                iter:    0           residual = " << norms[0] << std::endl;
      for (int i = 0; i < maxIts; ++i) {
        H->Iterate(*B, *X);
        norms = Utilities::ResidualNorm(*A, *X, *B);
        std::cout << "                iter:    " << i + 1 << "           residual = " << norms[0] << std::endl;
      }

    } else if (amgAsPrecond) {
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
        fancyout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      // Belos parameter list
      Teuchos::ParameterList belosList;
      belosList.set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
      belosList.set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
      belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList.set("Output Frequency", output);
      belosList.set("Output Style", Belos::Brief);
      // belosList.set("Orthogonalization",     "ICGS");
      if (convType == "none") {
        belosList.set("Explicit Residual Scaling", "None");
        belosList.set("Implicit Residual Scaling", "None");
      }

      // Create an iterative solver manager
      RCP<Belos::SolverManager<SC, MV, OP> > solver;
      if (krylovMethod == "cg") {
        solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
      } else if (krylovMethod == "gmres") {
        solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Invalid Krylov method.  Options are \"cg\" or \" gmres\".");
      }

      // Perform solve
      Belos::ReturnType ret = Belos::Unconverged;
      try {
        ret = solver->solve();

        // Get the number of iterations for this solve.
        fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

      } catch (...) {
        fancyout << std::endl
                 << "ERROR:  Belos threw an error! " << std::endl;
      }

      // Check convergence
      if (ret != Belos::Converged)
        fancyout << std::endl
                 << "ERROR:  Belos did not converge! " << std::endl;
      else
        fancyout << std::endl
                 << "SUCCESS:  Belos converged!" << std::endl;
#endif  // ifdef HAVE_MUELU_BELOS
    }
    comm->barrier();
    tm                = Teuchos::null;
    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
      MueLu::MutuallyExclusiveTime<MueLu::BaseClass>::PrintParentChildPairs();
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}

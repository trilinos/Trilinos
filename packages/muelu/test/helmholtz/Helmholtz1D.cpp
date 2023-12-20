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
#include <complex>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra and Galeri
#include <Xpetra_MultiVectorFactory.hpp>
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory_Helmholtz.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

// MueLu
#include "MueLu.hpp"
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_UseDefaultTypesComplex.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

// Belos
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  typedef Tpetra::Vector<SC, LO, GO, NO> TVEC;
  typedef Tpetra::MultiVector<SC, LO, GO, NO> TMV;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> TCRS;
  typedef Xpetra::CrsMatrix<SC, LO, GO, NO> XCRS;
  typedef Xpetra::TpetraCrsMatrix<SC, LO, GO, NO> XTCRS;
  typedef Xpetra::Matrix<SC, LO, GO, NO> XMAT;
  typedef Xpetra::CrsMatrixWrap<SC, LO, GO, NO> XWRAP;
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  typedef Belos::OperatorT<TMV> TOP;
  typedef Belos::OperatorTraits<SC, TMV, TOP> TOPT;
  typedef Belos::MultiVecTraits<SC, TMV> TMVT;
  typedef Belos::LinearProblem<SC, TMV, TOP> TProblem;
  typedef Belos::SolverManager<SC, TMV, TOP> TBelosSolver;
  typedef Belos::BlockGmresSolMgr<SC, TMV, TOP> TBelosGMRES;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::CommandLineProcessor clp(false);

    GO nx, ny, nz;
    nx = 100;
    ny = 100;
    nz = 100;
    double stretchx, stretchy, stretchz, h, delta;
    stretchx = 1.0;
    stretchy = 1.0;
    stretchz = 1.0;
    h        = 0.01;
    delta    = 2.0;
    int PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR;
    PMLXL = 10;
    PMLXR = 10;
    PMLYL = 10;
    PMLYR = 10;
    PMLZL = 10;
    PMLZR = 10;
    double omega, shift;
    omega = 20.0 * M_PI;
    shift = 0.5;

    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Helmholtz1D", 0, stretchx, stretchy, stretchz,
                                                    h, delta, PMLXL, PMLXR, PMLYL, PMLYR, PMLZL, PMLZR, omega, shift);
    Xpetra::Parameters xpetraParameters(clp);

    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
    RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

    Teuchos::ParameterList pl = matrixParameters.GetParameterList();
    RCP<RealValuedMultiVector> coordinates;
    Teuchos::ParameterList galeriList;
    galeriList.set("nx", pl.get("nx", nx));
    galeriList.set("ny", pl.get("ny", ny));
    galeriList.set("nz", pl.get("nz", nz));
    RCP<const Map> map;

    map         = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", map, matrixParameters.GetParameterList());

    RCP<const Tpetra::Map<LO, GO, NO> > tmap = Xpetra::toTpetra(map);

    Teuchos::ParameterList matrixParams = matrixParameters.GetParameterList();

    // Build problem
    RCP<Galeri::Xpetra::Problem_Helmholtz<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem_Helmholtz<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParams);
    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
    nullspace->putScalar((SC)1.0);

    comm->barrier();

    tm = Teuchos::null;

    // Construct a multigrid preconditioner
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

    // Multigrid Hierarchy
    RCP<Hierarchy> H = rcp(new Hierarchy(A));
    H->GetLevel(0)->Set("Nullspace", nullspace);
    FactoryManager Manager;
    H->Setup(Manager, 0, 5);
    // H->Write(-1,-1);

    tm = Teuchos::null;

    // Solve Ax = b
    tm          = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - LHS and RHS initialization")));
    RCP<TVEC> X = Tpetra::createVector<SC, LO, GO, NO>(tmap);
    RCP<TVEC> B = Tpetra::createVector<SC, LO, GO, NO>(tmap);
    X->putScalar((SC)0.0);
    B->putScalar((SC)0.0);
    if (comm->getRank() == 0) {
      B->replaceGlobalValue(0, (SC)1.0);
    }

    tm = Teuchos::null;

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - Belos Solve")));

    // Define Operator and Preconditioner
    RCP<TOP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));  // Turns a Xpetra::Matrix object into a Belos operator
    RCP<TOP> belosPrec = rcp(new Belos::MueLuOp<SC, LO, GO, NO>(H));   // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP<TProblem> belosProblem = rcp(new TProblem(belosOp, X, B));
    belosProblem->setRightPrec(belosPrec);
    bool set = belosProblem->setProblem();
    if (set == false) {
      if (comm->getRank() == 0)
        std::cout << std::endl
                  << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    int maxIts = 100;
    double tol = 1e-6;
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
    belosList.set("Flexible Gmres", false);       // set flexible GMRES on/off
    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency", 1);
    belosList.set("Output Style", Belos::Brief);

    // Create solver manager
    RCP<TBelosSolver> solver = rcp(new TBelosGMRES(belosProblem, rcp(&belosList, false)));

    // Perform solve
    Belos::ReturnType ret = Belos::Unconverged;
    ret                   = solver->solve();
    if (comm->getRank() == 0)
      std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

    // Check convergence
    if (ret != Belos::Converged) {
      if (comm->getRank() == 0) std::cout << std::endl
                                          << "ERROR:  Belos did not converge! " << std::endl;
    } else {
      if (comm->getRank() == 0) std::cout << std::endl
                                          << "SUCCESS:  Belos converged!" << std::endl;
    }

    // Get the number of iterations for this solve.
    if (comm->getRank() == 0)
      std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

    tm = Teuchos::null;

    globalTimeMonitor = Teuchos::null;

    TimeMonitor::summarize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main

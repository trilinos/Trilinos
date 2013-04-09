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
#include <complex>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// MueLu main header: include most common header files in one line
#include "MueLu.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_SmootherFactory.hpp"

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp

// Define default template types
typedef std::complex<double>                         SC;
typedef int                                          LO;
typedef int                                          GO;
typedef Kokkos::DefaultNode::DefaultNodeType         NO;
typedef Kokkos::DefaultKernels<SC,LO,NO>::SparseOps  LMO;

typedef Tpetra::Vector<SC,LO,GO,NO>                  TVEC;
typedef Tpetra::MultiVector<SC,LO,GO,NO>             MV;
typedef Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>           TCRS;
typedef Xpetra::CrsMatrix<SC,LO,GO,NO,LMO>           XCRS;
typedef Xpetra::TpetraCrsMatrix<SC,LO,GO,NO,LMO>     XTCRS;
typedef Xpetra::Matrix<SC,LO,GO,NO,LMO>              XMAT;
typedef Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO>       XWRAP;

typedef MueLu::Level                                 Level;
typedef MueLu::Hierarchy<SC,LO,GO,NO,LMO>            Hierarchy;
typedef MueLu::FactoryManager<SC,LO,GO>              FactoryManager;
typedef MueLu::TentativePFactory<SC,LO,GO,NO,LMO>    TPFactory;
typedef MueLu::SaPFactory<SC,LO,GO,NO,LMO>           SaPFactory;
typedef MueLu::GenericRFactory<SC,LO,GO,NO,LMO>      RFactory;
typedef MueLu::RAPFactory<SC,LO,GO,NO,LMO>           RAPFactory;
typedef MueLu::SmootherPrototype<SC,LO,GO,NO,LMO>    SmootherPrototype;
typedef MueLu::Ifpack2Smoother<SC,LO,GO,NO,LMO>      Ifpack2Smoother;
typedef MueLu::SmootherFactory<SC,LO,GO,NO,LMO>      SmootherFactory;
typedef MueLu::DirectSolver<SC,LO,GO,NO,LMO>         DirectSolver;

typedef Belos::OperatorT<MV>                         OP;
typedef Belos::OperatorTraits<SC,MV,OP>              OPT;
typedef Belos::MultiVecTraits<SC,MV>                 MVT;
typedef Belos::LinearProblem<SC,MV,OP>               Problem;
typedef Belos::SolverManager<SC,MV,OP>               BelosSolver;
typedef Belos::BlockGmresSolMgr<SC,MV,OP>            BelosGMRES;

int main(int argc, char *argv[]) {

  // RCPs
  using Teuchos::RCP;
  using Teuchos::rcp;

  // MPI initialization using Teuchos
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Parameters and constants
  GO numGlobalElements = 1000;
  SC h(1.0/((double) (numGlobalElements-1)),0.0);
  SC kh(0.625,0.0);
  SC kh2=kh*kh;
  SC k=kh/h;
  double alpha=1.0;
  double beta=0.0;
  SC shift(alpha,beta);
  SC one(1.0,0.0);
  SC two(2.0,0.0);
  SC complexone(0.0,1.0);

  // Construct a Map that puts approximately the same number of equations on each processor
  RCP<const Tpetra::Map<LO, GO, NO> > map = Tpetra::createUniformContigMap<LO, GO>(numGlobalElements, comm);

  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GO> myGlobalElements = map->getNodeElementList();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<TCRS> L = rcp(new TCRS(map, 3));
  RCP<TCRS> S = rcp(new TCRS(map, 3));
  RCP<TCRS> A = rcp(new TCRS(map, 3));

  // Laplacian
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      L->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<GO>(myGlobalElements[i], myGlobalElements[i] +1),
                            Teuchos::tuple<SC> (two, -one));
    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      L->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<GO>(myGlobalElements[i] -1, myGlobalElements[i]),
                            Teuchos::tuple<SC> (-one, two));
    }
    else {
      L->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<GO>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1),
                            Teuchos::tuple<SC> (-one, two, -one));
    }
  }
  L->fillComplete();

  // Shifted Helmholtz Operator
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      S->insertGlobalValues(myGlobalElements[i],
			    Teuchos::tuple<GO> (myGlobalElements[i]),
			    Teuchos::tuple<SC> (one));
    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      S->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<GO>(myGlobalElements[i] -1, myGlobalElements[i]),
                            Teuchos::tuple<SC> (-one, one-complexone*kh));
    }
    else {
      S->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<GO>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1),
                            Teuchos::tuple<SC> (-one, two-shift*kh2, -one));
    }
  }
  S->fillComplete();

  // Original Helmholtz Operator
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues(myGlobalElements[i],
			    Teuchos::tuple<GO> (myGlobalElements[i]),
			    Teuchos::tuple<SC> (one));

    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<GO>(myGlobalElements[i] -1, myGlobalElements[i]),
                            Teuchos::tuple<SC> (-one,one-complexone*kh) );
    }
    else {
      A->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<GO>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1),
                            Teuchos::tuple<SC> (-one, two-kh2, -one));
    }
  }
  A->fillComplete();

  // Turn Tpetra::CrsMatrix into MueLu::Matrix
  RCP<XCRS> mueluL_ = rcp(new XTCRS(L));
  RCP<XMAT> mueluL  = rcp(new XWRAP(mueluL_));
  RCP<XCRS> mueluS_ = rcp(new XTCRS(S));
  RCP<XMAT> mueluS  = rcp(new XWRAP(mueluS_));
  RCP<XCRS> mueluA_ = rcp(new XTCRS(A));
  RCP<XMAT> mueluA  = rcp(new XWRAP(mueluA_));

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(mueluL));
  FactoryManager M;

  // Prolongation/Restriction
  RCP<TPFactory>  TentPFact = rcp( new TPFactory  );
  RCP<SaPFactory> Pfact     = rcp( new SaPFactory );
  RCP<RFactory>   Rfact     = rcp( new RFactory   );
  RCP<RAPFactory> Acfact    = rcp( new RAPFactory );

  // Smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpack2Type;
  Teuchos::ParameterList ifpack2List;
  /*ifpack2Type = "KRYLOV";
  ifpack2List.set("krylov: iteration type",1);
  ifpack2List.set("krylov: number of iterations",4);
  ifpack2List.set("krylov: residual tolerance",1e-6);
  ifpack2List.set("krylov: block size",1);
  ifpack2List.set("krylov: zero starting solution",true);
  ifpack2List.set("krylov: preconditioner type",3);*/
  // Additive Schwarz smoother
  ifpack2Type = "SCHWARZ";
  ifpack2List.set("schwarz: compute condest", false);
  ifpack2List.set("schwarz: overlap level", 0);
  // ILUT smoother
  //ifpack2Type = "ILUT";
  //ifpack2List.set("fact: ilut level-of-fill", (double)1.0);
  //ifpack2List.set("fact: absolute threshold", (double)0.0);
  //ifpack2List.set("fact: relative threshold", (double)1.0);
  //ifpack2List.set("fact: relax value", (double)0.0);
  //smooProto = Teuchos::rcp( new MueLu::Ifpack2Smoother<SC, LO, GO, NO, LMO>("ILUT",ifpack2List) );
  // Gauss-Seidel smoother
  //ifpack2Type = "RELAXATION";
  //ifpack2List.set("relaxation: sweeps", (LO) 4);
  //ifpack2List.set("relaxation: damping factor", 1.0); // 0.7
  //ifpack2List.set("relaxation: type", "Gauss-Seidel");

  smooProto = Teuchos::rcp( new Ifpack2Smoother(ifpack2Type,ifpack2List) );
  RCP<SmootherFactory> SmooFact;
  LO maxLevels = 6;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  // create coarsest smoother
  RCP<SmootherPrototype> coarsestSmooProto;
  // Direct Solver
  std::string type = "";
  Teuchos::ParameterList coarsestSmooList;
#if defined(HAVE_AMESOS_SUPERLU)
  coarsestSmooProto = rcp( new DirectSolver("Superlu",coarsestSmooList) );
#else
  coarsestSmooProto = rcp( new DirectSolver("Klu",coarsestSmooList) );
#endif
  RCP<SmootherFactory> coarsestSmooFact = rcp(new SmootherFactory(coarsestSmooProto, Teuchos::null));

  // Setup R's and P's
  M.SetFactory("P", Pfact);
  M.SetFactory("R", Rfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Ptent", TentPFact);
  H->Keep("P", Pfact.get());
  H->Keep("R", Rfact.get());
  H->Keep("Ptent", TentPFact.get());
  H->Setup(M, 0, maxLevels);
  H->print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);

  // Setup coarse grid operators and smoothers
  RCP<Level> finestLevel = H->GetLevel();
  finestLevel->Set("A", mueluS);
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarsestSmooFact);
  H->Setup(M, 0, H->GetNumLevels());

  // right hand side and left hand side vectors
  RCP<TVEC> X = Tpetra::createVector<SC,LO,GO,NO>(map);
  RCP<TVEC> B = Tpetra::createVector<SC,LO,GO,NO>(map);
  X->putScalar((SC) 0.0);
  B->putScalar((SC) 0.0);
  if(comm->getRank()==0) {
    B->replaceGlobalValue(0, 1.0);
  }

  // Define Operator and Preconditioner
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC,LO,GO,NO,LMO>(mueluA));   // Turns a Xpetra::Matrix object into a Belos operator
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<SC,LO,GO,NO,LMO>(H));         // Turns a MueLu::Hierarchy object into a Belos operator

  // Construct a Belos LinearProblem object
  RCP<Problem> belosProblem = rcp(new Problem(belosOp,X,B));
  belosProblem->setRightPrec(belosPrec);
  bool set = belosProblem->setProblem();
  if (set == false) {
    if(comm->getRank()==0)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;

    return EXIT_FAILURE;
  }

  // Belos parameter list
  int maxIts = 100;
  double tol = 1e-6;
  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
  belosList.set("Flexible Gmres", false);          // set flexible GMRES on

  // Create a FGMRES solver manager
  RCP<BelosSolver> solver = rcp( new BelosGMRES(belosProblem, rcp(&belosList, false)) );

  // Perform solve
  Belos::ReturnType ret = solver->solve();

  // Get the number of iterations for this solve.
  if(comm->getRank()==0)
    std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

  // Compute actual residuals.
  int numrhs=1;
  bool badRes = false;
  std::vector<double> actual_resids(numrhs);
  std::vector<double> rhs_norm(numrhs);
  RCP<MV> resid = Tpetra::createMultiVector<SC,LO,GO,NO>(map, numrhs);
  OPT::Apply(*belosOp, *X, *resid);
  MVT::MvAddMv(-1.0, *resid, 1.0, *B, *resid);
  MVT::MvNorm(*resid, actual_resids);
  MVT::MvNorm(*B, rhs_norm);
  if(comm->getRank()==0) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  }
  for (int i = 0; i < numrhs; i++) {
    double actRes = abs(actual_resids[i])/rhs_norm[i];
    if(comm->getRank()==0) {
      std::cout <<"Problem " << i << " : \t" << actRes <<std::endl;
    }
    if (actRes > tol) { badRes = true; }
  }

  // Check convergence
  if (ret != Belos::Converged || badRes) {
    if(comm->getRank()==0) {
      std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    }
    return EXIT_FAILURE;
  }
  if(comm->getRank()==0) {
    std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
  }
  // print solution entries
  // using Teuchos::VERB_EXTREME;
  // Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream( Teuchos::rcpFromRef(std::cerr) );
  // X->describe(*out,VERB_EXTREME);

  return EXIT_SUCCESS;
}

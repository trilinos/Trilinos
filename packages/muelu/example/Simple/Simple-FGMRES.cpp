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
//#include "MueLu_PgPFactory.hpp"
//#include "MueLu_EminPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_SmootherFactory.hpp"

// Header files defining default types for template parameters.
// These headers must be included after other MueLu/Xpetra headers.
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=int, GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>    // => typedef MueLu::FooClass<Scalar, LocalOrdinal, ...> Foo

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp

// Define default template types
typedef std::complex<double> ScalarC;
//typedef int    LocalOrdinal;
//typedef int    GlobalOrdinal;

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;

  // MPI initialization using Teuchos

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Parameters

  // problem size
  GlobalOrdinal numGlobalElements = 10000;
  ScalarC h(1.0/((double) (numGlobalElements-1)),0.0);
  ScalarC kh(0.625,0.0);
  ScalarC kh2=kh*kh;
  ScalarC k=kh/h;
  double alpha=1.0;
  double beta=0.0;
  ScalarC shift(alpha,beta);
  ScalarC one(1.0,0.0);
  ScalarC two(2.0,0.0);
  ScalarC complexone(0.0,1.0);

  // Construct the problem

  // Construct a Map that puts approximately the same number of equations on each processor
  RCP<const Tpetra::Map<LO, GO, NO> > map = Tpetra::createUniformContigMap<LO, GO>(numGlobalElements, comm);

  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getNodeElementList();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<Tpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO> > L = rcp(new Tpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO>(map, 3));
  RCP<Tpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO> > S = rcp(new Tpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO>(map, 3));
  RCP<Tpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO> > A = rcp(new Tpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO>(map, 3));

  // Laplacian
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      L->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<ScalarC> (two, -one));
    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      L->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]), 
                            Teuchos::tuple<ScalarC> (-one, two));
    }
    else {
      L->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<ScalarC> (-one, two, -one));
    }
  }
  L->fillComplete();

  // Shifted Helmholtz Operator
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      S->insertGlobalValues(myGlobalElements[i],
			    Teuchos::tuple<GlobalOrdinal> (myGlobalElements[i]),
			    Teuchos::tuple<ScalarC> (one));
    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      S->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]), 
                            Teuchos::tuple<ScalarC> (-one, one-complexone*kh));
    }
    else {
      S->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<ScalarC> (-one, two-shift*kh2, -one));
    }
  }
  S->fillComplete();

  // Original Helmholtz Operator
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues(myGlobalElements[i],
			    Teuchos::tuple<GlobalOrdinal> (myGlobalElements[i]),
			    Teuchos::tuple<ScalarC> (one));

    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]), 
                            Teuchos::tuple<ScalarC> (-one,one-complexone*kh) );
    }
    else {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<ScalarC> (-one, two-kh2, -one));
    }
  }
  A->fillComplete();

  // Construct a multigrid preconditioner

  // Turn Tpetra::CrsMatrix into MueLu::Matrix
  RCP<Xpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO> > mueluL_ = rcp(new Xpetra::TpetraCrsMatrix<ScalarC, LO, GO, NO, LMO>(L)); 
  RCP<Xpetra::Matrix <ScalarC, LO, GO, NO, LMO> > mueluL  = rcp(new Xpetra::CrsMatrixWrap<ScalarC, LO, GO, NO, LMO>(mueluL_));
  RCP<Xpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO> > mueluS_ = rcp(new Xpetra::TpetraCrsMatrix<ScalarC, LO, GO, NO, LMO>(S)); 
  RCP<Xpetra::Matrix <ScalarC, LO, GO, NO, LMO> > mueluS  = rcp(new Xpetra::CrsMatrixWrap<ScalarC, LO, GO, NO, LMO>(mueluS_));
  RCP<Xpetra::CrsMatrix<ScalarC, LO, GO, NO, LMO> > mueluA_ = rcp(new Xpetra::TpetraCrsMatrix<ScalarC, LO, GO, NO, LMO>(A)); 
  RCP<Xpetra::Matrix <ScalarC, LO, GO, NO, LMO> > mueluA  = rcp(new Xpetra::CrsMatrixWrap<ScalarC, LO, GO, NO, LMO>(mueluA_));

  // Multigrid Hierarchy
  RCP< MueLu::Hierarchy<ScalarC, LO, GO, NO, LMO> > H = rcp(new MueLu::Hierarchy<ScalarC, LO, GO, NO, LMO>(mueluL));
  //H->setVerbLevel(Teuchos::VERB_HIGH);
  MueLu::FactoryManager<ScalarC, LocalOrdinal, GlobalOrdinal> M;

  // Prolongation/Restriction
  RCP< MueLu::TentativePFactory<ScalarC, LO, GO, NO, LMO> > TentPFact = rcp( new MueLu::TentativePFactory<ScalarC, LO, GO, NO, LMO> );
  RCP< MueLu::SaPFactory<ScalarC, LO, GO, NO, LMO> >        Pfact     = rcp( new MueLu::SaPFactory<ScalarC, LO, GO, NO, LMO> );
  RCP< MueLu::GenericRFactory<ScalarC, LO, GO, NO, LMO> >   Rfact     = rcp( new MueLu::GenericRFactory<ScalarC, LO, GO, NO, LMO> );
  RCP< MueLu::RAPFactory<ScalarC, LO, GO, NO, LMO> >        Acfact    = rcp( new MueLu::RAPFactory<ScalarC, LO, GO, NO, LMO> );

  // Smoothers
  RCP< MueLu::SmootherPrototype<ScalarC, LO, GO, NO, LMO> > smooProto;
  std::string ifpack2Type;
  Teuchos::ParameterList ifpack2List;
  ifpack2Type = "KRYLOV";
  ifpack2List.set("krylov: number of iterations",4);
  ifpack2List.set("krylov: residual tolerance",1e-6);
  ifpack2List.set("krylov: block size",1);
  ifpack2List.set("krylov: zero starting solution",true);
  ifpack2List.set("krylov: preconditioner type",1);
  // ILUT smoother
  //ifpack2Type = "ILUT";
  //ifpack2List.set("fact: ilut level-of-fill", (double)1.0);
  //ifpack2List.set("fact: absolute threshold", (double)0.0);
  //ifpack2List.set("fact: relative threshold", (double)1.0);
  //ifpack2List.set("fact: relax value", (double)0.0);
  //smooProto = Teuchos::rcp( new MueLu::Ifpack2Smoother<ScalarC, LO, GO, NO, LMO>("ILUT",ifpack2List) );
  // Gauss-Seidel smoother
  //ifpack2Type = "RELAXATION";
  //ifpack2List.set("relaxation: sweeps", (LO) 1);
  //ifpack2List.set("relaxation: damping factor", (SC) 1.0); // 0.7
  //ifpack2List.set("relaxation: type", "Gauss-Seidel");
  smooProto = Teuchos::rcp( new MueLu::Ifpack2Smoother<ScalarC, LO, GO, NO, LMO>(ifpack2Type, ifpack2List) );
  RCP< MueLu::SmootherFactory<ScalarC, LO, GO, NO, LMO> > SmooFact;
  LO maxLevels = 6;
  if (maxLevels > 1)
    SmooFact = rcp( new MueLu::SmootherFactory<ScalarC, LO, GO, NO, LMO>(smooProto) );
  // create coarsest smoother
  RCP< MueLu::SmootherPrototype<ScalarC, LO, GO, NO, LMO> > coarsestSmooProto;
  // Direct Solver
  std::string type = "";
  Teuchos::ParameterList coarsestSmooList;
#if defined(HAVE_AMESOS_SUPERLU)
  coarsestSmooProto = Teuchos::rcp( new MueLu::DirectSolver<ScalarC, LO, GO, NO, LMO>("Superlu", coarsestSmooList) );
#else
  coarsestSmooProto = Teuchos::rcp( new MueLu::DirectSolver<ScalarC, LO, GO, NO, LMO>("Klu", coarsestSmooList) );
#endif
  RCP< MueLu::SmootherFactory<ScalarC, LO, GO, NO, LMO> > coarsestSmooFact = rcp(new MueLu::SmootherFactory<ScalarC, LO, GO, NO, LMO>(coarsestSmooProto, Teuchos::null));

  // Setup R's and P's
  M.SetFactory("P", Pfact);
  M.SetFactory("R", Rfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Ptent", TentPFact);
  H->Keep("P", Pfact.get());
  H->Keep("R", Rfact.get());
  H->Keep("Ptent", TentPFact.get());
  H->Setup(M, 0, maxLevels);
  Teuchos::ParameterList status = H->Summarize(MueLu::None);
  int numLevels = status.get("number of levels",-1);
  H->print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);
  // Setup coarse grid operators and smoothers
  RCP<Level> finestLevel = H->GetLevel();
  finestLevel->Set("A", mueluS);
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarsestSmooFact);
  H->Setup(M, 0, numLevels);

  RCP<Tpetra::Vector<ScalarC, LO, GO, NO> > X = Tpetra::createVector<ScalarC, LO, GO, NO>(map);
  RCP<Tpetra::Vector<ScalarC, LO, GO, NO> > B = Tpetra::createVector<ScalarC, LO, GO, NO>(map);  
  X->putScalar((ScalarC) 0.0);
  B->putScalar((ScalarC) 0.0);
  B->replaceGlobalValue(0, 1.0);

  // Matrix and Multivector type that will be used with Belos
  typedef Tpetra::MultiVector<ScalarC, LO, GO, NO> MV;
  typedef Belos::OperatorT<MV>                OP;

  // Define Operator and Preconditioner
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<ScalarC, LO, GO, NO, LMO>(mueluA));   // Turns a Xpetra::Matrix object into a Belos operator
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<ScalarC, LO, GO, NO, LMO>(H));         // Turns a MueLu::Hierarchy object into a Belos operator

  // Construct a Belos LinearProblem object
  RCP< Belos::LinearProblem<ScalarC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<ScalarC, MV, OP>(belosOp, X, B));
  belosProblem->setRightPrec(belosPrec);
    
  bool set = belosProblem->setProblem();
  if (set == false) {
    std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return EXIT_FAILURE;
  }
    
  // Belos parameter list
  int maxIts = 100;
  double tol = 1e-6;
  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
  belosList.set("Flexible Gmres", true);          // set flexible GMRES on

  // Create a FGMRES solver manager
  RCP< Belos::SolverManager<ScalarC, MV, OP> > solver = rcp(new Belos::BlockGmresSolMgr<ScalarC, MV, OP>(belosProblem, rcp(&belosList, false)));
    
  // Perform solve
  Belos::ReturnType ret = solver->solve();
  
  // Get the number of iterations for this solve.
  std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
  
  // Compute actual residuals.
  int numrhs=1;
  bool badRes = false;
  std::vector<SC> actual_resids(numrhs);
  std::vector<SC> rhs_norm(numrhs);
  RCP<Tpetra::MultiVector<ScalarC, LO, GO, NO> > resid = Tpetra::createMultiVector<ScalarC, LO, GO, NO>(map, numrhs); 

  typedef Belos::OperatorTraits<ScalarC, MV, OP> OPT;
  typedef Belos::MultiVecTraits<ScalarC, MV>     MVT;
    
  OPT::Apply(*belosOp, *X, *resid);
  MVT::MvAddMv(-1.0, *resid, 1.0, *B, *resid);
  MVT::MvNorm(*resid, actual_resids);
  MVT::MvNorm(*B, rhs_norm);
  std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  for (int i = 0; i < numrhs; i++) {
    SC actRes = abs(actual_resids[i])/rhs_norm[i];
    std::cout <<"Problem " << i << " : \t" << actRes <<std::endl;
    if (actRes > tol) { badRes = true; }
  }

  // Check convergence
  if (ret != Belos::Converged || badRes) {
    std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;

  // print solution entries
  // using Teuchos::VERB_EXTREME;
  // Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream( Teuchos::rcpFromRef(std::cerr) );
  // X->describe(*out,VERB_EXTREME);  

  return EXIT_SUCCESS;
}

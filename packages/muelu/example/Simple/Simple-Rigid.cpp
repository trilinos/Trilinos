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

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_MultiVectorFactory.hpp>

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
#include "MueLu_RigidBodyModeFactory.hpp"
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp

// Define default template types
typedef Tpetra::Vector<SC,LO,GO,NO>                  TVEC;
typedef Tpetra::MultiVector<SC,LO,GO,NO>             MV;
typedef Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>           TCRS;

typedef MueLu::Level                                 Level;
typedef MueLu::Hierarchy<SC,LO,GO,NO,LMO>            Hierarchy;
typedef MueLu::FactoryManager<SC,LO,GO>              FactoryManager;
typedef MueLu::TentativePFactory<SC,LO,GO,NO,LMO>    TPFactory;
typedef MueLu::SaPFactory<SC,LO,GO,NO,LMO>           SaPFactory;
typedef MueLu::GenericRFactory<SC,LO,GO,NO,LMO>      GRFactory;
typedef MueLu::RAPFactory<SC,LO,GO,NO,LMO>           RAPFactory;
typedef MueLu::RigidBodyModeFactory<SC,LO,GO,NO,LMO> RBMFactory;
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

  // Figuring out the parallel distribution
  int nProcs = comm->getSize();
  LO  nTotalDOFs   = 243;
  int nDOFsPerNode = 3;
  int nTotalNodes  = nTotalDOFs/nDOFsPerNode;
  int nLocalNodes  = nTotalNodes/nProcs;
  int leftover     = nTotalNodes - (nLocalNodes*(nProcs-1));
  if(comm->getRank() == nProcs-1) {
    nLocalNodes = leftover;
  }
  int nnzeros = 10085;

  // Construct a Map that puts approximately the same number of mesh nodes per processor
  RCP<const Tpetra::Map<LO, GO, NO> > map = Tpetra::createUniformContigMap<LO, GO>(nTotalNodes, comm);
  // Tpetra map into Xpetra map
  RCP<const Map> xmap = Xpetra::toXpetra(map);
  // Map takes constant number of DOFs per node
  xmap = MapFactory::Build(xmap,nDOFsPerNode);
  map = Xpetra::toTpetra(xmap);

  // Create a CrsMatrix using the map
  RCP<CrsMatrix> A = CrsMatrixFactory::Build(xmap,50);

  // Read sample matrix from .txt file and put values into A
  std::ifstream matfile;
  matfile.open("stiff.txt");
  for (int i = 0; i < nnzeros; i++) {
    int current_row, current_column;
    SC current_value;
    matfile >> current_row >> current_column >> current_value ;
    if(map->isNodeGlobalElement(current_row)==true) {
      A->insertGlobalValues(current_row,
			    Teuchos::tuple<GO> (current_column),
			    Teuchos::tuple<SC> (current_value));
    }
  }
  A->fillComplete();

  // Turn Xpetra::CrsMatrix into Xpetra::Matrix
  RCP<Matrix> mueluA  = rcp(new CrsMatrixWrap(A));

  // MultiVector of coordinates
  // NOTE: must be of size equal to the number of DOFs!
  RCP<MultiVector> coordinates;
  coordinates = MultiVectorFactory::Build(xmap, 3);
  SC h=0.5;
  for(int k=0; k<9; k++) {
    for(int j=0; j<3; j++) {
      for(int i=0; i<3; i++) {
	int curidx = i+3*j+9*k;
        int curidx0 = curidx*3+0;
        int curidx1 = curidx*3+1;
        int curidx2 = curidx*3+2;
	if(xmap->isNodeGlobalElement(curidx0)==true) {
	  coordinates->replaceGlobalValue(curidx0,0,i*h);
	  coordinates->replaceGlobalValue(curidx1,0,i*h);
	  coordinates->replaceGlobalValue(curidx2,0,i*h);
	  coordinates->replaceGlobalValue(curidx0,1,j*h);
	  coordinates->replaceGlobalValue(curidx1,1,j*h);
	  coordinates->replaceGlobalValue(curidx2,1,j*h);
	  coordinates->replaceGlobalValue(curidx0,2,k*h);
	  coordinates->replaceGlobalValue(curidx1,2,k*h);
	  coordinates->replaceGlobalValue(curidx2,2,k*h);
	}
      }
    }
  }

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(mueluA));
  FactoryManager M;

  // Prolongation/Restriction
  RCP<TPFactory>  TentPFact = rcp( new TPFactory     );
  RCP<SaPFactory> Pfact     = rcp( new SaPFactory    );
  RCP<GRFactory>  Rfact     = rcp( new GRFactory     );
  RCP<RAPFactory> Acfact    = rcp( new RAPFactory    );
  RCP<RBMFactory> RBMfact   = rcp( new RBMFactory(3) );

  // Smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpack2Type;
  Teuchos::ParameterList ifpack2List;
  /*ifpack2Type = "KRYLOV";
  ifpack2List.set("krylov: number of iterations",4);
  ifpack2List.set("krylov: residual tolerance",1e-6);
  ifpack2List.set("krylov: block size",1);
  ifpack2List.set("krylov: zero starting solution",true);
  ifpack2List.set("krylov: preconditioner type",1);*/
  // ILUT smoother
  //ifpack2Type = "ILUT";
  //ifpack2List.set("fact: ilut level-of-fill", (double)1.0);
  //ifpack2List.set("fact: absolute threshold", (double)0.0);
  //ifpack2List.set("fact: relative threshold", (double)1.0);
  //ifpack2List.set("fact: relax value", (double)0.0);
  // Gauss-Seidel smoother
  ifpack2Type = "RELAXATION";
  ifpack2List.set("relaxation: sweeps", (LO) 1);
  ifpack2List.set("relaxation: damping factor", (SC) 1.0); // 0.7
  ifpack2List.set("relaxation: type", "Gauss-Seidel");
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
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarsestSmooFact);
  M.SetFactory("Nullspace", RBMfact );
  H->GetLevel(0)->Set("Coordinates",coordinates);
  H->Setup(M, 0, maxLevels);

  // right hand side and left hand side vectors
  RCP<TVEC> X = Tpetra::createVector<SC,LO,GO,NO>(map);
  RCP<TVEC> B = Tpetra::createVector<SC,LO,GO,NO>(map);  
  X->putScalar((SC) 0.0);
  B->randomize();
  //B->replaceGlobalValue(nTotalDOFs/2, 1.0);

  // Define Operator and Preconditioner
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC,LO,GO,NO,LMO>(mueluA));   // Turns a Xpetra::Matrix object into a Belos operator
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<SC,LO,GO,NO,LMO>(H));         // Turns a MueLu::Hierarchy object into a Belos operator

  // Construct a Belos LinearProblem object
  RCP<Problem> belosProblem = rcp(new Problem(belosOp,X,B));
  belosProblem->setRightPrec(belosPrec);    
  bool set = belosProblem->setProblem();
  if (set == false) {
    if(comm->getRank()==0) {
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    }
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
  RCP<BelosSolver> solver = rcp( new BelosGMRES(belosProblem, rcp(&belosList, false)) );
    
  // Perform solve
  Belos::ReturnType ret = solver->solve();
  
  // Get the number of iterations for this solve.
  if(comm->getRank()==0) {
    std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
  }  
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

  return EXIT_SUCCESS;
}

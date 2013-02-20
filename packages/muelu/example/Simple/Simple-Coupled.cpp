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
#include "MueLu_CoupledRBMFactory.hpp"

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
//typedef Kokkos::SerialNode                           NO;
//typedef Kokkos::AltSparseOps<void,int,
//			     Kokkos::SerialNode,
//			     Kokkos::details::AltSparseOpsDefaultAllocator
//			     <int, Kokkos::SerialNode> >     LMO;

typedef Tpetra::Vector<SC,LO,GO,NO>                  TVEC;
typedef Tpetra::MultiVector<SC,LO,GO,NO>             TMV;
typedef Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>           TCRS;
typedef Xpetra::MultiVector<SC,LO,GO,NO>             XMV;
typedef Xpetra::CrsMatrix<SC,LO,GO,NO,LMO>           XCRS;
typedef Xpetra::TpetraCrsMatrix<SC,LO,GO,NO,LMO>     XTCRS; 
typedef Xpetra::Matrix<SC,LO,GO,NO,LMO>              XMAT;
typedef Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO>       XWRAP;
typedef Xpetra::Map<LO,GO,NO>                        Map;
typedef Xpetra::MapFactory<LO,GO,NO>                 MapFactory;
typedef Xpetra::CrsMatrixFactory<SC,LO,GO,NO,LMO>    CrsMatrixFactory;
typedef Xpetra::MultiVectorFactory<SC,LO,GO,NO>      MultiVectorFactory;


typedef MueLu::Level                                 Level;
typedef MueLu::Hierarchy<SC,LO,GO,NO,LMO>            Hierarchy;
typedef MueLu::FactoryManager<SC,LO,GO>              FactoryManager;
typedef MueLu::TentativePFactory<SC,LO,GO,NO,LMO>    TPFactory;
typedef MueLu::SaPFactory<SC,LO,GO,NO,LMO>           SaPFactory;
typedef MueLu::GenericRFactory<SC,LO,GO,NO,LMO>      GRFactory;
typedef MueLu::RAPFactory<SC,LO,GO,NO,LMO>           RAPFactory;
typedef MueLu::CoupledRBMFactory<SC,LO,GO,NO,LMO>    RBMFactory;
typedef MueLu::SmootherPrototype<SC,LO,GO,NO,LMO>    SmootherPrototype;
typedef MueLu::Ifpack2Smoother<SC,LO,GO,NO,LMO>      Ifpack2Smoother;
typedef MueLu::SmootherFactory<SC,LO,GO,NO,LMO>      SmootherFactory;
typedef MueLu::DirectSolver<SC,LO,GO,NO,LMO>         DirectSolver;

typedef Belos::OperatorT<TMV>                        OP;
typedef Belos::OperatorTraits<SC,TMV,OP>             OPT;
typedef Belos::MultiVecTraits<SC,TMV>                MVT;
typedef Belos::LinearProblem<SC,TMV,OP>              Problem;
typedef Belos::SolverManager<SC,TMV,OP>              BelosSolver;
typedef Belos::BlockGmresSolMgr<SC,TMV,OP>           BelosGMRES;


int main(int argc, char *argv[]) {

  // RCPs
  using Teuchos::RCP;
  using Teuchos::rcp;

  // MPI initialization using Teuchos
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // predetermined DOF and parameter info
  // for case 1
  int nAcousticDOFs = 126;
  int nPaddedDOFs   = 3*nAcousticDOFs;
  int nElasticDOFs  = 243;
  int nTotalDOFs    = nPaddedDOFs+nElasticDOFs;
  int nDOFsPerNode  = 3;
  int nTotalNodes   = nTotalDOFs/nDOFsPerNode;
  int nnzeros_stiff = 12869;
  int nnzeros_damp  = 98;
  int nnzeros_mass  = 5635;
  int maxDOFsPerRow = 50;
  double freq = 11.0;
  double h=0.5;
  int nx=3;
  int ny=3;
  int nz=23;
  // for case 2
  /*int nAcousticDOFs = 600;
  int nPaddedDOFs   = 3*nAcousticDOFs;
  int nElasticDOFs  = 1425;
  int nTotalDOFs    = nPaddedDOFs+nElasticDOFs;
  int nDOFsPerNode  = 3;
  int nTotalNodes   = nTotalDOFs/nDOFsPerNode;
  int nnzeros_stiff = 94557;
  int nnzeros_damp  = 338;
  int nnzeros_mass  = 39715;
  int maxDOFsPerRow = 100;
  double freq = 22.0;
  double h=0.25;
  int nx=5;
  int ny=5;
  int nz=43;*/
  // for case 3
  /*int nAcousticDOFs = 3564;
  int nPaddedDOFs   = 3*nAcousticDOFs;
  int nElasticDOFs  = 9477;
  int nTotalDOFs    = nPaddedDOFs+nElasticDOFs;
  int nDOFsPerNode  = 3;
  int nTotalNodes   = nTotalDOFs/nDOFsPerNode;
  int nnzeros_stiff = 719287;
  int nnzeros_damp  = 1250;
  int nnzeros_mass  = 296875;
  int maxDOFsPerRow = 100;
  double freq = 44.0;
  double h=0.125;
  int nx=9;
  int ny=9;
  int nz=83;*/

  // Parameters
  double alpha  = 1.0;
  double beta   = 0.25;
  double omega  = 2.0*M_PI*freq;
  double omega2 = omega*omega;
  SC one(1.0,0.0);
  SC shift(alpha,beta);

  // Construct a Map that puts approximately the same number of mesh nodes per processor
  RCP<const Tpetra::Map<LO, GO, NO> > map = Tpetra::createUniformContigMap<LO, GO>(nTotalNodes, comm);
  // Tpetra map into Xpetra map
  RCP<const Map> xmap = Xpetra::toXpetra(map);
  // Map takes constant number of DOFs per node
  xmap = MapFactory::Build(xmap,nDOFsPerNode);
  map = Xpetra::toTpetra(xmap);

  // Create a CrsMatrix using the map
  // K - only stiffness matrix
  // A - original Helmholtz operator
  // S - shifted Helmholtz operator
  RCP<TCRS> K = rcp(new TCRS(map,maxDOFsPerRow));
  RCP<TCRS> A = rcp(new TCRS(map,maxDOFsPerRow));
  RCP<TCRS> S = rcp(new TCRS(map,maxDOFsPerRow));

  // Read data from .txt files and put into appropriate matrices
  // Note: must pad acoustic DOFs with identity matrix
  // stiffness matrix
  std::ifstream matfile_stiff;
  matfile_stiff.open("coupled_stiff.txt");
  for (int i = 0; i < nnzeros_stiff; i++) {
    int current_row, current_column;
    double current_value;
    matfile_stiff >> current_row >> current_column >> current_value ;
    SC cpx_current_value(current_value,0.0);
    if(current_row < nAcousticDOFs)    { current_row    = current_row*3;                             }
    else                               { current_row    = current_row-nAcousticDOFs+nPaddedDOFs;     }
    if(current_column < nAcousticDOFs) { current_column = current_column*3;                          }
    else                               { current_column = current_column-nAcousticDOFs+nPaddedDOFs;  }
    if(current_row >= nTotalDOFs || current_column >= nTotalDOFs) {
      std::cout<<"Error: matrix file contained corrupt values."<<std::endl;
      std::cout<<"current_row: "<<current_row<<" current_column: "<<current_column<<std::endl;
    }
    if(map->isNodeGlobalElement(current_row)==true) {
      K->insertGlobalValues(current_row,
			    Teuchos::tuple<GO> (current_column),
			    Teuchos::tuple<SC> (cpx_current_value));
      A->insertGlobalValues(current_row,
			    Teuchos::tuple<GO> (current_column),
			    Teuchos::tuple<SC> (cpx_current_value));
      S->insertGlobalValues(current_row,
			    Teuchos::tuple<GO> (current_column),
			    Teuchos::tuple<SC> (cpx_current_value));
    }
  }
  // pad with identity
  for(int i = 0; i < nAcousticDOFs; i++) {
    if(map->isNodeGlobalElement(3*i+1)==true) {
      K->insertGlobalValues(3*i+1,
			    Teuchos::tuple<GO> (3*i+1),
			    Teuchos::tuple<SC> (one));
      A->insertGlobalValues(3*i+1,
			    Teuchos::tuple<GO> (3*i+1),
			    Teuchos::tuple<SC> (one));
      S->insertGlobalValues(3*i+1,
			    Teuchos::tuple<GO> (3*i+1),
			    Teuchos::tuple<SC> (one));
    }
    if(map->isNodeGlobalElement(3*i+2)==true) {
      K->insertGlobalValues(3*i+2,
			    Teuchos::tuple<GO> (3*i+2),
			    Teuchos::tuple<SC> (one));
      A->insertGlobalValues(3*i+2,
			    Teuchos::tuple<GO> (3*i+2),
			    Teuchos::tuple<SC> (one));
      S->insertGlobalValues(3*i+2,
			    Teuchos::tuple<GO> (3*i+2),
			    Teuchos::tuple<SC> (one));
    }
    if(3*i+1 >= nTotalDOFs || 3*i+2 >= nTotalDOFs) {
      std::cout<<"Error: padding contained corrupt values."<<std::endl;
      std::cout<<"3*i+1: "<<3*i+1<<" 3*i+2: "<<3*i+2<<std::endl;
    }
  }
  // damping matrix
  std::ifstream matfile_damp;
  matfile_damp.open("coupled_damp.txt");
  for (int i = 0; i < nnzeros_damp; i++) {
    int current_row, current_column;
    double current_value;
    matfile_damp >> current_row >> current_column >> current_value ;
    SC cpx_current_value(0.0,omega*current_value);
    if(current_row < nAcousticDOFs)    { current_row    = current_row*3;                             }
    else                               { current_row    = current_row-nAcousticDOFs+nPaddedDOFs;     }
    if(current_column < nAcousticDOFs) { current_column = current_column*3;                          }
    else                               { current_column = current_column-nAcousticDOFs+nPaddedDOFs;  }
    if(current_row >= nTotalDOFs || current_column >= nTotalDOFs) {
      std::cout<<"Error: damping matrix file contained corrupt values."<<std::endl;
      std::cout<<"current_row: "<<current_row<<" current_column: "<<current_column<<std::endl;
    }
    if(map->isNodeGlobalElement(current_row)==true) {
      A->sumIntoGlobalValues(current_row,
			     Teuchos::tuple<GO> (current_column),
			     Teuchos::tuple<SC> (cpx_current_value));
      S->sumIntoGlobalValues(current_row,
			     Teuchos::tuple<GO> (current_column),
			     Teuchos::tuple<SC> (cpx_current_value));
    }
  }
  // mass matrix
  std::ifstream matfile_mass;
  matfile_mass.open("coupled_mass.txt");
  for (int i = 0; i < nnzeros_mass; i++) {
    int current_row, current_column;
    double current_value;
    matfile_mass >> current_row >> current_column >> current_value ;
    SC cpx_current_value(-omega2*current_value,0.0);
    SC cpx_current_value_shift=cpx_current_value*shift;
    if(current_row < nAcousticDOFs)    { current_row    = current_row*3;                             }
    else                               { current_row    = current_row-nAcousticDOFs+nPaddedDOFs;     }
    if(current_column < nAcousticDOFs) { current_column = current_column*3;                          }
    else                               { current_column = current_column-nAcousticDOFs+nPaddedDOFs;  }
    if(current_row >= nTotalDOFs || current_column >= nTotalDOFs) {
      std::cout<<"Error: mass matrix file contained corrupt values."<<std::endl;
      std::cout<<"current_row: "<<current_row<<" current_column: "<<current_column<<std::endl;
    }
    if(map->isNodeGlobalElement(current_row)==true) {
      A->sumIntoGlobalValues(current_row,
			     Teuchos::tuple<GO> (current_column),
			     Teuchos::tuple<SC> (cpx_current_value));
      S->sumIntoGlobalValues(current_row,
			     Teuchos::tuple<GO> (current_column),
			     Teuchos::tuple<SC> (cpx_current_value_shift));
    }
  }
  // Complete fill
  K->fillComplete();
  A->fillComplete();
  S->fillComplete();

  // Turn Tpetra::CrsMatrix into MueLu::Matrix
  RCP<XCRS> mueluK_ = rcp(new XTCRS(K)); 
  RCP<XMAT> mueluK  = rcp(new XWRAP(mueluK_));
  RCP<XCRS> mueluA_ = rcp(new XTCRS(A)); 
  RCP<XMAT> mueluA  = rcp(new XWRAP(mueluA_));
  RCP<XCRS> mueluS_ = rcp(new XTCRS(S)); 
  RCP<XMAT> mueluS  = rcp(new XWRAP(mueluS_));
  mueluK->SetFixedBlockSize(nDOFsPerNode);
  mueluA->SetFixedBlockSize(nDOFsPerNode);
  mueluS->SetFixedBlockSize(nDOFsPerNode);
  //xmap=mueluA->getDomainMap();
  //map=Xpetra::toTpetra(xmap);

  // MultiVector of coordinates
  // NOTE: must be of size equal to the number of DOFs!
  RCP<XMV> coordinates;
  coordinates = MultiVectorFactory::Build(xmap, 3);
  for(int k=0; k<nz; k++) {
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) {
	int curidx = i+nx*j+nx*ny*k;
        int curidx0 = curidx*3+0;
        int curidx1 = curidx*3+1;
        int curidx2 = curidx*3+2;
	SC ih = (SC) (i*h);
	SC jh = (SC) (j*h);
	SC kh = (SC) (k*h);
	if(xmap->isNodeGlobalElement(curidx0)==true) {
	  coordinates->replaceGlobalValue(curidx0,0,ih);
	  coordinates->replaceGlobalValue(curidx0,1,jh);
	  coordinates->replaceGlobalValue(curidx0,2,kh);
	}
	if(xmap->isNodeGlobalElement(curidx1)==true) {
	  coordinates->replaceGlobalValue(curidx1,0,ih);
	  coordinates->replaceGlobalValue(curidx1,1,jh);
	  coordinates->replaceGlobalValue(curidx1,2,kh);
	}
	if(xmap->isNodeGlobalElement(curidx2)==true) {
	  coordinates->replaceGlobalValue(curidx2,0,ih);
	  coordinates->replaceGlobalValue(curidx2,1,jh);
	  coordinates->replaceGlobalValue(curidx2,2,kh);
	}
      }
    }
  }

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(mueluK));
  FactoryManager M;

  // Prolongation/Restriction
  RCP<TPFactory>  TentPFact = rcp( new TPFactory     );
  RCP<SaPFactory> Pfact     = rcp( new SaPFactory    );
  RCP<RAPFactory> Acfact    = rcp( new RAPFactory    );
  RCP<RBMFactory> RBMfact   = rcp( new RBMFactory(3) );
  RBMfact->setLastAcousticDOF(nPaddedDOFs);

  // Smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpack2Type;
  Teuchos::ParameterList ifpack2List;
  // ifpack2Type = "KRYLOV";
  // ifpack2List.set("krylov: number of iterations",5);
  // ifpack2List.set("krylov: residual tolerance",1e-6);
  // ifpack2List.set("krylov: block size",1);
  // ifpack2List.set("krylov: zero starting solution",true);
  // ifpack2List.set("krylov: preconditioner type",1);
  // Additive Schwarz smoother
  ifpack2Type = "SCHWARZ";
  ifpack2List.set("schwarz: compute condest", false);
  ifpack2List.set("schwarz: combine mode", "Add"); // use string mode for this
  ifpack2List.set("schwarz: reordering type", "none");
  ifpack2List.set("schwarz: filter singletons", false);  
  ifpack2List.set("schwarz: overlap level", 0);
  // ILUT smoother
  //ifpack2Type = "ILUT";
  //ifpack2List.set("fact: ilut level-of-fill", (double)1.0);
  //ifpack2List.set("fact: absolute threshold", (double)0.0);
  //ifpack2List.set("fact: relative threshold", (double)1.0);
  //ifpack2List.set("fact: relax value", (double)0.0);
  // Gauss-Seidel smoother
  //ifpack2Type = "RELAXATION";
  //ifpack2List.set("relaxation: sweeps", (LO) 1);
  //ifpack2List.set("relaxation: damping factor", (SC) 1.0); // 0.7
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

  RCP<XMV> nullspace;
  RBMfact->BuildRBM(mueluA,coordinates,nullspace);

  // Setup R's and P's
  M.SetFactory("P", Pfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Ptent", TentPFact);
  H->GetLevel(0)->Set("Nullspace",nullspace);
  H->Keep("P", Pfact.get());
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
  if(comm->getRank()==0) 
    B->replaceGlobalValue(0, 1.0);

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
  belosList.set("Flexible Gmres", false);         // set flexible GMRES on

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
  RCP<TMV> resid = Tpetra::createMultiVector<SC,LO,GO,NO>(map, numrhs);     
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

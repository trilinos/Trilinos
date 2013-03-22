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
#include "Teuchos_ConfigDefs.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include "Xpetra_ConfigDefs.hpp"
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
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
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_RebalanceAcFactory.hpp"

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosXpetraAdapter.hpp" // this header defines Belos::XpetraOp()
#include "BelosMueLuAdapter.hpp"  // this header defines Belos::MueLuOp()
#endif

// Only run if we have long long
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef double                               Scalar;
typedef int                                  LocalOrdinal;
typedef long long                            GlobalOrdinal;
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar, LocalOrdinal, Node>::SparseOps LocalMatOps;
//
#include "MueLu_UseShortNames.hpp"

using Xpetra::global_size_t;

int main(int argc, char *argv[]) {
  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  int NumProcs = comm->getSize();
  int MyPID    = comm->getRank();

  if(!MyPID) printf("TwoBillion: Running Test\n");


  const long long FIRST_GID  = 3000000000L;
  //  const long long FIRST_GID  = 0L;
  
  const long long IndexBase  = 0L;
  //  const long long IndexBase  = 3000000000L;

  int num_per_proc = 10;
  
  global_size_t NumGlobalElements = NumProcs*num_per_proc;

  // Create Map w/ GIDs starting at > 2 billion
  RCP<const Map> map;
  RCP<CrsMatrix> Acrs;
  Teuchos::Array<GlobalOrdinal> mygids(num_per_proc);

  for(int i=0; i<num_per_proc; i++)
    mygids[i] = FIRST_GID + MyPID*num_per_proc + i;

  for(int i=0;i<num_per_proc;i++)
    printf("[%d] mygids[%d] = %lld\n",MyPID,i,mygids[i]);


  map = MapFactory::Build(Xpetra::UseTpetra, Teuchos::OrdinalTraits<global_size_t>::invalid(),mygids(),IndexBase,comm);

  //  RCP<Teuchos::FancyOStream> fox = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  //  fox->setOutputToRootOnly(-1);
  //  map->describe(*fox,Teuchos::VERB_EXTREME);

  // Create 1D Laplacian w/ GIDs starting at > 2 billion
  Teuchos::Array<Scalar> myvals(3);
  Teuchos::Array<GlobalOrdinal> mycols(3);
  Teuchos::ArrayView<Scalar> ValView;
  Teuchos::ArrayView<GlobalOrdinal> ColView;


  Acrs = CrsMatrixFactory::Build(map,3);
  for(int i=0; i<num_per_proc; i++) {
    if(mygids[i]==FIRST_GID ) { 
      mycols[0] = mygids[i];     myvals[0] = 2;
      mycols[1] = mygids[i]+1;   myvals[1] = -1;     
      ValView=myvals.view(0,2);
      ColView=mycols.view(0,2);
      //      printf("[%d %lld] cols %lld %lld\n",MyPID,mygids[i],mycols[0],mycols[1]);
    }
    else if(mygids[i] == FIRST_GID + (long long) NumGlobalElements - 1){
      mycols[0] = mygids[i]-1;   myvals[0] = -1;
      mycols[1] = mygids[i];     myvals[1] = 2;     
      ValView=myvals.view(0,2);
      ColView=mycols.view(0,2);
      //      printf("[%d %lld] cols %lld %lld\n",MyPID,mygids[i],mycols[0],mycols[1]);
    }
    else {
      mycols[0] = mygids[i]-1;   myvals[0] = -1;
      mycols[1] = mygids[i];     myvals[1] = -2;     
      mycols[2] = mygids[i]+1;   myvals[1] = -1;     
      ValView=myvals();
      ColView=mycols();
      //      printf("[%d %lld] cols %lld %lld %lld\n",MyPID,mygids[i],mycols[0],mycols[1],mycols[2]);
    }
    Acrs->insertGlobalValues(mygids[i],ColView,ValView);
  }
  Acrs->fillComplete();

  RCP<Matrix> A = rcp(new CrsMatrixWrap(Acrs));




  RCP<MueLu::Hierarchy<SC, LO, GO, NO, LMO> > H;
  //
  //
  // SETUP
  //
  //
  //
  // Hierarchy
  //
  
  H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  
  //
  // Finest level
  //
  
  RCP<Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",           A);

  //
  // FactoryManager
  //
  
  FactoryManager M;
  
  //
  //
  // Aggregation
  //

  RCP<CoupledAggregationFactory> AggregationFact = rcp(new CoupledAggregationFactory());
  M.SetFactory("Aggregates", AggregationFact);

  //
  // Non rebalanced factories
  //
  RCP<SaPFactory> PFact = rcp(new SaPFactory());
  RCP<Factory>    RFact = rcp(new TransPFactory());
  RCP<RAPFactory> AFact = rcp(new RAPFactory());
  AFact->setVerbLevel(Teuchos::VERB_HIGH);
  
  M.SetFactory("P", PFact);
  M.SetFactory("R", RFact);
  M.SetFactory("A", AFact);

  //
  // Smoothers
  //
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 2);
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(ifpackType, ifpackList));
  M.SetFactory("Smoother", rcp(new SmootherFactory(smootherPrototype)));

  //
  // Coarse Solver
  //
  Teuchos::ParameterList ifpackList2;
  ifpackList2.set("relaxation: sweeps", (LO) 6);
  ifpackType = "RELAXATION";
  ifpackList2.set("relaxation: type", "Symmetric Gauss-Seidel");
  RCP<SmootherPrototype> smootherPrototype2 = rcp(new TrilinosSmoother(ifpackType, ifpackList2));
  M.SetFactory("CoarseSolver", rcp(new SmootherFactory(smootherPrototype2)));

  //
  // Setup preconditioner
  //
  int startLevel = 0;
  H->Setup(M, startLevel, 10);

  //
  //
  // SOLVE
  //
  //

  // Define X, B
  RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
  RCP<MultiVector> B = MultiVectorFactory::Build(map, 1);
  Teuchos::Array<ST::magnitudeType> norms(1);

  X->setSeed(846930886);
  X->randomize();
  A->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);
  B->norm2(norms);
  B->scale(1.0/norms[0]);

  //
  // Use AMG as a preconditioner in Belos
  //
#ifdef HAVE_MUELU_BELOS
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
  double optTol = 1e-8;
  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  belosList.set("Convergence Tolerance", optTol);    // Relative convergence tolerance requested
  //belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
  belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList.set("Output Frequency", 1);
  belosList.set("Output Style", Belos::Brief);
  
  // Create an iterative solver manager
  RCP< Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
  
  // Perform solve
  Belos::ReturnType ret = Belos::Unconverged;
  try {
    ret = solver->solve();
    
    // Get the number of iterations for this solve.
    if (comm->getRank() == 0)
      std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
    
    // Compute actual residuals.
    int numrhs = 1;
    std::vector<double> actual_resids( numrhs ); //TODO: double?
    std::vector<double> rhs_norm( numrhs );
    RCP<MultiVector> resid = MultiVectorFactory::Build(map, numrhs);
    
    typedef Belos::OperatorTraits<SC, MV, OP>  OPT;
    typedef Belos::MultiVecTraits<SC, MV>     MVT;
    
    OPT::Apply( *belosOp, *X, *resid );
    MVT::MvAddMv( -1.0, *resid, 1.0, *B, *resid );
    MVT::MvNorm( *resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    *out<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    for ( int i = 0; i<numrhs; i++) {
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
#endif // HAVE_MUELU_BELOS

  return EXIT_SUCCESS;
}


#else
// if we don't have long longs...
int main(int argc, char *argv[]) {
  
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int MyPID = comm->getRank();
  if(!MyPID) printf("TwoBillion: Long long aren't compiled in.\n");
  return EXIT_SUCCESS;
}
#endif

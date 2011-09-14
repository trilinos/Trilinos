#include <unistd.h>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuPrecOp()
#endif

//
typedef double Scalar;
typedef int    LocalOrdinal;
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int GlobalOrdinal;
#else
typedef int GlobalOrdinal;
#endif
//
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
//
#include "MueLu_UseShortNames.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);
  Teuchos::Array<Teuchos::RCP<Teuchos::Time> > mtime;

  //out->setOutputToRootOnly(-1);
  //out->precision(12);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  MueLu::Gallery::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  // custom parameters
  LO maxLevels = 10;
  LO its=10;
  std::string smooType="sgs";
  int pauseForDebugger=0;
  int amgAsSolver=1;
  int amgAsPrecond=1;
  int useExplicitR=0;
  int sweeps=1;
  int maxCoarseSize=50;  //FIXME clp doesn't like long long int
  Scalar SADampingFactor=4./3;
  double tol = 1e-7;
  std::string aggOrdering = "natural";
  int minPerAgg=2;
  int maxNbrAlreadySelected=0;

  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  clp.setOption("its",&its,"number of multigrid cycles");
  clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");
  clp.setOption("fixPoint",&amgAsSolver,"apply multigrid as solver");
  clp.setOption("precond",&amgAsPrecond,"apply multigrid as preconditioner");
  clp.setOption("saDamping",&SADampingFactor,"prolongator damping factor");
  clp.setOption("explicitR",&useExplicitR,"restriction will be explicitly stored as transpose of prolongator");
  clp.setOption("sweeps",&sweeps,"sweeps to be used in SGS (or Chebyshev degree)");
  clp.setOption("maxCoarseSize",&maxCoarseSize,"maximum #dofs in coarse operator");
  clp.setOption("tol",&tol,"stopping tolerance for Krylov method");
  clp.setOption("aggOrdering",&aggOrdering,"aggregation ordering strategy (natural,random,graph)");
  clp.setOption("minPerAgg",&minPerAgg,"minimum #DOFs per aggregate");
  clp.setOption("maxNbrSel",&maxNbrAlreadySelected,"maximum # of nbrs allowed to be in other aggregates");
  clp.setOption("smooType",&smooType,"smoother type ('sgs 'or 'cheby')");
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  if (pauseForDebugger) {
    Utils::PauseForDebugger();
  }
  
  matrixParameters.check();
  xpetraParameters.check();
  // TODO: check custom parameters
  std::transform(smooType.begin(), smooType.end(), smooType.begin(), ::tolower);
  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  if (comm->getRank() == 0) {
    matrixParameters.print();
    xpetraParameters.print();
    // TODO: print custom parameters
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  mtime.push_back(M.getNewTimer("Matrix Build"));
  (mtime.back())->start();
  const RCP<const Map> map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  mtime.back()->stop();
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);

  mtime.push_back(M.getNewTimer("MueLu Setup"));
  mtime.back()->start();
  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H = rcp( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level> Finest = rcp( new MueLu::Level() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  H->SetLevel(Finest);

  Finest->Request("A");
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                //FIXME is implemented
  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  UCAggFact->SetMinNodesPerAggregate(minPerAgg);  //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
       *out << "aggregate ordering :                    NATURAL" << std::endl;
       UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
       *out << "aggregate ordering :                    RANDOM" << std::endl;
       UCAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
       *out << "aggregate ordering :                    GRAPH" << std::endl;
       UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  UCAggFact->SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));

  RCP<SaPFactory>       Pfact = rcp( new SaPFactory(TentPFact) );
  Pfact->SetDampingFactor(SADampingFactor);
  RCP<RAPFactory>       Acfact = rcp( new RAPFactory() );
  RCP<RFactory>         Rfact;
  RCP<GenericPRFactory> PRfact;
  if (useExplicitR) {
    Rfact = rcp( new TransPFactory() );
    PRfact = rcp( new GenericPRFactory(Pfact,Rfact));
  } else {
    PRfact = rcp( new GenericPRFactory(Pfact));
    H->SetImplicitTranspose(true);
    Acfact->SetImplicitTranspose(true);
    if (comm->getRank() == 0) std::cout << "\n\n* ***** USING IMPLICIT RESTRICTION OPERATOR ***** *\n" << std::endl;
  }
  PRfact->SetMaxCoarseSize((GO) maxCoarseSize);

  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) sweeps);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  if (smooType == "sgs") {
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  } else if (smooType == "cheby") {
    ifpackType = "CHEBYSHEV";
    ifpackList.set("chebyshev: degree", (LO) sweeps);
    ifpackList.set("chebyshev: ratio eigenvalue", (SC) 20);
    ifpackList.set("chebyshev: max eigenvalue", (double) -1.0);
    ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
    ifpackList.set("chebyshev: zero starting solution", true);
  }

  smooProto = rcp( new TrilinosSmoother(lib, ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1) 
    SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  status = H->FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  mtime.back()->stop();
  *out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(*out,Teuchos::ParameterList::PrintOptions().indent(2));

  Teuchos::ParameterList amesosList;
  amesosList.set("PrintTiming",true);
  RCP<DirectSolver> coarseProto = rcp( new DirectSolver(lib, "", amesosList) );
  SmootherFactory coarseSolveFact(coarseProto);

  //SmootherFactory coarseSolveFact(smooProto);    //JJH lazy man's way to have a one-level method with smoother
  H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  // Define RHS
  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
  RHS->norm2(norms);
  RHS->scale(1.0/norms[0]);
  
#define AMG_SOLVER
#ifdef AMG_SOLVER
  // Use AMG directly as an iterative method
  if (amgAsSolver) {
    //*out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  
    {
      X->putScalar( (SC) 0.0);

      H->PrintResidualHistory(true);
      mtime.push_back(M.getNewTimer("Fixed Point Solve"));
      mtime.back()->start();
      H->Iterate(*RHS,its,*X);
      mtime.back()->stop();
  
      //X->norm2(norms);
      //*out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
    }
  } //if (fixedPt)
#endif //ifdef AMG_SOLVER

  // Use AMG as a preconditioner in Belos
  if (amgAsPrecond && lib == Xpetra::UseTpetra)
  {
#if defined(HAVE_MUELU_BELOS) && defined(HAVE_MUELU_TPETRA)
#define BELOS_SOLVER
#endif

#ifdef BELOS_SOLVER

    X->putScalar( (SC) 0.0);

    int numrhs=1;  
    RCP<MultiVector> resid = MultiVectorFactory::Build(map,numrhs); 

    typedef ST::magnitudeType                 MT;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>  MV;
    typedef Belos::OperatorT<MV>              OP;

    // Vectors  
    RCP<MV> belosX     = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstTpetraMV(X);
    RCP<MV> belosRHS   = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstTpetraMV(RHS);
    RCP<MV> belosResid = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstTpetraMV(resid);

    // Construct a Belos LinearProblem object
    RCP<OP> belosOp   = rcp (new Belos::MueLuOp<SC,LO,GO,NO,LMO>(Op) );    // Turns a Xpetra::Operator object into a Belos 'OP'
    H->PrintResidualHistory(false);
    RCP<OP> belosPrec = rcp( new Belos::MueLuPrecOp<SC,LO,GO,NO,LMO>(H) ); // Turns a MueLu::Hierarchy  object into a Belos 'OP'

    RCP<Belos::LinearProblem<double,MV,OP> > problem = rcp( new Belos::LinearProblem<double,MV,OP>( belosOp, belosX, belosRHS ) );

    problem->setLeftPrec( belosPrec );
    
    bool set = problem->setProblem();
    if (set == false) {
      *out << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    
    // Create an iterative solver manager.

    // Belos parameter list
    int maxiters = 100;
    Teuchos::ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set("Output Frequency",1);
    belosList.set("Output Style",Belos::Brief);
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
    //belosList.set( "Verbosity", Belos::TimingDetails + Belos::StatusTestDetails);

    RCP< Belos::SolverManager<SC,MV,OP> > solver = rcp( new Belos::BlockCGSolMgr<SC,MV,OP>(problem, rcp(&belosList,false)) );
    
    Belos::ReturnType ret;
    bool badRes = false;

    try{
      // Perform solve
      mtime.push_back(M.getNewTimer("Belos Solve"));
      mtime.back()->start();
      ret = solver->solve();
      mtime.back()->stop();

      // Get the number of iterations for this solve.
      int numIters = solver->getNumIters();
      *out << "Number of iterations performed for this solve: " << numIters << std::endl;
    
      // Compute actual residuals.
      std::vector<double> actual_resids( numrhs ); //TODO: double?
      std::vector<double> rhs_norm( numrhs );

      typedef Belos::OperatorTraits<SC,MV,OP>  OPT;
      typedef Belos::MultiVecTraits<SC,MV>     MVT;
      
      OPT::Apply( *belosOp, *belosX, *belosResid );
      MVT::MvAddMv( -1.0, *belosResid, 1.0, *belosRHS, *belosResid );
      MVT::MvNorm( *belosResid, actual_resids );
      MVT::MvNorm( *belosRHS, rhs_norm );
      *out<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        *out<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) { badRes = true; }
      }
    } //try
    catch(...) {
      *out << std::endl << "ERROR:  Belos threw an error! " << std::endl;
    }

    // Check convergence
    if (ret!=Belos::Converged || badRes) {
      *out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    } else
      *out << std::endl << "SUCCESS:  Belos converged!" << std::endl;
#endif 
  } //if (amgAsPrecond)

  // Final summaries - this eats memory like a hot dog eating contest
  // M.summarize();
    
  int ntimers=mtime.size();
  Teuchos::ArrayRCP<double> lTime(ntimers);
  Teuchos::ArrayRCP<double> gTime(ntimers);

  for(int i=0;i<ntimers;i++) lTime[i]=mtime[i]->totalElapsedTime();

  // Allreduce is my friend.
#ifdef HAVE_MPI
  MPI_Allreduce(&*lTime,&*gTime,ntimers,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
  for(int i=0;i<ntimers;i++) gTime[i] = lTime[i];
#endif

  for(int i=0;i<ntimers;i++) *out<<mtime[i]->name()<<": \t"<<gTime[i]<<endl;

  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  return EXIT_SUCCESS;
}

// TODO: add warning if:
// DEBUG_MODE, LONG_LONG or KLU

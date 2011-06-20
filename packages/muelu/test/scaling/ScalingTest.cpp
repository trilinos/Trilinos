#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaLevel.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
//#include "MueLu_GaussSeidel.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_AggregationOptions.hpp"

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu_Vector.hpp>
#include <Cthulhu_VectorFactory.hpp>
#include <Cthulhu_MultiVectorFactory.hpp>
#include <Cthulhu_Parameters.hpp>

// Gallery
#define CTHULHU_ENABLED // == Gallery have to be build with the support of Cthulhu matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

//#include "MueLu_UseDefaultTypes.hpp"
typedef double Scalar;
typedef int    LocalOrdinal;
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int    GlobalOrdinal;
#else
typedef int GlobalOrdinal;
#warning Teuchos support for long long not enabled.
#endif
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

#include "MueLu_UseShortNames.hpp"
#include <unistd.h>
/**********************************************************************************/

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuPrecOp()

// How many timers do we need?
#define ntimers 3


int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;
 
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);
  int ctime=0;
  Teuchos::RCP<Teuchos::Time> mtime[ntimers];

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

  // custom parameters
  LO maxLevels = 5;
  int nx=100,ny=0,nz=0;
  std::string matrixType("Laplace1D");

  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  clp.setOption("matrixType",&matrixType,"matrix type: Laplace1D, Laplace2D, Laplace3D");
  clp.setOption("nx",&nx,"mesh points in x-direction.");
  clp.setOption("ny",&ny,"mesh points in y-direction.");
  clp.setOption("nz",&nz,"mesh points in z-direction.");

  GlobalOrdinal numGlobalElements = nx; 
  if (matrixType == "Laplace2D") numGlobalElements*=ny;
  else if (matrixType == "Laplace2D") numGlobalElements*=ny*nz;


  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  MueLu::Gallery::Parameters<GO> matrixParameters(clp, numGlobalElements); // manage parameters of the test case
  Cthulhu::Parameters cthulhuParameters(clp);             // manage parameters of cthulhu

  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  cthulhuParameters.check();
  // TODO: check custom parameters

  if (comm->getRank() == 0) {
    matrixParameters.print();
    cthulhuParameters.print();
    // TODO: print custom parameters
  }

/*
  if (cthulhuParameters.GetLib() != Cthulhu::UseTpetra) {
    *out << "This example is Tpetra only" << std::endl;
    return EXIT_FAILURE;
  }
*/

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  mtime[ctime]=M.getNewTimer("Matrix Build");
  mtime[ctime]->start();
  const RCP<const Map> map = MapFactory::Build(cthulhuParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  mtime[ctime]->stop();
  ctime++;

  //  return EXIT_SUCCESS;
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);

  mtime[ctime]=M.getNewTimer("MueLu Setup");
  mtime[ctime]->start();
  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H = rcp( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level<SC,LO,GO,NO,LMO> > Finest = rcp( new MueLu::Level<SC,LO,GO,NO,LMO>() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->SetA(Op);
  Finest->Save("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                //FIXME is implemented

  Finest->Save("NullSpace",nullSpace);
  H->SetLevel(Finest);

  MueLu::AggregationOptions aggOptions;
  aggOptions.SetPrintFlag(6);
  aggOptions.SetMinNodesPerAggregate(3);  //TODO should increase if run anything other than 1D
  aggOptions.SetMaxNeighAlreadySelected(0);
  aggOptions.SetOrdering(MueLu::AggOptions::NATURAL);
  aggOptions.SetPhase3AggCreation(0.5);
  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory(aggOptions));
  RCP<CoalesceDropFactory> cdFact;
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(cdFact,UCAggFact));

  RCP<SaPFactory>       Pfact = rcp( new SaPFactory(TentPFact) );
  RCP<GenericPRFactory> PRfact = rcp( new GenericPRFactory(Pfact));
  RCP<RAPFactory>       Acfact = rcp( new RAPFactory() );

  RCP<SmootherPrototype> smooProto;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  if (cthulhuParameters.GetLib() == Cthulhu::UseEpetra) {
#ifdef HAVE_MUELU_IFPACK
    smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
#endif
  } else if (cthulhuParameters.GetLib() == Cthulhu::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
    smooProto = rcp( new Ifpack2Smoother("RELAXATION",ifpackList) );
#endif
  }
  if (smooProto == Teuchos::null) {
    throw(MueLu::Exceptions::RuntimeError("main: smoother error"));
  }

  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
  //  RCP<SmootherFactory> SmooFact;
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);


  RCP<SmootherPrototype> coarseProto;
  SmootherFactory coarseSolveFact(smooProto);  // **** smoother for the coarse solver!
  H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  Teuchos::ParameterList status;
  status = H->FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  mtime[ctime]->stop();
  ctime++;
  *out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(*out,Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown

/*
  if (cthulhuParameters.GetLib() == Cthulhu::UseEpetra) {
#ifdef HAVE_MUELU_AMESOS
    Teuchos::ParameterList amesosList;
    amesosList.set("PrintTiming",true);
    coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
    //#elif HAVE_MUELU_IFPACK...
#endif
  } else if (cthulhuParameters.GetLib() == Cthulhu::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
    Teuchos::ParameterList ifpack2List;
    ifpack2List.set("fact: ilut level-of-fill",99); // TODO ??
    ifpack2List.set("fact: drop tolerance", 0);
    ifpack2List.set("fact: absolute threshold", 0);
    ifpack2List.set("fact: relative threshold", 0);
    coarseProto = rcp( new Ifpack2Smoother("ILUT",ifpack2List) );
#endif
  }
  if (coarseProto == Teuchos::null) {
    throw(MueLu::Exceptions::RuntimeError("main: coarse smoother error"));
  }
*/



  // Define RHS
  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  X->norm2(norms);
  
  RHS->setSeed(8675309);
  RHS->randomize();

#ifdef AMG_SOLVER
  *out << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  // Use AMG directly as an iterative method
  {
    X->putScalar( (SC) 0.0);

    H->PrintResidualHistory(true);
    H->Iterate(*RHS,its,*X);

    X->norm2(norms);
    *out << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  }

#endif

#define BELOS_SOLVER
#ifdef BELOS_SOLVER
  // Use AMG as a preconditioner in Belos
  {
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
    RCP<OP> belosOp   = rcp (new Belos::MueLuOp<SC,LO,GO,NO,LMO>(Op) );    // Turns a Cthulhu::Operator object into a Belos 'OP'
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
    double tol = 1e-7;
    Teuchos::ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set("Output Frequency",10);
    belosList.set("Output Style",Belos::Brief);
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
    //belosList.set( "Verbosity", Belos::TimingDetails + Belos::StatusTestDetails);

    RCP< Belos::SolverManager<SC,MV,OP> > solver = rcp( new Belos::BlockCGSolMgr<SC,MV,OP>(problem, rcp(&belosList,false)) );
    
    // Perform solve
    mtime[ctime]=M.getNewTimer("Belos Solve");
    mtime[ctime]->start();
    Belos::ReturnType ret = solver->solve();
    mtime[ctime]->stop();
    ctime++;

    // Get the number of iterations for this solve.
    int numIters = solver->getNumIters();
    *out << "Number of iterations performed for this solve: " << numIters << std::endl;
  
    // Compute actual residuals.
    bool badRes = false;
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

    // Final summaries - this eats memory like a hot dog eating contest
    // M.summarize();
    
    double lTime[ntimers];
    double gTime[ntimers];

    for(int i=0;i<ntimers;i++) lTime[i]=mtime[i]->totalElapsedTime();

    // Allreduce is my friend.
    MPI_Allreduce(&lTime,&gTime,ntimers,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    for(int i=0;i<ntimers;i++) *out<<mtime[i]->name()<<": \t"<<gTime[i]<<endl;

    // Check convergence
    if (ret!=Belos::Converged || badRes) {
      *out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
      return EXIT_FAILURE;
    }
    *out << std::endl << "SUCCESS:  Belos converged!" << std::endl;

  } // end of Belos
#endif // JG_TODO

  return EXIT_SUCCESS;

}

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

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"
#include <unistd.h>
/**********************************************************************************/

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuPrecOp()

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;
 
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);
  
  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  MueLu::Gallery::Parameters matrixParameters(clp, 8748); // manage parameters of the test case
  Cthulhu::Parameters cthulhuParameters(clp);             // manage parameters of cthulhu

  // custom parameters
  LO maxLevels = 3;
  LO its=10;
  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  clp.setOption("its",&its,"number of multigrid cycles");
  
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


#ifdef FOR_PARALLEL_DEBUGGING
  //Utils::BreakForDebugger(*comm);

  LO mypid = comm->getRank();

  if (mypid  == 0) std::cout << "Host and Process Ids for tasks" << std::endl;
  for (LO i = 0; i <comm->getSize() ; i++) {
    if (i == mypid ) {
      char buf[80];
      char hostname[80];
      gethostname(hostname, sizeof(hostname));
      LO pid = getpid();
      sprintf(buf, "Host: %s\tMPI rank: %d,\tPID: %d\n\tattach %d\n\tcontinue\n",
              hostname, mypid, pid, pid);
      printf("%s\n",buf);
      fflush(stdout);
      sleep(1);
    }
  }

  if (mypid == 0) {
    printf( "** Enter a character to continue > "); fflush(stdout);
    char go = ' ';
    scanf("%c",&go);
  }
  comm->barrier();
#endif

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(cthulhuParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);
  std::cout << "||NS|| = " << norms[0] << std::endl;

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
  aggOptions.SetMinNodesPerAggregate(3);
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
  ifpackList.set("relaxation: type", "Gauss-Seidel");
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
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  status = H->FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  std::cout  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(std::cout,Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown

  RCP<SmootherPrototype> coarseProto;
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

  SmootherFactory coarseSolveFact(coarseProto);
  H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  // Define RHS
  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  X->norm2(norms);
  std::cout << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  // Use AMG directly as an iterative method
  {
    X->putScalar( (SC) 0.0);

    H->PrintResidualHistory(true);
    H->Iterate(*RHS,its,*X);

    X->norm2(norms);
    std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  }
  
  // Use AMG as a preconditioner in Belos
  {
    X->putScalar( (SC) 0.0);
  
    typedef ST::magnitudeType                 MT;
    typedef Cthulhu::MultiVector<SC>          MV;
    typedef Belos::OperatorT<SC,LO,GO,NO,LMO> OP;
  
    // Construct a Belos LinearProblem object
    RCP<OP> belosOp   = rcp (new Belos::MueLuOp<SC,LO,GO,NO,LMO>(Op) );    // Turns a Cthulhu::Operator object into a Belos operator
    RCP<OP> belosPrec = rcp( new Belos::MueLuPrecOp<SC,LO,GO,NO,LMO>(H) ); // Turns a MueLu::Hierarchy  object into a Belos operator

    RCP<Belos::LinearProblem<double,MV,OP> > problem = rcp( new Belos::LinearProblem<double,MV,OP>( belosOp, X, RHS ) );
    problem->setLeftPrec( belosPrec );
    
    bool set = problem->setProblem();
    if (set == false) {
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    
    // Create an iterative solver manager.

    // Belos parameter list
    int maxiters = 10;
    double tol = 1e-4;
    Teuchos::ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);

    RCP< Belos::SolverManager<double,MV,OP> > solver = rcp( new Belos::BlockCGSolMgr<double,MV,OP>(problem, rcp(&belosList,false)) );
    
    // Perform solve
    Belos::ReturnType ret = solver->solve();
    
    // Get the number of iterations for this solve.
    int numIters = solver->getNumIters();
    std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;
  
    // Compute actual residuals.
    int numrhs=1;
    bool badRes = false;
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    RCP<MultiVector> resid = MultiVectorFactory::Build(map,numrhs); 

    typedef Belos::OperatorTraits<SC,MV,OP>  OPT;
    typedef Belos::MultiVecTraits<SC,MV>     MVT;
    
    OPT::Apply( *belosOp, *X, *resid );
    MVT::MvAddMv( -1.0, *resid, 1.0, *RHS, *resid );
    MVT::MvNorm( *resid, actual_resids );
    MVT::MvNorm( *RHS, rhs_norm );
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    for ( int i=0; i<numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      if (actRes > tol) { badRes = true; }
    }

    // Check convergence
    if (ret!=Belos::Converged || badRes) {
      std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;

  } // end of Belos

  return EXIT_SUCCESS;

}

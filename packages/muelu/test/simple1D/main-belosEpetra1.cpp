#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
//#include "MueLu_GaussSeidel.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Utilities.hpp"

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"
#include <unistd.h>
/**********************************************************************************/

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "Epetra_CrsMatrix.h"

#include "BelosEpetraAdapter.hpp" // this header defines Belos::MueLuPrecOp()
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuPrecOp()

int main(int argc, char *argv[]) {
#ifdef HAVE_MUELU_AMESOS

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
  MueLu::Gallery::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  // custom parameters
  LO maxLevels = 2;
  LO its=10;
  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  clp.setOption("its",&its,"number of multigrid cycles");
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  xpetraParameters.check();
  // TODO: check custom parameters

  if (comm->getRank() == 0) {
    matrixParameters.print();
    xpetraParameters.print();
    // TODO: print custom parameters
  }

  if (xpetraParameters.GetLib() != Xpetra::UseEpetra) {
    std::cout << "This example is Epetra only" << std::endl;
    return EXIT_FAILURE;
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
  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);
  std::cout << "||NS|| = " << norms[0] << std::endl;

  RCP<Hierarchy> H = rcp( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level> Finest = rcp( new MueLu::Level() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                //FIXME is implemented

  Finest->Set("NullSpace",nullSpace);
  H->SetLevel(Finest);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetMinNodesPerAggregate(3);
  UCAggFact->SetMaxNeighAlreadySelected(0);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);

  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));

  RCP<SaPFactory>       Pfact = rcp( new SaPFactory(TentPFact) );
  RCP<RFactory>         Rfact = rcp( new TransPFactory() );
  RCP<GenericPRFactory> PRfact = rcp( new GenericPRFactory(Pfact,Rfact));
  RCP<RAPFactory>       Acfact = rcp( new RAPFactory() );

  RCP<SmootherPrototype> smooProto;


  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  if (xpetraParameters.GetLib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_IFPACK
    ifpackList.set("relaxation: type", "symmetric Gauss-Seidel");
    smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
#endif
  } else if (xpetraParameters.GetLib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
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
  if (xpetraParameters.GetLib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_AMESOS
    Teuchos::ParameterList amesosList;
    amesosList.set("PrintTiming",true);
    coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
    //#elif HAVE_MUELU_IFPACK...
#else
#error ERROR
#endif


  } else if (xpetraParameters.GetLib() == Xpetra::UseTpetra) {
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

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();

  Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  X->norm2(norms);
  std::cout << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  X->putScalar( (SC) 0.0);

  H->PrintResidualHistory(true);
  H->Iterate(*RHS,its,*X);

  X->norm2(norms);
  std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  
  if (xpetraParameters.GetLib() == Xpetra::UseEpetra) {
    std::cout << "- - - - - - - - - - :" << std::endl;
    std::cout << "Epetra Belos run:" << std::endl;

    X->putScalar( (SC) 0.0);

    typedef double                            ST;
    typedef Teuchos::ScalarTraits<ST>        SCT;
    typedef SCT::magnitudeType                MT;
    typedef Belos::MultiVec<double>           MV;
    typedef Belos::Operator<double>           OP;
    typedef Belos::MultiVecTraits<ST,MV>     MVT;
    typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
    
    using Teuchos::ParameterList;

    RCP<Epetra_CrsMatrix> eA = MueLu::Utils<SC,LO,GO,NO,LMO>::Op2NonConstEpetraCrs(Op);
    RCP<Epetra_MultiVector> eX = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstEpetraMV(X);
    RCP<Epetra_MultiVector> eB = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstEpetraMV(RHS);

    RCP<Belos::MueLuEpetraPrecOp> belosPrec = rcp( new Belos::MueLuEpetraPrecOp( H ) );
    
    //
    // *******Construct a preconditioned linear problem********
    //
    RCP<Belos::EpetraOp>      belosOp = rcp (new Belos::EpetraOp(eA));
    RCP<Belos::EpetraMultiVec> belosX = rcp (new Belos::EpetraMultiVec(*eX));
    RCP<Belos::EpetraMultiVec> belosB = rcp (new Belos::EpetraMultiVec(*eB));

    RCP<Belos::LinearProblem<double,MV,OP> > problem
      = rcp( new Belos::LinearProblem<double,MV,OP>( belosOp, belosX, belosB ) );
    problem->setLeftPrec( belosPrec );
    
    bool set = problem->setProblem();
    if (set == false) {
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    
    // Create an iterative solver manager.
    // *****Create parameter list for the belos solver manager*****
    //
    int maxiters = 100;
    double tol = 1e-7;
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + 
		   Belos::TimingDetails + Belos::StatusTestDetails);

    RCP< Belos::SolverManager<double,MV,OP> > solver
      = rcp( new Belos::BlockCGSolMgr<double,MV,OP>(problem, rcp(&belosList,false)) );
    
    // Perform solve
    //
    Belos::ReturnType ret = solver->solve();
    
    //
    // Get the number of iterations for this solve.
    //
    int numIters = solver->getNumIters();
    std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;
  
    //
    // Compute actual residuals.
    //
    int numrhs=1;
    bool badRes = false;
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    RCP<MultiVector> cResid = MultiVectorFactory::Build(map,numrhs); 
    RCP<Epetra_MultiVector> resid = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstEpetraMV(cResid);
    Belos::EpetraMultiVec belosResid(*resid); 

    OPT::Apply( *belosOp, *belosX, belosResid );
    MVT::MvAddMv( -1.0, belosResid, 1.0, *belosB, belosResid );
    MVT::MvNorm( belosResid, actual_resids );
    MVT::MvNorm( *belosB, rhs_norm );

    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    for ( int i=0; i<numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      if (actRes > tol) { badRes = true; }
      }

    if (ret!=Belos::Converged || badRes) {
      std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    } else {
      std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
    }

  }

#endif // #ifdef HAVE_MUELU_AMESOS
  return EXIT_SUCCESS;

}

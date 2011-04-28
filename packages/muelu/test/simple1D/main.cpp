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

int main(int argc, char *argv[]) {
 
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

  MueLu::Hierarchy<SC,LO,GO,NO,LMO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level<SC,LO,GO,NO,LMO> > Finest = rcp( new MueLu::Level<SC,LO,GO,NO,LMO>() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->SetA(Op);
  Finest->Save("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                //FIXME is implemented

  Finest->Save("NullSpace",nullSpace);
  H.SetLevel(Finest);

  RCP<SaPFactory>       Pfact = rcp( new SaPFactory() );
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
  status = H.FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
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
  Teuchos::ParameterList ifpackList;
  ifpackList.set("fact: ilut level-of-fill",99); // TODO ??
  ifpackList.set("fact: drop tolerance", 0);
  ifpackList.set("fact: absolute threshold", 0);
  ifpackList.set("fact: relative threshold", 0);
  coarseProto = rcp( new Ifpack2Smoother("ILUT",ifpackList) );
#endif
  }
  if (coarseProto == Teuchos::null) {
    throw(MueLu::Exceptions::RuntimeError("main: coarse smoother error"));
  }

  SmootherFactory coarseSolveFact(coarseProto);
  H.SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  X->norm2(norms);
  std::cout << "||X_true|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  X->putScalar( (SC) 0.0);

  H.PrintResidualHistory(true);
  H.Iterate(*RHS,its,*X);

  X->norm2(norms);
  std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  
  return EXIT_SUCCESS;

}

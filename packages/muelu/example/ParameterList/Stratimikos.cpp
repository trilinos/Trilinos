#include <iostream>

// Teuchos includes
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

// Epetra includes
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

// EpetraExt includes
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// MueLu includes
#include <MueLu.hpp> // TODO Usefull?
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>  

int main(int argc,char * argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  // MPI initialization
  //
  
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // build global communicator TODO: convert from Teuchos::Comm
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false); 

  MueLu::Gallery::Parameters<int> matrixParameters(clp, 256); // manage parameters of the test case
  // Xpetra::Parameters              xpetraParameters(clp);   // manage parameters of xpetra
  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra; // Epetra only for the moment

  std::string xmlFileName = "stratimikos_ParameterList.xml"; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'stratimikos_ParameterList.xml'.");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  // Read in parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName); 

  //
  // Construct the problem
  //

  //   // Read in the matrix, store pointer as an RCP
  //   Epetra_CrsMatrix * ptrA = 0;
  //   EpetraExt::MatrixMarketFileToCrsMatrix("../data/nsjac_test.mm",Comm,ptrA);
  //   RCP<Epetra_CrsMatrix> A = rcp(ptrA);
  //
  //   // read in the RHS vector
  //   Epetra_Vector * ptrb = 0;
  //   EpetraExt::MatrixMarketFileToVector("../data/nsrhs_test.mm",A->OperatorRangeMap(),ptrb);
  //   RCP<const Epetra_Vector> b = rcp(ptrb);
  
  RCP<const Map> map = MapFactory::createUniformContigMap(lib, matrixParameters.GetNumGlobalElements(), comm);
  RCP<const Epetra_CrsMatrix> A = MueLu::Gallery::CreateCrsMatrix<double, int, int, Map,  Xpetra::EpetraCrsMatrix>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList())->getEpetra_CrsMatrix();

  //
  // Allocate vectors
  //

  RCP<Epetra_Vector> X = rcp(new Epetra_Vector(A->DomainMap()));
  X->PutScalar(0.0);

  RCP<Epetra_Vector> B = rcp(new Epetra_Vector(A->DomainMap()));
  B->SetSeed(846930886); B->Random();
  
  //
  // Build Thyra linear algebra objects
  //

  RCP<const Thyra::LinearOpBase<double> > thyraA = Thyra::epetraLinearOp(A);
  RCP<const Thyra::VectorBase<double> >   thyraB = Thyra::create_Vector(B, thyraA->range());
  RCP<Thyra::VectorBase<double> >         thyraX = Thyra::create_Vector(X, thyraA->domain());

  //
  // Build Stratimikos solver
  //

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;  // This is the Stratimikos main class (= factory of solver factory).
  Thyra::addMueLuToStratimikosBuilder(linearSolverBuilder);     // Register MueLu as a Stratimikos preconditioner strategy.
  linearSolverBuilder.setParameterList(paramList);              // Setup solver parameters using a Stratimikos parameter list.

  // Build a new "solver factory" according to the previously specified parameter list.
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder); 

  // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver.
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > thyraInverseA = Thyra::linearOpWithSolve(*solverFactory, thyraA);

  //
  // Solve Ax = b.
  //

  Thyra::SolveStatus<double> status = Thyra::solve<double>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());
  std::cout << status << std::endl;

  return (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED) ? EXIT_SUCCESS : EXIT_FAILURE;
}

// Thyra::assign(thyraX.ptr(), 0.0);

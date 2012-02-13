#include <iostream>

#include <MueLu_ConfigDefs.hpp>

#include <Teuchos_XMLParameterListHelpers.hpp> // getParametersFromXmlFile()
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
#include <Epetra_CrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MLInterpreter.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

// Default problem is Laplace1D with nx = 8748. Use --help to list available options.

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false); // Note: 

  MueLu::Gallery::Parameters<GO> matrixParameters(clp, 256); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);      // manage parameters of xpetra

  std::string xmlFileName; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default an hard-coded parameter list.");
  int muelu = true;        clp.setOption("muelu", &muelu,       "use muelu"); //TODO: bool instead of int
  int ml    = false;       
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
  clp.setOption("ml",    &ml,          "use ml");
#endif

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  // TODO: -ml and --linAlgebra

  if (comm->getRank() == 0) { std::cout << xpetraParameters << matrixParameters; }
  if (ml && xpetraParameters.GetLib() == Xpetra::UseTpetra) {
    ml = 0;
    std::cout << "ML preconditionner can only be built if --linAlgebra=0 (Epetra). Option --ml ignored" << std::endl;
  }

  //
  // Construct the problem
  //

  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator>  A   = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  //
  // Preconditionner configuration
  //
  
  // ML parameter list
  RCP<Teuchos::ParameterList> params;
  if (xmlFileName != "") {

    std::cout << "Reading " << xmlFileName << " ..." << std::endl;
    params = Teuchos::getParametersFromXmlFile(xmlFileName);

  } else {
    std::cout << "Hard-coded parameter list" << std::endl;
      params = rcp(new Teuchos::ParameterList());
      
      params->set("max levels", 2);
      
      Teuchos::ParameterList & l0 = params->sublist("smoother: list (level 0)");
      l0.set("smoother: damping factor", 0.9);
      l0.set("smoother: sweeps", 1);
      l0.set("smoother: pre or post", "both");
      l0.set("smoother: type", "symmetric Gauss-Seidel");
  }
    
  if (muelu) {

    //
    // Construct a multigrid preconditioner
    //
    
    // Multigrid Hierarchy
    MLInterpreter mueLuFactory(*params);
    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->GetLevel(0)->Set("A", A);
    
    mueLuFactory.SetupHierarchy(*H);

    //
    // Solve Ax = b
    //
    
    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);
    
    X->putScalar((Scalar) 0.0);
    B->setSeed(846930886); B->randomize();
    
    // Use AMG directly as an iterative solver (not as a preconditionner)
    int nIts = 9;
    
    H->Iterate(*B, nIts, *X);
    
    // Print relative residual norm
    ST::magnitudeType residualNorms = Utils::ResidualNorm(*A, *X, *B)[0];
    if (comm->getRank() == 0)
      std::cout << "||Residual|| = " << residualNorms << std::endl;
  }

#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
  if (ml) {
    
    //
    // Construct a multigrid preconditioner
    //

    // Multigrid Hierarchy
    RCP<CrsOperator>      crsOp         = Teuchos::rcp_dynamic_cast<CrsOperator>(A, true);
    RCP<CrsMatrix>        crsMtx        = crsOp->getCrsMatrix();
    RCP<EpetraCrsMatrix>  epetraCrsMtx  = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(crsMtx, true);
    RCP<const Epetra_CrsMatrix> epetra_CrsMtx = epetraCrsMtx->getEpetra_CrsMatrix();

    RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*epetra_CrsMtx, *params));
    
    //
    // Solve Ax = b
    //
  }
#endif

  return EXIT_SUCCESS;
}

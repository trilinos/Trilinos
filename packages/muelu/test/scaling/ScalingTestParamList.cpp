#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>
#include <MueLu_GalleryUtils.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>  

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false); // Note: 

  MueLu::Gallery::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);      // manage parameters of xpetra

  std::string xmlFileName = "scalingTest.xml"; clp.setOption("xml",   &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'muelu_ParameterList.xml'");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  if (comm->getRank() == 0) { std::cout << xpetraParameters << matrixParameters; }

  //
  // Construct the problem
  //

  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator>  A   = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  //
  // Construct a multigrid preconditioner
  //

  // Multigrid Hierarchy
  ParameterListInterpreter mueLuFactory(xmlFileName);
  RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();

  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  H->GetLevel(0)->Set("A", A);

  {
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(map,1);
    nullspace->putScalar( (SC) 1.0);
    Teuchos::Array<ST::magnitudeType> norms(1);
    H->GetLevel(0)->Set("Nullspace", nullspace);
  }
  
  {
    RCP<MultiVector> coordinates;
    
    if (matrixParameters.GetMatrixType() == "Laplace1D") {
      coordinates = MueLu::GalleryUtils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",map,matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace2D") {
      coordinates = MueLu::GalleryUtils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D",map,matrixParameters.GetParameterList());
    }
    else if (matrixParameters.GetMatrixType() == "Laplace3D") {
      coordinates = MueLu::GalleryUtils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D",map,matrixParameters.GetParameterList());
    }

    H->GetLevel(0)->Set("Coordinates", coordinates);
  }

  mueLuFactory.SetupHierarchy(*H);

  //
  // Solve Ax = b
  //

  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);
  
  {
    X->setSeed(846930886);
    X->randomize();
    A->apply(*X,*B,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Teuchos::Array<ST::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(1.0/norms[0]);
    X->putScalar( (SC) 0.0);
  }

  // Use AMG directly as an iterative solver (not as a preconditionner)
  int nIts = 10;

  H->IsPreconditioner(true);
  H->Iterate(*B,nIts,*X);
  H->IsPreconditioner(false);
  
  return EXIT_SUCCESS;
}

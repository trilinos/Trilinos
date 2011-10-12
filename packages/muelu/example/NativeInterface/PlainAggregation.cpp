// Note: use --help to list available options.

#include <iostream>

// MueLu
#include "MueLu.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"

// Define template parameters
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  //
  // MPI initialization
  //

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Process command line arguments
  //

  Teuchos::CommandLineProcessor  clp(false);
  MueLu::Gallery::Parameters<GO> matrixParameters(clp, 81); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);     // manage parameters of xpetra
  
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
  default:;
  }
  
  if (comm->getRank() == 0) std::cout << xpetraParameters << matrixParameters;

  //
  // Setup test case (Ax = b)
  //

  // Linear Algebra Library
  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  // Distribution
  RCP<const Map> map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);

  // Matrix
  RCP<Operator> A = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  // User defined nullspace
  RCP<MultiVector> nullSpace = VectorFactory::Build(map,1); nullSpace->putScalar((SC) 1.0);

  // Define B
  RCP<Vector> X = VectorFactory::Build(map,1);
  RCP<Vector> B = VectorFactory::Build(map,1);
  X->setSeed(846930886);
  X->randomize();
  A->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

  // X = 0
  X->putScalar((SC) 0.0);

  //
  // Create a multigrid configuration
  //

  // Transfer operators
  RCP<TentativePFactory> pFact = rcp( new TentativePFactory() );

  // Aggregation
  RCP<UCAggregationFactory> aggregationFact = rcp( new UCAggregationFactory() );
  aggregationFact->SetMinNodesPerAggregate(3);

  // Smoothers
  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LO) 1);
  smootherParamList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype> smootherPrototype     = rcp( new TrilinosSmoother(lib, "RELAXATION", smootherParamList) );
  RCP<SmootherFactory>   smootherFact          = rcp( new SmootherFactory(smootherPrototype) );

  // Coarse grid correction
  RCP<SmootherPrototype> coarseSolverPrototype = rcp( new DirectSolver(lib) );
  RCP<SmootherFactory>   coarseSolverFact      = rcp( new SmootherFactory(coarseSolverPrototype, Teuchos::null) );

  // 
  FactoryManager M;
  M.SetFactory("A",            rcp(new RAPFactory())); //TODO: to be remove, but will require some work
  M.SetFactory("P",            pFact);
  M.SetFactory("Aggregation",  aggregationFact);
  M.SetFactory("Smoother",     smootherFact);
  M.SetFactory("CoarseSolver", coarseSolverFact);

  //
  // Multigrid setup phase
  //  

  Hierarchy H;

  RCP<Level> finestLevel = H.GetLevel();
  finestLevel->Set("A", A);
  finestLevel->Set("Nullspace", nullSpace);

  H.Setup(M);

  //
  // Solve Ax = B
  //

  LO nIts = 9;
  H.Iterate(*B, nIts, *X);

  //
  // Print relative residual norm
  //

  ST::magnitudeType residualNorms = Utils::ResidualNorm(*A, *X, *B)[0];
  if (comm->getRank() == 0)
    std::cout << "||Residual|| = " << std::setiosflags(ios::fixed) << std::setprecision(20) << residualNorms << std::endl;

  return EXIT_SUCCESS;

}

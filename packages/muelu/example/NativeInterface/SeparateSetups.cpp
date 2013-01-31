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
// Note: use --help to list available options.

/*
  This example demonstrates:
    1) How to create the hierarchy in a separate step from the smoothers.
    2) How grid transfers can be generated using a different matrix than the smoothers.
    3) How grid transfers P & R can be generated using matrix A, and coarse grid operators
       can be generated via R*B*P.
*/

#include <iostream>

// MueLu
#include "MueLu.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_RAPFactory.hpp"

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

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
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 81); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);     // manage parameters of xpetra

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
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
  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
  RCP<Matrix> A = Pr->BuildMatrix();

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
  RCP<TentativePFactory> TentativePFact = rcp( new TentativePFactory() );
  RCP<SaPFactory>        SaPFact        = rcp( new SaPFactory() );
  RCP<TransPFactory>     RFact          = rcp( new TransPFactory());

  FactoryManager M;
  M.SetFactory("Ptent", TentativePFact);
  M.SetFactory("P",     SaPFact);
  M.SetFactory("R",     RFact);

  M.SetFactory("Smoother", Teuchos::null);      //skips smoother setup
  M.SetFactory("CoarseSolver", Teuchos::null);  //skips coarsest solve setup

  //
  // Multigrid setup phase
  //

  int startLevel = 0;
  int maxLevels = 10;

  std::cout << "=============== Setup transfers only ====================" << std::endl;

  Hierarchy H;
  H.SetDefaultVerbLevel(MueLu::High);

  RCP<Level> finestLevel = H.GetLevel();
  finestLevel->Set("A", A);
  finestLevel->Set("Nullspace", nullSpace);

  // Indicate which Hierarchy operators we want to keep
  H.Keep("P", SaPFact.get());  //SaPFact is the generating factory for P.
  H.Keep("R", RFact.get());    //RFact is the generating factory for R.
  H.Keep("Ptent", TentativePFact.get());  //SaPFact is the generating factory for P.

  H.Setup(M,startLevel,maxLevels);
  Teuchos::ParameterList status = H.Summarize(MueLu::None);
  int numLevels = status.get("number of levels",-1);

  std::cout << std::endl << std::endl;

  std::cout << "=========================================\n"
            << "Status of the preconditioner between runs\n"
            << "=========================================" << std::endl;
  std::cout << "# levels = " << numLevels << std::endl;
  H.print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);

  std::cout << "=============== Setup smoothers only ====================" << std::endl;

  // Create a new A.
  RCP<Matrix> newA = Pr->BuildMatrix();
  finestLevel->Set("A", newA);

  // Create Gauss-Seidel smoother.
  std::string ifpackType = "RELAXATION";
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 3);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(ifpackType, ifpackList));

  M.SetFactory("Smoother", rcp(new SmootherFactory(smootherPrototype)));

  // Create coarsest solver.
  RCP<SmootherPrototype> coarseSolverPrototype = rcp( new DirectSolver() );
  RCP<SmootherFactory>   coarseSolverFact      = rcp( new SmootherFactory(coarseSolverPrototype, Teuchos::null) );
  M.SetFactory("CoarseSolver", coarseSolverFact);

  // Note that we pass the number of levels back in.
  H.Setup(M,startLevel, numLevels);

  std::cout << "=============== Solve ====================" << std::endl;

  //
  // Solve Ax = B
  //

  LO nIts = 9;
  H.Iterate(*B, nIts, *X);

  status = H.Summarize(MueLu::High);


  //
  // Print relative residual norm
  //

  ST::magnitudeType residualNorms = Utils::ResidualNorm(*A, *X, *B)[0];
  if (comm->getRank() == 0)
    std::cout << "||Residual|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorms << std::endl;

  return EXIT_SUCCESS;

}

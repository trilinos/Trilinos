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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Factory.hpp>

#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_NoFactory.hpp>

#include <MueLu_UseDefaultTypes.hpp>

// This example demonstrates how to reuse some parts of a multigrid setup between runs.
//
// In this example, we suppose that the pattern of the matrix does not change between runs so that:
// - Aggregates can be reused
// - The tentative prolongator of Smoothed-Aggregation does not change (as it derived directly from the aggregate information).
// - The pattern of coarse grid A can be reused during its computation
//
// The resulting preconditioners are identical to multigrid preconditioners built without recycling the parts described above.
// This can be verified by using the --no-recycling option.

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  typedef Teuchos::ScalarTraits<SC> ST;
  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);
  Xpetra::Parameters xpetraParameters(clp);

  bool optRecycling           = true;  clp.setOption("recycling",             "no-recycling",             &optRecycling,           "Enable recycling of the multigrid preconditioner");


  switch (clp.parse(argc, argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  //
  // Construct the problems
  //

  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  Teuchos::RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
    Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
  RCP<Matrix> A1 = Pr->BuildMatrix();

  RCP<Matrix> A2 = Pr->BuildMatrix();
  A2->resumeFill();
  A2->scale(2.0);
  A2->fillComplete();

  RCP<Matrix> A3 = Pr->BuildMatrix();
  A3->resumeFill();
  A3->scale(3.0);
  A3->fillComplete();

  //
  // First solve
  //

  FactoryManager M;

  Hierarchy H(A1);

  if (optRecycling) {
    // Configuring "Keep" options before running the first setup.
    // PTENT:
    RCP<Factory> PtentFact = rcp(new TentativePFactory());
    M.SetFactory("Ptent", PtentFact);
    H.Keep("P",           PtentFact.get());
  }

  RCP<Factory> AcFact = rcp(new RAPFactory());
  M.SetFactory("A", AcFact);

  //

  H.Setup(M,0,3);

  {
    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);

    X->putScalar((Scalar) 0.0);
    B->setSeed(846930886); B->randomize();

    int nIts = 9;
    H.Iterate(*B, *X, nIts);

    ST::magnitudeType residualNorms = Utils::ResidualNorm(*A1, *X, *B)[0];
    if (comm->getRank() == 0)
      std::cout << "||Residual|| = " << residualNorms << std::endl;
  }

  //
  // Second solve
  //

  cout << "Status of the preconditioner between runs:" << std::endl;
  H.print(*getFancyOStream(Teuchos::rcpFromRef(cout)), MueLu::High);

  // Change the problem
  RCP<Level> finestLevel = H.GetLevel(0);
  finestLevel->Set("A", A2);

  if (optRecycling) {
    // Optional: this makes sure that the aggregates are never requested (and built) during the second run.
    // If someone request the aggregates, an exception will be thrown.
    M.SetFactory("Aggregates", MueLu::NoFactory::getRCP());
  }

  // Redo the setup
  H.Setup(M,0,3);

  {
    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);

    X->putScalar((Scalar) 0.0);
    B->setSeed(846930886); B->randomize();

    int nIts = 9;
    H.Iterate(*B, *X, nIts);

    ST::magnitudeType residualNorms = Utils::ResidualNorm(*A2, *X, *B)[0];
    if (comm->getRank() == 0)
      std::cout << "||Residual|| = " << residualNorms << std::endl;
  }

  //
  // Third solve
  //

  cout << "Status of the preconditioner between runs:" << std::endl;
  H.print(*getFancyOStream(Teuchos::rcpFromRef(cout)), MueLu::High);

  // Change the problem
  finestLevel = H.GetLevel(0);
  finestLevel->Set("A", A3);

  if (optRecycling) {
    // Optional: this makes sure that the aggregates are never requested (and built) during the second run.
    // If someone request the aggregates, an exception will be thrown.
    M.SetFactory("Aggregates", MueLu::NoFactory::getRCP());
  }

  // Redo the setup
  H.Setup(M,0,3);

  {
    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);

    X->putScalar((Scalar) 0.0);
    B->setSeed(846930886); B->randomize();

    int nIts = 9;
    H.Iterate(*B, *X, nIts);

    ST::magnitudeType residualNorms = Utils::ResidualNorm(*A2, *X, *B)[0];
    if (comm->getRank() == 0)
      std::cout << "||Residual|| = " << residualNorms << std::endl;
  }


  //
  // Clean-up
  //

  // Remove kept data from the preconditioner. This will force recomputation on future runs. "Keep" flags are also removed.

  if (optRecycling) {
    //if aggregates explicitly kept: H.Delete("Aggregates", M.GetFactory("Aggregates").get());
    H.Delete("P",           M.GetFactory("Ptent").get());
  }

  cout << "Status of the preconditioner at the end:" << std::endl;
  H.print(*getFancyOStream(Teuchos::rcpFromRef(cout)), MueLu::High);

  return EXIT_SUCCESS;
}

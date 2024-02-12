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

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include "MueLu_Exceptions.hpp"

#include "MueLu_Hierarchy.hpp"
#include "MueLu_EminPFactory.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_NullspacePresmoothFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#include "MueLu_ProjectorSmoother.hpp"

#include "MueLu_MergedSmoother.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_VerbosityLevel.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include <unistd.h>
/**********************************************************************************/

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosMueLuAdapter.hpp"  // this header defines Belos::MueLuOp()

#define NEUMANN
#define EMIN

using Teuchos::ArrayRCP;
using Teuchos::RCP;
using Teuchos::rcp;

namespace MueLuTests {

#include "MueLu_UseShortNames.hpp"

RCP<SmootherPrototype> gimmeGaussSeidelProto(Xpetra::UnderlyingLib lib) {
  RCP<SmootherPrototype> smooProto;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO)1);
  ifpackList.set("relaxation: damping factor", (SC)1.0);

  if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
#if defined(HAVE_MUELU_SERIAL)
    ifpackList.set("relaxation: type", "symmetric Gauss-Seidel");
    smooProto = MueLu::GetIfpackSmoother<SC, LO, GO, NO>("point relaxation stand-alone", ifpackList);
#else
    throw(MueLu::Exceptions::RuntimeError("gimmeGaussSeidelProto: IfpackSmoother only available with SerialNode."));
#endif
#endif
  } else if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_MUELU_IFPACK2)
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = rcp(new Ifpack2Smoother("RELAXATION", ifpackList));
#endif
  }
  if (smooProto == Teuchos::null) {
    throw(MueLu::Exceptions::RuntimeError("gimmeGaussSeidelSmoother: smoother error"));
  }

  return smooProto;
}

RCP<SmootherPrototype> gimmeCoarseProto(Xpetra::UnderlyingLib lib, const std::string& coarseSolver, int rank) {
  RCP<SmootherPrototype> coarseProto;
  if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)
#if defined(HAVE_MUELU_SERIAL)
    if (rank == 0) std::cout << "CoarseGrid: AMESOS" << std::endl;
    Teuchos::ParameterList amesosList;
    amesosList.set("PrintTiming", true);
    coarseProto = MueLu::GetAmesosSmoother<SC, LO, GO, NO>("Amesos_Klu", amesosList);
#else
    throw(MueLu::Exceptions::RuntimeError("gimmeGaussSeidelProto: AmesosSmoother only available with SerialNode."));
#endif
#endif
  } else if (lib == Xpetra::UseTpetra) {
    if (coarseSolver == "amesos2") {
#if defined(HAVE_MUELU_AMESOS2)
      if (rank == 0) std::cout << "CoarseGrid: AMESOS2" << std::endl;
      Teuchos::ParameterList paramList;  // unused
      coarseProto = rcp(new Amesos2Smoother("Superlu", paramList));
#else
      std::cout << "AMESOS2 not available (try --coarseSolver=ifpack2)" << std::endl;
      return Teuchos::null;  // TODO test for exception //EXIT_FAILURE;
#endif  // HAVE_MUELU_AMESOS2
    } else if (coarseSolver == "ifpack2") {
#if defined(HAVE_MUELU_IFPACK2)
      if (rank == 0) std::cout << "CoarseGrid: IFPACK2" << std::endl;
      Teuchos::ParameterList ifpack2List;
      ifpack2List.set("fact: ilut level-of-fill", (double)99);  // TODO ??
      ifpack2List.set("fact: drop tolerance", 0.0);
      ifpack2List.set("fact: absolute threshold", 0.0);
      ifpack2List.set("fact: relative threshold", 1.0);
      coarseProto = rcp(new Ifpack2Smoother("ILUT", ifpack2List));
#else
      std::cout << "IFPACK2 not available (try --coarseSolver=amesos2)" << std::endl;
      // TODO        TEUCHOS_TEST_FOR_EXCEPTION
      return Teuchos::null;
#endif
    } else {
      std::cout << "Unknow coarse grid solver (try  --coarseSolver=ifpack2 or --coarseSolver=amesos2)" << std::endl;
      return Teuchos::null;
    }
  }
  if (coarseProto == Teuchos::null) {
    throw(MueLu::Exceptions::RuntimeError("main: coarse smoother error"));
  }

  return coarseProto;
}

RCP<SmootherPrototype> gimmeMergedSmoother(int nSmoothers, Xpetra::UnderlyingLib lib, const std::string& coarseSolver, int rank) {
  ArrayRCP<RCP<SmootherPrototype> > smootherList(nSmoothers);

  for (int i = 0; i < nSmoothers; i++)
    smootherList[i] = gimmeGaussSeidelProto(lib);

  return rcp(new MergedSmoother(smootherList));
  // verbose mode: return rcp (new MergedSmoother(smootherList, true));
}

}  // namespace MueLuTests

int main(int argc, char* argv[]) {
#include "MueLu_UseShortNames.hpp"

  using namespace MueLuTests;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    /**********************************************************************************/
    /* SET TEST PARAMETERS                                                            */
    /**********************************************************************************/
    // Note: use --help to list available options.
    Teuchos::CommandLineProcessor clp(false);

    // Default is Laplace1D with nx = 8748.
    // It's a nice size for 1D and perfect aggregation. (6561=3^8)
    // Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                    // manage parameters of xpetra

    // custom parameters
    int nSmoothers           = 2;
    LO maxLevels             = 3;
    LO its                   = 10;
    std::string coarseSolver = "ifpack2";
    // std::string coarseSolver="amesos2";
    clp.setOption("nSmoothers", &nSmoothers, "number of Gauss-Seidel smoothers in the MergedSmoothers");
    clp.setOption("maxLevels", &maxLevels, "maximum number of levels allowed. If 1, then a MergedSmoother is used on the coarse grid");
    clp.setOption("its", &its, "number of multigrid cycles");
    clp.setOption("coarseSolver", &coarseSolver, "amesos2 or ifpack2 (Tpetra specific. Ignored for Epetra)");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    matrixParameters.check();
    xpetraParameters.check();
    // TODO: check custom parameters

    if (comm->getRank() == 0) {
      // matrixParameters.print();
      // xpetraParameters.print();
      // TODO: print custom parameters
    }

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/
    const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());  // TODO: Matrix vs. CrsMatrixWrap
    RCP<Matrix> Op = Pr->BuildMatrix();

#ifdef NEUMANN
    // Tranform matrix to Neumann b.c.
    // Essentially, we need to update two diagonal elements

    // TODO: calls to getLocalRowView not really needed

    Op->resumeFill();

    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const SC> values;
    Teuchos::Array<SC> newValues(2, 0.0);

    size_t myRank = Op->getRowMap()->getComm()->getRank();
    size_t nCpus  = Op->getRowMap()->getComm()->getSize();
    if (myRank == 0) {  // JG TODO: can we use rowMap->isNodeLocalElement(0) instead for more genericity?
      // LO firstRow = 0;
      newValues[0] = 1.0;
      newValues[1] = -1.0;
      Op->getLocalRowView(0, indices, values);
      Op->replaceLocalValues(0, indices, newValues);
    }
    if (myRank == nCpus - 1) {  // JG TODO: can we use rowMap->isNodeLocalElement(lastRow) instead for more genericity?
      LO lastRow   = Op->getLocalNumRows() - 1;
      newValues[0] = -1.0;
      newValues[1] = 1.0;
      Op->getLocalRowView(lastRow, indices, values);
      Op->replaceLocalValues(lastRow, indices, newValues);
    }

    Op->fillComplete();
#endif  // NEUMANN

    /**********************************************************************************/
    /*                                                                                */
    /**********************************************************************************/

    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
    nullSpace->putScalar((SC)1.0);
    Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
    nullSpace->norm1(norms);
    if (comm->getRank() == 0)
      std::cout << "||NS|| = " << norms[0] << std::endl;

    RCP<MueLu::Hierarchy<SC, LO, GO, NO> > H = rcp(new Hierarchy());
    H->SetDefaultVerbLevel(MueLu::Extreme);
    H->IsPreconditioner(false);

    RCP<MueLu::Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A", Op);
    Finest->Set("Nullspace", nullSpace);

    FactoryManager M;

    M.SetFactory("Aggregates", rcp(new UncoupledAggregationFactory()));
    M.SetFactory("Ptent", rcp(new TentativePFactory()));
    M.SetFactory("P", rcp(new SaPFactory()));

#ifdef EMIN
    // Energy-minimization
    RCP<PatternFactory> PatternFact = rcp(new PatternFactory());
#if 0
    PatternFact->SetFactory("P", M.GetFactory("Ptent"));
#else
    PatternFact->SetFactory("P", M.GetFactory("P"));
#endif
    M.SetFactory("Ppattern", PatternFact);

    RCP<EminPFactory> EminPFact = rcp(new EminPFactory());
    EminPFact->SetFactory("P", M.GetFactory("Ptent"));
    M.SetFactory("P", EminPFact);

    RCP<NullspacePresmoothFactory> NullPreFact = rcp(new NullspacePresmoothFactory());
    NullPreFact->SetFactory("Nullspace", M.GetFactory("Nullspace"));
    M.SetFactory("Nullspace", NullPreFact);
#endif

    RCP<SmootherPrototype> smooProto = gimmeMergedSmoother(nSmoothers, xpetraParameters.GetLib(), coarseSolver, comm->getRank());
    M.SetFactory("Smoother", rcp(new SmootherFactory(smooProto)));

    Teuchos::ParameterList status;

    RCP<SmootherPrototype> coarseProto;
    if (maxLevels != 1)
      coarseProto = gimmeCoarseProto(xpetraParameters.GetLib(), coarseSolver, comm->getRank());
    else
      coarseProto = gimmeMergedSmoother(nSmoothers, xpetraParameters.GetLib(), coarseSolver, comm->getRank());

    if (coarseProto == Teuchos::null)
      return EXIT_FAILURE;

#ifdef NEUMANN
    // Use coarse level projection solver
    RCP<SmootherPrototype> projectedSolver = rcp(new ProjectorSmoother(coarseProto));
    RCP<SmootherFactory> coarseSolveFact   = rcp(new SmootherFactory(projectedSolver, Teuchos::null));
#else
    RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(coarseProto, Teuchos::null));
#endif
    M.SetFactory("CoarseSolver", coarseSolveFact);

    H->EnableGraphDumping("graph.dot", 2);

    H->Setup(M, 0, maxLevels);
    // if (comm->getRank() == 0) {
    //   std::cout  << "======================\n Multigrid statistics \n======================" << std::endl;
    //   status.print(std::cout,Teuchos::ParameterList::PrintOptions().indent(2));
    // }

    // Define RHS
    RCP<MultiVector> X   = MultiVectorFactory::Build(map, 1);
    RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

    X->setSeed(846930886);
    X->randomize();
    X->norm2(norms);
    if (comm->getRank() == 0)
      std::cout << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

    Op->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

    // Use AMG directly as an iterative method
    {
      X->putScalar((SC)0.0);

      H->Iterate(*RHS, *X, its);

      X->norm2(norms);
      if (comm->getRank() == 0)
        std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

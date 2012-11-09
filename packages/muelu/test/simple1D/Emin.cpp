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
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

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

#include "MueLu_MergedSmoother.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_VerbosityLevel.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"
#include <unistd.h>
/**********************************************************************************/

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuOp()

// #define NEUMANN

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;

RCP<SmootherPrototype> gimmeGaussSeidelProto(Xpetra::UnderlyingLib lib) {

  RCP<SmootherPrototype> smooProto;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);

  if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_IFPACK
    ifpackList.set("relaxation: type", "symmetric Gauss-Seidel");
    smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
#endif
  } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = rcp( new Ifpack2Smoother("RELAXATION", ifpackList) );
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
#ifdef HAVE_MUELU_AMESOS
    if (rank == 0) std::cout << "CoarseGrid: AMESOS" << std::endl;
    Teuchos::ParameterList amesosList;
    amesosList.set("PrintTiming",true);
    coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
    //#elif HAVE_MUELU_IFPACK...
#endif
  } else if (lib == Xpetra::UseTpetra) {
    if (coarseSolver=="amesos2") {
#ifdef HAVE_MUELU_AMESOS2
      if (rank == 0) std::cout << "CoarseGrid: AMESOS2" << std::endl;
      Teuchos::ParameterList paramList; //unused
      coarseProto = rcp( new Amesos2Smoother("Superlu", paramList) );
#else
      std::cout  << "AMESOS2 not available (try --coarseSolver=ifpack2)" << std::endl;
      return Teuchos::null; // TODO test for exception //EXIT_FAILURE;
#endif // HAVE_MUELU_AMESOS2
    } else if(coarseSolver=="ifpack2") {
#if defined(HAVE_MUELU_IFPACK2)
      if (rank == 0) std::cout << "CoarseGrid: IFPACK2" << std::endl;
      Teuchos::ParameterList ifpack2List;
      ifpack2List.set("fact: ilut level-of-fill",99); // TODO ??
      ifpack2List.set("fact: drop tolerance", 0);
      ifpack2List.set("fact: absolute threshold", 0);
      ifpack2List.set("fact: relative threshold", 0);
      coarseProto = rcp( new Ifpack2Smoother("ILUT",ifpack2List) );
#else
      std::cout  << "IFPACK2 not available (try --coarseSolver=amesos2)" << std::endl;
      //TODO        TEUCHOS_TEST_FOR_EXCEPTION
      return Teuchos::null;
#endif
    } else {
      std::cout  << "Unknow coarse grid solver (try  --coarseSolver=ifpack2 or --coarseSolver=amesos2)" << std::endl;
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

  for (int i=0; i<nSmoothers; i++)
    smootherList[i] = gimmeGaussSeidelProto(lib);

  return rcp (new MergedSmoother(smootherList));
  //verbose mode: return rcp (new MergedSmoother(smootherList, true));
}


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
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  // custom parameters
  int nSmoothers=2;
  LO maxLevels = 3;
  LO its=10;
  std::string coarseSolver="ifpack2";
  // std::string coarseSolver="amesos2";
  int pauseForDebugger=0;
  clp.setOption("nSmoothers",&nSmoothers,"number of Gauss-Seidel smoothers in the MergedSmootehrs");
  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed. If 1, then a MergedSmoother is used on the coarse grid");
  clp.setOption("its",&its,"number of multigrid cycles");
  clp.setOption("coarseSolver",&coarseSolver,"amesos2 or ifpack2 (Tpetra specific. Ignored for Epetra)");
  clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");

  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  matrixParameters.check();
  xpetraParameters.check();
  // TODO: check custom parameters

  if (comm->getRank() == 0) {
    // matrixParameters.print();
    // xpetraParameters.print();
    // TODO: print custom parameters
  }

  if (pauseForDebugger) {
    Utils::PauseForDebugger();
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Matrix> Op = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Matrix vs. CrsMatrixWrap

#ifdef NEUMANN
  // Tranform matrix to Neumann b.c.
  // Essentially, we need to update two diagonal elements
  bool useTpetra = (Op->getRowMap()->lib() == Xpetra::UseTpetra);

  std::vector<SC> newValuesVector(2, 0.0);
  Teuchos::ArrayView<SC> newValues(newValuesVector);

#ifdef HAVE_MUELU_TPETRA
  if (useTpetra)
    Utils::Op2NonConstTpetraCrs(Op)->resumeFill();
#endif
  Teuchos::ArrayView<const LO> indices;
  Teuchos::ArrayView<const SC> values;

  size_t myRank = Op->getRowMap()->getComm()->getRank();
  size_t nCpus  = Op->getRowMap()->getComm()->getSize();
  if (myRank == 0) {
    LO firstRow = 0;
    newValues[0] = 1.0; newValues[1] = -1.0;
    Op->getLocalRowView(0, indices, values);
#ifdef HAVE_MUELU_TPETRA
    if (useTpetra)
      Utils::Op2NonConstTpetraCrs(Op)->replaceLocalValues(0, indices, newValues);
#endif
#ifdef HAVE_MUELU_EPETRA
    if (!useTpetra)
      Utils::Op2NonConstEpetraCrs(Op)->ReplaceMyValues(0, 2, newValues.getRawPtr(), indices.getRawPtr());
#endif
  }
  if (myRank == nCpus-1) {
    LO lastRow = Op->getNodeNumRows()-1;
    newValues[0] = -1.0; newValues[1] = 1.0;
    Op->getLocalRowView(lastRow, indices, values);
#ifdef HAVE_MUELU_TPETRA
    if (useTpetra)
      Utils::Op2NonConstTpetraCrs(Op)->replaceLocalValues(lastRow, indices, newValues);
#endif
#ifdef HAVE_MUELU_EPETRA
    if (!useTpetra)
      Utils::Op2NonConstEpetraCrs(Op)->ReplaceMyValues(lastRow, 2, newValues.getRawPtr(), indices.getRawPtr());
#endif
  }
#ifdef HAVE_MUELU_TPETRA
  Op->fillComplete();
#endif
#endif

  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0)
    std::cout << "||NS|| = " << norms[0] << std::endl;

  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H = rcp( new Hierarchy() );
  H->SetDefaultVerbLevel(MueLu::Extreme);
  H->IsPreconditioner(false);

  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A", Op);
  Finest->Set("Nullspace", nullSpace);

  RCP<FactoryManager>   myFactManager = rcp(new FactoryManager());

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  RCP<SaPFactory>           SaPFact   = rcp(new SaPFactory());
  RCP<TentativePFactory>    TentPFact = rcp(new TentativePFactory(UCAggFact));
#if 1
  // Energy-minimization
#if 0
  RCP<PatternFactory>       PatternFact = rcp(new PatternFactory(TentPFact));
#else
  RCP<PatternFactory>       PatternFact = rcp(new PatternFactory(SaPFact));
#endif
  RCP<ConstraintFactory>    Cfact = rcp(new ConstraintFactory(PatternFact));
  RCP<EminPFactory>         Pfact = rcp(new EminPFactory(TentPFact, Cfact));

  myFactManager->SetFactory("Ppattern", PatternFact);
  myFactManager->SetFactory("Constraint", Cfact);
#else
  RCP<SaPFactory>       Pfact = SaPFact;
#endif

  myFactManager->SetFactory("P", Pfact);
  myFactManager->SetFactory("A", rcp(new RAPFactory()));
  myFactManager->SetFactory("R", rcp(new TransPFactory()));

#if 1
  RCP<NullspacePresmoothFactory> NPFact = rcp(new NullspacePresmoothFactory(rcp(new NullspaceFactory("Nullspace", TentPFact))));
  myFactManager->SetFactory("Nullspace", NPFact);
#endif

  RCP<SmootherPrototype> smooProto = gimmeMergedSmoother(nSmoothers, xpetraParameters.GetLib(), coarseSolver, comm->getRank());
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smooProto) );
  myFactManager->SetFactory("Smoother", SmooFact);

  Teuchos::ParameterList status;

  status = H->Setup(*myFactManager, 0, maxLevels);
  if (comm->getRank() == 0) {
    std::cout  << "======================\n Multigrid statistics \n======================" << std::endl;
    status.print(std::cout,Teuchos::ParameterList::PrintOptions().indent(2));
  }

  RCP<SmootherPrototype> coarseProto;
  if (maxLevels != 1) {
    coarseProto = gimmeCoarseProto(xpetraParameters.GetLib(), coarseSolver, comm->getRank());
  } else {
    coarseProto = gimmeMergedSmoother(nSmoothers, xpetraParameters.GetLib(), coarseSolver, comm->getRank());
  }
  if (coarseProto == Teuchos::null) return EXIT_FAILURE;
  SmootherFactory coarseSolveFact(coarseProto);
  H->SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  // Define RHS
  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  X->setSeed(846930886);
  X->randomize();
  X->norm2(norms);
  if (comm->getRank() == 0)
    std::cout << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  // Use AMG directly as an iterative method
  {
    X->putScalar( (SC) 0.0);

    H->Iterate(*RHS,its,*X);

#ifdef NEUMANN
    // Project X
    X->norm1(norms);
    size_t numElements = X->getGlobalLength();
    SC alpha = norms[0]/numElements;
    for (LO i = 0; i < numElements; i++)
      X->getDataNonConst(0)[i] -= alpha;
#endif

    X->norm2(norms);
    if (comm->getRank() == 0)
      std::cout << "||X_" << std::setprecision(2) << its << "|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
  }

  return EXIT_SUCCESS;

}

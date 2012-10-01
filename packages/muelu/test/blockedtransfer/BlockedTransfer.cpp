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
/*
 * BlockedTransfer.cpp
 *
 *  Created on: 01.01.2012
 *      Author: tobias
 */

#include <unistd.h>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_GaussSeidelSmoother.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"

// Belos
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosXpetraAdapter.hpp" // this header defines Belos::XpetraOp()
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuOp()
#endif

//
typedef double Scalar;
typedef int    LocalOrdinal;
// FIXME
// #ifdef HAVE_TEUCHOS_LONG_LONG_INT
// typedef long long int GlobalOrdinal;
// #else
typedef int GlobalOrdinal;
// #endif

typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

#include "MueLu_UseShortNames.hpp"

/////////////////////////
// helper function

Teuchos::RCP<CrsMatrixWrap> GenerateProblemMatrix(const Teuchos::RCP<const Map> map, Scalar a = 2.0, Scalar b = -1.0, Scalar c = -1.0) {

  Teuchos::RCP<CrsMatrixWrap> mtx = Galeri::Xpetra::MatrixTraits<Map,CrsMatrixWrap>::Build(map, 3);

  LocalOrdinal NumMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
  GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
  GlobalOrdinal nIndexBase = map->getIndexBase();

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz=2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  for (LocalOrdinal i = 0; i < NumMyElements; ++i)
  {
    if (MyGlobalElements[i] == nIndexBase)
    {
      // off-diagonal for first row
      Indices[0] = nIndexBase;
      NumEntries = 1;
      Values[0] = c;
    }
    else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1)
    {
      // off-diagonal for last row
      Indices[0] = nIndexBase + NumGlobalElements - 2;
      NumEntries = 1;
      Values[0] = b;
    }
    else
    {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1] = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0] = c;
      NumEntries = 2;
    }

    // put the off-diagonal entries
    // Xpetra wants ArrayViews (sigh)
    Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
    Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
    mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

    // Put in the diagonal entry
    mtx->insertGlobalValues(MyGlobalElements[i],
        Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
        Teuchos::tuple<Scalar>(a) );

  } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)


  mtx->fillComplete(map,map);

  return mtx;
}

/////////////////
// MAIN

int main(int argc, char *argv[]) {
  using Teuchos::RCP; using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  //#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
  //#endif

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  //RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));

  xpetraParameters.check();

  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  if (comm->getRank() == 0) {
    std::cout << xpetraParameters;
    // TODO: print custom parameters // Or use paramList::print()!
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  RCP<const Map> bigMap;
  RCP<const Map> map1;
  RCP<const Map> map2;
  GO numElements = 500;
  GO numElements1 = 400;
  GO numElements2 = 100;

  //bigMap = MapFactory::Build(Xpetra::UseEpetra, numElements,  0, comm); // ok this is the problem :-)
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);
  map1   = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm, -1);
  map2   = StridedMapFactory::Build(lib, numElements2, numElements1, stridingInfo, comm, -1);

  std::vector<GlobalOrdinal> localGids; // vector with all local GIDs on cur proc
  Teuchos::ArrayView< const GlobalOrdinal > map1eleList = map1->getNodeElementList(); // append all local gids from map1 and map2
  localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
  Teuchos::ArrayView< const GlobalOrdinal > map2eleList = map2->getNodeElementList();
  localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
  Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0],localGids.size());
  bigMap = MapFactory::Build(lib, numElements, eleList, 0, comm); // create full big map (concatenation of map1 and map2)

  std::vector<Teuchos::RCP<const Map> > maps;
  maps.push_back(map1); maps.push_back(map2);

  Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar,LO,GO,Node>::Build(bigMap, maps);

  RCP<CrsMatrixWrap> Op11 = GenerateProblemMatrix(map1,2,-1,-1);
  RCP<CrsMatrixWrap> Op22 = GenerateProblemMatrix(map2,3,-2,-1);

  /*Op11->describe(*out,Teuchos::VERB_EXTREME);
  Op22->describe(*out,Teuchos::VERB_EXTREME);*/

  // build blocked operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(mapExtractor,mapExtractor,10));

  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat11 = Op11->getCrsMatrix();
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat22 = Op22->getCrsMatrix();
  bOp->setMatrix(0,0,crsMat11);
  bOp->setMatrix(1,1,crsMat22);

  bOp->fillComplete();

  // build hierarchy
  Hierarchy H;
  RCP<Level> levelOne = H.GetLevel();
  levelOne->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp)); // set blocked operator

  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

  RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());
  RCP<TransPFactory> R11Fact = rcp(new TransPFactory());

  RCP<TentativePFactory> P22TentFact = rcp(new TentativePFactory());
  RCP<PgPFactory> P22Fact = rcp(new PgPFactory());
  RCP<GenericRFactory> R22Fact = rcp(new GenericRFactory());

  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 5);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  RCP<SmootherPrototype> smoProto11     = rcp( new TrilinosSmoother(ifpackType, ifpackList, 0, A11Fact) );
  RCP<SmootherPrototype> smoProto22     = rcp( new TrilinosSmoother(ifpackType, ifpackList, 0, A22Fact) );

  //RCP<SmootherPrototype> smoProto11     = rcp( new DirectSolver("", Teuchos::ParameterList(), A11Fact) );
  //RCP<SmootherPrototype> smoProto22     = rcp( new DirectSolver("", Teuchos::ParameterList(), A22Fact) );

  RCP<SmootherFactory> Smoo11Fact = rcp( new SmootherFactory(smoProto11) );
  RCP<SmootherFactory> Smoo22Fact = rcp( new SmootherFactory(smoProto22) );

  RCP<FactoryManager> M11 = rcp(new FactoryManager());
  M11->SetFactory("A", A11Fact);
  M11->SetFactory("P", P11Fact);
  M11->SetFactory("Ptent", P11Fact); //for Nullspace
  M11->SetFactory("R", R11Fact);
  M11->SetFactory("Smoother", Smoo11Fact);
  M11->SetIgnoreUserData(true);

  RCP<FactoryManager> M22 = rcp(new FactoryManager());
  M22->SetFactory("A", A22Fact);
  M22->SetFactory("P", P22Fact);
  M22->SetFactory("R", R22Fact);
  M22->SetFactory("Ptent", P22TentFact); //for both P22 and Nullspace
  M22->SetFactory("Smoother", Smoo22Fact);
  M22->SetIgnoreUserData(true);

  RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory(Teuchos::null/*AFact necessary for row map index base*/));
  PFact->AddFactoryManager(M11);
  PFact->AddFactoryManager(M22);

  RCP<GenericRFactory> RFact = rcp(new GenericRFactory());

  RCP<RAPFactory> AcFact = rcp(new RAPFactory());

  // Smoothers
  //RCP<SmootherPrototype> smootherPrototype     = rcp( new GaussSeidelSmoother(1, 1.0) );
  RCP<BlockedGaussSeidelSmoother> smootherPrototype     = rcp( new BlockedGaussSeidelSmoother(2,1.0) );
  smootherPrototype->AddFactoryManager(M11);
  smootherPrototype->AddFactoryManager(M22);
  RCP<SmootherFactory>   smootherFact          = rcp( new SmootherFactory(smootherPrototype) );

  // Coarse grid correction
  //RCP<SmootherPrototype> coarseSolverPrototype = rcp( new DirectSolver() );
  RCP<BlockedGaussSeidelSmoother> coarseSolverPrototype = rcp( new BlockedGaussSeidelSmoother() );
  coarseSolverPrototype->AddFactoryManager(M11);
  coarseSolverPrototype->AddFactoryManager(M22);
  RCP<SmootherFactory>   coarseSolverFact      = rcp( new SmootherFactory(coarseSolverPrototype, Teuchos::null) );

  // main factory manager
  FactoryManager M;
  M.SetFactory("A",            AcFact);
  M.SetFactory("P",            PFact);
  M.SetFactory("R",            RFact);
  M.SetFactory("Smoother",     smootherFact); // TODO fix me
  M.SetFactory("CoarseSolver", coarseSolverFact);

  H.Setup(M);

  RCP<Level> l0 = H.GetLevel(0);
  RCP<Level> l1 = H.GetLevel(1);
  RCP<Level> l2 = H.GetLevel(2);

  l0->print(*out,Teuchos::VERB_EXTREME);
  l1->print(*out,Teuchos::VERB_EXTREME);
  l2->print(*out,Teuchos::VERB_EXTREME);

  // Define B
  RCP<Vector> X = VectorFactory::Build(bigMap,1);
  RCP<Vector> B = VectorFactory::Build(bigMap,1);
  X->setSeed(846930886);
  X->randomize();
  bOp->apply(*X, *B, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

  // X = 0
  X->putScalar((SC) 0.0);

  LO nIts = 9;
  H.Iterate(*B, nIts, *X); // we have no working smoother for blocked operators


  return EXIT_SUCCESS;
}

// TODO: add warning if:
// DEBUG_MODE, LONG_LONG or KLU



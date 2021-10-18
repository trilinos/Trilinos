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
/*
 * RebalanceAcFactory.cpp
 *
 *  Created on: 20.09.2011
 *      Author: tobias
 */

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_PgPFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_RepartitionHeuristicFactory.hpp>
//#include <MueLu_SaPFactory.hpp>
//#include <MueLu_RebalanceAcFactory.hpp>
//#include <MueLu_TrilinosSmoother.hpp>
//#include <MueLu_CoupledAggregationFactory.hpp>
//#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_RebalanceAcFactory.hpp>
#include <MueLu_RepartitionInterface.hpp>
#include <MueLu_IsorropiaInterface.hpp>
#include <MueLu_RebalanceTransferFactory.hpp>
#include <MueLu_RepartitionFactory.hpp>

namespace MueLuTests {

  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RebalanceAcFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    #include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    auto rAcFactory = rcp(new RebalanceAcFactory());
    TEST_ASSERT(!rAcFactory.is_null());
    TEST_EQUALITY(rAcFactory->NumRebalanceFactories() == 0, true);

    // Add Factory
    rAcFactory->AddRebalanceFactory(rcp(new RebalanceTransferFactory()));
    rAcFactory->AddRebalanceFactory(rcp(new RebalanceTransferFactory()));

    TEST_EQUALITY(rAcFactory->NumRebalanceFactories() == 2, true);

    auto paramList = rAcFactory->GetValidParameterList();
    TEST_ASSERT(!paramList.is_null());

  }

  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RebalanceAcFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    #include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;


    // create coarsest smoother
    RCP<SmootherPrototype> coarsestSmooProto;
    std::string type = "";
    Teuchos::ParameterList coarsestSmooList;
#if defined(HAVE_AMESOS_SUPERLU)
    coarsestSmooProto = Teuchos::rcp( new DirectSolver("Superlu", coarsestSmooList) );
#else
    coarsestSmooProto = Teuchos::rcp( new DirectSolver("Klu", coarsestSmooList) );
#endif
    RCP<SmootherFactory> coarsestSmooFact = rcp(new SmootherFactory(coarsestSmooProto, Teuchos::null));


    // Configure FactoryManager
    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("CoarseSolver", coarsestSmooFact);

    RCP<Hierarchy> H = rcp ( new Hierarchy() );
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    GO maxCoarseSize=1;
    H->SetMaxCoarseSize(maxCoarseSize);

    Teuchos::ParameterList matrixParams;
    matrixParams.set("matrixType","Laplace2D");
    double factor=10;
    matrixParams.set("stretchx",1.0);
    matrixParams.set("stretchy",factor);
//    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixParams,TestHelpers::Parameters::getLib());

    Teuchos::CommandLineProcessor clp(false); // Note:

    // - Levels
    LO  optMaxLevels     = 4;              clp.setOption("maxLevels",      &optMaxLevels,          "maximum number of levels allowed");
    int optMaxCoarseSize = 1;              clp.setOption("maxCoarseSize",  &optMaxCoarseSize,      "maximum #dofs in coarse operator"); //FIXME clp doesn't like long long int

    // - Repartitioning
    #if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
    int    optRepartition    = 1;             clp.setOption("repartition",    &optRepartition,        "enable repartitioning");
    LO     optMinRowsPerProc = 50;            clp.setOption("minRowsPerProc", &optMinRowsPerProc,     "min #rows allowable per proc before repartitioning occurs");
    double optNnzImbalance   = 1.2;           clp.setOption("nnzImbalance",   &optNnzImbalance,       "max allowable nonzero imbalance before repartitioning occurs");
    #else
    int optRepartition = 0;
    #endif

    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 256); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);      // manage parameters of xpetra

    // TUTORIALSPLIT ===========================================================
    RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<Matrix>  A = Pr->BuildMatrix();

    // build default null space
    LocalOrdinal numPDEs = 1;
    if(A->IsView("stridedMaps")==true) {
      Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
      numPDEs = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
      oldView = A->SwitchToView(oldView);
    }
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getDomainMap(), numPDEs);
    nullspace->putScalar( (SC) 1.0);
    Level *level = H->GetLevel().get();
    level->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    level->Set("A", A);
    level->Set("Nullspace",   nullspace);




    // Add Factory
    // build non-rebalanced transfer operators
    RCP<PgPFactory> Pfact = rcp( new PgPFactory() );
    RCP<Factory> Rfact  = rcp( new GenericRFactory());
    //RCP<SaPFactory> Pfact  = rcp( new SaPFactory() );
    //RCP<Factory>   Rfact  = rcp( new TransPFactory() );
    RCP<RAPFactory> Acfact = rcp( new RAPFactory() );
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);
    Rfact->SetFactory("P", Pfact);
    Acfact->SetFactory("P", Pfact);
    Acfact->SetFactory("R", Rfact);

    // define rebalancing factory for coarse block matrix A(1,1)
    RCP<AmalgamationFactory> rebAmalgFact = rcp(new AmalgamationFactory());
    rebAmalgFact->SetFactory("A", Acfact);

    // Repartitioning (decides how many partitions are built)
    RCP<Factory> RepartitionHeuristicFact = rcp(new RepartitionHeuristicFactory());
    {
      Teuchos::ParameterList paramList;
      paramList.set("repartition: min rows per proc", optMinRowsPerProc);
      paramList.set("repartition: max imbalance", optNnzImbalance);
      RepartitionHeuristicFact->SetParameterList(paramList);
    }
    RepartitionHeuristicFact->SetFactory("A", Acfact);

    // create amalgamated "Partition"
    RCP<MueLu::IsorropiaInterface<LO, GO, NO> > isoInterface = rcp(new MueLu::IsorropiaInterface<LO, GO, NO>());
    isoInterface->SetFactory("A", Acfact);
    isoInterface->SetFactory("number of partitions", RepartitionHeuristicFact);
    isoInterface->SetFactory("UnAmalgamationInfo", rebAmalgFact);


    // create "Partition" by unamalgamtion
    RCP<MueLu::RepartitionInterface<LO, GO, NO> > repInterface = rcp(new MueLu::RepartitionInterface<LO, GO, NO>());
    repInterface->SetFactory("A", Acfact);
    repInterface->SetFactory("number of partitions", RepartitionHeuristicFact);
    repInterface->SetFactory("AmalgamatedPartition", isoInterface);

    // Repartitioning (creates "Importer" from "Partition")
    RCP<Factory> RepartitionFact = rcp(new RepartitionFactory());
    RepartitionFact->SetFactory("A", Acfact);
    RepartitionFact->SetFactory("number of partitions", RepartitionHeuristicFact);
    RepartitionFact->SetFactory("Partition", repInterface);


    // Reordering of the transfer operators
    RCP<Factory> RebalancedPFact = rcp(new RebalanceTransferFactory());
    RebalancedPFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Interpolation")));
    RebalancedPFact->SetFactory("P", Pfact);
    RebalancedPFact->SetFactory("Nullspace", M.GetFactory("Ptent")); // TODO

    RCP<Factory> RebalancedRFact = rcp(new RebalanceTransferFactory());
    RebalancedRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
    RebalancedRFact->SetFactory("R", Rfact);

    auto RebalancedAFact = rcp(new RebalanceAcFactory());
    RebalancedAFact->SetFactory("A", Acfact);

//    RebalancedAFact->AddRebalanceFactory(RebalancedPFact);
//    RebalancedAFact->AddRebalanceFactory(RebalancedRFact);


    M.SetFactory("A", RebalancedAFact);
    M.SetFactory("P", RebalancedPFact);
    M.SetFactory("R", RebalancedRFact);
    M.SetFactory("Nullspace",   RebalancedPFact);
//    M.SetFactory("Coordinates", RebalancedPFact);
    M.SetFactory("Importer",    RepartitionFact);

    int startLevel = 0;
    H->Setup(M, startLevel, optMaxLevels);

    RebalancedAFact->Build(*level, *level);


  }

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RebalanceAcFactory, Constructor, Scalar, LO, GO, Node) \
//  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RebalanceAcFactory, Build, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>


}//namespace MueLuTests



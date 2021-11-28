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
#include <MueLu_SaPFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
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
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"

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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RebalanceAcFactory, BuildWithImporter, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    #include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    if (TestHelpers::Parameters::getLib() == Xpetra::UseEpetra) {
      out << "skipping test for linAlgebra==UseEpetra" << std::endl;
      return;
    }

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Teuchos::CommandLineProcessor clp(false);
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case

    Level aLevel;
    Level corseLevel;
    RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

    TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(aLevel, corseLevel);
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2); //can be an empty operator
    corseLevel.Set("A", A);
    RCP<const Import> importer = ImportFactory::Build(A->getRowMap(), map);
    corseLevel.Set("Importer", importer);

    RCP<RebalanceAcFactory> RebalancedAFact = rcp(new RebalanceAcFactory());
    RebalancedAFact->SetDefaultVerbLevel(MueLu::Extreme);
    RebalancedAFact->Build(aLevel, corseLevel);

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RebalanceAcFactory, BuildWithoutImporter, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    #include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    Level aLevel;
    Level corseLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(aLevel, corseLevel);
    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2); //can be an empty operator
    corseLevel.Set("A", A);
    RCP<const Import> importer  = Teuchos::null;
    corseLevel.Set("Importer", importer);

    RCP<RebalanceAcFactory> RebalancedAFact = rcp(new RebalanceAcFactory());
    RebalancedAFact->SetDefaultVerbLevel(MueLu::Extreme);
    RebalancedAFact->Build(aLevel, corseLevel);

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RebalanceAcFactory, BuildWithCompleteFactoryManager, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    #include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;

    // - Smoothed-Aggregation
    int optMaxCoarseSize = 50;
    LO optMinRowsPerProc = 2000;
    double optNnzImbalance = 1.2;
    Teuchos::CommandLineProcessor clp(false);
    LO  optMaxLevels     = 10;
    // Default is Laplace1D with nx = 8748.
    // It's a nice size for 1D and perfect aggregation. (6561 = 3^8)
    //Nice size for 1D and perfect aggregation on small numbers of processors. (8748 = 4*3^7)
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                   // manage parameters of xpetra

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

    RCP<MultiVector> coordinates;

    // USER GUIDE // define communicator
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<Matrix>  A = Pr->BuildMatrix();
    int nDofsPerNode = 3;
    A->SetFixedBlockSize(nDofsPerNode);   // 2 velocity dofs and 1 pressure dof per node.


    // Laplace 1D
    coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, MultiVector>("1D", map, matrixParameters.GetParameterList());

    // USER GUIDE   // define near null space
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getDomainMap(), 1);
    nullspace->putScalar( (SC) 1.0);
    // USER GUIDE   //

    // USER GUIDE   // create new hierarchy
    RCP<MueLu::Hierarchy<SC, LO, GO, NO> > H;
    // USER GUIDE   //


    // USER GUIDE     // instantiate new Hierarchy object
    H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    H->SetMaxCoarseSize((GO) optMaxCoarseSize);

    // USER GUIDE     // create a fine level object
    RCP<Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A",           A);
    Finest->Set("Nullspace",   nullspace);
    Finest->Set("Coordinates", coordinates); //FIXME: XCoordinates, YCoordinate

    RCP<FactoryManager> M = rcp(new FactoryManager());
    M->SetKokkosRefactor(false);
    M->SetFactory("CoarseSolver", coarsestSmooFact);

    // USER GUIDE       // declare some factories (potentially overwrite default factories)
    RCP<PgPFactory> PFact = rcp(new PgPFactory());

    RCP<Factory>    RFact = rcp(new GenericRFactory());

    RCP<RAPFactory> AFact = rcp(new RAPFactory());
    AFact->setVerbLevel(Teuchos::VERB_HIGH);

#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN)
    // Repartitioning

    // The Factory Manager will be configured to return the rebalanced versions of P, R, A by default.
    // Everytime we want to use the non-rebalanced versions, we need to explicitly define the generating factory.
    RFact->SetFactory("P", PFact);
    AFact->SetFactory("P", PFact);
    AFact->SetFactory("R", RFact);

    // Transfer coordinates
    RCP<CoordinatesTransferFactory> TransferCoordinatesFact = rcp(new CoordinatesTransferFactory());
    AFact->AddTransferFactory(TransferCoordinatesFact); // FIXME REMOVE

    // define rebalancing factory for coarse block matrix A(1,1)
    RCP<AmalgamationFactory> rebAmalgFact = rcp(new AmalgamationFactory());
    rebAmalgFact->SetFactory("A", AFact);

    // Repartitioning (decides how many partitions are built)
    RCP<Factory> RepartitionHeuristicFact = rcp(new RepartitionHeuristicFactory());
    {
      Teuchos::ParameterList paramList;
      paramList.set("repartition: min rows per proc", optMinRowsPerProc);
      paramList.set("repartition: max imbalance", optNnzImbalance);
      RepartitionHeuristicFact->SetParameterList(paramList);
    }
    RepartitionHeuristicFact->SetFactory("A", AFact);

    // create amalgamated "Partition"
    RCP<MueLu::IsorropiaInterface<LO, GO, NO> > isoInterface = rcp(new MueLu::IsorropiaInterface<LO, GO, NO>());
    isoInterface->SetFactory("A", AFact);
    isoInterface->SetFactory("number of partitions", RepartitionHeuristicFact);
    isoInterface->SetFactory("UnAmalgamationInfo", rebAmalgFact);

    // create "Partition" by unamalgamtion
    RCP<MueLu::RepartitionInterface<LO, GO, NO> > repInterface = rcp(new MueLu::RepartitionInterface<LO, GO, NO>());
    repInterface->SetFactory("A", AFact);
    repInterface->SetFactory("number of partitions", RepartitionHeuristicFact);
    repInterface->SetFactory("AmalgamatedPartition", isoInterface);

    // Repartitioning (creates "Importer" from "Partition")
    RCP<Factory> RepartitionFact = rcp(new RepartitionFactory());
    RepartitionFact->SetFactory("A", AFact);
    RepartitionFact->SetFactory("number of partitions", RepartitionHeuristicFact);
    RepartitionFact->SetFactory("Partition", repInterface);

    // Reordering of the transfer operators
    RCP<Factory> RebalancedPFact = rcp(new RebalanceTransferFactory());
    RebalancedPFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Interpolation")));
    RebalancedPFact->SetFactory("P", PFact);
    RebalancedPFact->SetFactory("Coordinates", TransferCoordinatesFact);
    RebalancedPFact->SetFactory("Nullspace", M->GetFactory("Ptent")); // TODO

    RCP<Factory> RebalancedRFact = rcp(new RebalanceTransferFactory());
    RebalancedRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
    RebalancedRFact->SetFactory("R", RFact);

    // Compute Ac from rebalanced P and R
    RCP<RebalanceAcFactory> RebalancedAFact = rcp(new RebalanceAcFactory());
    RebalancedAFact->SetFactory("A", AFact);
    // Add Factory
    RebalancedAFact->AddRebalanceFactory(rcp(new RebalanceTransferFactory()));
    RebalancedAFact->AddRebalanceFactory(rcp(new RebalanceTransferFactory()));
    RebalancedAFact->SetDefaultVerbLevel(MueLu::Extreme);

    // Configure FactoryManager
    M->SetFactory("P", RebalancedPFact);
    M->SetFactory("R", RebalancedRFact);
    M->SetFactory("A", RebalancedAFact);
    M->SetFactory("Nullspace",   RebalancedPFact);
    M->SetFactory("Coordinates", RebalancedPFact);
    M->SetFactory("Importer",    RepartitionFact);

  #else
          TEUCHOS_TEST_FOR_EXCEPT(true);
  #endif
    //
    // Setup preconditioner
    //

    // USER GUIDE     // setup multigrid hierarchy
    int startLevel = 0;
    H->Setup(*M.get(), startLevel, optMaxLevels);
    RCP<Level> l0 = H->GetLevel(0);
    l0->SetFactoryManager(M);
    TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);

  }

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RebalanceAcFactory, Constructor, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RebalanceAcFactory, BuildWithoutImporter, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RebalanceAcFactory, BuildWithImporter, Scalar, LO, GO, Node) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RebalanceAcFactory, BuildWithCompleteFactoryManager, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>


}//namespace MueLuTests



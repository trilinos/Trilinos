// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_StridedMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_BraessSarazinSmoother.hpp"

#include "MueLu_CoarseMapFactory.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"

#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_RebalanceBlockAcFactory.hpp"
#include "MueLu_RebalanceBlockInterpolationFactory.hpp"
#include "MueLu_RebalanceBlockRestrictionFactory.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_CloneRepartitionInterface.hpp"
#endif

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include "Navier2D_Helpers.h"

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 *
 *  3 level multigrid with Braess-Sarazin smoother
 *  Reuse aggregates of block 0 for aggregates in block 1
 *
 */

int main(int argc, char* argv[]) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using namespace MueLuTests;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

    // Timing
    Teuchos::Time myTime("global");
    Teuchos::TimeMonitor MM(myTime);

    // read in some command line parameters
    Teuchos::CommandLineProcessor clp(false);

    int rebalanceBlocks = 1;
    clp.setOption("rebalanceBlocks", &rebalanceBlocks, "rebalance blocks (1=yes, else=no)");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
#ifndef HAVE_XPETRA_INT_LONG_LONG
    *out << "Warning: scaling test was not compiled with long long int support" << std::endl;

    // custom parameters
    LocalOrdinal maxLevels = 3;

    GlobalOrdinal maxCoarseSize = 1;  // FIXME clp doesn't like long long int

    int globalNumDofs = 1500;  // used for the maps
    int nDofsPerNode  = 3;     // used for generating the fine level null-space

    // build strided maps
    // striding information: 2 velocity dofs and 1 pressure dof = 3 dofs per node
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(2);
    stridingInfo.push_back(1);

    /////////////////////////////////////// build strided maps
    // build strided maps:
    // xstridedfullmap: full map (velocity and pressure dof gids), continous
    // xstridedvelmap: only velocity dof gid maps (i.e. 0,1,3,4,6,7...)
    // xstridedpremap: only pressure dof gid maps (i.e. 2,5,8,...)
    Xpetra::UnderlyingLib lib             = Xpetra::UseEpetra;
    RCP<const StridedMap> xstridedfullmap = StridedMapFactory::Build(lib, globalNumDofs, 0, stridingInfo, comm, -1);
    RCP<const StridedMap> xstridedvelmap  = StridedMapFactory::Build(xstridedfullmap, 0);
    RCP<const StridedMap> xstridedpremap  = StridedMapFactory::Build(xstridedfullmap, 1);

    /////////////////////////////////////// transform Xpetra::Map objects to Epetra
    // this is needed for AztecOO
    const RCP<const Epetra_Map> fullmap = rcpFromRef(Xpetra::toEpetra(*xstridedfullmap));
    RCP<const Epetra_Map> velmap        = rcpFromRef(Xpetra::toEpetra(*xstridedvelmap));
    RCP<const Epetra_Map> premap        = rcpFromRef(Xpetra::toEpetra(*xstridedpremap));

    /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

    // read in problem
    Epetra_CrsMatrix* ptrA    = 0;
    Epetra_Vector* ptrf       = 0;
    Epetra_MultiVector* ptrNS = 0;

    *out << "Reading matrix market file" << std::endl;

    EpetraExt::MatrixMarketFileToCrsMatrix("A_re1000_5932.txt", *fullmap, *fullmap, *fullmap, ptrA);
    EpetraExt::MatrixMarketFileToVector("b_re1000_5932.txt", *fullmap, ptrf);
    // EpetraExt::MatrixMarketFileToCrsMatrix("/home/tobias/promotion/trilinos/fc17-dyn/packages/muelu/test/navierstokes/A_re1000_5932.txt",*fullmap,*fullmap,*fullmap,ptrA);
    // EpetraExt::MatrixMarketFileToVector("/home/tobias/promotion/trilinos/fc17-dyn/packages/muelu/test/navierstokes/b_re1000_5932.txt",*fullmap,ptrf);

    RCP<Epetra_CrsMatrix> epA    = Teuchos::rcp(ptrA);
    RCP<Epetra_Vector> epv       = Teuchos::rcp(ptrf);
    RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);

    /////////////////////////////////////// split system into 2x2 block system

    *out << "Split matrix into 2x2 block matrix" << std::endl;

    // split fullA into A11,..., A22
    Teuchos::RCP<Epetra_CrsMatrix> A11;
    Teuchos::RCP<Epetra_CrsMatrix> A12;
    Teuchos::RCP<Epetra_CrsMatrix> A21;
    Teuchos::RCP<Epetra_CrsMatrix> A22;

    if (SplitMatrix2x2(epA, *velmap, *premap, A11, A12, A21, A22) == false)
      *out << "Problem with splitting matrix" << std::endl;

    /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(A11));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(A12));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(A21));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(A22));

    /////////////////////////////////////// generate MapExtractor object

    std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > xmaps;
    xmaps.push_back(xstridedvelmap);
    xmaps.push_back(xstridedpremap);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xstridedfullmap, xmaps);

    /////////////////////////////////////// build blocked transfer operator
    // using the map extractor
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map_extractor, map_extractor, 10));
    bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA11)));
    bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA12)));
    bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA21)));
    bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA22)));

    bOp->fillComplete();

    //////////////////////////////////////////////////// create Hierarchy
    RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    // H->setDefaultVerbLevel(Teuchos::VERB_NONE);
    H->SetMaxCoarseSize(maxCoarseSize);

    //////////////////////////////////////////////////////// finest Level
    RCP<MueLu::Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));

    ////////////////////////////////////////// prepare null space for A11
    RCP<MultiVector> nullspace11 = MultiVectorFactory::Build(xstridedvelmap, 2);  // this is a 2D standard null space

    for (int i = 0; i < nDofsPerNode - 1; ++i) {
      Teuchos::ArrayRCP<Scalar> nsValues = nullspace11->getDataNonConst(i);
      int numBlocks                      = nsValues.size() / (nDofsPerNode - 1);
      for (int j = 0; j < numBlocks; ++j) {
        nsValues[j * (nDofsPerNode - 1) + i] = 1.0;
      }
    }

    Finest->Set("Nullspace1", nullspace11);

    ////////////////////////////////////////// prepare null space for A22
    RCP<MultiVector> nullspace22         = MultiVectorFactory::Build(xstridedpremap, 1);  // this is a 2D standard null space
    Teuchos::ArrayRCP<Scalar> nsValues22 = nullspace22->getDataNonConst(0);
    for (int j = 0; j < nsValues22.size(); ++j) {
      nsValues22[j] = 1.0;
    }

    Finest->Set("Nullspace2", nullspace22);

    /////////////////////////////////////////// define rebalanced block AC factory
    // This is the main factory for "A" and defines the input for
    //   - the SubBlockAFactory objects
    //   - the rebalanced block Ac factory
    RCP<RebalanceBlockAcFactory> RebalancedAcFact = rcp(new RebalanceBlockAcFactory());

    /////////////////////////////////////////// define non-rebalanced blocked transfer ops
    RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());  // use row map index base from bOp
    RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
    RFact->SetFactory("P", PFact);

    // non-rebalanced block coarse matrix factory
    // output is non-rebalanced coarse block matrix Ac
    // used as input for rebalanced block coarse factory RebalancedAcFact
    RCP<Factory> AcFact = rcp(new BlockedRAPFactory());
    AcFact->SetFactory("A", MueLu::NoFactory::getRCP());
    AcFact->SetFactory("P", PFact);  // use non-rebalanced block prolongator as input
    AcFact->SetFactory("R", RFact);  // use non-rebalanced block restrictor as input

    // Repartitioning (decides how many partitions are built)
    RCP<Factory> RepartitionHeuristicFact = rcp(new RepartitionHeuristicFactory());
    {
      Teuchos::ParameterList paramList;
      paramList.set("repartition: min rows per proc", 200);
      paramList.set("repartition: max imbalance", 1.3);
      if (rebalanceBlocks == 1)
        paramList.set("repartition: start level", 1);
      else
        paramList.set("repartition: start level", 10);  // supress rebalancing
      RepartitionHeuristicFact->SetParameterList(paramList);
    }
    RepartitionHeuristicFact->SetFactory("A", AcFact);

    // define matrix sub-blocks of possibly rebalanced block matrix A
    // These are used as input for
    //   - the sub blocks of the transfer operators
    RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
    A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
    A11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
    /*A11Fact->SetParameter("Range map: Striding info",Teuchos::ParameterEntry(std::string("{ 2 1 }")));
    A11Fact->SetParameter("Range map: Strided block id",Teuchos::ParameterEntry(0));
    A11Fact->SetParameter("Domain map: Striding info",Teuchos::ParameterEntry(std::string("{ 2 1 }")));
    A11Fact->SetParameter("Domain map: Strided block id",Teuchos::ParameterEntry(0));*/

    RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
    A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
    A22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));
    /*A22Fact->SetParameter("Range map: Striding info",Teuchos::ParameterEntry(std::string("{ 2 1 }")));
    A22Fact->SetParameter("Range map: Strided block id",Teuchos::ParameterEntry(1));
    A22Fact->SetParameter("Domain map: Striding info",Teuchos::ParameterEntry(std::string("{ 2 1 }")));
    A22Fact->SetParameter("Domain map: Strided block id",Teuchos::ParameterEntry(1));*/

    /////////////////////////////////////////// define rebalancing factories
    // define sub blocks of the coarse non-rebalanced block matrix Ac
    // input is the block operator generated by AcFact
    RCP<SubBlockAFactory> rebA11Fact = Teuchos::rcp(new SubBlockAFactory());
    rebA11Fact->SetFactory("A", AcFact);
    rebA11Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
    rebA11Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
    rebA11Fact->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(std::string("{ 2, 1 }")));
    rebA11Fact->SetParameter("Range map: Strided block id", Teuchos::ParameterEntry(0));
    rebA11Fact->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(std::string("{ 2, 1 }")));
    rebA11Fact->SetParameter("Domain map: Strided block id", Teuchos::ParameterEntry(0));

    RCP<SubBlockAFactory> rebA22Fact = Teuchos::rcp(new SubBlockAFactory());
    rebA22Fact->SetFactory("A", AcFact);
    rebA22Fact->SetParameter("block row", Teuchos::ParameterEntry(1));
    rebA22Fact->SetParameter("block col", Teuchos::ParameterEntry(1));
    rebA22Fact->SetParameter("Range map: Striding info", Teuchos::ParameterEntry(std::string("{ 2, 1 }")));
    rebA22Fact->SetParameter("Range map: Strided block id", Teuchos::ParameterEntry(1));
    rebA22Fact->SetParameter("Domain map: Striding info", Teuchos::ParameterEntry(std::string("{ 2, 1 }")));
    rebA22Fact->SetParameter("Domain map: Strided block id", Teuchos::ParameterEntry(1));

    // define rebalancing factory for coarse block matrix A(1,1)
    RCP<AmalgamationFactory> rebAmalgFact11 = rcp(new AmalgamationFactory());
    rebAmalgFact11->SetFactory("A", rebA11Fact);
    rebAmalgFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);

    RCP<MueLu::IsorropiaInterface<LO, GO, NO> > isoInterface1 = rcp(new MueLu::IsorropiaInterface<LO, GO, NO>());
    isoInterface1->SetFactory("A", rebA11Fact);
    isoInterface1->SetFactory("number of partitions", RepartitionHeuristicFact);
    isoInterface1->SetFactory("UnAmalgamationInfo", rebAmalgFact11);

    RCP<MueLu::RepartitionInterface<LO, GO, NO> > repInterface1 = rcp(new MueLu::RepartitionInterface<LO, GO, NO>());
    repInterface1->SetFactory("A", rebA11Fact);
    repInterface1->SetFactory("number of partitions", RepartitionHeuristicFact);
    repInterface1->SetFactory("AmalgamatedPartition", isoInterface1);

    // Repartitioning (creates "Importer" from "Partition")
    RCP<Factory> RepartitionFact = rcp(new RepartitionFactory());
    RepartitionFact->SetFactory("A", rebA11Fact);
    RepartitionFact->SetFactory("number of partitions", RepartitionHeuristicFact);
    RepartitionFact->SetFactory("Partition", repInterface1);
    RepartitionFact->SetParameter("repartition: print partition distribution", Teuchos::ParameterEntry(true));
    RepartitionFact->SetParameter("repartition: remap parts", Teuchos::ParameterEntry(true));

    // define rebalancing factory for coarse block matrix A(1,1)
    RCP<AmalgamationFactory> rebAmalgFact22 = rcp(new AmalgamationFactory());
    rebAmalgFact22->SetFactory("A", rebA22Fact);
    rebAmalgFact22->setDefaultVerbLevel(Teuchos::VERB_EXTREME);

    RCP<MueLu::CloneRepartitionInterface<SC, LO, GO, NO> > repInterface2 = rcp(new MueLu::CloneRepartitionInterface<SC, LO, GO, NO>());
    repInterface2->SetFactory("A", rebA22Fact);
    repInterface2->SetFactory("number of partitions", RepartitionHeuristicFact);
    repInterface2->SetFactory("Partition", repInterface1);

    // second repartition factory
    RCP<Factory> RepartitionFact2 = rcp(new RepartitionFactory());
    RepartitionFact2->SetFactory("A", rebA22Fact);
    RepartitionFact2->SetFactory("number of partitions", RepartitionHeuristicFact);
    RepartitionFact2->SetFactory("Partition", repInterface2);  // this is not valid
    RepartitionFact2->SetParameter("repartition: print partition distribution", Teuchos::ParameterEntry(true));
    RepartitionFact2->SetParameter("repartition: remap parts", Teuchos::ParameterEntry(false)); /* do not remap! */

    ////////////////////////////////////////// build non-rebalanced matrix blocks
    // build factories for transfer operator P(1,1) and R(1,1)
    RCP<AmalgamationFactory> amalgFact11 = rcp(new AmalgamationFactory());
    amalgFact11->SetFactory("A", A11Fact);
    amalgFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);

    RCP<CoalesceDropFactory> dropFact11 = rcp(new CoalesceDropFactory());
    dropFact11->SetFactory("A", A11Fact);
    dropFact11->SetFactory("UnAmalgamationInfo", amalgFact11);
    dropFact11->setDefaultVerbLevel(Teuchos::VERB_EXTREME);

    RCP<UncoupledAggregationFactory> UncoupledAggFact11 = rcp(new UncoupledAggregationFactory());
    UncoupledAggFact11->SetFactory("Graph", dropFact11);
    UncoupledAggFact11->SetMinNodesPerAggregate(9);
    UncoupledAggFact11->SetMaxNeighAlreadySelected(2);
    UncoupledAggFact11->SetOrdering("natural");

    RCP<CoarseMapFactory> coarseMapFact11 = Teuchos::rcp(new CoarseMapFactory());
    coarseMapFact11->setStridingData(stridingInfo);
    coarseMapFact11->setStridedBlockId(0);

    RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());
    RCP<TransPFactory> R11Fact     = rcp(new TransPFactory());

    Teuchos::RCP<NullspaceFactory> nspFact11 = Teuchos::rcp(new NullspaceFactory("Nullspace1"));
    nspFact11->SetFactory("Nullspace1", P11Fact);  // pick "Nullspace1" from Finest level

    //////////////////////////////// define factory manager for (1,1) block
    RCP<FactoryManager> M11 = rcp(new FactoryManager());
    M11->SetFactory("A", A11Fact);  // rebalanced fine-level block operator
    M11->SetFactory("P", P11Fact);  // non-rebalanced transfer operator block P(1,1)
    M11->SetFactory("R", R11Fact);  // non-rebalanced transfer operator block R(1,1)
    M11->SetFactory("Aggregates", UncoupledAggFact11);
    M11->SetFactory("Graph", dropFact11);
    M11->SetFactory("DofsPerNode", dropFact11);
    M11->SetFactory("UnAmalgamationInfo", amalgFact11);
    M11->SetFactory("Nullspace", nspFact11);  // TODO check me?
    M11->SetFactory("CoarseMap", coarseMapFact11);
    M11->SetIgnoreUserData(true);  // always use data from factories defined in factory manager

    ////////////////////////////////////////// build non-rebalanced matrix blocks
    // build factories for transfer operator P(2,2) and R(2,2)
    RCP<AmalgamationFactory> amalgFact22 = rcp(new AmalgamationFactory());
    RCP<TentativePFactory> P22Fact       = rcp(new TentativePFactory());
    RCP<TransPFactory> R22Fact           = rcp(new TransPFactory());

    // connect null space and tentative PFactory
    Teuchos::RCP<NullspaceFactory> nspFact22 = Teuchos::rcp(new NullspaceFactory("Nullspace2"));
    nspFact22->SetFactory("Nullspace2", P22Fact);  // define null space generated by P22Fact as null space for coarse level (non-rebalanced)

    RCP<CoarseMapFactory> coarseMapFact22 = Teuchos::rcp(new CoarseMapFactory());
    coarseMapFact22->setStridingData(stridingInfo);
    coarseMapFact22->setStridedBlockId(1);

    //////////////////////////////// define factory manager for (2,2) block
    RCP<FactoryManager> M22 = rcp(new FactoryManager());
    M22->SetFactory("A", A22Fact);                      // rebalanced fine-level block operator
    M22->SetFactory("P", P22Fact);                      // non-rebalanced transfer operator P(2,2)
    M22->SetFactory("R", R22Fact);                      // non-rebalanced transfer operator R(2,2)
    M22->SetFactory("Aggregates", UncoupledAggFact11);  // aggregates from block (1,1)
    M22->SetFactory("Nullspace", nspFact22);
    M22->SetFactory("UnAmalgamationInfo", amalgFact22);
    M22->SetFactory("Ptent", P22Fact);
    M22->SetFactory("CoarseMap", coarseMapFact22);
    M22->SetIgnoreUserData(true);

    /////////////////////////////////////////// define rebalanced blocked transfer ops
    //////////////////////////////// define factory manager for (1,1) block
    RCP<FactoryManager> rebM11 = rcp(new FactoryManager());
    rebM11->SetFactory("A", AcFact);  // important: must be a 2x2 block A Factory
    rebM11->SetFactory("Importer", RepartitionFact);
    rebM11->SetFactory("number of partitions", RepartitionHeuristicFact);
    rebM11->SetFactory("Nullspace", nspFact11);
    // rebM11->SetIgnoreUserData(true);

    RCP<FactoryManager> rebM22 = rcp(new FactoryManager());
    rebM22->SetFactory("A", AcFact);                   // important: must be a 2x2 block A Factory
    rebM22->SetFactory("Importer", RepartitionFact2);  // use dummy repartitioning factory
    rebM22->SetFactory("number of partitions", RepartitionHeuristicFact);
    rebM22->SetFactory("Nullspace", nspFact22);

    // Reordering of the transfer operators
    RCP<RebalanceBlockInterpolationFactory> RebalancedBlockPFact = rcp(new RebalanceBlockInterpolationFactory());
    RebalancedBlockPFact->SetFactory("P", PFact);  // use non-rebalanced block P operator as input
    RebalancedBlockPFact->AddFactoryManager(rebM11);
    RebalancedBlockPFact->AddFactoryManager(rebM22);

    RCP<RebalanceBlockRestrictionFactory> RebalancedBlockRFact = rcp(new RebalanceBlockRestrictionFactory());
    // RebalancedBlockRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
    RebalancedBlockRFact->SetFactory("R", RFact);  // non-rebalanced block P operator
    RebalancedBlockRFact->SetParameter("repartition: use subcommunicators", Teuchos::ParameterEntry(true));
    RebalancedBlockRFact->AddFactoryManager(rebM11);
    RebalancedBlockRFact->AddFactoryManager(rebM22);

    ///////////////////////////////////////// initialize non-rebalanced block transfer operators
    // output are the non-rebalanced block transfer operators used as input in AcFact to build
    // the non-rebalanced coarse level block matrix Ac
    PFact->AddFactoryManager(M11);  // use non-rebalanced information from sub block factory manager M11
    PFact->AddFactoryManager(M22);  // use non-rebalanced information from sub block factory manager M22

    ///////////////////////////////////////// initialize rebalanced coarse block AC factory
    RebalancedAcFact->SetFactory("A", AcFact);  // use non-rebalanced block operator as input
    RebalancedAcFact->SetParameter("repartition: use subcommunicators", Teuchos::ParameterEntry(true));
    RebalancedAcFact->AddFactoryManager(rebM11);
    RebalancedAcFact->AddFactoryManager(rebM22);

    //////////////////////////////////////////////////////////////////////
    // Smoothers

    // Another factory manager for braes sarazin smoother
    // Schur Complement Factory, using the factory to generate AcFact
    SC omega                          = 1.3;
    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(omega));
    SFact->SetFactory("A", MueLu::NoFactory::getRCP());  // this finally be the rebalanced block operator!

    // Smoother Factory, using SFact as a factory for A
    std::string ifpackSCType;
    Teuchos::ParameterList ifpackSCList;
    ifpackSCList.set("relaxation: sweeps", (LocalOrdinal)3);
    ifpackSCList.set("relaxation: damping factor", (Scalar)1.0);
    ifpackSCType = "RELAXATION";
    ifpackSCList.set("relaxation: type", "Gauss-Seidel");
    RCP<SmootherPrototype> smoProtoSC = rcp(new TrilinosSmoother(ifpackSCType, ifpackSCList, 0));
    smoProtoSC->SetFactory("A", SFact);
    RCP<SmootherFactory> SmooSCFact = rcp(new SmootherFactory(smoProtoSC));

    RCP<BraessSarazinSmoother> smootherPrototype = rcp(new BraessSarazinSmoother());
    smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(3));
    smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));
    smootherPrototype->SetFactory("A", MueLu::NoFactory::getRCP());
    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(smootherPrototype));

    RCP<BraessSarazinSmoother> coarseSolverPrototype = rcp(new BraessSarazinSmoother());
    coarseSolverPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(3));
    coarseSolverPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));
    coarseSolverPrototype->SetFactory("A", MueLu::NoFactory::getRCP());
    RCP<SmootherFactory> coarseSolverFact = rcp(new SmootherFactory(coarseSolverPrototype, Teuchos::null));

    RCP<FactoryManager> MB = rcp(new FactoryManager());
    MB->SetFactory("A", SFact);
    MB->SetFactory("Smoother", SmooSCFact);
    MB->SetIgnoreUserData(true);  // always use data from factories defined in factory manager
    smootherPrototype->AddFactoryManager(MB, 0);
    coarseSolverPrototype->AddFactoryManager(MB, 0);

    ////////////////////////////////////////// define main factory manager
    FactoryManager M;
    M.SetFactory("A", RebalancedAcFact);      // rebalance block AC Factory using importer
    M.SetFactory("P", RebalancedBlockPFact);  // rebalance prolongator using non-balanced Ac
    M.SetFactory("R", RebalancedBlockRFact);  // rebalance restrictor and null space using non-balanced Ac
    M.SetFactory("Smoother", smootherFact);
    M.SetFactory("PreSmoother", smootherFact);
    M.SetFactory("PostSmoother", smootherFact);
    M.SetFactory("CoarseSolver", coarseSolverFact);

    H->Setup(M, 0, maxLevels);

    /**out << std::endl;
     *out << "print content of multigrid levels:" << std::endl;

     Finest->print(*out);

     RCP<Level> coarseLevel = H->GetLevel(1);
     coarseLevel->print(*out);

     RCP<Level> coarseLevel2 = H->GetLevel(2);
     coarseLevel2->print(*out);*/

    RCP<MultiVector> xLsg = MultiVectorFactory::Build(xstridedfullmap, 1);

    // Use AMG directly as an iterative method
#if 0
    {
      xLsg->putScalar( (SC) 0.0);

      // Epetra_Vector -> Xpetra::Vector
      RCP<Vector> xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

      // calculate initial (absolute) residual
      Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
      xRhs->norm2(norms);
      *out << "||x_0|| = " << norms[0] << std::endl;

      // apply ten multigrid iterations
      H->Iterate(*xRhs,*xLsg,100);


      // calculate and print residual
      RCP<MultiVector> xTmp = MultiVectorFactory::Build(xstridedfullmap,1);
      bOp->apply(*xLsg,*xTmp,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
      xRhs->update((SC)-1.0,*xTmp,(SC)1.0);
      xRhs->norm2(norms);
      *out << "||x|| = " << norms[0] << std::endl;
    }
#endif

    //
    // Solve Ax = b using AMG as a preconditioner in AztecOO
    //
    {
      RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
      X->PutScalar(0.0);
      Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

      AztecOO aztecSolver(epetraProblem);
      aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

      MueLu::EpetraOperator aztecPrec(H);
      aztecSolver.SetPrecOperator(&aztecPrec);

      int maxIts = 50;
      double tol = 1e-8;

      aztecSolver.Iterate(maxIts, tol);
    }

#endif  // end ifndef HAVE_LONG_LONG_INT
#endif  // #if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  std::cout << "Epetra (and/or EpetraExt) are not available. Skip test." << std::endl;
  return EXIT_SUCCESS;
#endif
}

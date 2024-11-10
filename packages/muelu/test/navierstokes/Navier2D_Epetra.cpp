// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Navier2D_epetra.cpp
 *
 *  Created on: Mar 26, 2011
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
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
#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#endif

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

/*!
 *  2d Navier Stokes example (for Epetra)
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

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

    // Timing
    Teuchos::Time myTime("global");
    Teuchos::TimeMonitor m(myTime);

    //
    // SET TEST PARAMETERS
    //
    // Note: use --help to list available options.
    Teuchos::CommandLineProcessor clp(false);

    // - Levels
    LO optMaxLevels = 4;
    clp.setOption("maxLevels", &optMaxLevels, "maximum number of levels allowed");
    int optMaxCoarseSize = 1;
    clp.setOption("maxCoarseSize", &optMaxCoarseSize, "maximum #dofs in coarse operator");  // FIXME clp doesn't like long long int

    // - Repartitioning
#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
    int optRepartition = 1;
    clp.setOption("repartition", &optRepartition, "enable repartitioning");
    LO optMinRowsPerProc = 50;
    clp.setOption("minRowsPerProc", &optMinRowsPerProc, "min #rows allowable per proc before repartitioning occurs");
    double optNnzImbalance = 1.2;
    clp.setOption("nnzImbalance", &optNnzImbalance, "max allowable nonzero imbalance before repartitioning occurs");
#else
    int optRepartition = 0;
#endif

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    /////////////////////////////////////////////
    // custom parameters
    LO maxLevels = optMaxLevels;

    GO maxCoarseSize          = optMaxCoarseSize;
    std::string aggOrdering   = "natural";
    int minPerAgg             = 3;
    int maxNbrAlreadySelected = 0;

    int globalNumDofs = 1500;
    int nProcs        = comm->getSize();
    int nDofsPerNode  = 3;

    int nLocalDofs     = (int)globalNumDofs / nProcs;
    nLocalDofs         = nLocalDofs - (nLocalDofs % nDofsPerNode);
    int nCumulatedDofs = 0;
    MueLu_sumAll(comm, nLocalDofs, nCumulatedDofs);
    // Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM, 1, nLocalDofs, &nCumulatedDofs );

    if (comm->getRank() == nProcs - 1) {
      nLocalDofs += globalNumDofs - nCumulatedDofs;
    }

    std::cout << "PROC: " << comm->getRank() << " numLocalDofs=" << nLocalDofs << std::endl;

    // read in problem
    Epetra_Map emap(globalNumDofs, nLocalDofs, 0, *Xpetra::toEpetra(comm));
    Epetra_CrsMatrix* ptrA    = 0;
    Epetra_Vector* ptrf       = 0;
    Epetra_MultiVector* ptrNS = 0;

    std::cout << "Reading matrix market file" << std::endl;
    EpetraExt::MatrixMarketFileToCrsMatrix("A_re1000_5932.txt", emap, emap, emap, ptrA);
    EpetraExt::MatrixMarketFileToVector("b_re1000_5932.txt", emap, ptrf);
    RCP<Epetra_CrsMatrix> epA    = Teuchos::rcp(ptrA);
    RCP<Epetra_Vector> epv       = Teuchos::rcp(ptrf);
    RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);

    // Epetra_CrsMatrix -> Xpetra::Matrix
    RCP<CrsMatrix> exA       = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO, NO>(epA));
    RCP<CrsMatrixWrap> crsOp = Teuchos::rcp(new CrsMatrixWrap(exA));
    RCP<Matrix> Op           = Teuchos::rcp_dynamic_cast<Matrix>(crsOp);

    Op->SetFixedBlockSize(nDofsPerNode);  // 2 velocity dofs and 1 pressure dof per node.

    // Epetra_Vector -> Xpetra::Vector
    RCP<Vector> xRhs = Teuchos::rcp(new Xpetra::EpetraVectorT<int, Node>(epv));

    RCP<MultiVector> xNS = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<int, Node>(epNS));

    // Epetra_Map -> Xpetra::Map
    const RCP<const Map> map = Xpetra::toXpetra<int, Node>(emap);

    RCP<Hierarchy> H = rcp(new Hierarchy());
    H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    H->SetMaxCoarseSize(maxCoarseSize);

    // build finest Level
    RCP<MueLu::Level> Finest = H->GetLevel();
    Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
    Finest->Set("A", Op);
    // Finest->Set("Nullspace",xNS);

    if (optRepartition == 1) {
      // create null space

      RCP<MultiVector> nullspace;
      // determine numPDEs
      LocalOrdinal numPDEs = 1;
      if (Op->IsView("stridedMaps") == true) {
        Xpetra::viewLabel_t oldView = Op->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
        // TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const StridedMap>(Op->getRowMap()) == Teuchos::null, Exceptions::BadCast, "cast to strided row map failed.");
        numPDEs = Teuchos::rcp_dynamic_cast<const StridedMap>(Op->getRowMap())->getFixedBlockSize();
        oldView = Op->SwitchToView(oldView);
      }

      // GetOStream(Runtime1, 0) << "Generating canonical nullspace: dimension = " << numPDEs << std::endl;
      nullspace = MultiVectorFactory::Build(Op->getDomainMap(), numPDEs);

      for (int i = 0; i < numPDEs; ++i) {
        Teuchos::ArrayRCP<Scalar> nsValues = nullspace->getDataNonConst(i);
        int numBlocks                      = nsValues.size() / numPDEs;
        for (int j = 0; j < numBlocks; ++j) {
          nsValues[j * numPDEs + i] = 1.0;
        }
      }
      Finest->Set("Nullspace", nullspace);
    }

    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    dropFact->SetParameter("lightweight wrap", Teuchos::ParameterEntry(false));
    dropFact->SetVerbLevel(MueLu::Extreme);

    // RCP<PreDropFunctionConstVal> predrop = rcp(new PreDropFunctionConstVal(0.00001));
    // dropFact->SetPreDropFunction(predrop);
    RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
    UncoupledAggFact->SetFactory("Graph", dropFact);
    *out << "========================= Aggregate option summary =========================" << std::endl;
    *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
    *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
    UncoupledAggFact->SetMinNodesPerAggregate(minPerAgg);  // TODO should increase if run anything other than 1D
    UncoupledAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
    std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
    if (aggOrdering == "natural" || aggOrdering == "random" || aggOrdering == "graph") {
      *out << "aggregate ordering :                    " << aggOrdering << std::endl;
      UncoupledAggFact->SetOrdering(aggOrdering);
    } else {
      std::string msg =
          "main: bad aggregation option "
          "" +
          aggOrdering +
          ""
          ".";
      throw(MueLu::Exceptions::RuntimeError(msg));
    }
    *out << "=============================================================================" << std::endl;

    // build non-rebalanced transfer operators
    RCP<PgPFactory> Pfact = rcp(new PgPFactory());
    RCP<Factory> Rfact    = rcp(new GenericRFactory());
    // RCP<SaPFactory> Pfact  = rcp( new SaPFactory() );
    // RCP<Factory>   Rfact  = rcp( new TransPFactory() );
    RCP<RAPFactory> Acfact = rcp(new RAPFactory());
    Acfact->setVerbLevel(Teuchos::VERB_HIGH);

    // build level smoothers
    RCP<SmootherPrototype> smooProto;
    std::string ifpackType;
    Teuchos::ParameterList ifpackList;
    ifpackList.set("relaxation: sweeps", (LO)3);
    ifpackList.set("relaxation: damping factor", (SC)0.6);  // 0.7
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");

    smooProto = Teuchos::rcp(new TrilinosSmoother(ifpackType, ifpackList));
    RCP<SmootherFactory> SmooFact;
    if (maxLevels > 1)
      SmooFact = rcp(new SmootherFactory(smooProto));

    // create coarsest smoother
    RCP<SmootherPrototype> coarsestSmooProto;
    std::string type = "";
    Teuchos::ParameterList coarsestSmooList;
#if defined(HAVE_AMESOS_SUPERLU)
    coarsestSmooProto = Teuchos::rcp(new DirectSolver("Superlu", coarsestSmooList));
#else
    coarsestSmooProto  = Teuchos::rcp(new DirectSolver("Klu", coarsestSmooList));
#endif
    RCP<SmootherFactory> coarsestSmooFact = rcp(new SmootherFactory(coarsestSmooProto, Teuchos::null));

    FactoryManager M;
    M.SetKokkosRefactor(false);
    M.SetFactory("Graph", dropFact);
    // M.SetFactory("UnAmalgamationInfo", amalgFact);
    M.SetFactory("Aggregates", UncoupledAggFact);
    M.SetFactory("Smoother", SmooFact);
    M.SetFactory("CoarseSolver", coarsestSmooFact);

    if (optRepartition == 0) {
      // no rebalancing
      M.SetFactory("P", Pfact);
      M.SetFactory("R", Rfact);
      M.SetFactory("A", Acfact);
    } else {
#if defined(HAVE_MPI) && defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MUELU_ISORROPIA)
      // The Factory Manager will be configured to return the rebalanced versions of P, R, A by default.
      // Everytime we want to use the non-rebalanced versions, we need to explicitly define the generating factory.
      Rfact->SetFactory("P", Pfact);
      //
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
      RebalancedPFact->SetFactory("Nullspace", M.GetFactory("Ptent"));

      RCP<Factory> RebalancedRFact = rcp(new RebalanceTransferFactory());
      RebalancedRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
      RebalancedRFact->SetFactory("R", Rfact);
      // RebalancedRFact->SetFactory("Coordinates", TransferCoordinatesFact);

      // Compute Ac from rebalanced P and R
      RCP<Factory> RebalancedAFact = rcp(new RebalanceAcFactory());
      RebalancedAFact->SetFactory("A", Acfact);

      // Configure FactoryManager
      M.SetFactory("A", RebalancedAFact);
      M.SetFactory("P", RebalancedPFact);
      M.SetFactory("R", RebalancedRFact);
      M.SetFactory("Nullspace", RebalancedPFact);
      M.SetFactory("Importer", RepartitionFact);
#else
      // no re-balancing available
      M.SetFactory("P", Pfact);
      M.SetFactory("R", Rfact);
      M.SetFactory("A", Acfact);
#endif
    }

    H->Setup(M, 0, maxLevels);

    {  // some debug output
      // print out content of levels
      std::cout << "FINAL CONTENT of multigrid levels" << std::endl;
      for (LO l = 0; l < H->GetNumLevels(); l++) {
        RCP<Level> coarseLevel = H->GetLevel(l);
        coarseLevel->print(*out);
      }
      std::cout << "END FINAL CONTENT of multigrid levels" << std::endl;
    }  // end debug output

    RCP<MultiVector> xLsg = MultiVectorFactory::Build(map, 1);

    // Use AMG directly as an iterative method
    {
      xLsg->putScalar((SC)0.0);

      // calculate initial (absolute) residual
      Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
      xRhs->norm2(norms);
      *out << "||x_0|| = " << norms[0] << std::endl;

      // apply ten multigrid iterations
      H->Iterate(*xRhs, *xLsg, 10);

      // calculate and print residual
      RCP<MultiVector> xTmp = MultiVectorFactory::Build(map, 1);
      Op->apply(*xLsg, *xTmp, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);
      xRhs->update((SC)-1.0, *xTmp, (SC)1.0);
      xRhs->norm2(norms);
      *out << "||x|| = " << norms[0] << std::endl;
    }

    //
    // Solve Ax = b using AMG as a preconditioner in AztecOO
    //
    {
      RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
      X->PutScalar(0.0);
      Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

      AztecOO aztecSolver(epetraProblem);
      aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

#if 0
      // TODO TAW: 4/8/2016
      // temporarely deactivate this due to runtime error on perseus:
      // Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed
      // if SERIAL=OFF, OPENMP=OFF, PTHREAD=ON, CUDA=OFF
      // probably a fix necessary in EpetraOperator (which only supports
      // SERIAL or OPENMP, but not PTHREAD of course).
      MueLu::EpetraOperator aztecPrec(H);
      aztecSolver.SetPrecOperator(&aztecPrec);

      int maxIts = 50;
      double tol = 1e-8;

      aztecSolver.Iterate(maxIts, tol);
#endif
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  std::cout << "Epetra (and/or EpetraExt) are not available. Skip test." << std::endl;
  return EXIT_SUCCESS;
#endif
}

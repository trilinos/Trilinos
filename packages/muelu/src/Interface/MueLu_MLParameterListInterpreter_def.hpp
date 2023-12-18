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
#ifndef MUELU_MLPARAMETERLISTINTERPRETER_DEF_HPP
#define MUELU_MLPARAMETERLISTINTERPRETER_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ML)
#include <ml_ValidateParameters.h>
#endif

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Operator.hpp>

#include "MueLu_MLParameterListInterpreter_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_ParameterListUtils.hpp"

#include "MueLu_CoalesceDropFactory_kokkos.hpp"
// #include "MueLu_CoordinatesTransferFactory_kokkos.hpp"
// #include "MueLu_NullspaceFactory_kokkos.hpp"
#include "MueLu_SaPFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"

#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)
#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
//#include "MueLu_RebalanceMapFactory.hpp"
#endif

// Note: do not add options that are only recognized by MueLu.

// TODO: this parameter list interpreter should force MueLu to use default ML parameters
// - Ex: smoother sweep=2 by default for ML

// Read a parameter value from a parameter list and store it into a variable named 'varName'
#define MUELU_READ_PARAM(paramList, paramStr, varType, defaultValue, varName) \
  varType varName = defaultValue;                                             \
  if (paramList.isParameter(paramStr)) varName = paramList.get<varType>(paramStr);

// Read a parameter value from a paraeter list and copy it into a new parameter list (with another parameter name)
#define MUELU_COPY_PARAM(paramList, paramStr, varType, defaultValue, outParamList, outParamStr) \
  if (paramList.isParameter(paramStr))                                                          \
    outParamList.set(outParamStr, paramList.get<varType>(paramStr));                            \
  else                                                                                          \
    outParamList.set(outParamStr, static_cast<varType>(defaultValue));

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MLParameterListInterpreter(Teuchos::ParameterList& paramList, Teuchos::RCP<const Teuchos::Comm<int> > comm, std::vector<RCP<FactoryBase> > factoryList)
  : nullspace_(NULL)
  , xcoord_(NULL)
  , ycoord_(NULL)
  , zcoord_(NULL)
  , TransferFacts_(factoryList)
  , blksize_(1) {
  if (paramList.isParameter("xml parameter file")) {
    std::string filename = paramList.get("xml parameter file", "");
    if (filename.length() != 0) {
      TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), Exceptions::RuntimeError, "xml parameter file requires a valid comm");
      Teuchos::ParameterList paramList2 = paramList;
      Teuchos::updateParametersFromXmlFileAndBroadcast(filename, Teuchos::Ptr<Teuchos::ParameterList>(&paramList2), *comm);
      paramList2.remove("xml parameter file");
      SetParameterList(paramList2);
    } else
      SetParameterList(paramList);
  } else
    SetParameterList(paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MLParameterListInterpreter(const std::string& xmlFileName, std::vector<RCP<FactoryBase> > factoryList)
  : nullspace_(NULL)
  , TransferFacts_(factoryList)
  , blksize_(1) {
  Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
  SetParameterList(*paramList);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList_in) {
  Teuchos::ParameterList paramList = paramList_in;

  //
  // Read top-level of the parameter list
  //

  // hard-coded default values == ML defaults according to the manual
  MUELU_READ_PARAM(paramList, "ML output", int, 0, verbosityLevel);
  MUELU_READ_PARAM(paramList, "max levels", int, 10, maxLevels);
  MUELU_READ_PARAM(paramList, "PDE equations", int, 1, nDofsPerNode);

  MUELU_READ_PARAM(paramList, "coarse: max size", int, 128, maxCoarseSize);

  MUELU_READ_PARAM(paramList, "aggregation: type", std::string, "Uncoupled", agg_type);
  // MUELU_READ_PARAM(paramList, "aggregation: threshold",                double,                 0.0,       agg_threshold);
  MUELU_READ_PARAM(paramList, "aggregation: damping factor", double, (double)4 / (double)3, agg_damping);
  // MUELU_READ_PARAM(paramList, "aggregation: smoothing sweeps",            int,                   1,       agg_smoothingsweeps);
  MUELU_READ_PARAM(paramList, "aggregation: nodes per aggregate", int, 1, minPerAgg);
  MUELU_READ_PARAM(paramList, "aggregation: keep Dirichlet bcs", bool, false, bKeepDirichletBcs);                // This is a MueLu specific extension that does not exist in ML
  MUELU_READ_PARAM(paramList, "aggregation: max neighbours already aggregated", int, 0, maxNbrAlreadySelected);  // This is a MueLu specific extension that does not exist in M
  MUELU_READ_PARAM(paramList, "aggregation: aux: enable", bool, false, agg_use_aux);
  MUELU_READ_PARAM(paramList, "aggregation: aux: threshold", double, false, agg_aux_thresh);

  MUELU_READ_PARAM(paramList, "null space: type", std::string, "default vectors", nullspaceType);
  MUELU_READ_PARAM(paramList, "null space: dimension", int, -1, nullspaceDim);      // TODO: ML default not in documentation
  MUELU_READ_PARAM(paramList, "null space: vectors", double*, NULL, nullspaceVec);  // TODO: ML default not in documentation

  MUELU_READ_PARAM(paramList, "energy minimization: enable", bool, false, bEnergyMinimization);

  MUELU_READ_PARAM(paramList, "RAP: fix diagonal", bool, false, bFixDiagonal);  // This is a MueLu specific extension that does not exist in ML

  MUELU_READ_PARAM(paramList, "x-coordinates", double*, NULL, xcoord);
  MUELU_READ_PARAM(paramList, "y-coordinates", double*, NULL, ycoord);
  MUELU_READ_PARAM(paramList, "z-coordinates", double*, NULL, zcoord);

  //
  // Move smoothers/aggregation/coarse parameters to sublists
  //

  // ML allows to have level-specific smoothers/aggregation/coarse parameters at the top level of the list or/and defined in sublists:
  // See also: ML Guide section 6.4.1, MueLu::CreateSublists, ML_CreateSublists
  ParameterList paramListWithSubList;
  MueLu::CreateSublists(paramList, paramListWithSubList);
  paramList = paramListWithSubList;  // swap

  // pull out "use kokkos refactor"
  bool setKokkosRefactor = false;
  bool useKokkosRefactor = !Node::is_serial;
  if (paramList.isType<bool>("use kokkos refactor")) {
    useKokkosRefactor = paramList.get<bool>("use kokkos refactor");
    setKokkosRefactor = true;
    paramList.remove("use kokkos refactor");
  }

  //
  // Validate parameter list
  //

  {
    bool validate = paramList.get("ML validate parameter list", true); /* true = default in ML */
    if (validate) {
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
      // Validate parameter list using ML validator
      int depth = paramList.get("ML validate depth", 5); /* 5 = default in ML */
      TEUCHOS_TEST_FOR_EXCEPTION(!ML_Epetra::ValidateMLPParameters(paramList, depth), Exceptions::RuntimeError,
                                 "ERROR: ML's Teuchos::ParameterList contains incorrect parameter!");
#else
      // If no validator available: issue a warning and set parameter value to false in the output list
      this->GetOStream(Warnings0) << "Warning: MueLu_ENABLE_ML=OFF. The parameter list cannot be validated." << std::endl;
      paramList.set("ML validate parameter list", false);

#endif  // HAVE_MUELU_ML
    }   // if(validate)
  }     // scope

  // Matrix option
  blksize_ = nDofsPerNode;

  // Translate verbosity parameter

  // Translate verbosity parameter
  MsgType eVerbLevel = None;
  if (verbosityLevel == 0) eVerbLevel = None;
  if (verbosityLevel >= 1) eVerbLevel = Low;
  if (verbosityLevel >= 5) eVerbLevel = Medium;
  if (verbosityLevel >= 10) eVerbLevel = High;
  if (verbosityLevel >= 11) eVerbLevel = Extreme;
  if (verbosityLevel >= 42) eVerbLevel = Test;
  if (verbosityLevel >= 43) eVerbLevel = InterfaceTest;
  this->verbosity_ = eVerbLevel;

  TEUCHOS_TEST_FOR_EXCEPTION(agg_type != "Uncoupled", Exceptions::RuntimeError,
                             "MueLu::MLParameterListInterpreter::SetParameterList(): parameter \"aggregation: type\": only 'Uncoupled' aggregation is supported.");

  // Create MueLu factories
  RCP<Factory> dropFact;
  if (useKokkosRefactor)
    dropFact = rcp(new CoalesceDropFactory_kokkos());
  else
    dropFact = rcp(new CoalesceDropFactory());

  if (agg_use_aux) {
    dropFact->SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
    dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(agg_aux_thresh));
  }

  // Uncoupled aggregation
  RCP<Factory> AggFact = Teuchos::null;
  if (useKokkosRefactor) {
    AggFact = rcp(new UncoupledAggregationFactory_kokkos());
  } else
    AggFact = rcp(new UncoupledAggregationFactory());

  AggFact->SetFactory("Graph", dropFact);
  AggFact->SetFactory("DofsPerNode", dropFact);
  AggFact->SetParameter("aggregation: preserve Dirichlet points", Teuchos::ParameterEntry(bKeepDirichletBcs));
  AggFact->SetParameter("aggregation: ordering", Teuchos::ParameterEntry(std::string("natural")));
  AggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(maxNbrAlreadySelected));
  AggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(minPerAgg));

  if (verbosityLevel > 3) {
    std::ostringstream oss;
    oss << "========================= Aggregate option summary  =========================" << std::endl;
    oss << "min Nodes per aggregate :               " << minPerAgg << std::endl;
    oss << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
    oss << "aggregate ordering :                    natural" << std::endl;
    oss << "=============================================================================" << std::endl;
    this->GetOStream(Runtime1) << oss.str();
  }

  RCP<Factory> PFact;
  RCP<Factory> RFact;
  RCP<Factory> PtentFact;
  if (useKokkosRefactor)
    PtentFact = rcp(new TentativePFactory_kokkos());
  else
    PtentFact = rcp(new TentativePFactory());
  if (agg_damping == 0.0 && bEnergyMinimization == false) {
    // tentative prolongation operator (PA-AMG)
    PFact = PtentFact;
    RFact = rcp(new TransPFactory());
  } else if (agg_damping != 0.0 && bEnergyMinimization == false) {
    // smoothed aggregation (SA-AMG)
    RCP<Factory> SaPFact;
    if (useKokkosRefactor)
      SaPFact = rcp(new SaPFactory_kokkos());
    else
      SaPFact = rcp(new SaPFactory());
    SaPFact->SetParameter("sa: damping factor", ParameterEntry(agg_damping));
    PFact = SaPFact;
    RFact = rcp(new TransPFactory());
  } else if (bEnergyMinimization == true) {
    // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
    PFact = rcp(new PgPFactory());
    RFact = rcp(new GenericRFactory());
  }

  RCP<RAPFactory> AcFact = rcp(new RAPFactory());
  AcFact->SetParameter("RepairMainDiagonal", Teuchos::ParameterEntry(bFixDiagonal));
  for (size_t i = 0; i < TransferFacts_.size(); i++) {
    AcFact->AddTransferFactory(TransferFacts_[i]);
  }

  //
  // introduce rebalancing
  //
#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)
  Teuchos::RCP<Factory> RebalancedPFact            = Teuchos::null;
  Teuchos::RCP<Factory> RebalancedRFact            = Teuchos::null;
  Teuchos::RCP<Factory> RepartitionFact            = Teuchos::null;
  Teuchos::RCP<RebalanceAcFactory> RebalancedAFact = Teuchos::null;

  MUELU_READ_PARAM(paramList, "repartition: enable", int, 0, bDoRepartition);
  if (bDoRepartition == 1) {
    // The Factory Manager will be configured to return the rebalanced versions of P, R, A by default.
    // Everytime we want to use the non-rebalanced versions, we need to explicitly define the generating factory.
    RFact->SetFactory("P", PFact);
    //
    AcFact->SetFactory("P", PFact);
    AcFact->SetFactory("R", RFact);

    // define rebalancing factory for coarse matrix
    Teuchos::RCP<MueLu::AmalgamationFactory<SC, LO, GO, NO> > rebAmalgFact = Teuchos::rcp(new MueLu::AmalgamationFactory<SC, LO, GO, NO>());
    rebAmalgFact->SetFactory("A", AcFact);

    MUELU_READ_PARAM(paramList, "repartition: max min ratio", double, 1.3, maxminratio);
    MUELU_READ_PARAM(paramList, "repartition: min per proc", int, 512, minperproc);

    // Repartitioning heuristic
    RCP<RepartitionHeuristicFactory> RepartitionHeuristicFact = Teuchos::rcp(new RepartitionHeuristicFactory());
    {
      Teuchos::ParameterList paramListRepFact;
      paramListRepFact.set("repartition: min rows per proc", minperproc);
      paramListRepFact.set("repartition: max imbalance", maxminratio);
      RepartitionHeuristicFact->SetParameterList(paramListRepFact);
    }
    RepartitionHeuristicFact->SetFactory("A", AcFact);

    // create "Partition"
    Teuchos::RCP<MueLu::IsorropiaInterface<LO, GO, NO> > isoInterface = Teuchos::rcp(new MueLu::IsorropiaInterface<LO, GO, NO>());
    isoInterface->SetFactory("A", AcFact);
    isoInterface->SetFactory("number of partitions", RepartitionHeuristicFact);
    isoInterface->SetFactory("UnAmalgamationInfo", rebAmalgFact);

    // create "Partition" by unamalgamtion
    Teuchos::RCP<MueLu::RepartitionInterface<LO, GO, NO> > repInterface = Teuchos::rcp(new MueLu::RepartitionInterface<LO, GO, NO>());
    repInterface->SetFactory("A", AcFact);
    repInterface->SetFactory("number of partitions", RepartitionHeuristicFact);
    repInterface->SetFactory("AmalgamatedPartition", isoInterface);
    // repInterface->SetFactory("UnAmalgamationInfo", rebAmalgFact); // not necessary?

    // Repartitioning (creates "Importer" from "Partition")
    RepartitionFact = Teuchos::rcp(new RepartitionFactory());
    RepartitionFact->SetFactory("A", AcFact);
    RepartitionFact->SetFactory("number of partitions", RepartitionHeuristicFact);
    RepartitionFact->SetFactory("Partition", repInterface);

    // Reordering of the transfer operators
    RebalancedPFact = Teuchos::rcp(new RebalanceTransferFactory());
    RebalancedPFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Interpolation")));
    RebalancedPFact->SetFactory("P", PFact);
    RebalancedPFact->SetFactory("Nullspace", PtentFact);
    RebalancedPFact->SetFactory("Importer", RepartitionFact);

    RebalancedRFact = Teuchos::rcp(new RebalanceTransferFactory());
    RebalancedRFact->SetParameter("type", Teuchos::ParameterEntry(std::string("Restriction")));
    RebalancedRFact->SetFactory("R", RFact);
    RebalancedRFact->SetFactory("Importer", RepartitionFact);

    // Compute Ac from rebalanced P and R
    RebalancedAFact = Teuchos::rcp(new RebalanceAcFactory());
    RebalancedAFact->SetFactory("A", AcFact);
  }
#else  // #ifdef HAVE_MUELU_ISORROPIA
  // Get rid of [-Wunused] warnings
  //(void)
  //
  // ^^^ FIXME (mfh 17 Nov 2013) That definitely doesn't compile.
#endif

  //
  // Nullspace factory
  //

  // Set fine level nullspace
  // extract pre-computed nullspace from ML parameter list
  // store it in nullspace_ and nullspaceDim_
  if (nullspaceType != "default vectors") {
    TEUCHOS_TEST_FOR_EXCEPTION(nullspaceType != "pre-computed", Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no valid nullspace (no pre-computed null space). error.");
    TEUCHOS_TEST_FOR_EXCEPTION(nullspaceDim == -1, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no valid nullspace (nullspace dim == -1). error.");
    TEUCHOS_TEST_FOR_EXCEPTION(nullspaceVec == NULL, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no valid nullspace (nullspace == NULL). You have to provide a valid fine-level nullspace in \'null space: vectors\'");

    nullspaceDim_ = nullspaceDim;
    nullspace_    = nullspaceVec;
  }

  Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory("Nullspace"));
  nspFact->SetFactory("Nullspace", PtentFact);

  // Stash coordinates
  xcoord_ = xcoord;
  ycoord_ = ycoord;
  zcoord_ = zcoord;

  //
  // Hierarchy + FactoryManager
  //

  // Hierarchy options
  this->numDesiredLevel_ = maxLevels;
  this->maxCoarseSize_   = maxCoarseSize;

  //
  // Coarse Smoother
  //
  ParameterList& coarseList = paramList.sublist("coarse: list");
  // check whether coarse solver is set properly. If not, set default coarse solver.
  if (!coarseList.isParameter("smoother: type"))
    coarseList.set("smoother: type", "Amesos-KLU");  // set default coarse solver according to ML 5.0 guide
  RCP<SmootherFactory> coarseFact = GetSmootherFactory(coarseList, Teuchos::null);

  // Smoothers Top Level Parameters

  RCP<ParameterList> topLevelSmootherParam = ExtractSetOfParameters(paramList, "smoother");

  //

  // Prepare factory managers
  // TODO: smootherFact can be reuse accross level if same parameters/no specific parameterList

  for (int levelID = 0; levelID < maxLevels; levelID++) {
    //
    // Level FactoryManager
    //

    RCP<FactoryManager> manager = rcp(new FactoryManager());
    if (setKokkosRefactor)
      manager->SetKokkosRefactor(useKokkosRefactor);

    //
    // Smoothers
    //

    {
      // Merge level-specific parameters with global parameters. level-specific parameters takes precedence.
      // TODO: unit-test this part alone

      ParameterList levelSmootherParam = GetMLSubList(paramList, "smoother", levelID);  // copy
      MergeParameterList(*topLevelSmootherParam, levelSmootherParam, false);            /* false = do no overwrite levelSmootherParam parameters by topLevelSmootherParam parameters */
      // std::cout << std::endl << "Merged List for level  " << levelID << std::endl;
      // std::cout << levelSmootherParam << std::endl;

      RCP<SmootherFactory> smootherFact = GetSmootherFactory(levelSmootherParam, Teuchos::null);  // TODO: missing AFact input arg.

      manager->SetFactory("Smoother", smootherFact);
    }

    //
    // Misc
    //

    manager->SetFactory("CoarseSolver", coarseFact);  // TODO: should not be done in the loop
    manager->SetFactory("Graph", dropFact);
    manager->SetFactory("Aggregates", AggFact);
    manager->SetFactory("DofsPerNode", dropFact);
    manager->SetFactory("Ptent", PtentFact);

#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)
    if (bDoRepartition == 1) {
      manager->SetFactory("A", RebalancedAFact);
      manager->SetFactory("P", RebalancedPFact);
      manager->SetFactory("R", RebalancedRFact);
      manager->SetFactory("Nullspace", RebalancedPFact);
      manager->SetFactory("Importer", RepartitionFact);
    } else {
#endif                                            // #ifdef HAVE_MUELU_ISORROPIA
      manager->SetFactory("Nullspace", nspFact);  // use same nullspace factory throughout all multigrid levels
      manager->SetFactory("A", AcFact);           // same RAP factory for all levels
      manager->SetFactory("P", PFact);            // same prolongator and restrictor factories for all levels
      manager->SetFactory("R", RFact);            // same prolongator and restrictor factories for all levels
#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)
    }
#endif

    this->AddFactoryManager(levelID, 1, manager);
  }  // for (level loop)
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupHierarchy(Hierarchy& H) const {
  // if nullspace_ has already been extracted from ML parameter list
  // make nullspace available for MueLu
  if (nullspace_ != NULL) {
    RCP<Level> fineLevel = H.GetLevel(0);
    RCP<Operator> Op     = fineLevel->Get<RCP<Operator> >("A");
    RCP<Matrix> A        = rcp_dynamic_cast<Matrix>(Op);
    if (!A.is_null()) {
      const RCP<const Map> rowMap = fineLevel->Get<RCP<Matrix> >("A")->getRowMap();
      RCP<MultiVector> nullspace  = MultiVectorFactory::Build(rowMap, nullspaceDim_, true);

      for (size_t i = 0; i < Teuchos::as<size_t>(nullspaceDim_); i++) {
        Teuchos::ArrayRCP<Scalar> nullspacei = nullspace->getDataNonConst(i);
        const size_t myLength                = nullspace->getLocalLength();

        for (size_t j = 0; j < myLength; j++) {
          nullspacei[j] = nullspace_[i * myLength + j];
        }
      }

      fineLevel->Set("Nullspace", nullspace);
    }
  }

  // Do the same for coordinates
  size_t num_coords = 0;
  double* coordPTR[3];
  if (xcoord_) {
    coordPTR[0] = xcoord_;
    num_coords++;
    if (ycoord_) {
      coordPTR[1] = ycoord_;
      num_coords++;
      if (zcoord_) {
        coordPTR[2] = zcoord_;
        num_coords++;
      }
    }
  }
  if (num_coords) {
    Teuchos::RCP<Level> fineLevel = H.GetLevel(0);
    Teuchos::RCP<Operator> Op     = fineLevel->Get<RCP<Operator> >("A");
    Teuchos::RCP<Matrix> A        = rcp_dynamic_cast<Matrix>(Op);
    if (!A.is_null()) {
      const Teuchos::RCP<const Map> rowMap  = fineLevel->Get<RCP<Matrix> >("A")->getRowMap();
      Teuchos::RCP<MultiVector> coordinates = MultiVectorFactory::Build(rowMap, num_coords, true);

      for (size_t i = 0; i < num_coords; i++) {
        Teuchos::ArrayRCP<Scalar> coordsi = coordinates->getDataNonConst(i);
        const size_t myLength             = coordinates->getLocalLength();
        for (size_t j = 0; j < myLength; j++) {
          coordsi[j] = coordPTR[i][j];
        }
      }
      fineLevel->Set("Coordinates", coordinates);
    }
  }

  HierarchyManager::SetupHierarchy(H);
}

// TODO: code factorization with MueLu_ParameterListInterpreter.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetSmootherFactory(const Teuchos::ParameterList& paramList,
                       const RCP<FactoryBase>& AFact) {
  typedef Teuchos::ScalarTraits<Scalar> STS;
  SC one = STS::one();

  std::string type = "symmetric Gauss-Seidel";  // default

  //
  // Get 'type'
  //

  //     //TODO: fix defaults!!

  //     // Default coarse grid smoother
  //     std::string type;
  //     if ("smoother" == "coarse") {
  // #if (defined(HAVE_MUELU_EPETRA) && defined( HAVE_MUELU_AMESOS)) || (defined(HAVE_MUELU_AMESOS2)) // FIXME: test is wrong (ex: compiled with Epetra&&Tpetra&&Amesos2 but without Amesos => error running Epetra problem)
  //       type = ""; // use default defined by AmesosSmoother or Amesos2Smoother
  // #else
  //       type = "symmetric Gauss-Seidel"; // use a sym Gauss-Seidel (with no damping) as fallback "coarse solver" (TODO: needs Ifpack(2))
  // #endif
  //     } else {
  //       // TODO: default smoother?
  //       type = "";
  //     }

  if (paramList.isParameter("smoother: type")) type = paramList.get<std::string>("smoother: type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no \"smoother: type\" in the smoother parameter list" << std::endl
                                                                                                                                                               << paramList);

  //
  // Create the smoother prototype
  //

  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList smootherParamList;

  if (type == "Jacobi" || type == "Gauss-Seidel" || type == "symmetric Gauss-Seidel") {
    if (type == "symmetric Gauss-Seidel") type = "Symmetric Gauss-Seidel";  // FIXME

    ifpackType = "RELAXATION";
    smootherParamList.set("relaxation: type", type);

    MUELU_COPY_PARAM(paramList, "smoother: sweeps", int, 2, smootherParamList, "relaxation: sweeps");
    MUELU_COPY_PARAM(paramList, "smoother: damping factor", Scalar, one, smootherParamList, "relaxation: damping factor");

    smooProto = rcp(new TrilinosSmoother(ifpackType, smootherParamList, 0));
    smooProto->SetFactory("A", AFact);

  } else if (type == "Chebyshev" || type == "MLS") {
    ifpackType = "CHEBYSHEV";

    MUELU_COPY_PARAM(paramList, "smoother: sweeps", int, 2, smootherParamList, "chebyshev: degree");
    if (paramList.isParameter("smoother: MLS alpha")) {
      MUELU_COPY_PARAM(paramList, "smoother: MLS alpha", double, 20, smootherParamList, "chebyshev: ratio eigenvalue");
    } else {
      MUELU_COPY_PARAM(paramList, "smoother: Chebyshev alpha", double, 20, smootherParamList, "chebyshev: ratio eigenvalue");
    }

    smooProto = rcp(new TrilinosSmoother(ifpackType, smootherParamList, 0));
    smooProto->SetFactory("A", AFact);

  } else if (type == "Hiptmair") {
    ifpackType                  = "HIPTMAIR";
    std::string subSmootherType = "Chebyshev";
    if (paramList.isParameter("subsmoother: type"))
      subSmootherType = paramList.get<std::string>("subsmoother: type");
    std::string subSmootherIfpackType;
    if (subSmootherType == "Chebyshev")
      subSmootherIfpackType = "CHEBYSHEV";
    else if (subSmootherType == "Jacobi" || subSmootherType == "Gauss-Seidel" || subSmootherType == "symmetric Gauss-Seidel") {
      if (subSmootherType == "symmetric Gauss-Seidel") subSmootherType = "Symmetric Gauss-Seidel";  // FIXME
      subSmootherIfpackType = "RELAXATION";
    } else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: unknown smoother type. '" << subSmootherType << "' not supported by MueLu.");

    smootherParamList.set("hiptmair: smoother type 1", subSmootherIfpackType);
    smootherParamList.set("hiptmair: smoother type 2", subSmootherIfpackType);

    auto smoother1ParamList = smootherParamList.sublist("hiptmair: smoother list 1");
    auto smoother2ParamList = smootherParamList.sublist("hiptmair: smoother list 2");

    if (subSmootherType == "Chebyshev") {
      MUELU_COPY_PARAM(paramList, "subsmoother: edge sweeps", int, 2, smoother1ParamList, "chebyshev: degree");
      MUELU_COPY_PARAM(paramList, "subsmoother: node sweeps", int, 2, smoother2ParamList, "chebyshev: degree");

      MUELU_COPY_PARAM(paramList, "subsmoother: Chebyshev", double, 20, smoother1ParamList, "chebyshev: ratio eigenvalue");
      MUELU_COPY_PARAM(paramList, "subsmoother: Chebyshev", double, 20, smoother2ParamList, "chebyshev: ratio eigenvalue");
    } else {
      MUELU_COPY_PARAM(paramList, "subsmoother: edge sweeps", int, 2, smoother1ParamList, "relaxation: sweeps");
      MUELU_COPY_PARAM(paramList, "subsmoother: node sweeps", int, 2, smoother2ParamList, "relaxation: sweeps");

      MUELU_COPY_PARAM(paramList, "subsmoother: SGS damping factor", double, 0.8, smoother2ParamList, "relaxation: damping factor");
    }

    smooProto = rcp(new TrilinosSmoother(ifpackType, smootherParamList, 0));
    smooProto->SetFactory("A", AFact);

  } else if (type == "IFPACK") {  // TODO: this option is not described in the ML Guide v5.0

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
    ifpackType = paramList.get<std::string>("smoother: ifpack type");

    if (ifpackType == "ILU") {
      // TODO fix this (type mismatch double vs. int)
      // MUELU_COPY_PARAM(paramList, "smoother: ifpack level-of-fill", double /*int*/, 0.0 /*2*/,  smootherParamList, "fact: level-of-fill");
      if (paramList.isParameter("smoother: ifpack level-of-fill"))
        smootherParamList.set("fact: level-of-fill", Teuchos::as<int>(paramList.get<double>("smoother: ifpack level-of-fill")));
      else
        smootherParamList.set("fact: level-of-fill", as<int>(0));

      MUELU_COPY_PARAM(paramList, "smoother: ifpack overlap", int, 2, smootherParamList, "partitioner: overlap");

      // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
      smooProto =
          MueLu::GetIfpackSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(ifpackType,
                                                                              smootherParamList,
                                                                              paramList.get<int>("smoother: ifpack overlap"));
      smooProto->SetFactory("A", AFact);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
    }
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: MueLu compiled without Ifpack support");
#endif

  } else if (type.length() > strlen("Amesos") && type.substr(0, strlen("Amesos")) == "Amesos") { /* catch Amesos-* */
    std::string solverType = type.substr(strlen("Amesos") + 1);                                  /* ("Amesos-KLU" -> "KLU") */

    // Validator: following upper/lower case is what is allowed by ML
    bool valid                           = false;
    const int validatorSize              = 5;
    std::string validator[validatorSize] = {"Superlu", "Superludist", "KLU", "UMFPACK", "MUMPS"}; /* TODO: should "" be allowed? */
    for (int i = 0; i < validatorSize; i++) {
      if (validator[i] == solverType) valid = true;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!valid, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: unknown smoother type. '" << type << "' not supported.");

    // FIXME: MueLu should accept any Upper/Lower case. Not the case for the moment
    std::transform(solverType.begin() + 1, solverType.end(), solverType.begin() + 1, ::tolower);

    smooProto = Teuchos::rcp(new DirectSolver(solverType, Teuchos::ParameterList()));
    smooProto->SetFactory("A", AFact);

  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: unknown smoother type. '" << type << "' not supported by MueLu.");
  }
  TEUCHOS_TEST_FOR_EXCEPTION(smooProto == Teuchos::null, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no smoother prototype. fatal error.");

  //
  // Create the smoother factory
  //

  RCP<SmootherFactory> SmooFact = rcp(new SmootherFactory());

  // Set parameters of the smoother factory
  MUELU_READ_PARAM(paramList, "smoother: pre or post", std::string, "both", preOrPost);
  if (preOrPost == "both") {
    SmooFact->SetSmootherPrototypes(smooProto, smooProto);
  } else if (preOrPost == "pre") {
    SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
  } else if (preOrPost == "post") {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
  }

  return SmooFact;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "Transfer factory is not derived from TwoLevelFactoryBase. Since transfer factories will be handled by the RAPFactory they have to be derived from TwoLevelFactoryBase!");
  TransferFacts_.push_back(factory);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NumTransferFactories() const {
  return TransferFacts_.size();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupOperator(Operator& Op) const {
  try {
    Matrix& A = dynamic_cast<Matrix&>(Op);
    if (A.IsFixedBlockSizeSet() && (A.GetFixedBlockSize() != blksize_))
      this->GetOStream(Warnings0) << "Setting matrix block size to " << blksize_ << " (value of the parameter in the list) "
                                  << "instead of " << A.GetFixedBlockSize() << " (provided matrix)." << std::endl;

    A.SetFixedBlockSize(blksize_);

#ifdef HAVE_MUELU_DEBUG
    MatrixUtils::checkLocalRowMapMatchesColMap(A);
#endif  // HAVE_MUELU_DEBUG

  } catch (std::bad_cast&) {
    this->GetOStream(Warnings0) << "Skipping setting block size as the operator is not a matrix" << std::endl;
  }
}

}  // namespace MueLu

#define MUELU_MLPARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_MLPARAMETERLISTINTERPRETER_DEF_HPP */

// TODO: see if it can be factorized with ML interpreter (ex: generation of Ifpack param list)

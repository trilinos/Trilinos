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
#ifndef MUELU_EASYPARAMETERLISTINTERPRETER_DEF_HPP
#define MUELU_EASYPARAMETERLISTINTERPRETER_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_EasyParameterListInterpreter_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_EminPFactory.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"

namespace MueLu {

  // This macro is tricky. The use case is when we do not have a level specific parameter, so we
  // need to take the default value from the general list for all levels, or set it to the default.
#define MUELU_READ_2LIST_PARAM(paramList, defaultList, paramStr, varType, defaultValue, varName) \
  varType varName; \
  if      (paramList.isParameter(paramStr))   varName = paramList.get<varType>(paramStr); \
  else if (defaultList.isParameter(paramStr)) varName = defaultList.get<varType>(paramStr); \
  else                                        varName = paramList.get<varType>(paramStr, defaultValue);

  // This macro check whether the variable is in the list.
  // If it is, it copies its value to the second list, possibly with a new name
  // Similar to the above macro, we all try to take a value from the default list
  // NOTE: this essentially converts UserAPI parameter names into MueLu internal ones
#define MUELU_TEST_AND_SET_PARAM(listWrite, varNameWrite, paramList, defaultList, varNameRead, T) \
  if      (paramList.isParameter(varNameRead))   listWrite.set(varNameWrite, paramList.get<T>(varNameRead)); \
  else if (defaultList.isParameter(varNameRead)) listWrite.set(varNameWrite, defaultList.get<T>(varNameRead));

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::EasyParameterListInterpreter(Teuchos::ParameterList& paramList) {
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::EasyParameterListInterpreter(const std::string& xmlFileName, const Teuchos::Comm<int>& comm) {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), comm);
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameterList(const Teuchos::ParameterList& paramList_) {
    // Create a non const copy of the parameter list
    // Working with a modifiable list is much much easier than with original one
    ParameterList paramList = paramList_;

    // Translate cycle type parameter
    Cycle_ = Hierarchy::GetDefaultCycle();
    if (paramList.isParameter("cycle type")) {
      std::map<std::string,CycleType> cycleMap;
      cycleMap["V"] = VCYCLE;
      cycleMap["W"] = WCYCLE;

      std::string cycleType = paramList.get<std::string>("cycle type");
      TEUCHOS_TEST_FOR_EXCEPTION(cycleMap.count(cycleType) == 0, Exceptions::RuntimeError, "Invalid cycle type: \"" << cycleType << "\"");
      Cycle_ = cycleMap[cycleType];
    }

    this->maxCoarseSize_    = paramList.get<int>("coarse: max size",    Hierarchy::GetDefaultMaxCoarseSize());
    this->numDesiredLevel_  = paramList.get<int>("max levels",          Hierarchy::GetDefaultMaxLevels());
    this->graphOutputLevel_ = paramList.get<int>("debug: graph level", -1);
    this->blockSize_        = paramList.get<int>("number of equations", 1);

    // Translate verbosity parameter
    this->verbosity_ = static_cast<MsgType>(Hierarchy::GetDefaultVerbLevel());      // cast int to enum
    if (paramList.isParameter("verbosity")) {
      std::map<std::string,MsgType> verbMap;
      verbMap["none"]    = None;
      verbMap["low"]     = Low;
      verbMap["medium"]  = Medium;
      verbMap["high"]    = High;
      verbMap["extreme"] = Extreme;

      std::string verbosityLevel = paramList.get<std::string>("verbosity");
      TEUCHOS_TEST_FOR_EXCEPTION(verbMap.count(verbosityLevel) == 0, Exceptions::RuntimeError, "Invalid verbosity level: \"" << verbosityLevel << "\"");
      this->verbosity_ = verbMap[verbosityLevel];
    }

    // Create default manager
    RCP<FactoryManager> defaultManager = rcp(new FactoryManager());
    UpdateFactoryManager(paramList, ParameterList(), *defaultManager, defaultManager);
    defaultManager->Print();

    for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
      RCP<FactoryManager> levelManager;

      try {
        // Some level specific parameters, update default manager
        bool mustAlreadyExist = true;
        ParameterList& levelList = paramList.sublist("level " + toString(levelID), mustAlreadyExist);
        UpdateFactoryManager(levelList, paramList, *defaultManager, levelManager);

      } catch(std::exception) {
        // No level specific parameter, use default manager
        levelManager = defaultManager;
      }

      this->AddFactoryManager(levelID, 1, levelManager);
    }
    // FIXME: parameters passed to packages, like Ifpack2, are not touched by us, resulting in "[unused]" flag
    // being displayed. On the other hand, we don't want to simply iterate through them touching. I don't know
    // what a good solution looks like
    this->GetOStream(Runtime1, 0) << paramList << std::endl;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupMatrix(Matrix& A) const {
    if (A.GetFixedBlockSize() != blockSize_)
      this->GetOStream(Warnings0,  0) << "Warning: setting matrix block size to " << blockSize_ << " (value of \"number of equations\" parameter in the list) "
          << "instead of " << A.GetFixedBlockSize() << " (provided matrix)." << std::endl;
    A.SetFixedBlockSize(blockSize_);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupHierarchy(Hierarchy& H) const {
    HierarchyManager::SetupHierarchy(H);
    H.SetCycle(Cycle_);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UpdateFactoryManager(Teuchos::ParameterList& paramList,
        const Teuchos::ParameterList& defaultList, const FactoryManager& managerIn, RCP<FactoryManager>& manager) {
    manager = rcp(new FactoryManager(managerIn));

    // NOTE: Factory::SetParameterList must be called prior to Factory::SetFactory, as
    // SetParameterList sets default values for non mentioned parameters, including factories

    // === Smoothing ===
    if (paramList.isParameter("smoother: type")) {
      // FIXME: get default values from the factory
      // NOTE: none of the smoothers at the moment use parameter validation framework, so we
      // cannot get the default values from it.
      manager->SetFactory("Smoother", rcp(new SmootherFactory(rcp(new TrilinosSmoother(paramList.get<std::string>("smoother: type",    ""),
                                                                                       paramList.sublist         ("smoother: params",  false),
                                                                                       paramList.get<int>        ("smoother: overlap", 0))))));
    }

    // === Aggregation ===
    // Aggregation graph
    RCP<CoalesceDropFactory> dropFactory = rcp(new CoalesceDropFactory());
    ParameterList dropParams = *(dropFactory->GetValidParameterList());
    dropParams.set                      ("lightweight wrap", true);
    MUELU_TEST_AND_SET_PARAM(dropParams, "algorithm",                     paramList, defaultList, "aggregation: drop scheme",         std::string);
    MUELU_TEST_AND_SET_PARAM(dropParams, "aggregation threshold",         paramList, defaultList, "aggregation: drop tol",            double);
    MUELU_TEST_AND_SET_PARAM(dropParams, "Dirichlet detection threshold", paramList, defaultList, "aggregation: Dirichlet threshold", double);
    if (paramList.isParameter("aggregation: Dirichlet detection"))
      dropParams.set("disable Dirichlet detection", !paramList.get<bool>("aggregation: Dirichlet detection"));

    dropFactory->SetParameterList(dropParams);
    manager->SetFactory("Graph", dropFactory);

    // Aggregation sheme
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "aggregation: type", std::string, "uncoupled", aggType);
    RCP<Factory> aggFactory;
    if      (aggType == "uncoupled") aggFactory = rcp(new UncoupledAggregationFactory());
    else if (aggType == "coupled")   aggFactory = rcp(new CoupledAggregationFactory());
    aggFactory->SetFactory("Graph",       manager->GetFactory("Graph"));
    aggFactory->SetFactory("DofsPerNode", manager->GetFactory("Graph"));
    manager->SetFactory("Aggregates", aggFactory);

    // Coarse map
    RCP<CoarseMapFactory> coarseMap = rcp(new CoarseMapFactory());
    coarseMap->SetFactory("Aggregates", manager->GetFactory("Aggregates"));
    manager->SetFactory("CoarseMap", coarseMap);

    // Tentative P
    RCP<TentativePFactory> Ptent = rcp(new TentativePFactory());
    Ptent->SetFactory("Aggregates", manager->GetFactory("Aggregates"));
    Ptent->SetFactory("CoarseMap",  manager->GetFactory("CoarseMap"));
    manager->SetFactory("Ptent",     Ptent);
    manager->SetFactory("Nullspace", Ptent);

    // Coordinates
    RCP<CoordinatesTransferFactory> coords = rcp(new CoordinatesTransferFactory());
    coords->SetFactory("Aggregates", manager->GetFactory("Aggregates"));
    coords->SetFactory("CoarseMap",  manager->GetFactory("CoarseMap"));
    manager->SetFactory("Coordinates", coords);

    // === Prolongation ===
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "multigrid algorithm", std::string, "sa", multigridAlgo);
    if (multigridAlgo == "sa") {
      // Smoothed aggregation
      RCP<SaPFactory> P = rcp(new SaPFactory());
      ParameterList Pparams = *(P->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(Pparams, "Damping factor", paramList, defaultList, "sa: damping factor", double);
      P->SetParameterList(Pparams);

      if (paramList.isParameter("sa: use filtered matrix") && paramList.get<bool>("sa: use filtered matrix")) {
        // Filtering
        RCP<FilteredAFactory> filterFactory = rcp(new FilteredAFactory());
        filterFactory->SetFactory("Graph", manager->GetFactory("Graph"));
        P->SetFactory("A", filterFactory);
      }

      P->SetFactory("P", manager->GetFactory("Ptent"));
      manager->SetFactory("P", P);

    } else if (multigridAlgo == "emin") {
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "emin: pattern", std::string, "AkPtent", patternType);
      TEUCHOS_TEST_FOR_EXCEPTION(patternType != "AkPtent", Exceptions::InvalidArgument, "Invalid pattern name: \"" << patternType << "\". Valid options: \"AkPtent\"");

      // Pattern
      RCP<PatternFactory> patternFactory = rcp(new PatternFactory());
      ParameterList patternParams = *(patternFactory->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(patternParams, "k", paramList, defaultList, "emin: pattern order", int);
      patternFactory->SetParameterList(patternParams);
      patternFactory->SetFactory("P", manager->GetFactory("Ptent"));
      manager->SetFactory("Ppattern", patternFactory);

      // Constraint
      RCP<ConstraintFactory> constraintFactory = rcp(new ConstraintFactory());
      constraintFactory->SetFactory("Ppattern",        manager->GetFactory("Ppattern"));
      constraintFactory->SetFactory("CoarseNullspace", manager->GetFactory("Ptent"));
      manager->SetFactory("Constraint", constraintFactory);

      // Energy minimization
      RCP<EminPFactory> P = rcp(new EminPFactory());
      ParameterList Pparams = *(P->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(Pparams, "Niterations", paramList, defaultList, "emin: num iterations", int);
      P->SetParameterList(Pparams);
      P->SetFactory("P",          manager->GetFactory("Ptent"));
      P->SetFactory("Constraint", manager->GetFactory("Constraint"));
      manager->SetFactory("P", P);

    } else if (multigridAlgo == "pg") {
      // Petrov-Galerkin
      RCP<PgPFactory> P = rcp(new PgPFactory());
      P->SetFactory("P", manager->GetFactory("Ptent"));
      manager->SetFactory("P", P);
    }

    // === Restriction ===
    RCP<TransPFactory> R = rcp(new TransPFactory());
    R->SetFactory("P", manager->GetFactory("P"));
    manager->SetFactory("R", R);

    // === RAP ===
    RCP<RAPFactory> RAP = rcp(new RAPFactory());
    RAP->SetFactory("P", manager->GetFactory("P"));
    RAP->SetFactory("R", manager->GetFactory("R"));
    RAP->AddTransferFactory(manager->GetFactory("Coordinates"));
    manager->SetFactory("A", RAP);

    // === Repartitioning ===
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "repartition: enable", bool, false, enableRepart);
    if (enableRepart) {
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "repartition: partitioner", std::string, "zoltan", partName);
      TEUCHOS_TEST_FOR_EXCEPTION(partName != "zoltan" && partName != "zoltan2", Exceptions::InvalidArgument,
                                 "Invalid partitioner name: \"" << partName << "\". Valid options: \"zoltan\", \"zoltan2\"");
      // Partitioner
      RCP<Factory> partitioner;
      if (partName == "zoltan") {
        partitioner = rcp(new ZoltanInterface());
        // NOTE: ZoltanInteface ("zoltan") does not support external parameters through ParameterList
      } else if (partName == "zoltan2") {
        partitioner = rcp(new Zoltan2Interface());
        ParameterList partParams = *(partitioner->GetValidParameterList());
        RCP<const ParameterList> partpartParams = rcp(new ParameterList(paramList.sublist("repartition: params", false)));
        partParams.set("ParameterList", partpartParams);
        partitioner->SetParameterList(partParams);
      }
      partitioner->SetFactory("A",           manager->GetFactory("A"));
      partitioner->SetFactory("Coordinates", manager->GetFactory("Coordinates"));
      manager->SetFactory("Partition", partitioner);

      // Repartitioner
      RCP<RepartitionFactory> repartFactory = rcp(new RepartitionFactory());
      ParameterList repartParams = *(repartFactory->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(repartParams, "startLevel",          paramList, defaultList, "repartition: start level",       int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "startLevel",          paramList, defaultList, "repartition: start level",       int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "minRowsPerProcessor", paramList, defaultList, "repartition: min rows per proc", int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "nonzeroImbalance",    paramList, defaultList, "repartition: max imbalance",     double);
      MUELU_TEST_AND_SET_PARAM(repartParams, "remapPartitions",     paramList, defaultList, "repartition: remap parts",       bool);
      repartFactory->SetParameterList(repartParams);
      repartFactory->SetFactory("A",         manager->GetFactory("A"));
      repartFactory->SetFactory("Partition", manager->GetFactory("Partition"));
      manager->SetFactory("Importer", repartFactory);

      // Rebalanced A
      RCP<RebalanceAcFactory>       newA = rcp(new RebalanceAcFactory());
      newA->SetFactory("A",         manager->GetFactory("A"));
      newA->SetFactory("Importer",  manager->GetFactory("Importer"));
      manager->SetFactory("A", newA);

      // Rebalanced P
      RCP<RebalanceTransferFactory> newP = rcp(new RebalanceTransferFactory());
      ParameterList newPparams;
      newPparams.set("type", "Interpolation");
      newP->SetParameterList(newPparams);
      newP->SetFactory("Importer",    manager->GetFactory("Importer"));
      newP->SetFactory("P",           manager->GetFactory("P"));
      manager->SetFactory("P", newP);

      // Rebalanced R
      RCP<RebalanceTransferFactory> newR = rcp(new RebalanceTransferFactory());
      ParameterList newRparams;
      newRparams.set("type", "Restriction");
      newR->SetParameterList(newRparams);
      newR->SetFactory("Importer",    manager->GetFactory("Importer"));
      newR->SetFactory("R",           manager->GetFactory("R"));
      newR->SetFactory("Nullspace",   manager->GetFactory("Nullspace"));
      newR->SetFactory("Coordinates", manager->GetFactory("Coordinates"));
      manager->SetFactory("R",           newR);
      manager->SetFactory("Coordinates", newR);
      manager->SetFactory("Nullspace",   newR);
    }
  }
#undef MUELU_READ_2LIST_PARAM
#undef MUELU_TEST_AND_SET_PARAM

} // namespace MueLu

#define MUELU_EASYPARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_EASYPARAMETERLISTINTERPRETER_DEF_HPP */

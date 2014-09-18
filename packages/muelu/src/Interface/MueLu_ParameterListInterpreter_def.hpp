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
#ifndef MUELU_PARAMETERLISTINTERPRETER_DEF_HPP
#define MUELU_PARAMETERLISTINTERPRETER_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_ParameterListInterpreter_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EminPFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryFactory.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_NullspaceFactory.hpp"
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ParameterListInterpreter(Teuchos::ParameterList& paramList) {
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ParameterListInterpreter(const std::string& xmlFileName, const Teuchos::Comm<int>& comm) {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), comm);
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
    Cycle_     = Hierarchy::GetDefaultCycle();
    blockSize_ = 1;

    if (paramList.isSublist("Hierarchy")) {
      SetFactoryParameterList(paramList);

    } else {
      Validate(paramList);
      SetEasyParameterList(paramList);
    }
  }

  // =====================================================================================================
  // ====================================== EASY interpreter =============================================
  // =====================================================================================================
  //! Helper functions to compare two paramter lists
  static inline bool areSame(const ParameterList& list1, const ParameterList& list2);

  // This macro is tricky. The use case is when we do not have a level specific parameter, so we
  // need to take the default value from the general list for all levels, or set it to the default.
#define MUELU_READ_2LIST_PARAM(paramList, defaultList, paramStr, varType, defaultValue, varName) \
  varType varName; \
  if      (paramList.isParameter(paramStr))   varName = paramList.get<varType>(paramStr); \
  else if (defaultList.isParameter(paramStr)) varName = defaultList.get<varType>(paramStr); \
  else                                        varName = defaultValue;

  // This macro check whether the variable is in the list.
  // If it is, it copies its value to the second list, possibly with a new name
  // Similar to the above macro, we all try to take a value from the default list
  // NOTE: this essentially converts UserAPI parameter names into MueLu internal ones
#define MUELU_TEST_AND_SET_PARAM(listWrite, varName, paramList, defaultList, T) \
  if      (paramList.isParameter(varName)) { \
    try { \
      listWrite.set(varName, paramList.get<T>(varName)); \
    } \
    catch(Teuchos::Exceptions::InvalidParameterType) { \
      TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true, Teuchos::Exceptions::InvalidParameterType, \
                                          "Error: parameter \"" << varName << "\" must be of type " << Teuchos::TypeNameTraits<T>::name()); \
    } \
  } \
  else if (defaultList.isParameter(varName)) listWrite.set(varName, defaultList.get<T>(varName));

#define MUELU_TEST_AND_SET_VAR(var, varName, paramList, T) \
  (paramList.isParameter(varName) ? var = paramList.get<T>(varName), true : false)

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetEasyParameterList(const Teuchos::ParameterList& constParamList) {
    // Create a non const copy of the parameter list
    // Working with a modifiable list is much much easier than with original one
    ParameterList paramList = constParamList;

    // Translate cycle type parameter
    if (paramList.isParameter("cycle type")) {
      std::map<std::string,CycleType> cycleMap;
      cycleMap["V"] = VCYCLE;
      cycleMap["W"] = WCYCLE;

      std::string cycleType = paramList.get<std::string>("cycle type");
      TEUCHOS_TEST_FOR_EXCEPTION(cycleMap.count(cycleType) == 0, Exceptions::RuntimeError, "Invalid cycle type: \"" << cycleType << "\"");
      Cycle_ = cycleMap[cycleType];
    }

    this->maxCoarseSize_    = paramList.get<int> ("coarse: max size",    MasterList::getDefault<int>("coarse: max size"));
    this->numDesiredLevel_  = paramList.get<int> ("max levels",          MasterList::getDefault<int>("max levels"));
    blockSize_              = paramList.get<int> ("number of equations", 1);

    (void)MUELU_TEST_AND_SET_VAR(this->graphOutputLevel_, "debug: graph level", paramList, int);

    // Save level data
    if (paramList.isSublist("export data")) {
      ParameterList printList = paramList.sublist("export data");

      if (printList.isParameter("A"))
        this->matricesToPrint_     = Teuchos::getArrayFromStringParameter<int>(printList, "A");
      if (printList.isParameter("P"))
        this->prolongatorsToPrint_ = Teuchos::getArrayFromStringParameter<int>(printList, "P");
      if (printList.isParameter("R"))
        this->restrictorsToPrint_  = Teuchos::getArrayFromStringParameter<int>(printList, "R");
    }

    // Translate verbosity parameter
    this->verbosity_ = static_cast<MsgType>(Hierarchy::GetDefaultVerbLevel());      // cast int to enum
    if (paramList.isParameter("verbosity")) {
      std::map<std::string,MsgType> verbMap;
      verbMap["none"]    = None;
      verbMap["low"]     = Low;
      verbMap["medium"]  = Medium;
      verbMap["high"]    = High;
      verbMap["extreme"] = Extreme;
      verbMap["test"]    = Test;

      std::string verbosityLevel = paramList.get<std::string>("verbosity");
      TEUCHOS_TEST_FOR_EXCEPTION(verbMap.count(verbosityLevel) == 0, Exceptions::RuntimeError, "Invalid verbosity level: \"" << verbosityLevel << "\"");
      this->verbosity_ = verbMap[verbosityLevel];
      this->SetVerbLevel(this->verbosity_);
    }

    // Detect if we need to transfer coordinates to coarse levels. We do that iff
    //  - we use "distance laplacian" dropping on some level, or
    //  - we use repartitioning on some level
    // This is not ideal, as we may have "repartition: enable" turned on by default
    // and not present in the list, but it is better than nothing.
    useCoordinates_ = false;
    if ((paramList.isParameter("repartition: enable")      && paramList.get<bool>("repartition: enable")             == true) ||
        (paramList.isParameter("aggregation: drop scheme") && paramList.get<std::string>("aggregation: drop scheme") == "distance laplacian")) {
      useCoordinates_ = true;

    } else {
      for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
        std::string levelStr = "level" + toString(levelID);

        if (paramList.isSublist(levelStr)) {
          const ParameterList& levelList = paramList.sublist(levelStr);

          if ((levelList.isParameter("repartition: enable")      && levelList.get<bool>("repartition: enable")             == true) ||
              (levelList.isParameter("aggregation: drop scheme") && levelList.get<std::string>("aggregation: drop scheme") == "distance laplacian")) {
            useCoordinates_ = true;
            break;
          }
        }
      }
    }

    // Detect if we do implicit P and R rebalance
    changedPRrebalance_ = false;
    if (paramList.isParameter("repartition: enable") && paramList.get<bool>("repartition: enable") == true)
      changedPRrebalance_ = MUELU_TEST_AND_SET_VAR(this->doPRrebalance_, "repartition: rebalance P and R", paramList, bool);

    // Detect if we use implicit transpose
    changedImplicitTranspose_ = MUELU_TEST_AND_SET_VAR(this->implicitTranspose_, "transpose: use implicit", paramList, bool);

    // Create default manager
    RCP<FactoryManager> defaultManager = rcp(new FactoryManager());
    defaultManager->SetVerbLevel(this->verbosity_);
    UpdateFactoryManager(paramList, ParameterList(), *defaultManager);
    defaultManager->Print();

    for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
      RCP<FactoryManager> levelManager;

      if (paramList.isSublist("level " + toString(levelID))) {
        // Some level specific parameters, update default manager
        bool mustAlreadyExist = true;
        ParameterList& levelList = paramList.sublist("level " + toString(levelID), mustAlreadyExist);

        levelManager = rcp(new FactoryManager(*defaultManager));
        levelManager->SetVerbLevel(defaultManager->GetVerbLevel());

        UpdateFactoryManager(levelList, paramList, *levelManager);

      } else {
        // No level specific parameter, use default manager
        levelManager = defaultManager;
      }

      this->AddFactoryManager(levelID, 1, levelManager);
    }

    // FIXME: parameters passed to packages, like Ifpack2, are not touched by us, resulting in "[unused]" flag
    // being displayed. On the other hand, we don't want to simply iterate through them touching. I don't know
    // what a good solution looks like
    if (!paramList.isParameter("print initial parameters") || paramList.get<bool>("print initial parameters") == true)
      this->GetOStream(static_cast<MsgType>(Runtime1), 0) << paramList << std::endl;

    if (!paramList.isParameter("print unused parameters") || paramList.get<bool>("print unused parameters") == true) {
      // Check unused parameters
      ParameterList unusedParamList;

      // Check for unused parameters that aren't lists
      for (ParameterList::ConstIterator it = paramList.begin(); it != paramList.end(); it++) {
        const ParameterEntry& entry = paramList.entry(it);

        if (!entry.isList() && !entry.isUsed())
          unusedParamList.setEntry(paramList.name(it), entry);
      }

      // Check for unused parameters in level-specific sublists
      for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
        std::string levelStr = "level " + toString(levelID);

        if (paramList.isSublist(levelStr)) {
          const ParameterList& levelList = paramList.sublist(levelStr);

          for (ParameterList::ConstIterator itr = levelList.begin(); itr != levelList.end(); ++itr) {
            const ParameterEntry& entry = levelList.entry(itr);

            if (!entry.isList() && !entry.isUsed())
              unusedParamList.sublist(levelStr).setEntry(levelList.name(itr), entry);
          }
        }
      }

      if (unusedParamList.numParams() > 0) {
        std::ostringstream unusedParamsStream;
        int indent = 4;
        unusedParamList.print(unusedParamsStream, indent);

        this->GetOStream(Warnings1) << "The following parameters were not used:\n" << unusedParamsStream.str() << std::endl;
      }
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::UpdateFactoryManager(Teuchos::ParameterList& paramList,
        const Teuchos::ParameterList& defaultList, FactoryManager& manager) const {
    // NOTE: Factory::SetParameterList must be called prior to Factory::SetFactory, as
    // SetParameterList sets default values for non mentioned parameters, including factories

    // === Smoothing ===
    bool isCustomSmoother =
        paramList.isParameter("smoother: pre or post") ||
        paramList.isParameter("smoother: type")    || paramList.isParameter("smoother: pre type")    || paramList.isParameter("smoother: post type")   ||
        paramList.isSublist  ("smoother: params")  || paramList.isSublist  ("smoother: pre params")  || paramList.isSublist  ("smoother: post params") ||
        paramList.isParameter("smoother: sweeps")  || paramList.isParameter("smoother: pre sweeps")  || paramList.isParameter("smoother: post sweeps") ||
        paramList.isParameter("smoother: overlap") || paramList.isParameter("smoother: pre overlap") || paramList.isParameter("smoother: post overlap");;
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "smoother: pre or post", std::string, "both", PreOrPost);
    if (PreOrPost == "none") {
      manager.SetFactory("Smoother", Teuchos::null);

    } else if (isCustomSmoother) {
      // FIXME: get default values from the factory
      // NOTE: none of the smoothers at the moment use parameter validation framework, so we
      // cannot get the default values from it.
#define TEST_MUTUALLY_EXCLUSIVE(arg1,arg2) \
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter(#arg1) && paramList.isParameter(#arg2), \
                                 Exceptions::InvalidArgument, "You cannot specify both \""#arg1"\" and \""#arg2"\"");
#define TEST_MUTUALLY_EXCLUSIVE_S(arg1,arg2) \
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist(#arg1) && paramList.isSublist(#arg2), \
                                 Exceptions::InvalidArgument, "You cannot specify both \""#arg1"\" and \""#arg2"\"");

      TEST_MUTUALLY_EXCLUSIVE  ("smoother: type",    "smoother: pre type");
      TEST_MUTUALLY_EXCLUSIVE  ("smoother: type",    "smoother: post type");
      TEST_MUTUALLY_EXCLUSIVE  ("smoother: sweeps",  "smoother: pre sweeps");
      TEST_MUTUALLY_EXCLUSIVE  ("smoother: sweeps",  "smoother: post sweeps");
      TEST_MUTUALLY_EXCLUSIVE  ("smoother: overlap", "smoother: pre overlap");
      TEST_MUTUALLY_EXCLUSIVE  ("smoother: overlap", "smoother: post overlap");
      TEST_MUTUALLY_EXCLUSIVE_S("smoother: params",  "smoother: pre params");
      TEST_MUTUALLY_EXCLUSIVE_S("smoother: params",  "smoother: post params");
      TEUCHOS_TEST_FOR_EXCEPTION(PreOrPost == "both" && (paramList.isParameter("smoother: pre type") != paramList.isParameter("smoother: post type")),
                                 Exceptions::InvalidArgument, "You must specify both \"smoother: pre type\" and \"smoother: post type\"");

      // Default values
      int overlap = 0;
      ParameterList defaultSmootherParams;
      defaultSmootherParams.set("relaxation: type",           "Symmetric Gauss-Seidel");
      defaultSmootherParams.set("relaxation: sweeps",         Teuchos::OrdinalTraits<LO>::one());
      defaultSmootherParams.set("relaxation: damping factor", Teuchos::ScalarTraits<Scalar>::one());

      RCP<SmootherPrototype> preSmoother = Teuchos::null, postSmoother = Teuchos::null;
      std::string            preSmootherType,             postSmootherType;
      ParameterList          preSmootherParams,           postSmootherParams;

      if (paramList.isParameter("smoother: overlap"))
        overlap = paramList.get<int>("smoother: overlap");

      if (PreOrPost == "pre" || PreOrPost == "both") {
        if (paramList.isParameter("smoother: pre type")) {
          preSmootherType = paramList.get<std::string>("smoother: pre type");
        } else {
          MUELU_READ_2LIST_PARAM(paramList, defaultList, "smoother: type", std::string, "RELAXATION", preSmootherTypeTmp);
          preSmootherType = preSmootherTypeTmp;
        }
        if (paramList.isParameter("smoother: pre overlap"))
          overlap = paramList.get<int>("smoother: pre overlap");

        if (paramList.isSublist("smoother: pre params"))
          preSmootherParams = paramList.sublist("smoother: pre params");
        else if (paramList.isSublist("smoother: params"))
          preSmootherParams = paramList.sublist("smoother: params");
        else if (defaultList.isSublist("smoother: params"))
          preSmootherParams = defaultList.sublist("smoother: params");
        else if (preSmootherType == "RELAXATION")
          preSmootherParams = defaultSmootherParams;

        preSmoother = rcp(new TrilinosSmoother(preSmootherType, preSmootherParams, overlap));
      }

      if (PreOrPost == "post" || PreOrPost == "both") {
        if (paramList.isParameter("smoother: post type"))
          postSmootherType = paramList.get<std::string>("smoother: post type");
        else {
          MUELU_READ_2LIST_PARAM(paramList, defaultList, "smoother: type", std::string, "RELAXATION", postSmootherTypeTmp);
          postSmootherType = postSmootherTypeTmp;
        }

        if (paramList.isSublist("smoother: post params"))
          postSmootherParams = paramList.sublist("smoother: post params");
        else if (paramList.isSublist("smoother: params"))
          postSmootherParams = paramList.sublist("smoother: params");
        else if (defaultList.isSublist("smoother: params"))
          postSmootherParams = defaultList.sublist("smoother: params");
        else if (postSmootherType == "RELAXATION")
          postSmootherParams = defaultSmootherParams;
        if (paramList.isParameter("smoother: post overlap"))
          overlap = paramList.get<int>("smoother: post overlap");

        if (postSmootherType == preSmootherType && areSame(preSmootherParams, postSmootherParams))
          postSmoother = preSmoother;
        else
          postSmoother = rcp(new TrilinosSmoother(postSmootherType, postSmootherParams, overlap));
      }

      manager.SetFactory("Smoother", rcp(new SmootherFactory(preSmoother, postSmoother)));
    }

    // === Coarse solver ===
    bool isCustomCoarseSolver =
        paramList.isParameter("coarse: type")   ||
        paramList.isParameter("coarse: params");
    if (paramList.isParameter("coarse: type") && paramList.get<std::string>("coarse: type") == "none") {
      manager.SetFactory("CoarseSolver", Teuchos::null);

    } else if (isCustomCoarseSolver) {
      // FIXME: get default values from the factory
      // NOTE: none of the smoothers at the moment use parameter validation framework, so we
      // cannot get the default values from it.
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "coarse: type", std::string, "", coarseType);

      int overlap = 0;
      if (paramList.isParameter("coarse: overlap"))
        overlap = paramList.get<int>("coarse: overlap");

      ParameterList coarseParams;
      if (paramList.isSublist("coarse: params"))
        coarseParams = paramList.sublist("coarse: params");
      else if (defaultList.isSublist("coarse: params"))
        coarseParams = defaultList.sublist("coarse: params");

      RCP<SmootherPrototype> coarseSmoother;
      // TODO: this is not a proper place to check. If we consider direct solver to be a special
      // case of smoother, we would like to unify Amesos and Ifpack2 smoothers in src/Smoothers, and
      // have a single factory responsible for those. Then, this check would belong there.
      if (coarseType == "RELAXATION" || coarseType == "CHEBYSHEV" ||
          coarseType == "ILUT" || coarseType == "ILU" || coarseType == "RILUK" || coarseType == "SCHWARZ" ||
          coarseType == "Amesos")
        coarseSmoother = rcp(new TrilinosSmoother(coarseType, coarseParams, overlap));
      else
        coarseSmoother = rcp(new DirectSolver(coarseType, coarseParams));

      manager.SetFactory("CoarseSolver", rcp(new SmootherFactory(coarseSmoother, Teuchos::null)));
    }

    // === Aggregation ===
    // Aggregation graph
    RCP<CoalesceDropFactory> dropFactory = rcp(new CoalesceDropFactory());
    ParameterList dropParams;
    dropParams.set("lightweight wrap", true);
    MUELU_TEST_AND_SET_PARAM(dropParams, "aggregation: drop scheme",         paramList, defaultList, std::string);
    // Rename classical to original
    if (dropParams.isParameter("aggregation: drop scheme") && dropParams.get<std::string>("aggregation: drop scheme") == "classical")
      dropParams.set("aggregation: drop scheme", "original");
    MUELU_TEST_AND_SET_PARAM(dropParams, "aggregation: drop tol",            paramList, defaultList, double);
    MUELU_TEST_AND_SET_PARAM(dropParams, "aggregation: Dirichlet threshold", paramList, defaultList, double);

    dropFactory->SetParameterList(dropParams);
    manager.SetFactory("Graph", dropFactory);

    // Aggregation sheme
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "aggregation: type", std::string, "uncoupled", aggType);
    RCP<Factory> aggFactory;
    if      (aggType == "uncoupled") {
      aggFactory = rcp(new UncoupledAggregationFactory());
      ParameterList aggParams;
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: mode",                      paramList, defaultList, std::string);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: ordering",                  paramList, defaultList, std::string);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: min agg size",              paramList, defaultList, int);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: max agg size",              paramList, defaultList, int);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: max selected neighbors",    paramList, defaultList, int);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: enable phase 1",            paramList, defaultList, bool);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: enable phase 2a",           paramList, defaultList, bool);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: enable phase 2b",           paramList, defaultList, bool);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: enable phase 3",            paramList, defaultList, bool);
      MUELU_TEST_AND_SET_PARAM(aggParams, "aggregation: preserve Dirichlet points", paramList, defaultList, bool);
      aggFactory->SetParameterList(aggParams);

    } else if (aggType == "coupled") {
      aggFactory = rcp(new CoupledAggregationFactory());
    }
    aggFactory->SetFactory("Graph",       manager.GetFactory("Graph"));
    aggFactory->SetFactory("DofsPerNode", manager.GetFactory("Graph"));
    manager.SetFactory("Aggregates", aggFactory);

    // Coarse map
    RCP<CoarseMapFactory> coarseMap = rcp(new CoarseMapFactory());
    coarseMap->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
    manager.SetFactory("CoarseMap", coarseMap);

    // Tentative P
    RCP<Factory> Ptent = rcp(new TentativePFactory());
    Ptent->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
    Ptent->SetFactory("CoarseMap",  manager.GetFactory("CoarseMap"));
    manager.SetFactory("Ptent",     Ptent);

    // Nullspace
    RCP<NullspaceFactory> nullSpace = rcp(new NullspaceFactory());
    nullSpace->SetFactory("Nullspace", manager.GetFactory("Ptent"));
    manager.SetFactory("Nullspace", nullSpace);

    // === Prolongation ===
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "multigrid algorithm", std::string, "sa", multigridAlgo);
    if (multigridAlgo == "unsmoothed") {
      manager.SetFactory("P", Ptent);

    } else if (multigridAlgo == "sa") {
      // Smoothed aggregation
      RCP<SaPFactory> P = rcp(new SaPFactory());
      ParameterList Pparams;
      MUELU_TEST_AND_SET_PARAM(Pparams, "sa: damping factor", paramList, defaultList, double);
      P->SetParameterList(Pparams);

      if (paramList.isParameter("sa: use filtered matrix") && paramList.get<bool>("sa: use filtered matrix")) {
        // Filtering
        RCP<FilteredAFactory> filterFactory = rcp(new FilteredAFactory());
        ParameterList fParams;
        MUELU_TEST_AND_SET_PARAM(fParams, "filtered matrix: use lumping", paramList, defaultList, bool);
        filterFactory->SetParameterList(fParams);
        filterFactory->SetFactory("Graph", manager.GetFactory("Graph"));
        P->SetFactory("A", filterFactory);
      }

      P->SetFactory("P", manager.GetFactory("Ptent"));
      manager.SetFactory("P", P);

    } else if (multigridAlgo == "emin") {
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "emin: pattern", std::string, "AkPtent", patternType);
      TEUCHOS_TEST_FOR_EXCEPTION(patternType != "AkPtent", Exceptions::InvalidArgument,
                                 "Invalid pattern name: \"" << patternType << "\". Valid options: \"AkPtent\"");
      // Pattern
      RCP<PatternFactory> patternFactory = rcp(new PatternFactory());
      ParameterList patternParams;
      MUELU_TEST_AND_SET_PARAM(patternParams, "emin: pattern order", paramList, defaultList, int);
      patternFactory->SetParameterList(patternParams);
      patternFactory->SetFactory("P", manager.GetFactory("Ptent"));
      manager.SetFactory("Ppattern", patternFactory);

      // Constraint
      RCP<ConstraintFactory> constraintFactory = rcp(new ConstraintFactory());
      constraintFactory->SetFactory("Ppattern",        manager.GetFactory("Ppattern"));
      constraintFactory->SetFactory("CoarseNullspace", manager.GetFactory("Ptent"));
      manager.SetFactory("Constraint", constraintFactory);

      // Energy minimization
      RCP<EminPFactory> P = rcp(new EminPFactory());
      ParameterList Pparams;
      MUELU_TEST_AND_SET_PARAM(Pparams, "emin: num iterations",   paramList, defaultList, int);
      MUELU_TEST_AND_SET_PARAM(Pparams, "emin: iterative method", paramList, defaultList, std::string);
      P->SetParameterList(Pparams);
      P->SetFactory("P",          manager.GetFactory("Ptent"));
      P->SetFactory("Constraint", manager.GetFactory("Constraint"));
      manager.SetFactory("P", P);

    } else if (multigridAlgo == "pg") {
      // Petrov-Galerkin
      RCP<PgPFactory> P = rcp(new PgPFactory());
      P->SetFactory("P", manager.GetFactory("Ptent"));
      manager.SetFactory("P", P);
    }

    // === Restriction ===
    if (!this->implicitTranspose_) {
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "problem: symmetric", bool, true, isSymmetric);
      if (isSymmetric == false && (multigridAlgo == "unsmoothed" || multigridAlgo == "emin")) {
        this->GetOStream(Warnings0) << "Switching to symmetric problem as multigrid algorithm \"" << multigridAlgo << "\" is restricted to symmetric case" << std::endl;
        isSymmetric = true;
      }

      RCP<Factory> R;
      if (isSymmetric)  R = rcp(new TransPFactory());
      else              R = rcp(new GenericRFactory());

      R->SetFactory("P", manager.GetFactory("P"));
      manager.SetFactory("R", R);

    } else {
      manager.SetFactory("R", Teuchos::null);
    }

    // === RAP ===
    RCP<RAPFactory> RAP = rcp(new RAPFactory());
    ParameterList RAPparams;
    MUELU_TEST_AND_SET_PARAM(RAPparams, "transpose: use implicit", paramList, defaultList, bool);
    RAP->SetParameterList(RAPparams);
    RAP->SetFactory("P", manager.GetFactory("P"));
    if (!this->implicitTranspose_)
      RAP->SetFactory("R", manager.GetFactory("R"));
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "aggregation: export visualization data", bool, false, visAgg);
    if (visAgg) {
      RCP<AggregationExportFactory> aggExport = rcp(new AggregationExportFactory());
      aggExport->SetFactory("DofsPerNode", manager.GetFactory("Graph"));
      RAP->AddTransferFactory(aggExport);
    }
    manager.SetFactory("A", RAP);

    // === Coordinates ===
    if (useCoordinates_) {
      RCP<CoordinatesTransferFactory> coords = rcp(new CoordinatesTransferFactory());
      coords->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
      coords->SetFactory("CoarseMap",  manager.GetFactory("CoarseMap"));
      manager.SetFactory("Coordinates", coords);

      RAP->AddTransferFactory(manager.GetFactory("Coordinates"));
    }

    // === Repartitioning ===
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "repartition: enable", bool, false, enableRepart);
    if (enableRepart) {
#ifdef HAVE_MPI
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "repartition: partitioner", std::string, "zoltan", partName);
      TEUCHOS_TEST_FOR_EXCEPTION(partName != "zoltan" && partName != "zoltan2", Exceptions::InvalidArgument,
                                 "Invalid partitioner name: \"" << partName << "\". Valid options: \"zoltan\", \"zoltan2\"");
      // Partitioner
      RCP<Factory> partitioner;
      if (partName == "zoltan") {
#ifdef HAVE_MUELU_ZOLTAN
        partitioner = rcp(new ZoltanInterface());
        // NOTE: ZoltanInteface ("zoltan") does not support external parameters through ParameterList
#else
        throw Exceptions::RuntimeError("Zoltan interface is not available");
#endif
      } else if (partName == "zoltan2") {
#ifdef HAVE_MUELU_ZOLTAN2
        partitioner = rcp(new Zoltan2Interface());
        ParameterList partParams;
        RCP<const ParameterList> partpartParams = rcp(new ParameterList(paramList.sublist("repartition: params", false)));
        partParams.set("ParameterList", partpartParams);
        partitioner->SetParameterList(partParams);
#else
        throw Exceptions::RuntimeError("Zoltan2 interface is not available");
#endif
      }
      partitioner->SetFactory("A",           manager.GetFactory("A"));
      partitioner->SetFactory("Coordinates", manager.GetFactory("Coordinates"));
      manager.SetFactory("Partition", partitioner);

      // Repartitioner
      RCP<RepartitionFactory> repartFactory = rcp(new RepartitionFactory());
      ParameterList repartParams;
      MUELU_TEST_AND_SET_PARAM(repartParams, "repartition: start level",       paramList, defaultList, int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "repartition: min rows per proc", paramList, defaultList, int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "repartition: max imbalance",     paramList, defaultList, double);
      MUELU_TEST_AND_SET_PARAM(repartParams, "repartition: keep proc 0",       paramList, defaultList, bool);
      MUELU_TEST_AND_SET_PARAM(repartParams, "repartition: print partition distribution", paramList, defaultList, bool);
      MUELU_TEST_AND_SET_PARAM(repartParams, "repartition: remap parts",       paramList, defaultList, bool);
      MUELU_TEST_AND_SET_PARAM(repartParams, "repartition: remap num values",  paramList, defaultList, int);
      repartFactory->SetParameterList(repartParams);
      repartFactory->SetFactory("A",         manager.GetFactory("A"));
      repartFactory->SetFactory("Partition", manager.GetFactory("Partition"));
      manager.SetFactory("Importer", repartFactory);

      // Rebalanced A
      RCP<RebalanceAcFactory> newA = rcp(new RebalanceAcFactory());
      newA->  SetFactory("A",         manager.GetFactory("A"));
      newA->  SetFactory("Importer",  manager.GetFactory("Importer"));
      manager.SetFactory("A",         newA);

      // Rebalanced P
      RCP<RebalanceTransferFactory> newP = rcp(new RebalanceTransferFactory());
      ParameterList newPparams;
      newPparams.set("type",                           "Interpolation");
      if (changedPRrebalance_)
        newPparams.set("repartition: rebalance P and R", this->doPRrebalance_);
      newP->  SetParameterList(newPparams);
      newP->  SetFactory("Importer",    manager.GetFactory("Importer"));
      newP->  SetFactory("P",           manager.GetFactory("P"));
      newP->  SetFactory("Nullspace",   manager.GetFactory("Ptent"));
      newP->  SetFactory("Coordinates", manager.GetFactory("Coordinates"));
      manager.SetFactory("P",           newP);
      manager.SetFactory("Coordinates", newP);

      // Rebalanced R
      RCP<RebalanceTransferFactory> newR = rcp(new RebalanceTransferFactory());
      ParameterList newRparams;
      newRparams.set("type",                           "Restriction");
      if (changedPRrebalance_)
        newRparams.set("repartition: rebalance P and R", this->doPRrebalance_);
      if (changedImplicitTranspose_)
        newRparams.set("transpose: use implicit",        this->implicitTranspose_);
      newR->  SetParameterList(newRparams);
      newR->  SetFactory("Importer",       manager.GetFactory("Importer"));
      if (!this->implicitTranspose_) {
        newR->SetFactory("R",              manager.GetFactory("R"));
        manager.SetFactory("R",            newR);
      }

      // NOTE: the role of NullspaceFactory is to provide nullspace on the finest
      // level if a user does not do that. For all other levels it simply passes
      // nullspace from a real factory to whoever needs it. If we don't use
      // repartitioning, that factory is "TentativePFactory"; if we do, it is
      // "RebalanceTransferFactory". But we still have to have NullspaceFactory as
      // the "Nullspace" of the manager
      nullSpace->SetFactory("Nullspace", newP);
#else
      throw Exceptions::RuntimeError("No repartitioning available for a serial run");
#endif
    }
  }
#undef MUELU_READ_2LIST_PARAM
#undef MUELU_TEST_AND_SET_PARAM

  int LevenshteinDistance(const char* s, size_t len_s, const char* t, size_t len_t);

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Validate(const Teuchos::ParameterList& constParamList) const {
    ParameterList paramList = constParamList;
    const ParameterList& validList = *MasterList::List();

    // Validate up to maxLevels level specific parameter sublists
    const int maxLevels = 100;

    // Extract level specific list
    std::vector<ParameterList> paramLists;
    for (int levelID = 0; levelID < maxLevels; levelID++) {
      std::string sublistName = "level " + toString(levelID);
      if (paramList.isSublist(sublistName)) {
        paramLists.push_back(paramList.sublist(sublistName));
        // paramLists.back().setName(sublistName);
        paramList.remove(sublistName);
      }
    }
    paramLists.push_back(paramList);
    // paramLists.back().setName("main");

    const int maxDepth = 0;
    for (size_t i = 0; i < paramLists.size(); i++) {
      // validate every sublist
      try {
        paramLists[i].validateParameters(validList, maxDepth);
      } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
        std::string eString = e.what();

        // Parse name from: <Error, the parameter {name="smoothe: type",...>
        size_t nameStart = eString.find_first_of('"') + 1;
        size_t nameEnd   = eString.find_first_of('"', nameStart);
        std::string name = eString.substr(nameStart, nameEnd - nameStart);

        int         bestScore = 100;
        std::string bestName  = "";
        for (ParameterList::ConstIterator it = validList.begin(); it != validList.end(); it++) {
          const std::string& pName = validList.name(it);
          this->GetOStream(Runtime1) << "| " << pName;

          int score = LevenshteinDistance(name.c_str(), name.length(), pName.c_str(), pName.length());
          this->GetOStream(Runtime1) << " -> " << score << std::endl;
          if (score < bestScore) {
            bestScore = score;
            bestName  = pName;
          }
        }

        if (bestScore < 10 && bestName != "") {
          TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameterName, eString << "The parameter name \"" + name + "\" is not valid. Did you mean \"" + bestName << "\"?\n");
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameterName, eString << "The parameter name \"" + name + "\" is not valid.\n");
        }

      }
    }
  }

  // =====================================================================================================
  // ==================================== FACTORY interpreter ============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetFactoryParameterList(const Teuchos::ParameterList& constParamList) {
    // Create a non const copy of the parameter list
    // Working with a modifiable list is much much easier than with original one
    ParameterList paramList = constParamList;

    // Parameter List Parsing:
    // ---------
    //   <ParameterList name="MueLu">
    //     <ParameterList name="Matrix">
    //   </ParameterList>
    if (paramList.isSublist("Matrix"))
      blockSize_ = paramList.sublist("Matrix").get<int>("number of equations", 1);

    // Parameter List Parsing:
    // ---------
    //   <ParameterList name="MueLu">
    //     <ParameterList name="Factories"> <== call BuildFactoryMap() on this parameter list
    //     ...
    //     </ParameterList>
    //   </ParameterList>
    FactoryMap factoryMap;
    FactoryManagerMap factoryManagers;
    if (paramList.isSublist("Factories"))
      this->BuildFactoryMap(paramList.sublist("Factories"), factoryMap, factoryMap, factoryManagers);

    // Parameter List Parsing:
    // ---------
    //   <ParameterList name="MueLu">
    //     <ParameterList name="Hierarchy">
    //       <Parameter name="verbose"  type="string" value="Warnings"/> <== get
    //       <Parameter name="numDesiredLevel" type="int" value="10"/>   <== get
    //
    //       <ParameterList name="firstLevel">                           <== parse first args and call BuildFactoryMap() on the rest of this parameter list
    //         ...
    //       </ParameterList>
    //     </ParameterList>
    //   </ParameterList>
    if (paramList.isSublist("Hierarchy")) {
      ParameterList hieraList = paramList.sublist("Hierarchy"); // copy because list temporally modified (remove 'id')

      // Get hierarchy options
      if (hieraList.isParameter("max levels")) {
        this->numDesiredLevel_ = hieraList.get<int>("max levels");
        hieraList.remove("max levels");
      }

      if (hieraList.isParameter("coarse: max size")) {
        this->maxCoarseSize_ = hieraList.get<int>("coarse: max size");
        hieraList.remove("coarse: max size");
      }

      if (hieraList.isParameter("repartition: rebalance P and R")) {
        this->doPRrebalance_ = hieraList.get<bool>("repartition: rebalance P and R");
        hieraList.remove("repartition: rebalance P and R");
      }

      //TODO Move this its own class or MueLu::Utils?
      std::map<std::string,MsgType> verbMap;
      //for developers
      verbMap["Errors"]         = Errors;
      verbMap["Warnings0"]      = Warnings0;
      verbMap["Warnings00"]     = Warnings00;
      verbMap["Warnings1"]      = Warnings1;
      verbMap["PerfWarnings"]   = PerfWarnings;
      verbMap["Runtime0"]       = Runtime0;
      verbMap["Runtime1"]       = Runtime1;
      verbMap["RuntimeTimings"] = RuntimeTimings;
      verbMap["NoTimeReport"]   = NoTimeReport;
      verbMap["Parameters0"]    = Parameters0;
      verbMap["Parameters1"]    = Parameters1;
      verbMap["Statistics0"]    = Statistics0;
      verbMap["Statistics1"]    = Statistics1;
      verbMap["Timings0"]       = Timings0;
      verbMap["Timings1"]       = Timings1;
      verbMap["TimingsByLevel"] = TimingsByLevel;
      verbMap["External"]       = External;
      verbMap["Debug"]          = Debug;
      verbMap["Test"]           = Test;
      //for users and developers
      verbMap["None"]           = None;
      verbMap["Low"]            = Low;
      verbMap["Medium"]         = Medium;
      verbMap["High"]           = High;
      verbMap["Extreme"]        = Extreme;
      if (hieraList.isParameter("verbosity")) {
        std::string vl = hieraList.get<std::string>("verbosity");
        hieraList.remove("verbosity");
        //TODO Move this to its own class or MueLu::Utils?
        if (verbMap.find(vl) != verbMap.end())
          this->verbosity_ = verbMap[vl];
        else
          TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::ParameterListInterpreter():: invalid verbosity level");
      }

      if (hieraList.isParameter("dependencyOutputLevel"))
        this->graphOutputLevel_ = hieraList.get<int>("dependencyOutputLevel");

      // Check for the reuse case
      if (hieraList.isParameter("reuse"))
        Factory::DisableMultipleCheckGlobally();

      if (hieraList.isSublist("DataToWrite")) {
        //TODO We should be able to specify any data.  If it exists, write it.
        //TODO This would requires something like std::set<dataName,Array<int> >
        Teuchos::ParameterList foo = hieraList.sublist("DataToWrite");
        std::string dataName = "Matrices";
        if (foo.isParameter(dataName))
          this->matricesToPrint_ = Teuchos::getArrayFromStringParameter<int>(foo,dataName);
        dataName = "Prolongators";
        if (foo.isParameter(dataName))
          this->prolongatorsToPrint_ = Teuchos::getArrayFromStringParameter<int>(foo,dataName);
        dataName = "Restrictors";
        if (foo.isParameter(dataName))
          this->restrictorsToPrint_ = Teuchos::getArrayFromStringParameter<int>(foo,dataName);
      }

      // Get level configuration
      for (ParameterList::ConstIterator param = hieraList.begin(); param != hieraList.end(); ++param) {
        const std::string & paramName  = hieraList.name(param);

        if (paramName != "DataToWrite" && hieraList.isSublist(paramName)) {
          ParameterList levelList = hieraList.sublist(paramName); // copy because list temporally modified (remove 'id')

          int startLevel = 0;       if(levelList.isParameter("startLevel"))      { startLevel      = levelList.get<int>("startLevel");      levelList.remove("startLevel"); }
          int numDesiredLevel = 1;  if(levelList.isParameter("numDesiredLevel")) { numDesiredLevel = levelList.get<int>("numDesiredLevel"); levelList.remove("numDesiredLevel"); }

          // Parameter List Parsing:
          // ---------
          //   <ParameterList name="firstLevel">
          //      <Parameter name="startLevel"       type="int" value="0"/>
          //      <Parameter name="numDesiredLevel"  type="int" value="1"/>
          //      <Parameter name="verbose"          type="string" value="Warnings"/>
          //
          //      [] <== call BuildFactoryMap() on the rest of the parameter list
          //
          //  </ParameterList>
          FactoryMap levelFactoryMap;
          BuildFactoryMap(levelList, factoryMap, levelFactoryMap, factoryManagers);

          RCP<FactoryManagerBase> m = rcp(new FactoryManager(levelFactoryMap));

          if (startLevel >= 0)
            this->AddFactoryManager(startLevel, numDesiredLevel, m);
          else
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::ParameterListInterpreter():: invalid level id");
        } /* TODO: else { } */
      }
    }
  }

  // Parameter List Parsing:
  // Create an entry in factoryMap for each parameter of the list paramList
  // ---------
  //   <ParameterList name="...">
  //     <Parameter name="smootherFact0" type="string" value="TrilinosSmoother"/>
  //
  //     <ParameterList name="smootherFact1">
  //       <Parameter name="type" type="string" value="TrilinosSmoother"/>
  //       ...
  //     </ParameterList>
  //    </ParameterList>
  //
  //TODO: static?
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildFactoryMap(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, FactoryMap& factoryMapOut, FactoryManagerMap& factoryManagers) const {
    for (Teuchos::ParameterList::ConstIterator param = paramList.begin(); param != paramList.end(); ++param) {
      const std::string             & paramName  = paramList.name(param);
      const Teuchos::ParameterEntry & paramValue = paramList.entry(param);

      //TODO: do not allow name of existing MueLu classes (can be tested using FactoryFactory)

      // TODO: add support for "factory groups" which are stored in a map.
      // A factory group has a name and a list of factories

      if (paramValue.isList()) {
        Teuchos::ParameterList paramList1 = Teuchos::getValue<Teuchos::ParameterList>(paramValue);
        if (paramList1.isParameter("factory")) { // default: just a factory definition
          factoryMapOut[paramName] = FactoryFactory().BuildFactory(paramValue, factoryMapIn, factoryManagers);

        } else if (paramList1.isParameter("group")) { // definitiion of a factory group (for a factory manager)
          std::string groupType = paramList1.get<std::string>("group");
          TEUCHOS_TEST_FOR_EXCEPTION(groupType!="FactoryManager", Exceptions::RuntimeError, "group must be of type \"FactoryManager\".");

          Teuchos::ParameterList groupList = paramList1; // copy because list temporally modified (remove 'id')
          groupList.remove("group");

          FactoryMap groupFactoryMap;
          BuildFactoryMap(groupList, factoryMapIn, groupFactoryMap, factoryManagers);

          // do not store groupFactoryMap in factoryMapOut
          // Create a factory manager object from groupFactoryMap
          RCP<FactoryManagerBase> m = rcp(new FactoryManager(groupFactoryMap));

          factoryManagers[paramName] = m;

        } else {
          this->GetOStream(Warnings0) << "Could not interpret parameter list " << paramList1 << std::endl;
          TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "XML Parameter list must either be of type \"factory\" or of type \"group\".");
        }
      } else {
        // default: just a factory (no parameter list)
        factoryMapOut[paramName] = FactoryFactory().BuildFactory(paramValue, factoryMapIn, factoryManagers);
      }
    }
  }

  // =====================================================================================================
  // ======================================= MISC functions ==============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupMatrix(Matrix& A) const {
    if (A.GetFixedBlockSize() != blockSize_)
      this->GetOStream(Warnings0) << "Setting matrix block size to " << blockSize_ << " (value of the parameter in the list) "
          << "instead of " << A.GetFixedBlockSize() << " (provided matrix)." << std::endl;

    A.SetFixedBlockSize(blockSize_);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupHierarchy(Hierarchy& H) const {
    HierarchyManager::SetupHierarchy(H);
    H.SetCycle(Cycle_);
  }


  static bool compare(const ParameterList& list1, const ParameterList& list2) {
    // First loop through and validate the parameters at this level.
    // In addition, we generate a list of sublists that we will search next
    for (ParameterList::ConstIterator it = list1.begin(); it != list1.end(); it++) {
      const std::string&             name   = it->first;
      const Teuchos::ParameterEntry& entry1 = it->second;

      const Teuchos::ParameterEntry *entry2 = list2.getEntryPtr(name);
      if (!entry2)                                           // entry is not present in the second list
        return false;
      if (entry1.isList() && entry2->isList()) {             // sublist check
        compare(Teuchos::getValue<ParameterList>(entry1), Teuchos::getValue<ParameterList>(*entry2));
        continue;
      }
      if (entry1.getAny(false) != entry2->getAny(false))     // entries have different types or different values
        return false;
    }

    return true;
  }
  static inline bool areSame(const ParameterList& list1, const ParameterList& list2) {
    return compare(list1, list2) && compare(list2, list1);
  }

} // namespace MueLu

#define MUELU_PARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_PARAMETERLISTINTERPRETER_DEF_HPP */

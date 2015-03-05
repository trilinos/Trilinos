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

#include <Xpetra_Matrix.hpp>

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

// These code chunks should only be enabled once Tpetra supports proper graph
// reuse in MMM. At the moment, only Epetra does, while Tpetra throws
// #define REUSE_MATRIX_GRAPHS

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ParameterListInterpreter(Teuchos::ParameterList& paramList,Teuchos::RCP<FactoryFactory> factFact) : factFact_(factFact) {
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ParameterListInterpreter(const std::string& xmlFileName, const Teuchos::Comm<int>& comm,Teuchos::RCP<FactoryFactory> factFact) : factFact_(factFact) {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), comm);
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
    Cycle_     = Hierarchy::GetDefaultCycle();
    blockSize_ = 1;
    dofOffset_ = 0;

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

  // Get value from one of the lists, or set it to default
  // Use case: check for a parameter value in a level-specific sublist, then in a root level list;
  // if it is absent from both, set it to default
#define MUELU_SET_VAR_2LIST(paramList, defaultList, paramName, paramType, varName) \
  paramType varName; \
  if      (paramList.isParameter(paramName))   varName = paramList.get<paramType>(paramName); \
  else if (defaultList.isParameter(paramName)) varName = defaultList.get<paramType>(paramName); \
  else                                         varName = MasterList::getDefault<paramType>(paramName);

#define MUELU_TEST_AND_SET_VAR(paramList, paramName, paramType, varName) \
  (paramList.isParameter(paramName) ? varName = paramList.get<paramType>(paramName), true : false)

  // Set parameter in a list if it is present in any of two lists
  // User case: set factory specific parameter, first checking for a level-specific value, then cheking root level value
#define MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, paramName, paramType, listWrite) \
  try { \
      if      (paramList  .isParameter(paramName)) listWrite.set(paramName, paramList  .get<paramType>(paramName)); \
      else if (defaultList.isParameter(paramName)) listWrite.set(paramName, defaultList.get<paramType>(paramName)); \
  } \
  catch(Teuchos::Exceptions::InvalidParameterType) { \
    TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true, Teuchos::Exceptions::InvalidParameterType, \
        "Error: parameter \"" << paramName << "\" must be of type " << Teuchos::TypeNameTraits<paramType>::name()); \
  } \

#define MUELU_TEST_PARAM_2LIST(paramList, defaultList, paramName, paramType, cmpValue) \
  (cmpValue == ( \
    paramList.isParameter(paramName)   ? paramList  .get<paramType>(paramName) : ( \
    defaultList.isParameter(paramName) ? defaultList.get<paramType>(paramName) : \
                                         MasterList::getDefault<paramType>(paramName) ) ) )

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetEasyParameterList(const Teuchos::ParameterList& constParamList) {
    ParameterList paramList;
    MUELU_SET_VAR_2LIST(constParamList, constParamList, "problem: type", std::string, problemType);
    if (problemType != "unknown") {
      paramList = *MasterList::GetProblemSpecificList(problemType);
      paramList.setParameters(constParamList);

    } else {
      // Create a non const copy of the parameter list
      // Working with a modifiable list is much much easier than with original one
      paramList = constParamList;
    }

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
    blockSize_              = paramList.get<int> ("number of equations", MasterList::getDefault<int>("number of equations"));

    (void)MUELU_TEST_AND_SET_VAR(paramList, "debug: graph level", int, this->graphOutputLevel_);

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

    // Set verbosity parameter
    {
      std::map<std::string,MsgType> verbMap;
      verbMap["none"]    = None;
      verbMap["low"]     = Low;
      verbMap["medium"]  = Medium;
      verbMap["high"]    = High;
      verbMap["extreme"] = Extreme;
      verbMap["test"]    = Test;

      MUELU_SET_VAR_2LIST(paramList, paramList, "verbosity", std::string, verbosityLevel);

      TEUCHOS_TEST_FOR_EXCEPTION(verbMap.count(verbosityLevel) == 0, Exceptions::RuntimeError,
                                 "Invalid verbosity level: \"" << verbosityLevel << "\"");
      this->verbosity_ = verbMap[verbosityLevel];
      this->SetVerbLevel(this->verbosity_);
    }

    // Detect if we need to transfer coordinates to coarse levels. We do that iff
    //  - we use "distance laplacian" dropping on some level, or
    //  - we use repartitioning on some level
    // This is not ideal, as we may have "repartition: enable" turned on by default
    // and not present in the list, but it is better than nothing.
    useCoordinates_ = false;
    if (MUELU_TEST_PARAM_2LIST(paramList, paramList, "repartition: enable", bool, true) ||
        MUELU_TEST_PARAM_2LIST(paramList, paramList, "aggregation: drop scheme", std::string, "distance laplacian")) {
      useCoordinates_ = true;

    } else {
      for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
        std::string levelStr = "level" + toString(levelID);

        if (paramList.isSublist(levelStr)) {
          const ParameterList& levelList = paramList.sublist(levelStr);

          if (MUELU_TEST_PARAM_2LIST(levelList, paramList, "repartition: enable", bool, true) ||
              MUELU_TEST_PARAM_2LIST(levelList, paramList, "aggregation: drop scheme", std::string, "distance laplacian")) {
            useCoordinates_ = true;
            break;
          }
        }
      }
    }

    // Detect if we do implicit P and R rebalance
    changedPRrebalance_ = false;
    if (MUELU_TEST_PARAM_2LIST(paramList, paramList, "repartition: enable", bool, true))
      changedPRrebalance_ = MUELU_TEST_AND_SET_VAR(paramList, "repartition: rebalance P and R", bool, this->doPRrebalance_);

    // Detect if we use implicit transpose
    changedImplicitTranspose_ = MUELU_TEST_AND_SET_VAR(paramList, "transpose: use implicit", bool, this->implicitTranspose_);

    // Create default manager
    // FIXME: should it be here, or higher up
    RCP<FactoryManager> defaultManager = rcp(new FactoryManager());
    defaultManager->SetVerbLevel(this->verbosity_);

    std::vector<keep_pair> keeps0;
    UpdateFactoryManager(paramList, ParameterList(), *defaultManager, 0, keeps0);

    // Create level specific factory managers
    for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
      // Note, that originally if there were no level specific parameters, we
      // simply copied the defaultManager However, with the introduction of
      // levelID to UpdateFactoryManager (required for reuse), we can no longer
      // guarantee that the kept variables are the same for each level even if
      // dependency structure does not change.
      RCP<FactoryManager> levelManager = rcp(new FactoryManager(*defaultManager));
      levelManager->SetVerbLevel(defaultManager->GetVerbLevel());

      ParameterList levelList;
      if (paramList.isSublist("level " + toString(levelID)))
        levelList = paramList.sublist("level " + toString(levelID), true/*mustAlreadyExist*/);

      std::vector<keep_pair> keeps;
      UpdateFactoryManager(levelList, paramList, *levelManager, levelID, keeps);
      this->keep_[levelID] = keeps;

      this->AddFactoryManager(levelID, 1, levelManager);
    }

    // FIXME: parameters passed to packages, like Ifpack2, are not touched by us, resulting in "[unused]" flag
    // being displayed. On the other hand, we don't want to simply iterate through them touching. I don't know
    // what a good solution looks like
    if (MUELU_TEST_PARAM_2LIST(paramList, paramList, "print initial parameters", bool, true))
      this->GetOStream(static_cast<MsgType>(Runtime1), 0) << paramList << std::endl;

    if (MUELU_TEST_PARAM_2LIST(paramList, paramList, "print unused parameters", bool, true)) {
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
        const Teuchos::ParameterList& defaultList, FactoryManager& manager, int levelID, std::vector<keep_pair>& keeps) const {
    // NOTE: Factory::SetParameterList must be called prior to Factory::SetFactory, as
    // SetParameterList sets default values for non mentioned parameters, including factories

    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type", std::string, reuseType);
    TEUCHOS_TEST_FOR_EXCEPTION(reuseType != "none" && reuseType != "tP" && reuseType != "RP" && reuseType != "emin" && reuseType != "RAP" && reuseType != "full",
                               Exceptions::RuntimeError, "Unknown \"reuse: type\" value: \"" << reuseType << "\". Please consult User's Guide.");

    MUELU_SET_VAR_2LIST(paramList, defaultList, "multigrid algorithm", std::string, multigridAlgo);
    TEUCHOS_TEST_FOR_EXCEPTION(multigridAlgo != "unsmoothed" && multigridAlgo != "sa" && multigridAlgo != "pg" && multigridAlgo != "emin",
                               Exceptions::RuntimeError, "Unknown \"multigrid algorithm\" value: \"" << multigridAlgo << "\". Please consult User's Guide.");

    // Only some combinations of reuse and multigrid algorithms are tested, all
    // other are considered invalid at the moment
    if (reuseType == "none" || reuseType == "RP" || reuseType == "RAP") {
      // This works for all kinds of multigrid algorithms

    } else if (reuseType == "tP" && (multigridAlgo != "sa" && multigridAlgo != "unsmoothed")) {
      reuseType = "none";
      this->GetOStream(Warnings0) << "Ignoring \"tP\" reuse option as it is only compatible with \"sa\", or \"unsmoothed\" multigrid algorithms" << std::endl;

    } else if (reuseType == "emin" && multigridAlgo != "emin") {
      reuseType = "none";
      this->GetOStream(Warnings0) << "Ignoring \"emin\" reuse option it is only compatible with \"emin\" multigrid algorithm" << std::endl;
    }

    // === Smoothing ===
    // FIXME: should custom smoother check default list too?
    bool isCustomSmoother =
        paramList.isParameter("smoother: pre or post") ||
        paramList.isParameter("smoother: type")    || paramList.isParameter("smoother: pre type")    || paramList.isParameter("smoother: post type")   ||
        paramList.isSublist  ("smoother: params")  || paramList.isSublist  ("smoother: pre params")  || paramList.isSublist  ("smoother: post params") ||
        paramList.isParameter("smoother: sweeps")  || paramList.isParameter("smoother: pre sweeps")  || paramList.isParameter("smoother: post sweeps") ||
        paramList.isParameter("smoother: overlap") || paramList.isParameter("smoother: pre overlap") || paramList.isParameter("smoother: post overlap");
    MUELU_SET_VAR_2LIST(paramList, defaultList, "smoother: pre or post", std::string, PreOrPost);
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

      RCP<SmootherFactory> preSmoother = Teuchos::null, postSmoother = Teuchos::null;
      std::string          preSmootherType,             postSmootherType;
      ParameterList        preSmootherParams,           postSmootherParams;

      if (paramList.isParameter("smoother: overlap"))
        overlap = paramList.get<int>("smoother: overlap");

      if (PreOrPost == "pre" || PreOrPost == "both") {
        if (paramList.isParameter("smoother: pre type")) {
          preSmootherType = paramList.get<std::string>("smoother: pre type");
        } else {
          MUELU_SET_VAR_2LIST(paramList, defaultList, "smoother: type", std::string, preSmootherTypeTmp);
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

        preSmoother = rcp(new SmootherFactory(rcp(new TrilinosSmoother(preSmootherType, preSmootherParams, overlap))));
      }

      if (PreOrPost == "post" || PreOrPost == "both") {
        if (paramList.isParameter("smoother: post type"))
          postSmootherType = paramList.get<std::string>("smoother: post type");
        else {
          MUELU_SET_VAR_2LIST(paramList, defaultList, "smoother: type", std::string, postSmootherTypeTmp);
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
          postSmoother = rcp(new SmootherFactory(rcp(new TrilinosSmoother(postSmootherType, postSmootherParams, overlap))));
      }

      if (preSmoother == postSmoother)
        manager.SetFactory("Smoother",     preSmoother);
      else {
        manager.SetFactory("PreSmoother",  preSmoother);
        manager.SetFactory("PostSmoother", postSmoother);
      }
    }
    if ((reuseType == "RAP" && levelID) || (reuseType == "full")) {
      // The difference between "RAP" and "full" is keeping smoothers. However,
      // as in both cases we keep coarse matrices, we do not need to update
      // coarse smoothers. On the other hand, if a user changes fine level
      // matrix, "RAP" would update the fine level smoother, while "full" would
      // not
      keeps.push_back(keep_pair("PreSmoother",  manager.GetFactory("PreSmoother") .get()));
      keeps.push_back(keep_pair("PostSmoother", manager.GetFactory("PostSmoother").get()));
    }

    // === Coarse solver ===
    // FIXME: should custom coarse solver check default list too?
    bool isCustomCoarseSolver =
        paramList.isParameter("coarse: type")   ||
        paramList.isParameter("coarse: params");
    if (MUELU_TEST_PARAM_2LIST(paramList, defaultList, "coarse: type", std::string, "none")) {
      manager.SetFactory("CoarseSolver", Teuchos::null);

    } else if (isCustomCoarseSolver) {
      // FIXME: get default values from the factory
      // NOTE: none of the smoothers at the moment use parameter validation framework, so we
      // cannot get the default values from it.
      MUELU_SET_VAR_2LIST(paramList, defaultList, "coarse: type", std::string, coarseType);

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

      manager.SetFactory("CoarseSolver", rcp(new SmootherFactory(coarseSmoother)));
    }
    if ((reuseType == "RAP" && levelID) || (reuseType == "full")) {
      // We do keep_pair("PreSmoother", // manager.GetFactory("CoarseSolver").get())
      // as the coarse solver factory is in fact a smoothing factory, so the
      // only pieces of data it generates are PreSmoother and PostSmoother
      keeps.push_back(keep_pair("PreSmoother", manager.GetFactory("CoarseSolver").get()));
    }

    // === Aggregation ===
    // Aggregation graph
    RCP<CoalesceDropFactory> dropFactory = rcp(new CoalesceDropFactory());
    ParameterList dropParams;
    dropParams.set("lightweight wrap", true);
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: drop scheme",     std::string, dropParams);
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: drop scheme",     std::string, dropParams);
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: drop tol",             double, dropParams);
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: Dirichlet threshold",  double, dropParams);
    dropFactory->SetParameterList(dropParams);
    manager.SetFactory("Graph",       dropFactory);

    // Aggregation sheme
    MUELU_SET_VAR_2LIST(paramList, defaultList, "aggregation: type", std::string, aggType);
    TEUCHOS_TEST_FOR_EXCEPTION(aggType != "uncoupled" && aggType != "coupled", Exceptions::RuntimeError,
                               "Unknown aggregation algorithm: \"" << aggType << "\". Please consult User's Guide.");
    RCP<Factory> aggFactory;
    if      (aggType == "uncoupled") {
      aggFactory = rcp(new UncoupledAggregationFactory());
      ParameterList aggParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: mode",               std::string, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: ordering",           std::string, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: min agg size",               int, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: max agg size",               int, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: max selected neighbors",     int, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 1",            bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 2a",           bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 2b",           bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 3",            bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: preserve Dirichlet points", bool, aggParams);
      aggFactory->SetParameterList(aggParams);
      // I'm not sure why we need this line, but without it we construct CoalesceDropFactory twice
      // SaPFactory
      //   FilteredAFactory
      //     CoalesceDropFactory
      //   TentativePFactory
      //     UncoupledAggregationFactory
      //       CoalesceDropFactory
      aggFactory->SetFactory("DofsPerNode", manager.GetFactory("Graph"));

    } else if (aggType == "coupled") {
      aggFactory = rcp(new CoupledAggregationFactory());
    }
    aggFactory->SetFactory("Graph", manager.GetFactory("Graph"));
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

    if (reuseType == "tP") {
      keeps.push_back(keep_pair("Nullspace", manager.GetFactory("Ptent").get()));
      keeps.push_back(keep_pair("P",         manager.GetFactory("Ptent").get()));
    }

    // Nullspace
    RCP<NullspaceFactory> nullSpace = rcp(new NullspaceFactory());
    nullSpace->SetFactory("Nullspace", manager.GetFactory("Ptent"));
    manager.SetFactory("Nullspace", nullSpace);

    // === Prolongation ===
    if (multigridAlgo == "unsmoothed") {
      manager.SetFactory("P", Ptent);

    } else if (multigridAlgo == "sa") {
      // Smoothed aggregation
      RCP<SaPFactory> P = rcp(new SaPFactory());
      ParameterList Pparams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "sa: damping factor", double, Pparams);
#if REUSE_MATRIX_GRAPHS
      if (reuseType == "tP" && MUELU_TEST_PARAM_2LIST(paramList, defaultList, "sa: use filtered matrix", bool, false)) {
        // Pattern can only be reuse when we don't use filtered matrix, as
        // otherwise graph can potentially change
        Pparams.set("Keep AP Pattern", true);
      }
#endif
      P->SetParameterList(Pparams);

      // Filtering
      if (MUELU_TEST_PARAM_2LIST(paramList, defaultList, "sa: use filtered matrix", bool, true)) {
        RCP<FilteredAFactory> filterFactory = rcp(new FilteredAFactory());
        ParameterList fParams;
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "filtered matrix: use lumping", bool, fParams);
        filterFactory->SetParameterList(fParams);
        filterFactory->SetFactory("Graph",      manager.GetFactory("Graph"));
        // I'm not sure why we need this line. See comments for DofsPerNode for UncoupledAggregation above
        filterFactory->SetFactory("Filtering",  manager.GetFactory("Graph"));
        P->SetFactory("A", filterFactory);
      }

      P->SetFactory("P", manager.GetFactory("Ptent"));
      manager.SetFactory("P", P);

    } else if (multigridAlgo == "emin") {
      MUELU_SET_VAR_2LIST(paramList, defaultList, "emin: pattern", std::string, patternType);
      TEUCHOS_TEST_FOR_EXCEPTION(patternType != "AkPtent", Exceptions::InvalidArgument,
                                 "Invalid pattern name: \"" << patternType << "\". Valid options: \"AkPtent\"");
      // Pattern
      RCP<PatternFactory> patternFactory = rcp(new PatternFactory());
      ParameterList patternParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "emin: pattern order", int, patternParams);
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
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "emin: num iterations",           int, Pparams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "emin: iterative method", std::string, Pparams);
      if (reuseType == "emin") {
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "emin: num reuse iterations",   int, Pparams);
        Pparams.set("Keep P0",          true);
        Pparams.set("Keep Constraint0", true);
      }
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
      MUELU_SET_VAR_2LIST(paramList, defaultList, "problem: symmetric", bool, isSymmetric);
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
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "transpose: use implicit", bool, RAPparams);
#if REUSE_MATRIX_GRAPHS
    if (reuseType == "RP") {
      RAPparams.set("Keep AP Pattern",  true);
      RAPparams.set("Keep RAP Pattern", true);
    }
#endif
    RAP->SetParameterList(RAPparams);
    RAP->SetFactory("P", manager.GetFactory("P"));
    if (!this->implicitTranspose_)
      RAP->SetFactory("R", manager.GetFactory("R"));
    if (MUELU_TEST_PARAM_2LIST(paramList, defaultList, "aggregation: export visualization data", bool, true)) {
      RCP<AggregationExportFactory> aggExport = rcp(new AggregationExportFactory());
      aggExport->SetFactory("DofsPerNode", manager.GetFactory("DofsPerNode"));
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

    if (reuseType == "RP" || reuseType == "RAP" || reuseType == "full")
      keeps.push_back(keep_pair("Nullspace", manager.GetFactory("Nullspace").get()));

    if (reuseType == "RP") {
      keeps.push_back(keep_pair("P", manager.GetFactory("P").get()));
      if (!this->implicitTranspose_)
        keeps.push_back(keep_pair("R", manager.GetFactory("R").get()));
    }
    if ((reuseType == "tP" || reuseType == "RP" || reuseType == "emin") && useCoordinates_)
      keeps.push_back(keep_pair("Coordinates", manager.GetFactory("Coordinates").get()));

    // === Repartitioning ===
    MUELU_SET_VAR_2LIST(paramList, defaultList, "repartition: enable", bool, enableRepart);
    if (enableRepart) {
#ifdef HAVE_MPI
      // Short summary of the issue: RebalanceTransferFactory shares ownership
      // of "P" with SaPFactory, and therefore, changes the stored version.
      // That means that if SaPFactory generated P, and stored it on the level,
      // then after rebalancing the value in that storage changed. It goes
      // against the concept of factories (I think), that every factory is
      // responsible for its own objects, and they are immutable outside.
      //
      // In reuse, this is what happens: as we reuse Importer across setups,
      // the order of factories changes, and coupled with shared ownership
      // leads to problems.
      // *First setup*
      //    SaP               builds     P [and stores it]
      //    TransP            builds     R [and stores it]
      //    RAP               builds     A [and stores it]
      //    RebalanceTransfer rebalances P [and changes the P stored by SaP]   (*)
      //    RebalanceTransfer rebalances R
      //    RebalanceAc       rebalances A
      // *Second setup* ("RP" reuse)
      //    RebalanceTransfer rebalances P [which is incorrect due to (*)]
      //    RebalanceTransfer rebalances R
      //    RAP               builds     A [which is incorrect due to (*)]
      //    RebalanceAc       rebalances A [which throws due to map inconsistency]
      //    ...
      // *Second setup* ("tP" reuse)
      //    SaP               builds     P [and stores it]
      //    RebalanceTransfer rebalances P [and changes the P stored by SaP]   (**)
      //    TransP            builds     R [which is incorrect due to (**)]
      //    RebalanceTransfer rebalances R
      //    ...
      //
      // Couple solutions to this:
      //    1. [implemented] Requre "tP" and "PR" reuse to only be used with
      //       implicit rebalancing.
      //    2. Do deep copy of P, and changed domain map and importer there.
      //       Need to investigate how expensive this is.
      TEUCHOS_TEST_FOR_EXCEPTION(this->doPRrebalance_ && (reuseType == "tP" || reuseType == "RP"), Exceptions::InvalidArgument,
                                 "Reuse types \"tP\" and \"PR\" require \"repartition: rebalance P and R\" set to \"false\"");

      MUELU_SET_VAR_2LIST(paramList, defaultList, "repartition: partitioner", std::string, partName);
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
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: start level",                   int, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: min rows per proc",             int, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: max imbalance",              double, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: keep proc 0",                  bool, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: print partition distribution", bool, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: remap parts",                  bool, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: remap num values",              int, repartParams);
      repartFactory->SetParameterList(repartParams);
      repartFactory->SetFactory("A",         manager.GetFactory("A"));
      repartFactory->SetFactory("Partition", manager.GetFactory("Partition"));
      manager.SetFactory("Importer", repartFactory);
      if (reuseType != "none")
        keeps.push_back(keep_pair("Importer", manager.GetFactory("Importer").get()));

      // Rebalanced A
      RCP<RebalanceAcFactory> newA = rcp(new RebalanceAcFactory());
      ParameterList rebAcParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: use subcommunicators", bool, rebAcParams);
      newA->SetParameterList(rebAcParams);
      newA->SetFactory("A",         manager.GetFactory("A"));
      newA->SetFactory("Importer",  manager.GetFactory("Importer"));
      manager.SetFactory("A", newA);

      // Rebalanced P
      RCP<RebalanceTransferFactory> newP = rcp(new RebalanceTransferFactory());
      ParameterList newPparams;
      newPparams.set("type",                           "Interpolation");
      if (changedPRrebalance_)
        newPparams.set("repartition: rebalance P and R", this->doPRrebalance_);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: use subcommunicators", bool, newPparams);
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
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: use subcommunicators", bool, newRparams);
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
    if (reuseType == "RAP" || reuseType == "full") {
      keeps.push_back(keep_pair("P", manager.GetFactory("P").get()));
      if (!this->implicitTranspose_)
        keeps.push_back(keep_pair("R", manager.GetFactory("R").get()));
      keeps.push_back(keep_pair("A", manager.GetFactory("A").get()));
    }
  }
#undef MUELU_SET_VAR_2LIST
#undef MUELU_TEST_AND_SET_VAR
#undef MUELU_TEST_AND_SET_PARAM_2LIST
#undef MUELU_TEST_PARAM_2LIST

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
    if (paramList.isSublist("Matrix")) {
      blockSize_ = paramList.sublist("Matrix").get<int>("number of equations", MasterList::getDefault<int>("number of equations"));
      dofOffset_ = paramList.sublist("Matrix").get<GlobalOrdinal>("DOF offset", 0); // undocumented parameter allowing to define a DOF offset of the global dofs of an operator (defaul = 0)
    }

    // create new FactoryFactory object if necessary
    if (factFact_ == Teuchos::null)
      factFact_ = Teuchos::rcp(new FactoryFactory());

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

      if (hieraList.isParameter("transpose: use implicit")) {
        this->implicitTranspose_ = hieraList.get<bool>("transpose: use implicit");
        hieraList.remove("transpose: use implicit");
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
          factoryMapOut[paramName] = factFact_->BuildFactory(paramValue, factoryMapIn, factoryManagers);

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
        factoryMapOut[paramName] = factFact_->BuildFactory(paramValue, factoryMapIn, factoryManagers);
      }
    }
  }

  // =====================================================================================================
  // ======================================= MISC functions ==============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupOperator(Operator& Op) const {
    try {
      Matrix& A = dynamic_cast<Matrix&>(Op);
      if (A.GetFixedBlockSize() != blockSize_)
        this->GetOStream(Warnings0) << "Setting matrix block size to " << blockSize_ << " (value of the parameter in the list) "
            << "instead of " << A.GetFixedBlockSize() << " (provided matrix)." << std::endl;

      A.SetFixedBlockSize(blockSize_, dofOffset_);

    } catch (std::bad_cast& e) {
      this->GetOStream(Warnings0) << "Skipping setting block size as the operator is not a matrix" << std::endl;
    }
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

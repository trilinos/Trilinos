
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
#include "MueLu_BrickAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EminPFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_FacadeClassFactory.hpp"
#include "MueLu_FactoryFactory.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_LineDetectionFactory.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_RAPShiftFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SemiCoarsenPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TogglePFactory.hpp"
#include "MueLu_ToggleCoordinatesTransferFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_HybridAggregationFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory_kokkos.hpp"
#include "MueLu_CoordinatesTransferFactory_kokkos.hpp"
#include "MueLu_NullspaceFactory_kokkos.hpp"
#include "MueLu_SaPFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"
#endif

#ifdef HAVE_MUELU_MATLAB
#include "../matlab/src/MueLu_MatlabSmoother_decl.hpp"
#include "../matlab/src/MueLu_MatlabSmoother_def.hpp"
#include "../matlab/src/MueLu_TwoLevelMatlabFactory_decl.hpp"
#include "../matlab/src/MueLu_TwoLevelMatlabFactory_def.hpp"
#include "../matlab/src/MueLu_SingleLevelMatlabFactory_decl.hpp"
#include "../matlab/src/MueLu_SingleLevelMatlabFactory_def.hpp"
#endif

#ifdef HAVE_MUELU_INTREPID2
#include "MueLu_IntrepidPCoarsenFactory.hpp"
#endif

#include <unordered_set>

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ParameterListInterpreter(ParameterList& paramList, Teuchos::RCP<const Teuchos::Comm<int> > comm, Teuchos::RCP<FactoryFactory> factFact, Teuchos::RCP<FacadeClassFactory> facadeFact) : factFact_(factFact) {
    RCP<Teuchos::TimeMonitor> tM = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("MueLu: ParameterListInterpreter (ParameterList)"))));
    if(facadeFact == Teuchos::null)
      facadeFact_ = Teuchos::rcp(new FacadeClassFactory());
    else
      facadeFact_ = facadeFact;

    if (paramList.isParameter("xml parameter file")) {
      std::string filename = paramList.get("xml parameter file", "");
      if (filename.length() != 0) {
        TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), Exceptions::RuntimeError, "xml parameter file requires a valid comm");

        ParameterList paramList2 = paramList;
        Teuchos::updateParametersFromXmlFileAndBroadcast(filename, Teuchos::Ptr<Teuchos::ParameterList>(&paramList2), *comm);
        SetParameterList(paramList2);

      } else {
        SetParameterList(paramList);
      }

    } else {
      SetParameterList(paramList);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ParameterListInterpreter(const std::string& xmlFileName, const Teuchos::Comm<int>& comm, Teuchos::RCP<FactoryFactory> factFact, Teuchos::RCP<FacadeClassFactory> facadeFact) : factFact_(factFact) {
    RCP<Teuchos::TimeMonitor> tM = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("MueLu: ParameterListInterpreter (XML)"))));
    if(facadeFact == Teuchos::null)
      facadeFact_ = Teuchos::rcp(new FacadeClassFactory());
    else
      facadeFact_ = facadeFact;

    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), comm);
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const ParameterList& paramList) {
    Cycle_     = Hierarchy::GetDefaultCycle();
    scalingFactor_= Teuchos::ScalarTraits<double>::one();
    blockSize_ = 1;
    dofOffset_ = 0;

    if (paramList.isSublist("Hierarchy")) {
      SetFactoryParameterList(paramList);

    } else if (paramList.isParameter("MueLu preconditioner") == true) {
      this->GetOStream(Runtime0) << "Use facade class: " << paramList.get<std::string>("MueLu preconditioner")  << std::endl;
      Teuchos::RCP<ParameterList> pp = facadeFact_->SetParameterList(paramList);
      SetFactoryParameterList(*pp);

    }  else {
      // The validator doesn't work correctly for non-serializable data (Hint: template parameters), so strip it out
      ParameterList serialList, nonSerialList;

      ExtractNonSerializableData(paramList, serialList, nonSerialList);
      Validate(serialList);
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

#ifndef HAVE_MUELU_KOKKOS_REFACTOR
#define MUELU_KOKKOS_FACTORY(varName, oldFactory, newFactory) \
  RCP<Factory> varName = rcp(new oldFactory());
#define MUELU_KOKKOS_FACTORY_NO_DECL(varName, oldFactory, newFactory) \
  varName = rcp(new oldFactory());
#else
#define MUELU_KOKKOS_FACTORY(varName, oldFactory, newFactory) \
  RCP<Factory> varName; \
  if (!useKokkos_) varName = rcp(new oldFactory()); \
  else             varName = rcp(new newFactory());
#define MUELU_KOKKOS_FACTORY_NO_DECL(varName, oldFactory, newFactory) \
  if (!useKokkos_) varName = rcp(new oldFactory()); \
  else             varName = rcp(new newFactory());
#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  SetEasyParameterList(const ParameterList& constParamList) {
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

    // Check for Kokkos
#if !defined(HAVE_MUELU_KOKKOS_REFACTOR)
    useKokkos_ = false;
#elif defined(HAVE_MUELU_KOKKOS_REFACTOR_USE_BY_DEFAULT)
    ParameterList tempList("tempList");
    tempList.set("use kokkos refactor",true);
    MUELU_SET_VAR_2LIST(paramList, tempList, "use kokkos refactor", bool, useKokkos);
    useKokkos_ = useKokkos;
#else
    MUELU_SET_VAR_2LIST(paramList, paramList, "use kokkos refactor", bool, useKokkos);
    useKokkos_ = useKokkos;
#endif

    // Check for timer synchronization
    MUELU_SET_VAR_2LIST(paramList, paramList, "synchronize factory timers", bool, syncTimers);
    if (syncTimers)
        Factory::EnableTimerSync();

    // Translate cycle type parameter
    if (paramList.isParameter("cycle type")) {
      std::map<std::string, CycleType> cycleMap;
      cycleMap["V"] = VCYCLE;
      cycleMap["W"] = WCYCLE;

      auto cycleType = paramList.get<std::string>("cycle type");
      TEUCHOS_TEST_FOR_EXCEPTION(cycleMap.count(cycleType) == 0, Exceptions::RuntimeError,
                                 "Invalid cycle type: \"" << cycleType << "\"");
      Cycle_ = cycleMap[cycleType];
    }

    if (paramList.isParameter("coarse grid correction scaling factor"))
      scalingFactor_ = paramList.get<double>("coarse grid correction scaling factor");

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
      if (printList.isParameter("Nullspace"))
        this->nullspaceToPrint_  = Teuchos::getArrayFromStringParameter<int>(printList, "Nullspace");
      if (printList.isParameter("Coordinates"))
        this->coordinatesToPrint_  = Teuchos::getArrayFromStringParameter<int>(printList, "Coordinates");
      if (printList.isParameter("pcoarsen: element to node map"))
        this->elementToNodeMapsToPrint_  = Teuchos::getArrayFromStringParameter<int>(printList, "pcoarsen: element to node map");
    }

    // Set verbosity parameter
    VerbLevel oldVerbLevel = VerboseObject::GetDefaultVerbLevel();
    {
      std::map<std::string, MsgType> verbMap;
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
      VerboseObject::SetDefaultVerbLevel(this->verbosity_);
    }

    // Detect if we need to transfer coordinates to coarse levels. We do that iff
    //  - we use "distance laplacian" dropping on some level, or
    //  - we use repartitioning on some level
    //  - we use brick aggregation
    //  - we use Ifpack2 line partitioner
    // This is not ideal, as we may have "repartition: enable" turned on by default
    // and not present in the list, but it is better than nothing.
    useCoordinates_ = false;
    if (MUELU_TEST_PARAM_2LIST(paramList, paramList, "repartition: enable",      bool,        true) ||
        MUELU_TEST_PARAM_2LIST(paramList, paramList, "aggregation: drop scheme", std::string, "distance laplacian") ||
        MUELU_TEST_PARAM_2LIST(paramList, paramList, "aggregation: type",        std::string, "brick") ||
        MUELU_TEST_PARAM_2LIST(paramList, paramList, "aggregation: export visualization data", bool, true)) {
      useCoordinates_ = true;
    } else if(paramList.isSublist("smoother: params")) {
      const auto smooParamList = paramList.sublist("smoother: params");
      if(smooParamList.isParameter("partitioner: type") &&
         (smooParamList.get<std::string>("partitioner: type") == "line")) {
        useCoordinates_ = true;
      }
    } else {
      for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
        std::string levelStr = "level " + toString(levelID);

        if (paramList.isSublist(levelStr)) {
          const ParameterList& levelList = paramList.sublist(levelStr);

          if (MUELU_TEST_PARAM_2LIST(levelList, paramList, "repartition: enable",      bool,        true) ||
              MUELU_TEST_PARAM_2LIST(levelList, paramList, "aggregation: drop scheme", std::string, "distance laplacian") ||
              MUELU_TEST_PARAM_2LIST(levelList, paramList, "aggregation: type",        std::string, "brick") ||
              MUELU_TEST_PARAM_2LIST(levelList, paramList, "aggregation: export visualization data", bool, true)) {
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

    // We will ignore keeps0
    std::vector<keep_pair> keeps0;
    UpdateFactoryManager(paramList, ParameterList(), *defaultManager, 0/*levelID*/, keeps0);

    // Create level specific factory managers
    for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
      // Note, that originally if there were no level specific parameters, we
      // simply copied the defaultManager However, with the introduction of
      // levelID to UpdateFactoryManager (required for reuse), we can no longer
      // guarantee that the kept variables are the same for each level even if
      // dependency structure does not change.
      RCP<FactoryManager> levelManager = rcp(new FactoryManager(*defaultManager));
      levelManager->SetVerbLevel(defaultManager->GetVerbLevel());

      std::vector<keep_pair> keeps;
      if (paramList.isSublist("level " + toString(levelID))) {
        // We do this so the parameters on the level get flagged correctly as "used"
        ParameterList& levelList = paramList.sublist("level " + toString(levelID), true/*mustAlreadyExist*/);
        UpdateFactoryManager(levelList, paramList, *levelManager, levelID, keeps);

      } else {
        ParameterList levelList;
        UpdateFactoryManager(levelList, paramList, *levelManager, levelID, keeps);
      }

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

    VerboseObject::SetDefaultVerbLevel(oldVerbLevel);
  }


  // =====================================================================================================
  // ==================================== UpdateFactoryManager ===========================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                       int levelID, std::vector<keep_pair>& keeps) const
  {
    // NOTE: Factory::SetParameterList must be called prior to Factory::SetFactory, as
    // SetParameterList sets default values for non mentioned parameters, including factories

    using strings = std::unordered_set<std::string>;

    // shortcut
    if (paramList.numParams() == 0 && defaultList.numParams() > 0)
      paramList = ParameterList(defaultList);

    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type", std::string, reuseType);
    TEUCHOS_TEST_FOR_EXCEPTION(strings({"none", "tP", "RP", "emin", "RAP", "full", "S"}).count(reuseType) == 0,
        Exceptions::RuntimeError, "Unknown \"reuse: type\" value: \"" << reuseType << "\". Please consult User's Guide.");

    MUELU_SET_VAR_2LIST(paramList, defaultList, "multigrid algorithm", std::string, multigridAlgo);
    TEUCHOS_TEST_FOR_EXCEPTION(strings({"unsmoothed", "sa", "pg", "emin", "matlab", "pcoarsen"}).count(multigridAlgo) == 0,
        Exceptions::RuntimeError, "Unknown \"multigrid algorithm\" value: \"" << multigridAlgo << "\". Please consult User's Guide.");
#ifndef HAVE_MUELU_MATLAB
    TEUCHOS_TEST_FOR_EXCEPTION(multigridAlgo == "matlab", Exceptions::RuntimeError,
        "Cannot use matlab for multigrid algorithm - MueLu was not configured with MATLAB support.");
#endif
#ifndef HAVE_MUELU_INTREPID2
    TEUCHOS_TEST_FOR_EXCEPTION(multigridAlgo == "pcoarsen", Exceptions::RuntimeError,
        "Cannot use IntrepidPCoarsen prolongator factory - MueLu was not configured with Intrepid support.");
#endif

    // Only some combinations of reuse and multigrid algorithms are tested, all
    // other are considered invalid at the moment
    if (reuseType == "none" || reuseType == "S" || reuseType == "RP" || reuseType == "RAP") {
      // This works for all kinds of multigrid algorithms

    } else if (reuseType == "tP" && (multigridAlgo != "sa" && multigridAlgo != "unsmoothed")) {
      reuseType = "none";
      this->GetOStream(Warnings0) << "Ignoring \"tP\" reuse option as it is only compatible with \"sa\", "
          "or \"unsmoothed\" multigrid algorithms" << std::endl;

    } else if (reuseType == "emin" && multigridAlgo != "emin") {
      reuseType = "none";
      this->GetOStream(Warnings0) << "Ignoring \"emin\" reuse option it is only compatible with "
          "\"emin\" multigrid algorithm" << std::endl;
    }

    // == Non-serializable data ===
    // Check both the parameter and the type
    bool have_userP = false;
    if (paramList.isParameter("P") && !paramList.get<RCP<Matrix> >("P").is_null())
      have_userP  = true;

    // == Smoothers ==
    UpdateFactoryManager_Smoothers(paramList, defaultList, manager, levelID, keeps);

    // === Coarse solver ===
    UpdateFactoryManager_CoarseSolvers(paramList, defaultList, manager, levelID, keeps);

    // === Aggregation ===
    UpdateFactoryManager_Aggregation_TentativeP(paramList, defaultList, manager, levelID, keeps);

    // === Nullspace ===
    RCP<Factory> nullSpaceFactory; // Cache this guy for the combination of semi-coarsening & repartitioning
    UpdateFactoryManager_Nullspace(paramList, defaultList, manager, levelID, keeps, nullSpaceFactory);

    // === Prolongation ===
    // NOTE: None of the UpdateFactoryManager routines called here check the
    // multigridAlgo. This is intentional, to allow for reuse of components
    // underneath. Thus, the multigridAlgo was checked in the beginning of the
    // function.
    if (have_userP) {
      // User prolongator
      manager.SetFactory("P", NoFactory::getRCP());

    } else if (multigridAlgo == "unsmoothed") {
      // Unsmoothed aggregation
      manager.SetFactory("P", manager.GetFactory("Ptent"));

    } else if (multigridAlgo == "sa") {
      // Smoothed aggregation
      UpdateFactoryManager_SA(paramList, defaultList, manager, levelID, keeps);

    } else if (multigridAlgo == "emin") {
      // Energy minimization
      UpdateFactoryManager_Emin(paramList, defaultList, manager, levelID, keeps);

    } else if (multigridAlgo == "pg") {
      // Petrov-Galerkin
      UpdateFactoryManager_PG(paramList, defaultList, manager, levelID, keeps);

    } else if (multigridAlgo == "matlab") {
      // Matlab Coarsneing
      UpdateFactoryManager_Matlab(paramList, defaultList, manager, levelID, keeps);

    } else if (multigridAlgo == "pcoarsen") {
      // P-Coarsening
      UpdateFactoryManager_PCoarsen(paramList, defaultList, manager, levelID, keeps);
    }

    // === Semi-coarsening ===
    UpdateFactoryManager_SemiCoarsen(paramList, defaultList, manager, levelID, keeps);

    // === Restriction ===
    UpdateFactoryManager_Restriction(paramList, defaultList, manager, levelID, keeps);

    // === RAP ===
    UpdateFactoryManager_RAP(paramList, defaultList, manager, levelID, keeps);

    // === Coordinates ===
    UpdateFactoryManager_Coordinates(paramList, defaultList, manager, levelID, keeps);

    // === Pre-Repartition Keeps for Reuse ===
    if ((reuseType == "RP" || reuseType == "RAP" || reuseType == "full") && levelID)
      keeps.push_back(keep_pair("Nullspace", manager.GetFactory("Nullspace").get()));

    if (reuseType == "RP" && levelID) {
      keeps.push_back(keep_pair("P", manager.GetFactory("P").get()));
      if (!this->implicitTranspose_)
        keeps.push_back(keep_pair("R", manager.GetFactory("R").get()));
    }
    if ((reuseType == "tP" || reuseType == "RP" || reuseType == "emin") && useCoordinates_ && levelID)
      keeps.push_back(keep_pair("Coordinates", manager.GetFactory("Coordinates").get()));

    // === Repartitioning ===
    UpdateFactoryManager_Repartition(paramList, defaultList, manager, levelID, keeps, nullSpaceFactory);

    // === Final Keeps for Reuse ===
    if ((reuseType == "RAP" || reuseType == "full") && levelID) {
      keeps.push_back(keep_pair("P", manager.GetFactory("P").get()));
      if (!this->implicitTranspose_)
        keeps.push_back(keep_pair("R", manager.GetFactory("R").get()));
      keeps.push_back(keep_pair("A", manager.GetFactory("A").get()));
    }
  }

  // =====================================================================================================
  // ========================================= Smoothers =================================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Smoothers(ParameterList& paramList, const ParameterList& defaultList,
                                 FactoryManager& manager, int levelID, std::vector<keep_pair>& keeps) const
  {
    MUELU_SET_VAR_2LIST(paramList, defaultList, "multigrid algorithm", std::string, multigridAlgo);
    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type", std::string, reuseType);

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

#ifdef HAVE_MUELU_INTREPID2
      // Propagate P-coarsening for Topo smoothing
      if (multigridAlgo == "pcoarsen" && preSmootherType == "TOPOLOGICAL" &&
          defaultList.isParameter("pcoarsen: schedule") && defaultList.isParameter("pcoarsen: element")) {
        // P-Coarsening by schedule (new interface)
        // NOTE: levelID represents the *coarse* level in this case
        auto pcoarsen_schedule = Teuchos::getArrayFromStringParameter<int>(defaultList, "pcoarsen: schedule");
        auto pcoarsen_element  = defaultList.get<std::string>("pcoarsen: element");

        if (levelID < (int)pcoarsen_schedule.size()) {
          // Topo info for P-Coarsening
          auto lo = pcoarsen_element + std::to_string(pcoarsen_schedule[levelID]);
          preSmootherParams.set("pcoarsen: hi basis", lo);
        }
      }
#endif

#ifdef HAVE_MUELU_MATLAB
        if (preSmootherType == "matlab")
          preSmoother = rcp(new SmootherFactory(rcp(new MatlabSmoother(preSmootherParams))));
        else
#endif
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
        else {
#ifdef HAVE_MUELU_INTREPID2
          // Propagate P-coarsening for Topo smoothing
          if (multigridAlgo == "pcoarsen" && preSmootherType == "TOPOLOGICAL" &&
              defaultList.isParameter("pcoarsen: schedule") && defaultList.isParameter("pcoarsen: element")) {
            // P-Coarsening by schedule (new interface)
            // NOTE: levelID represents the *coarse* level in this case
            auto pcoarsen_schedule = Teuchos::getArrayFromStringParameter<int>(defaultList,"pcoarsen: schedule");
            auto pcoarsen_element = defaultList.get<std::string>("pcoarsen: element");

            if (levelID < (int)pcoarsen_schedule.size()) {
              // Topo info for P-Coarsening
              auto lo = pcoarsen_element + std::to_string(pcoarsen_schedule[levelID]);
              postSmootherParams.set("pcoarsen: hi basis", lo);
            }
          }
#endif

#ifdef HAVE_MUELU_MATLAB
          if (postSmootherType == "matlab")
            postSmoother = rcp(new SmootherFactory(rcp(new MatlabSmoother(postSmootherParams))));
          else
#endif
            postSmoother = rcp(new SmootherFactory(rcp(new TrilinosSmoother(postSmootherType, postSmootherParams, overlap))));
        }
      }

      if (preSmoother == postSmoother)
        manager.SetFactory("Smoother",     preSmoother);
      else {
        manager.SetFactory("PreSmoother",  preSmoother);
        manager.SetFactory("PostSmoother", postSmoother);
      }
    }

    // The first clause is not necessary, but it is here for clarity Smoothers
    // are reused if smoother explicitly said to reuse them, or if any other
    // reuse option is enabled
    bool reuseSmoothers = (reuseType == "S" || reuseType != "none");
    if (reuseSmoothers) {
      auto preSmootherFactory = rcp_const_cast<Factory>(rcp_dynamic_cast<const Factory>(manager.GetFactory("PreSmoother")));

      if (preSmootherFactory != Teuchos::null) {
        ParameterList postSmootherFactoryParams;
        postSmootherFactoryParams.set("keep smoother data", true);
        preSmootherFactory->SetParameterList(postSmootherFactoryParams);

        keeps.push_back(keep_pair("PreSmoother data", preSmootherFactory.get()));
      }

      auto postSmootherFactory = rcp_const_cast<Factory>(rcp_dynamic_cast<const Factory>(manager.GetFactory("PostSmoother")));
      if (postSmootherFactory != Teuchos::null) {
        ParameterList postSmootherFactoryParams;
        postSmootherFactoryParams.set("keep smoother data", true);
        postSmootherFactory->SetParameterList(postSmootherFactoryParams);

        keeps.push_back(keep_pair("PostSmoother data", postSmootherFactory.get()));
      }

      auto coarseFactory = rcp_const_cast<Factory>(rcp_dynamic_cast<const Factory>(manager.GetFactory("CoarseSolver")));
      if (coarseFactory != Teuchos::null) {
        ParameterList coarseFactoryParams;
        coarseFactoryParams.set("keep smoother data", true);
        coarseFactory->SetParameterList(coarseFactoryParams);

        keeps.push_back(keep_pair("PreSmoother data", coarseFactory.get()));
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

      // We do keep_pair("PreSmoother", manager.GetFactory("CoarseSolver").get())
      // as the coarse solver factory is in fact a smoothing factory, so the
      // only pieces of data it generates are PreSmoother and PostSmoother
      keeps.push_back(keep_pair("PreSmoother", manager.GetFactory("CoarseSolver").get()));
    }
  }

  // =====================================================================================================
  // ====================================== Coarse Solvers ===============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_CoarseSolvers(ParameterList& paramList, const ParameterList& defaultList,
                                     FactoryManager& manager, int levelID, std::vector<keep_pair>& keeps) const
  {
    // FIXME: should custom coarse solver check default list too?
    bool isCustomCoarseSolver =
        paramList.isParameter("coarse: type")   ||
        paramList.isParameter("coarse: params");
    if (MUELU_TEST_PARAM_2LIST(paramList, defaultList, "coarse: type", std::string, "none")) {
      this->GetOStream(Warnings0) << "No coarse grid solver" << std::endl;
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

      using strings = std::unordered_set<std::string>;

      RCP<SmootherPrototype> coarseSmoother;
      // TODO: this is not a proper place to check. If we consider direct solver to be a special
      // case of smoother, we would like to unify Amesos and Ifpack2 smoothers in src/Smoothers, and
      // have a single factory responsible for those. Then, this check would belong there.
      if (strings({"RELAXATION", "CHEBYSHEV", "ILUT", "ILU", "RILUK", "SCHWARZ", "Amesos",
        "BLOCK RELAXATION", "BLOCK_RELAXATION", "BLOCKRELAXATION" ,
        "SPARSE BLOCK RELAXATION", "SPARSE_BLOCK_RELAXATION", "SPARSEBLOCKRELAXATION",
        "LINESMOOTHING_BANDEDRELAXATION", "LINESMOOTHING_BANDED_RELAXATION", "LINESMOOTHING_BANDED RELAXATION",
        "LINESMOOTHING_TRIDIRELAXATION", "LINESMOOTHING_TRIDI_RELAXATION", "LINESMOOTHING_TRIDI RELAXATION",
        "LINESMOOTHING_TRIDIAGONALRELAXATION", "LINESMOOTHING_TRIDIAGONAL_RELAXATION", "LINESMOOTHING_TRIDIAGONAL RELAXATION",             
        "TOPOLOGICAL", "FAST_ILU", "FAST_IC", "FAST_ILDL"}).count(coarseType)) {
        coarseSmoother = rcp(new TrilinosSmoother(coarseType, coarseParams, overlap));
      } else {
#ifdef HAVE_MUELU_MATLAB
        if (coarseType == "matlab")
          coarseSmoother = rcp(new MatlabSmoother(coarseParams));
        else
#endif
        coarseSmoother = rcp(new DirectSolver(coarseType, coarseParams));
      }

      manager.SetFactory("CoarseSolver", rcp(new SmootherFactory(coarseSmoother)));
    }
  }

  // =====================================================================================================
  // ========================================= Smoothers =================================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Aggregation_TentativeP(ParameterList& paramList, const ParameterList& defaultList,
                                              FactoryManager& manager, int levelID, std::vector<keep_pair>& keeps) const
  {
    using strings = std::unordered_set<std::string>;

    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type", std::string, reuseType);

    // Aggregation graph
    RCP<Factory> dropFactory;

    if (MUELU_TEST_PARAM_2LIST(paramList, paramList, "aggregation: drop scheme", std::string, "matlab")) {
#ifdef HAVE_MUELU_MATLAB
      dropFactory = rcp(new SingleLevelMatlabFactory());
      ParameterList socParams  = paramList.sublist("strength-of-connection: params");
      dropFactory->SetParameterList(socParams);
#else
      throw std::runtime_error("Cannot use MATLAB evolutionary strength-of-connection - MueLu was not configured with MATLAB support.");
#endif
    } else {
      MUELU_KOKKOS_FACTORY_NO_DECL(dropFactory, CoalesceDropFactory, CoalesceDropFactory_kokkos);
      ParameterList dropParams;
      dropParams.set("lightweight wrap", true);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: drop scheme",     std::string, dropParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: drop tol",             double, dropParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: Dirichlet threshold",  double, dropParams);
      if (useKokkos_) {
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "filtered matrix: use lumping",      bool, dropParams);
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "filtered matrix: reuse graph",      bool, dropParams);
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "filtered matrix: reuse eigenvalue", bool, dropParams);
      }
      dropFactory->SetParameterList(dropParams);
    }
    manager.SetFactory("Graph", dropFactory);

    // Aggregation scheme
    MUELU_SET_VAR_2LIST(paramList, defaultList, "aggregation: type", std::string, aggType);
    TEUCHOS_TEST_FOR_EXCEPTION(!strings({"uncoupled", "coupled", "brick", "matlab"}).count(aggType),
        Exceptions::RuntimeError, "Unknown aggregation algorithm: \"" << aggType << "\". Please consult User's Guide.");
    #ifndef HAVE_MUELU_MATLAB
    if (aggType == "matlab")
        throw std::runtime_error("Cannot use MATLAB aggregation - MueLu was not configured with MATLAB support.");
    #endif
    RCP<Factory> aggFactory;
    if (aggType == "uncoupled") {
      MUELU_KOKKOS_FACTORY_NO_DECL(aggFactory, UncoupledAggregationFactory, UncoupledAggregationFactory_kokkos);
      ParameterList aggParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: mode",               std::string, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: ordering",           std::string, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: min agg size",               int, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: max agg size",               int, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: max selected neighbors",     int, aggParams);
      if(useKokkos_) {
        //if not using kokkos refactor Uncoupled, there is no algorithm option (always Serial)
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: phase 1 algorithm",  std::string, aggParams);
      }
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 1",            bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 2a",           bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 2b",           bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: enable phase 3",            bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: preserve Dirichlet points", bool, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: error on nodes with no on-rank neighbors", bool, aggParams);
      aggFactory->SetParameterList(aggParams);
      // make sure that the aggregation factory has all necessary data
      aggFactory->SetFactory("DofsPerNode", manager.GetFactory("Graph"));
      aggFactory->SetFactory("Graph", manager.GetFactory("Graph"));

    } else if (aggType == "coupled") {
      aggFactory = rcp(new CoupledAggregationFactory());
      aggFactory->SetFactory("Graph", manager.GetFactory("Graph"));

    } else if (aggType == "brick") {
      aggFactory = rcp(new BrickAggregationFactory());
      ParameterList aggParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: brick x size", int, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: brick y size", int, aggParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: brick z size", int, aggParams);
      aggFactory->SetParameterList(aggParams);

      if (levelID > 1) {
        // We check for levelID > 0, as in the interpreter aggFactory for
        // levelID really corresponds to level 0. Managers are clunky, as they
        // contain factories for two different levels
        aggFactory->SetFactory("Coordinates", this->GetFactoryManager(levelID-1)->GetFactory("Coordinates"));
      }
    }
#ifdef HAVE_MUELU_MATLAB
    else if(aggType == "matlab") {
      ParameterList aggParams = paramList.sublist("aggregation: params");
      aggFactory = rcp(new SingleLevelMatlabFactory());
      aggFactory->SetParameterList(aggParams);
    }
#endif
    manager.SetFactory("Aggregates", aggFactory);

    // Coarse map
    MUELU_KOKKOS_FACTORY(coarseMap, CoarseMapFactory, CoarseMapFactory_kokkos);
    coarseMap->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
    manager.SetFactory("CoarseMap", coarseMap);

    // Tentative P
    MUELU_KOKKOS_FACTORY(Ptent, TentativePFactory, TentativePFactory_kokkos);
    ParameterList ptentParams;
    if (paramList.isSublist("matrixmatrix: kernel params"))
      ptentParams.sublist("matrixmatrix: kernel params", false) = paramList.sublist("matrixmatrix: kernel params");
    if (defaultList.isSublist("matrixmatrix: kernel params"))
      ptentParams.sublist("matrixmatrix: kernel params", false) = defaultList.sublist("matrixmatrix: kernel params");
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "tentative: calculate qr", bool, ptentParams);
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "tentative: build coarse coordinates", bool, ptentParams);
    Ptent->SetParameterList(ptentParams);
    Ptent->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
    Ptent->SetFactory("CoarseMap",  manager.GetFactory("CoarseMap"));
    manager.SetFactory("Ptent",     Ptent);

    if (reuseType == "tP" && levelID) {
      keeps.push_back(keep_pair("Nullspace", Ptent.get()));
      keeps.push_back(keep_pair("P",         Ptent.get()));
    }
  }

  // =====================================================================================================
  // ============================================ RAP ====================================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_RAP(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                           int levelID, std::vector<keep_pair>& keeps) const
  {
    if (paramList.isParameter("A") && !paramList.get<RCP<Matrix> >("A").is_null()) {
      // We have user matrix A
      manager.SetFactory("A", NoFactory::getRCP());
      return;
    }

    ParameterList RAPparams;

    RCP<RAPFactory>      RAP;
    RCP<RAPShiftFactory> RAPs;
    // Allow for Galerkin or shifted RAP
    // FIXME: Should this not be some form of MUELU_SET_VAR_2LIST?
    std::string alg = paramList.get("rap: algorithm", "galerkin");
    if (alg == "shift" || alg == "non-galerkin") {
      RAPs = rcp(new RAPShiftFactory());
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "rap: shift", double, RAPparams);

    } else {
      RAP = rcp(new RAPFactory());
    }

    if (paramList.isSublist("matrixmatrix: kernel params"))
      RAPparams.sublist("matrixmatrix: kernel params", false) = paramList.sublist("matrixmatrix: kernel params");
    if (defaultList.isSublist("matrixmatrix: kernel params"))
      RAPparams.sublist("matrixmatrix: kernel params", false) = defaultList.sublist("matrixmatrix: kernel params");
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "transpose: use implicit", bool, RAPparams);
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "rap: fix zero diagonals", bool, RAPparams);
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "rap: triple product", bool, RAPparams);

    try {
      if (paramList.isParameter("aggregation: allow empty prolongator columns")) {
        RAPparams.set("CheckMainDiagonal",  paramList.get<bool>("aggregation: allow empty prolongator columns"));
        RAPparams.set("RepairMainDiagonal", paramList.get<bool>("aggregation: allow empty prolongator columns"));
      }
      else if (defaultList.isParameter("aggregation: allow empty prolongator columns")) {
        RAPparams.set("CheckMainDiagonal",  defaultList.get<bool>("aggregation: allow empty prolongator columns"));
        RAPparams.set("RepairMainDiagonal", defaultList.get<bool>("aggregation: allow empty prolongator columns"));
      }

    } catch (Teuchos::Exceptions::InvalidParameterType) {
      TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true, Teuchos::Exceptions::InvalidParameterType,
          "Error: parameter \"aggregation: allow empty prolongator columns\" must be of type " << Teuchos::TypeNameTraits<bool>::name());
    }

    if (!RAP.is_null()) {
      RAP->SetParameterList(RAPparams);
      RAP->SetFactory("P", manager.GetFactory("P"));
    } else {
      RAPs->SetParameterList(RAPparams);
      RAPs->SetFactory("P", manager.GetFactory("P"));
    }

    if (!this->implicitTranspose_) {
       if (!RAP.is_null())
         RAP->SetFactory("R", manager.GetFactory("R"));
       else
         RAPs->SetFactory("R", manager.GetFactory("R"));
    }

    if (MUELU_TEST_PARAM_2LIST(paramList, defaultList, "aggregation: export visualization data", bool, true)) {
      RCP<AggregationExportFactory> aggExport = rcp(new AggregationExportFactory());
      ParameterList aggExportParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: output filename",             std::string, aggExportParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: output file: agg style",      std::string, aggExportParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: output file: iter",                   int, aggExportParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: output file: time step",              int, aggExportParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: output file: fine graph edges",      bool, aggExportParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: output file: coarse graph edges",    bool, aggExportParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "aggregation: output file: build colormap",        bool, aggExportParams);
      aggExport->SetParameterList(aggExportParams);
      aggExport->SetFactory("DofsPerNode", manager.GetFactory("DofsPerNode"));

      if (!RAP.is_null())
        RAP->AddTransferFactory(aggExport);
      else
        RAPs->AddTransferFactory(aggExport);
    }
    if (!RAP.is_null())
      manager.SetFactory("A", RAP);
    else
      manager.SetFactory("A", RAPs);

    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type",      std::string, reuseType);
    MUELU_SET_VAR_2LIST(paramList, defaultList, "sa: use filtered matrix", bool, useFiltering);
    bool filteringChangesMatrix = useFiltering && !MUELU_TEST_PARAM_2LIST(paramList, defaultList, "aggregation: drop tol", double, 0);

    if (reuseType == "RP" || (reuseType == "tP" && !filteringChangesMatrix)) {
      if (!RAP.is_null()) {
        keeps.push_back(keep_pair("AP reuse data",  RAP.get()));
        keeps.push_back(keep_pair("RAP reuse data", RAP.get()));

      } else {
        keeps.push_back(keep_pair("AP reuse data",  RAPs.get()));
        keeps.push_back(keep_pair("RAP reuse data", RAPs.get()));
      }
    }
  }

  // =====================================================================================================
  // ======================================= Restriction =================================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Coordinates(ParameterList& paramList, const ParameterList& defaultList,
                                   FactoryManager& manager, int levelID, std::vector<keep_pair>& keeps) const
  {
    bool have_userCO = false;
    if (paramList.isParameter("Coordinates") && !paramList.get<RCP<MultiVector> >("Coordinates").is_null())
      have_userCO = true;

    if (useCoordinates_) {
      if (have_userCO) {
        manager.SetFactory("Coordinates", NoFactory::getRCP());

      } else {
        MUELU_KOKKOS_FACTORY(coords, CoordinatesTransferFactory, CoordinatesTransferFactory_kokkos);
        coords->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
        coords->SetFactory("CoarseMap",  manager.GetFactory("CoarseMap"));
        manager.SetFactory("Coordinates", coords);

        auto RAP = rcp_const_cast<RAPFactory>(rcp_dynamic_cast<const RAPFactory>(manager.GetFactory("A")));
        if (!RAP.is_null()) {
          RAP->AddTransferFactory(manager.GetFactory("Coordinates"));
        } else {
          auto RAPs = rcp_const_cast<RAPShiftFactory>(rcp_dynamic_cast<const RAPShiftFactory>(manager.GetFactory("A")));
          RAPs->AddTransferFactory(manager.GetFactory("Coordinates"));
        }
      }
    }
  }

  // =====================================================================================================
  // ======================================= Restriction =================================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Restriction(ParameterList& paramList, const ParameterList& defaultList,
                                   FactoryManager& manager, int levelID, std::vector<keep_pair>& keeps) const
  {
    MUELU_SET_VAR_2LIST(paramList, defaultList, "multigrid algorithm", std::string, multigridAlgo);
    bool have_userR = false;
    if (paramList.isParameter("R") && !paramList.get<RCP<Matrix> >("R").is_null())
      have_userR = true;

    // === Restriction ===
    if (!this->implicitTranspose_) {
      MUELU_SET_VAR_2LIST(paramList, defaultList, "problem: symmetric", bool, isSymmetric);

      if (isSymmetric == false && (multigridAlgo == "unsmoothed" || multigridAlgo == "emin")) {
        this->GetOStream(Warnings0) <<
            "Switching \"problem: symmetric\" parameter to symmetric as multigrid algorithm. " <<
            multigridAlgo << " is primarily supposed to be used for symmetric problems.\n\n" <<
            "Please note: if you are using \"unsmoothed\" transfer operators the \"problem: symmetric\" parameter " <<
            "has no real mathematical meaning, i.e. you can use it for non-symmetric\n" <<
            "problems, too. With \"problem: symmetric\"=\"symmetric\" you can use implicit transpose for building " <<
            "the restriction operators which may drastically reduce the amount of consumed memory." << std::endl;
        isSymmetric = true;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(multigridAlgo == "pg" && isSymmetric == true, Exceptions::RuntimeError,
          "Petrov-Galerkin smoothed transfer operators are only allowed for non-symmetric problems: Set \"problem: symmetric\" to false!\n" \
          "While PG smoothed transfer operators generally would also work for symmetric problems this is an unusual use case. " \
          "You can use the factory-based xml interface though if you need PG-AMG for symmetric problems.");

        if (have_userR) {
          manager.SetFactory("R", NoFactory::getRCP());
        } else {
          RCP<Factory> R;
          if (isSymmetric)  R = rcp(new TransPFactory());
          else              R = rcp(new GenericRFactory());

          R->SetFactory("P", manager.GetFactory("P"));
          manager.SetFactory("R", R);
        }

      } else {
        manager.SetFactory("R", Teuchos::null);
      }
  }

  // =====================================================================================================
  // ========================================= Repartition ===============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Repartition(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                                   int levelID, std::vector<keep_pair>& keeps, RCP<Factory> & nullSpaceFactory) const
  {
    // === Repartitioning ===
    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type", std::string, reuseType);
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

      // TEUCHOS_TEST_FOR_EXCEPTION(aggType == "brick", Exceptions::InvalidArgument,
      //                            "Aggregation type \"brick\" requires \"repartition: enable\" set to \"false\"");

      MUELU_SET_VAR_2LIST(paramList, defaultList, "repartition: partitioner", std::string, partName);
      TEUCHOS_TEST_FOR_EXCEPTION(partName != "zoltan" && partName != "zoltan2", Exceptions::InvalidArgument,
                                 "Invalid partitioner name: \"" << partName << "\". Valid options: \"zoltan\", \"zoltan2\"");

      bool switched = false;
      (void)switched;
#ifndef HAVE_MUELU_ZOLTAN
      if (partName == "zoltan") {
        this->GetOStream(Warnings0) << "Zoltan interface is not available, trying to switch to Zoltan2" << std::endl;
        partName = "zoltan2";
        switched = true;
      }
#endif
#ifndef HAVE_MUELU_ZOLTAN2
      if (partName == "zoltan2" && !switched) {
        this->GetOStream(Warnings0) << "Zoltan2 interface is not available, trying to switch to Zoltan" << std::endl;
        partName = "zoltan";
      }
#endif

      // RepartitionHeuristic
      auto repartheurFactory = rcp(new RepartitionHeuristicFactory());
      ParameterList repartheurParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: start level",          int, repartheurParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: min rows per proc",    int, repartheurParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: target rows per proc", int, repartheurParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: max imbalance",     double, repartheurParams);
      repartheurFactory->SetParameterList(repartheurParams);
      repartheurFactory->SetFactory("A",         manager.GetFactory("A"));
      manager.SetFactory("number of partitions", repartheurFactory);

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
      partitioner->SetFactory("A",                    manager.GetFactory("A"));
      partitioner->SetFactory("number of partitions", manager.GetFactory("number of partitions"));
      partitioner->SetFactory("Coordinates",          manager.GetFactory("Coordinates"));
      manager.SetFactory("Partition", partitioner);

      // Repartitioner
      auto repartFactory = rcp(new RepartitionFactory());
      ParameterList repartParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: print partition distribution", bool, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: remap parts",                  bool, repartParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: remap num values",              int, repartParams);
      repartFactory->SetParameterList(repartParams);
      repartFactory->SetFactory("A",                    manager.GetFactory("A"));
      repartFactory->SetFactory("number of partitions", manager.GetFactory("number of partitions"));
      repartFactory->SetFactory("Partition",            manager.GetFactory("Partition"));
      manager.SetFactory("Importer", repartFactory);
      if (reuseType != "none" && reuseType != "S" && levelID)
        keeps.push_back(keep_pair("Importer", manager.GetFactory("Importer").get()));

      // Rebalanced A
      auto newA = rcp(new RebalanceAcFactory());
      ParameterList rebAcParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: use subcommunicators", bool, rebAcParams);
      newA->SetParameterList(rebAcParams);
      newA->SetFactory("A",         manager.GetFactory("A"));
      newA->SetFactory("Importer",  manager.GetFactory("Importer"));
      manager.SetFactory("A", newA);

      // Rebalanced P
      auto newP = rcp(new RebalanceTransferFactory());
      ParameterList newPparams;
      newPparams.set("type", "Interpolation");
      if (changedPRrebalance_)
        newPparams.set("repartition: rebalance P and R", this->doPRrebalance_);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: use subcommunicators", bool, newPparams);
      newP->  SetParameterList(newPparams);
      newP->  SetFactory("Importer",    manager.GetFactory("Importer"));
      newP->  SetFactory("P",           manager.GetFactory("P"));
      if (!paramList.isParameter("semicoarsen: number of levels"))
        newP->SetFactory("Nullspace",   manager.GetFactory("Ptent"));
      else
        newP->SetFactory("Nullspace",   manager.GetFactory("P")); // TogglePFactory
      newP->  SetFactory("Coordinates", manager.GetFactory("Coordinates"));
      manager.SetFactory("P",           newP);
      manager.SetFactory("Coordinates", newP);

      // Rebalanced R
      auto newR = rcp(new RebalanceTransferFactory());
      ParameterList newRparams;
      newRparams.set("type", "Restriction");
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "repartition: use subcommunicators", bool, newRparams);
      if (changedPRrebalance_)
        newRparams.set("repartition: rebalance P and R", this->doPRrebalance_);
      if (changedImplicitTranspose_)
        newRparams.set("transpose: use implicit",        this->implicitTranspose_);
      newR->  SetParameterList(newRparams);
      newR->  SetFactory("Importer", manager.GetFactory("Importer"));
      if (!this->implicitTranspose_) {
        newR->SetFactory("R",        manager.GetFactory("R"));
        manager.SetFactory("R",      newR);
      }

      // NOTE: the role of NullspaceFactory is to provide nullspace on the finest
      // level if a user does not do that. For all other levels it simply passes
      // nullspace from a real factory to whoever needs it. If we don't use
      // repartitioning, that factory is "TentativePFactory"; if we do, it is
      // "RebalanceTransferFactory". But we still have to have NullspaceFactory as
      // the "Nullspace" of the manager
      // NOTE: This really needs to be set on the *NullSpaceFactory*, not manager.get("Nullspace").
      nullSpaceFactory->SetFactory("Nullspace", newP);
#else
      throw Exceptions::RuntimeError("No repartitioning available for a serial run");
#endif
    }
  }

  // =====================================================================================================
  // =========================================== Nullspace ===============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Nullspace(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                                 int levelID, std::vector<keep_pair>& keeps, RCP<Factory> & nullSpaceFactory) const
  {
    // Nullspace
    MUELU_KOKKOS_FACTORY(nullSpace, NullspaceFactory, NullspaceFactory_kokkos);

    bool have_userNS = false;
    if (paramList.isParameter("Nullspace") && !paramList.get<RCP<MultiVector> >("Nullspace").is_null())
      have_userNS = true;

    if (!have_userNS) {
      nullSpace->SetFactory("Nullspace", manager.GetFactory("Ptent"));
      manager.SetFactory("Nullspace", nullSpace);
    }
    nullSpaceFactory = nullSpace;
  }

  // =====================================================================================================
  // ================================= Algorithm: SemiCoarsening =========================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_SemiCoarsen(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                                   int levelID, std::vector<keep_pair>& keeps) const
  {
    // === Semi-coarsening ===
    RCP<SemiCoarsenPFactory>  semicoarsenFactory = Teuchos::null;
    if (paramList.isParameter("semicoarsen: number of levels") &&
        paramList.get<int>("semicoarsen: number of levels") > 0) {

      ParameterList togglePParams;
      ParameterList semicoarsenPParams;
      ParameterList linedetectionParams;
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "semicoarsen: number of levels", int,         togglePParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "semicoarsen: coarsen rate",     int,         semicoarsenPParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "linedetection: orientation",    std::string, linedetectionParams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "linedetection: num layers",     int,         linedetectionParams);

      semicoarsenFactory                             = rcp(new SemiCoarsenPFactory());
      RCP<LineDetectionFactory> linedetectionFactory = rcp(new LineDetectionFactory());
      RCP<TogglePFactory>       togglePFactory       = rcp(new TogglePFactory());

      linedetectionFactory->SetParameterList(linedetectionParams);
      semicoarsenFactory  ->SetParameterList(semicoarsenPParams);
      togglePFactory      ->SetParameterList(togglePParams);

      togglePFactory->AddCoarseNullspaceFactory (semicoarsenFactory);
      togglePFactory->AddProlongatorFactory     (semicoarsenFactory);
      togglePFactory->AddPtentFactory           (semicoarsenFactory);
      togglePFactory->AddCoarseNullspaceFactory (manager.GetFactory("Ptent"));
      togglePFactory->AddProlongatorFactory     (manager.GetFactory("P"));
      togglePFactory->AddPtentFactory           (manager.GetFactory("Ptent"));

      manager.SetFactory("CoarseNumZLayers",          linedetectionFactory);
      manager.SetFactory("LineDetection_Layers",      linedetectionFactory);
      manager.SetFactory("LineDetection_VertLineIds", linedetectionFactory);

      manager.SetFactory("P",         togglePFactory);
      manager.SetFactory("Ptent",     togglePFactory);
      manager.SetFactory("Nullspace", togglePFactory);
    }


    if (paramList.isParameter("semicoarsen: number of levels")) {
      auto tf = rcp(new ToggleCoordinatesTransferFactory());
      tf->SetFactory("Chosen P", manager.GetFactory("P"));
      tf->AddCoordTransferFactory(semicoarsenFactory);

      MUELU_KOKKOS_FACTORY(coords, CoordinatesTransferFactory, CoordinatesTransferFactory_kokkos);
      coords->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
      coords->SetFactory("CoarseMap",  manager.GetFactory("CoarseMap"));
      tf->AddCoordTransferFactory(coords);
      manager.SetFactory("Coordinates", tf);
    }
  }


  // =====================================================================================================
  // ================================== Algorithm: P-Coarsening ==========================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_PCoarsen(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                                int levelID, std::vector<keep_pair>& keeps) const
  {
#ifdef HAVE_MUELU_INTREPID2
    // This only makes sense to invoke from the default list.
    if (defaultList.isParameter("pcoarsen: schedule") && defaultList.isParameter("pcoarsen: element")) {
      // P-Coarsening by schedule (new interface)
      // NOTE: levelID represents the *coarse* level in this case
      auto pcoarsen_schedule = Teuchos::getArrayFromStringParameter<int>(defaultList,"pcoarsen: schedule");
      auto pcoarsen_element  = defaultList.get<std::string>("pcoarsen: element");

      if (levelID >= (int)pcoarsen_schedule.size()) {
        // Past the p-coarsening levels, we do Smoothed Aggregation
        // NOTE: We should probably consider allowing other options past p-coarsening
        UpdateFactoryManager_SA(paramList, defaultList, manager, levelID, keeps);

      } else {
        // P-Coarsening
        ParameterList Pparams;
        auto P = rcp(new IntrepidPCoarsenFactory());
        std::string lo = pcoarsen_element + std::to_string(pcoarsen_schedule[levelID]);
        std::string hi = (levelID ? pcoarsen_element + std::to_string(pcoarsen_schedule[levelID-1]) : lo);
        Pparams.set("pcoarsen: hi basis", hi);
        Pparams.set("pcoarsen: lo basis", lo);
        P->SetParameterList(Pparams);
        manager.SetFactory("P", P);

        // Add special nullspace handling
        rcp_dynamic_cast<Factory>(manager.GetFactoryNonConst("Nullspace"))->SetFactory("Nullspace", manager.GetFactory("P"));
      }

    } else {
      // P-Coarsening by manual specification (old interface)
      ParameterList Pparams;
      auto P = rcp(new IntrepidPCoarsenFactory());
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "pcoarsen: hi basis", std::string, Pparams);
      MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "pcoarsen: lo basis", std::string, Pparams);
      P->SetParameterList(Pparams);
      manager.SetFactory("P", P);

      // Add special nullspace handling
      rcp_dynamic_cast<Factory>(manager.GetFactoryNonConst("Nullspace"))->SetFactory("Nullspace", manager.GetFactory("P"));
    }

#endif
  }

  // =====================================================================================================
  // ============================== Algorithm: Smoothed Aggregation ======================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_SA(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager, int levelID, std::vector<keep_pair>& keeps) const {
    // Smoothed aggregation
    MUELU_KOKKOS_FACTORY(P, SaPFactory, SaPFactory_kokkos);
    ParameterList Pparams;
    if (paramList.isSublist("matrixmatrix: kernel params"))
      Pparams.sublist("matrixmatrix: kernel params", false) = paramList.sublist("matrixmatrix: kernel params");
    if (defaultList.isSublist("matrixmatrix: kernel params"))
      Pparams.sublist("matrixmatrix: kernel params", false) = defaultList.sublist("matrixmatrix: kernel params");
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "sa: damping factor", double, Pparams);
    P->SetParameterList(Pparams);

    // Filtering
    MUELU_SET_VAR_2LIST(paramList, defaultList, "sa: use filtered matrix", bool, useFiltering);
    if (useFiltering) {
      // NOTE: Here, non-Kokkos and Kokkos versions diverge in the way the
      // dependency tree is setup. The Kokkos version has merged the the
      // FilteredAFactory into the CoalesceDropFactory.
      if (!useKokkos_) {
        RCP<Factory> filterFactory = rcp(new FilteredAFactory());

        ParameterList fParams;
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "filtered matrix: use lumping",      bool, fParams);
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "filtered matrix: reuse graph",      bool, fParams);
        MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "filtered matrix: reuse eigenvalue", bool, fParams);
        filterFactory->SetParameterList(fParams);
        filterFactory->SetFactory("Graph",      manager.GetFactory("Graph"));
        // I'm not sure why we need this line. See comments for DofsPerNode for UncoupledAggregation above
        filterFactory->SetFactory("Filtering",  manager.GetFactory("Graph"));

        P->SetFactory("A", filterFactory);

      } else {
        P->SetFactory("A", manager.GetFactory("Graph"));
      }
    }

    P->SetFactory("P", manager.GetFactory("Ptent"));
    manager.SetFactory("P", P);

    bool filteringChangesMatrix = useFiltering && !MUELU_TEST_PARAM_2LIST(paramList, defaultList, "aggregation: drop tol", double, 0);
    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type", std::string, reuseType);
    if (reuseType == "tP" && !filteringChangesMatrix)
      keeps.push_back(keep_pair("AP reuse data", P.get()));
  }

  // =====================================================================================================
  // =============================== Algorithm: Energy Minimization ======================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Emin(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                            int levelID, std::vector<keep_pair>& keeps) const
  {
    MUELU_SET_VAR_2LIST(paramList, defaultList, "emin: pattern", std::string, patternType);
    MUELU_SET_VAR_2LIST(paramList, defaultList, "reuse: type", std::string, reuseType);
    TEUCHOS_TEST_FOR_EXCEPTION(patternType != "AkPtent", Exceptions::InvalidArgument,
                               "Invalid pattern name: \"" << patternType << "\". Valid options: \"AkPtent\"");
    // Pattern
    auto patternFactory = rcp(new PatternFactory());
    ParameterList patternParams;
    MUELU_TEST_AND_SET_PARAM_2LIST(paramList, defaultList, "emin: pattern order", int, patternParams);
    patternFactory->SetParameterList(patternParams);
    patternFactory->SetFactory("P", manager.GetFactory("Ptent"));
    manager.SetFactory("Ppattern", patternFactory);

    // Constraint
    auto constraintFactory = rcp(new ConstraintFactory());
    constraintFactory->SetFactory("Ppattern",        manager.GetFactory("Ppattern"));
    constraintFactory->SetFactory("CoarseNullspace", manager.GetFactory("Ptent"));
    manager.SetFactory("Constraint", constraintFactory);

    // Energy minimization
    auto P = rcp(new EminPFactory());
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
  }

  // =====================================================================================================
  // ================================= Algorithm: Petrov-Galerkin ========================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_PG(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                          int levelID, std::vector<keep_pair>& keeps) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(this->implicitTranspose_, Exceptions::RuntimeError,
        "Implicit transpose not supported with Petrov-Galerkin smoothed transfer operators: Set \"transpose: use implicit\" to false!\n" \
        "Petrov-Galerkin transfer operator smoothing for non-symmetric problems requires a separate handling of the restriction operator which " \
        "does not allow the usage of implicit transpose easily.");

    // Petrov-Galerkin
    auto P = rcp(new PgPFactory());
    P->SetFactory("P", manager.GetFactory("Ptent"));
    manager.SetFactory("P", P);
  }


  // =====================================================================================================
  // ====================================== Algorithm: Matlab ============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  UpdateFactoryManager_Matlab(ParameterList& paramList, const ParameterList& defaultList, FactoryManager& manager,
                              int levelID, std::vector<keep_pair>& keeps) const {
#ifdef HAVE_MUELU_MATLAB
    ParameterList Pparams = paramList.sublist("transfer: params");
    auto P = rcp(new TwoLevelMatlabFactory());
    P->SetParameterList(Pparams);
    P->SetFactory("P", manager.GetFactory("Ptent"));
    manager.SetFactory("P", P);
#endif
  }

#undef MUELU_SET_VAR_2LIST
#undef MUELU_TEST_AND_SET_VAR
#undef MUELU_TEST_AND_SET_PARAM_2LIST
#undef MUELU_TEST_PARAM_2LIST
#undef MUELU_KOKKOS_FACTORY

  size_t LevenshteinDistance(const char* s, size_t len_s, const char* t, size_t len_t);

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Validate(const ParameterList& constParamList) const {
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
#ifdef HAVE_MUELU_MATLAB
    // If Muemex is supported, hide custom level variables from validator by removing them from paramList's sublists
    for (size_t i = 0; i < paramLists.size(); i++) {
      std::vector<std::string> customVars; // list of names (keys) to be removed from list

      for(Teuchos::ParameterList::ConstIterator it = paramLists[i].begin(); it != paramLists[i].end(); it++) {
        std::string paramName = paramLists[i].name(it);

        if (IsParamMuemexVariable(paramName))
          customVars.push_back(paramName);
      }

      // Remove the keys
      for (size_t j = 0; j < customVars.size(); j++)
        paramLists[i].remove(customVars[j], false);
    }
#endif

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

        size_t bestScore = 100;
        std::string bestName  = "";
        for (ParameterList::ConstIterator it = validList.begin(); it != validList.end(); it++) {
          const std::string& pName = validList.name(it);
          this->GetOStream(Runtime1) << "| " << pName;
          size_t score = LevenshteinDistance(name.c_str(), name.length(), pName.c_str(), pName.length());
          this->GetOStream(Runtime1) << " -> " << score << std::endl;
          if (score < bestScore) {
            bestScore = score;
            bestName  = pName;
          }
        }
        if (bestScore < 10 && bestName != "") {
          TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameterName,
            eString << "The parameter name \"" + name + "\" is not valid. Did you mean \"" + bestName << "\"?\n");

        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameterName,
            eString << "The parameter name \"" + name + "\" is not valid.\n");
        }
      }
    }
  }

  // =====================================================================================================
  // ==================================== FACTORY interpreter ============================================
  // =====================================================================================================
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  SetFactoryParameterList(const ParameterList& constParamList) {
    // Create a non const copy of the parameter list
    // Working with a modifiable list is much much easier than with original one
    ParameterList paramList = constParamList;

    // Parameter List Parsing:
    // ---------
    //   <ParameterList name="MueLu">
    //     <ParameterList name="Matrix">
    //   </ParameterList>
    if (paramList.isSublist("Matrix")) {
      blockSize_ = paramList.sublist("Matrix").get<int>("PDE equations", MasterList::getDefault<int>("number of equations"));
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

      if (hieraList.isParameter("coarse grid correction scaling factor")) {
        this->scalingFactor_ = hieraList.get<double>("coarse grid correction scaling factor");
        hieraList.remove("coarse grid correction scaling factor");
      }

      // Translate cycle type parameter
      if (hieraList.isParameter("cycle type")) {
        std::map<std::string, CycleType> cycleMap;
        cycleMap["V"] = VCYCLE;
        cycleMap["W"] = WCYCLE;

        std::string cycleType = hieraList.get<std::string>("cycle type");
        TEUCHOS_TEST_FOR_EXCEPTION(cycleMap.count(cycleType) == 0, Exceptions::RuntimeError, "Invalid cycle type: \"" << cycleType << "\"");
        this->Cycle_ = cycleMap[cycleType];
      }

      //TODO Move this its own class or MueLu::Utils?
      std::map<std::string, MsgType> verbMap;
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
        //TODO This would requires something like std::set<dataName, Array<int> >
        ParameterList foo = hieraList.sublist("DataToWrite");
        std::string dataName = "Matrices";
        if (foo.isParameter(dataName))
          this->matricesToPrint_ = Teuchos::getArrayFromStringParameter<int>(foo, dataName);
        dataName = "Prolongators";
        if (foo.isParameter(dataName))
          this->prolongatorsToPrint_ = Teuchos::getArrayFromStringParameter<int>(foo, dataName);
        dataName = "Restrictors";
        if (foo.isParameter(dataName))
          this->restrictorsToPrint_ = Teuchos::getArrayFromStringParameter<int>(foo, dataName);
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

          RCP<FactoryManager> m = rcp(new FactoryManager(levelFactoryMap));
          if (hieraList.isParameter("use kokkos refactor"))
            m->SetKokkosRefactor(hieraList.get<bool>("use kokkos refactor"));

          if (startLevel >= 0)
            this->AddFactoryManager(startLevel, numDesiredLevel, m);
          else
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::ParameterListInterpreter():: invalid level id");
        } /* TODO: else { } */
      }
    }
  }


  //TODO: static?
  /// \brief Interpret "Factories" sublist
  ///
  /// \param paramList [in]: "Factories" ParameterList
  /// \param factoryMapIn [in]: FactoryMap maps variable names to factories. This factory map is used to resolve data dependencies of previously defined factories.
  /// \param factoryMapOut [out]: FactoryMap maps variable names to factories. New factory entries are added to that FactoryMap. Usually, factoryMapIn and factoryMapOut should use the same object, such that new factories are added. We have to distinguish input and output if we build sub-factory managers, though.
  /// \param factoryManagers [in/out]: FacotryManagerMap maps group names to a FactoryManager object.
  ///
  /// Interpret "Factories" parameter list. For each "factory" entry, add a new entry in the factoryMapOut map or create a new FacotryManager
  ///
  /// Parameter List Parsing:
  /// Create an entry in factoryMapOut for each parameter of the list paramList
  /// ---------
  ///    <ParameterList name="...">
  ///     <Parameter name="smootherFact0" type="string" value="TrilinosSmoother"/>
  ///
  ///     <ParameterList name="smootherFact1">
  ///       <Parameter name="type" type="string" value="TrilinosSmoother"/>
  ///       ...
  ///     </ParameterList>
  ///    </ParameterList>
  ///
  /// ---------
  /// Group factories
  /// We can group factories using parameter sublists with the "group" parameter
  ///
  ///    <ParameterList name="myFirstGroup">
  ///      <Parameter name="group" type="string" value="FactoryManager"/>
  ///      <Parameter name="A" type="string" value="mySubBlockAFactory1"/>
  ///      <Parameter name="P" type="string" value="myTentativePFact1"/>
  ///      <Parameter name="Aggregates" type="string" value="myAggFact1"/>
  ///      <Parameter name="Nullspace" type="string" value="myNspFact1"/>
  ///      <Parameter name="CoarseMap" type="string" value="myCoarseMap1"/>
  ///    </ParameterList>

  ///    <ParameterList name="mySecondGroup">
  ///      <Parameter name="group" type="string" value="FactoryManager"/>
  ///      <Parameter name="A" type="string" value="mySubBlockAFactory2"/>
  ///      <Parameter name="P" type="string" value="myTentativePFact2"/>
  ///      <Parameter name="Aggregates" type="string" value="myAggFact1"/><!-- reuse aggs -->
  ///      <Parameter name="Nullspace" type="string" value="myNspFact2"/>
  ///      <Parameter name="CoarseMap" type="string" value="myCoarseMap2"/>
  ///    </ParameterList>
  ///
  ///  These factory groups can be used with factories for blocked operators (such as the BlockedPFactory)
  ///  to easily define the operations on the sub-blocks.
  ///
  ///    <ParameterList name="myBlockedPFact">
  ///      <Parameter name="factory" type="string" value="BlockedPFactory"/>
  ///      <!-- factory manager for block 1 -->
  ///      <ParameterList name="block1">
  ///        <Parameter name="group" type="string" value="myFirstGroup"/>
  ///      </ParameterList>
  ///      <!-- factory manager for block 2 -->
  ///      <ParameterList name="block2">
  ///        <Parameter name="group" type="string" value="mySecondGroup"/>
  ///      </ParameterList>
  ///    </ParameterList>
  ///
  ///
  ///  As an alternative one can also directly specify the factories in the sublists "block1", "block2", etc..., of course.
  ///  But using blocks has the advantage that one can reuse them in all blocked factories.
  ///
  ///    <ParameterList name="myBlockedPFact">
  ///      <Parameter name="factory" type="string" value="BlockedPFactory"/>
  ///      <!-- factory manager for block 1 -->
  ///      <ParameterList name="block1">
  ///        <Parameter name="A" type="string" value="mySubBlockAFactory1"/>
  ///        <Parameter name="P" type="string" value="myTentativePFact1"/>
  ///        <Parameter name="Aggregates" type="string" value="myAggFact1"/>
  ///        <Parameter name="Nullspace" type="string" value="myNspFact1"/>
  ///        <Parameter name="CoarseMap" type="string" value="myCoarseMap1"/>
  ///      </ParameterList>
  ///      <!-- factory manager for block 2 -->
  ///      <ParameterList name="block2">
  ///        <Parameter name="A" type="string" value="mySubBlockAFactory2"/>
  ///        <Parameter name="P" type="string" value="myTentativePFact2"/>
  ///        <Parameter name="Aggregates" type="string" value="myAggFact1"/><!-- reuse aggs -->
  ///        <Parameter name="Nullspace" type="string" value="myNspFact2"/>
  ///        <Parameter name="CoarseMap" type="string" value="myCoarseMap2"/>
  ///      </ParameterList>
  ///    </ParameterList>
  ///
  ///  As an alternative one can also directly specify the factories in the sublists "block1", "block2", etc..., of course.
  ///
  ///

  /// ---------
  /// add more dependencies (circular dependencies)
  ///
  ///  The NullspaceFactory needs to know which factory generates the null space on the coarse level (e.g., the TentativePFactory or the RebalancedPFactory).
  ///  However, we cannot set the information in this place in the xml file, since the tentative prolongator facotry is typically defined later.
  ///  We have to add that dependency later to the NullspaceFactory:
  ///
  ///    <ParameterList name="myNspFact">
  ///      <Parameter name="factory" type="string" value="NullspaceFactory"/>
  ///      <!--<Paramter name="Nullspace" type="string" value="myRebalanceProlongatorFact"/>-->
  ///    </ParameterList>
  ///
  ///    <ParameterList name="myTentativePFact">
  ///      <Parameter name="factory"   type="string" value="TentativePFactory"/>
  ///      <...>
  ///      <Parameter name="Nullspace" type="string" value="myNspFact"/>
  ///      <Parameter name="CoarseMap" type="string" value="myCoarseMap"/>
  ///    </ParameterList>
  ///
  ///    <ParameterList name="myRebalanceProlongatorFact">
  ///      <Parameter name="factory"   type="string" value="RebalanceTransferFactory"/>
  ///      <...>
  ///      <Parameter name="Nullspace" type="string" value="myTentativePFact"/>
  ///    </ParameterList>
  ///
  ///  After the definition of the generating factory for the nullspace (in this case myRebalanceProlongatorFact)
  ///  we add that dependency to the NullspaceFactory instance myNspFact
  ///
  ///    <ParameterList name="myNspFactDeps">
  ///      <Parameter name="dependency for" type="string" value="myNspFact"/>
  ///      <Parameter name="Nullspace" type="string" value="myRebalanceProlongatorFact"/>
  ///    </ParameterList>
  ///
  ///  We have to create a new block (with a different name than myNspFact). In the example we use "myNspFactDeps".
  ///  It should contain a parameter "dependency for" with the name of the factory that we want the dependencies to be addded to.
  ///  With above block we do not need the entry for the Nullspace in the global FactoryManager any more.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildFactoryMap(const ParameterList& paramList, const FactoryMap& factoryMapIn, FactoryMap& factoryMapOut, FactoryManagerMap& factoryManagers) const {
    for (ParameterList::ConstIterator param = paramList.begin(); param != paramList.end(); ++param) {
      const std::string             & paramName  = paramList.name(param);  //< paramName contains the user chosen factory name (e.g., "smootherFact1")
      const Teuchos::ParameterEntry & paramValue = paramList.entry(param); //< for factories, paramValue should be either a list or just a MueLu Factory (e.g., TrilinosSmoother)

      //TODO: do not allow name of existing MueLu classes (can be tested using FactoryFactory)

      if (paramValue.isList()) {
        ParameterList paramList1 = Teuchos::getValue<ParameterList>(paramValue);
        if (paramList1.isParameter("factory")) { // default: just a factory definition
          // New Factory is a sublist with internal parameters and/or data dependencies
          TEUCHOS_TEST_FOR_EXCEPTION(paramList1.isParameter("dependency for") == true, Exceptions::RuntimeError,
            "MueLu::ParameterListInterpreter(): It seems that in the parameter lists for defining " << paramName <<
            " there is both a 'factory' and 'dependency for' parameter. This is not allowed. Please remove the 'dependency for' parameter.");

          factoryMapOut[paramName] = factFact_->BuildFactory(paramValue, factoryMapIn, factoryManagers);

        } else if (paramList1.isParameter("dependency for")) { // add more data dependencies to existing factory
          TEUCHOS_TEST_FOR_EXCEPTION(paramList1.isParameter("factory") == true, Exceptions::RuntimeError,
            "MueLu::ParameterListInterpreter(): It seems that in the parameter lists for defining " << paramName <<
            " there is both a 'factory' and 'dependency for' parameter. This is not allowed.");

          std::string factoryName = paramList1.get<std::string>("dependency for");

          RCP<const FactoryBase> factbase = factoryMapIn.find(factoryName /*paramName*/)->second; // access previously defined factory
          TEUCHOS_TEST_FOR_EXCEPTION(factbase.is_null() == true, Exceptions::RuntimeError,
            "MueLu::ParameterListInterpreter(): could not find factory " + factoryName + " in factory map. Did you define it before?");

          RCP<const Factory> factoryconst = Teuchos::rcp_dynamic_cast<const Factory>(factbase);
          RCP<      Factory> factory      = Teuchos::rcp_const_cast<Factory>(factoryconst);

          // Read the RCP<Factory> parameters of the class T
          RCP<const ParameterList> validParamList = factory->GetValidParameterList();
          for (ParameterList::ConstIterator vparam = validParamList->begin(); vparam != validParamList->end(); ++vparam) {
            const std::string& pName = validParamList->name(vparam);

            if (!paramList1.isParameter(pName)) {
              // Ignore unknown parameters
              continue;
            }

            if (validParamList->isType< RCP<const FactoryBase> >(pName)) {
              // Generate or get factory described by pName and set dependency
              RCP<const FactoryBase> generatingFact = factFact_->BuildFactory(paramList1.getEntry(pName), factoryMapIn, factoryManagers);
              factory->SetFactory(pName, generatingFact.create_weak());

            } else if (validParamList->isType<RCP<const ParameterList> >(pName)) {
              if (pName == "ParameterList") {
                // NOTE: we cannot use
                //     subList = sublist(rcpFromRef(paramList), pName)
                // here as that would result in sublist also being a reference to a temporary object.
                // The resulting dereferencing in the corresponding factory would then segfault
                RCP<const ParameterList> subList = Teuchos::sublist(rcp(new ParameterList(paramList1)), pName);
                factory->SetParameter(pName, ParameterEntry(subList));
              }
            } else {
              factory->SetParameter(pName, paramList1.getEntry(pName));
            }
          }

        } else if (paramList1.isParameter("group")) { // definitiion of a factory group (for a factory manager)
          // Define a new (sub) FactoryManager
          std::string groupType = paramList1.get<std::string>("group");
          TEUCHOS_TEST_FOR_EXCEPTION(groupType!="FactoryManager", Exceptions::RuntimeError,
                                     "group must be of type \"FactoryManager\".");

          ParameterList groupList = paramList1; // copy because list temporally modified (remove 'id')
          groupList.remove("group");

          FactoryMap groupFactoryMap;
          BuildFactoryMap(groupList, factoryMapIn, groupFactoryMap, factoryManagers);

          // do not store groupFactoryMap in factoryMapOut
          // Create a factory manager object from groupFactoryMap
          RCP<FactoryManagerBase> m = rcp(new FactoryManager(groupFactoryMap));

          factoryManagers[paramName] = m;

        } else {
          this->GetOStream(Warnings0) << "Could not interpret parameter list " << paramList1 << std::endl;
          TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError,
                                     "XML Parameter list must either be of type \"factory\" or of type \"group\".");
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
            << "instead of " << A.GetFixedBlockSize() << " (provided matrix)." << std::endl
            << "You may want to check \"number of equations\" (or \"PDE equations\" for factory style list) parameter." << std::endl;

      A.SetFixedBlockSize(blockSize_, dofOffset_);

    } catch (std::bad_cast& e) {
      this->GetOStream(Warnings0) << "Skipping setting block size as the operator is not a matrix" << std::endl;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetupHierarchy(Hierarchy& H) const {
    H.SetCycle(Cycle_);
    H.SetProlongatorScalingFactor(scalingFactor_);
    HierarchyManager::SetupHierarchy(H);
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

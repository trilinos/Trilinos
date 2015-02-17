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
#ifndef MUELU_FACTORYFACTORY_DECL_HPP
#define MUELU_FACTORYFACTORY_DECL_HPP

#include <string>
#include <vector>

#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_Array.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryFactory_fwd.hpp"

#include "MueLu_HierarchyFactory.hpp"

#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Hierarchy_fwd.hpp"

#include "MueLu_Monitor.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_AggregationExportFactory.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#ifdef HAVE_MUELU_EXPERIMENTAL
#include "MueLu_BlockedCoarseMapFactory.hpp"
#include "MueLu_BlockedDirectSolver.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_BraessSarazinSmoother.hpp"
#endif
#include "MueLu_BrickAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EminPFactory.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#ifdef HAVE_MUELU_EXPERIMENTAL
#include "MueLu_IndefBlockedDiagonalSmoother.hpp"
#endif
#include "MueLu_IsorropiaInterface.hpp"
#include "MueLu_RepartitionInterface.hpp"
#include "MueLu_MapTransferFactory.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_NullspacePresmoothFactory.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#ifdef HAVE_MUELU_EXPERIMENTAL
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_SimpleSmoother.hpp"
#endif
#include "MueLu_SmootherFactory.hpp"
#ifdef HAVE_MUELU_EXPERIMENTAL
#include "MueLu_SubBlockAFactory.hpp"
#endif
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_UserAggregationFactory.hpp"
#include "MueLu_UserPFactory.hpp"
#include "MueLu_SemiCoarsenPFactory.hpp"
#ifdef HAVE_MUELU_EXPERIMENTAL
#include "MueLu_UzawaSmoother.hpp"
#endif
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"

namespace MueLu {

  /*! class FactoryFactory

  @brief Factory that can generate other factories from


  */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class FactoryFactory : public BaseClass {
#undef MUELU_FACTORYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef std::map<std::string, RCP<const FactoryBase>  > FactoryMap; // TODO: remove
    typedef std::map<std::string, RCP<FactoryManagerBase> > FactoryManagerMap;

  public:

    // Parameter List Parsing:
    // ---------
    //     <Parameter name="smootherFact0" type="string" value="TrilinosSmoother"/>
    //
    // or:
    //
    //     <ParameterList name="smootherFact1">
    //       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
    //       ...
    //     </ParameterList>
    //
    virtual RCP<const FactoryBase> BuildFactory(const Teuchos::ParameterEntry& param, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      // Find factory
      std::string            factoryName;
      Teuchos::ParameterList paramList;
      if (!param.isList()) {
        factoryName = Teuchos::getValue<std::string>(param);
      } else {
        paramList = Teuchos::getValue<Teuchos::ParameterList>(param);
        factoryName = paramList.get<std::string>("factory");
      }

      // TODO: see how Teko handles this (=> register factories).
      if (factoryName == "AggregationExportFactory")        return Build2<AggregationExportFactory>     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "AmalgamationFactory")             return Build2<AmalgamationFactory>          (paramList, factoryMapIn, factoryManagersIn);
#ifdef HAVE_MUELU_EXPERIMENTAL
      if (factoryName == "BlockedCoarseMapFactory")         return Build2<BlockedCoarseMapFactory>      (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedRAPFactory")               return BuildRAPFactory<BlockedRAPFactory>   (paramList, factoryMapIn, factoryManagersIn);
#endif
#if defined(HAVE_MPI)
      if (factoryName == "BrickAggregationFactory")         return Build2<BrickAggregationFactory>      (paramList, factoryMapIn, factoryManagersIn);
#endif
      if (factoryName == "CoarseMapFactory")                return Build2<CoarseMapFactory>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoalesceDropFactory")             return Build2<CoalesceDropFactory>          (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "ConstraintFactory")               return Build2<ConstraintFactory>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoupledAggregationFactory")       return BuildCoupledAggregationFactory       (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "CoordinatesTransferFactory")      return Build2<CoordinatesTransferFactory>   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "DirectSolver")                    return BuildDirectSolver                    (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "EminPFactory")                    return Build2<EminPFactory>                 (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "FilteredAFactory")                return Build2<FilteredAFactory>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "GenericRFactory")                 return Build2<GenericRFactory>              (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "MapTransferFactory")              return Build2<MapTransferFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "MultiVectorTransferFactory")      return Build2<MultiVectorTransferFactory>   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "NoFactory")                       return Teuchos::null;
      if (factoryName == "NullspaceFactory")                return Build2<NullspaceFactory>             (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "NullspacePresmoothFactory")       return Build2<NullspacePresmoothFactory>    (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "PatternFactory")                  return Build2<PatternFactory>               (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "PgPFactory")                      return Build2<PgPFactory>                   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SaPFactory")                      return Build2<SaPFactory>                   (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RAPFactory")                      return BuildRAPFactory<RAPFactory>          (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RebalanceAcFactory")              return Build2<RebalanceAcFactory>           (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RebalanceTransferFactory")        return Build2<RebalanceTransferFactory>     (paramList, factoryMapIn, factoryManagersIn);
#ifdef HAVE_MUELU_EXPERIMENTAL
      if (factoryName == "SubBlockAFactory")                return Build2<SubBlockAFactory>             (paramList, factoryMapIn, factoryManagersIn);
#endif
      if (factoryName == "TentativePFactory")               return Build2<TentativePFactory>            (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "TransPFactory")                   return Build2<TransPFactory>                (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "TrilinosSmoother")                return BuildTrilinosSmoother                (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UncoupledAggregationFactory")     return BuildUncoupledAggregationFactory     (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UserAggregationFactory")          return Build2<UserAggregationFactory>       (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UserPFactory")                    return Build2<UserPFactory>                 (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SemiCoarsenPFactory")             return Build2<SemiCoarsenPFactory>          (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "RepartitionInterface")            return Build2<RepartitionInterface>         (paramList, factoryMapIn, factoryManagersIn);

      if (factoryName == "ZoltanInterface") {
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
        return Build2<ZoltanInterface>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a ZoltanInterface object: Zoltan is disabled: HAVE_MUELU_ZOLTAN && HAVE_MPI == false.");
#endif // HAVE_MUELU_ZOLTAN && HAVE_MPI
      }
      if (factoryName == "Zoltan2Interface") {
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)
        return Build2<Zoltan2Interface>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a Zoltan2Interface object: Zoltan2 is disabled: HAVE_MUELU_ZOLTAN2 && HAVE_MPI == false.");
#endif // HAVE_MUELU_ZOLTAN2 && HAVE_MPI
      }
      if (factoryName == "IsorropiaInterface") {
#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)
        return Build2<IsorropiaInterface>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a IsorropiaInterface object: Isorropia is disabled: HAVE_MUELU_ISORROPIA && HAVE_MPI == false.");
#endif // HAVE_MUELU_ZOLTAN2 && HAVE_MPI
      }

      if (factoryName == "RepartitionFactory") {
#ifdef HAVE_MPI
        return Build2<RepartitionFactory>(paramList, factoryMapIn, factoryManagersIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a RepartitionFactory object: HAVE_MPI == false.");
#endif // HAVE_MPI
      }
      // Blocked factories
#ifdef HAVE_MUELU_EXPERIMENTAL
      if (factoryName == "BlockedDirectSolver")             return BuildBlockedDirectSolver(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedGaussSeidelSmoother")      return BuildBlockedSmoother<BlockedGaussSeidelSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BlockedPFactory")                 return BuildBlockedPFactory(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "BraessSarazinSmoother")           return BuildBlockedSmoother<BraessSarazinSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "IndefiniteBlockDiagonalSmoother") return BuildBlockedSmoother<IndefBlockedDiagonalSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SimpleSmoother")                  return BuildBlockedSmoother<SimpleSmoother>(paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "SchurComplementFactory")          return Build2<SchurComplementFactory> (paramList, factoryMapIn, factoryManagersIn);
      if (factoryName == "UzawaSmoother")                   return BuildBlockedSmoother<UzawaSmoother>(paramList, factoryMapIn, factoryManagersIn);
#endif
      // Use a user defined factories (in <Factories> node)
      if (factoryMapIn.find(factoryName) != factoryMapIn.end()) {
        TEUCHOS_TEST_FOR_EXCEPTION((param.isList() && (++paramList.begin() != paramList.end())), Exceptions::RuntimeError,
                                   "MueLu::FactoryFactory: Error during the parsing of: " << std::endl << paramList << std::endl
                                   << "'" << factoryName << "' is not a factory name but an existing instance of a factory." << std::endl
                                   << "Extra parameters cannot be specified after the creation of the object." << std::endl << std::endl
                                   << "Correct syntaxes includes:" << std::endl
                                   << " <Parameter name=\"...\" type=\"string\" value=\"" << factoryName << "\"/>" << std::endl
                                   << "or" << std::endl
                                   << " <ParameterList name=\"...\"><Parameter name=\"factory\" type=\"string\" value=\"" << factoryName << "\"/></ParameterList>" << std::endl
                                   );

        return factoryMapIn.find(factoryName)->second;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory: unknown factory name : " << factoryName);

      return Teuchos::null;
    }

    //
    //
    //

    // FOLLOWING FUNCTIONS SHOULD LIVE WITH THE CORRESPONDING CLASS

    //
    //
    //

#define arraysize(ar)  (sizeof(ar) / sizeof(ar[0]))

    template <class T> // T must implement the Factory interface
    RCP<T> Build(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> factory = rcp(new T());

      const char* strarray[] = {"A", "P", "R", "Graph", "UnAmalgamationInfo", "Aggregates", "Nullspace", "TransferFactory", "DofsPerNode"};
      std::vector<std::string> v(strarray, strarray + arraysize(strarray));
      for (size_t i = 0; i < v.size(); ++i)
        if (paramList.isParameter(v[i]))
          factory->SetFactory(v[i], BuildFactory(paramList.getEntry(v[i]), factoryMapIn, factoryManagersIn));

      return factory;
    }

    template <class T> // T must implement the Factory interface
    RCP<T> Build2(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> factory = rcp(new T());

      ParameterList paramListWithFactories;

      // Read the RCP<Factory> parameters of the class T
      RCP<const ParameterList> validParamList = factory->GetValidParameterList();
      for (ParameterList::ConstIterator param = validParamList->begin(); param != validParamList->end(); ++param) {
        const std::string& pName = validParamList->name(param);

        if (!paramList.isParameter(pName)) {
          // Ignore unknown parameters
          continue;
        }

        if (validParamList->isType< RCP<const FactoryBase> >(pName)) {
          // Generate or get factory described by param
          RCP<const FactoryBase> generatingFact = BuildFactory(paramList.getEntry(pName), factoryMapIn, factoryManagersIn);
          paramListWithFactories.set(pName, generatingFact);

        } else if (validParamList->isType<RCP<const ParameterList> >(pName)) {
          if (pName == "ParameterList") {
            // NOTE: we cannot use
            //     subList = sublist(rcpFromRef(paramList), pName)
            // here as that would result in sublist also being a reference to a temporary object.
            // The resulting dereferencing in the corresponding factory would then segfault
            RCP<const ParameterList> subList = Teuchos::sublist(rcp(new ParameterList(paramList)), pName);
            paramListWithFactories.set(pName, subList);
          }
        } else {
          paramListWithFactories.setEntry(pName, paramList.getEntry(pName));
        }
      }

      // Configure the factory
      factory->SetParameterList(paramListWithFactories);

      return factory;
    }

    template <class T> // T must implement the Factory interface
    RCP<T> BuildRAPFactory(const Teuchos::ParameterList & paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<T> factory;
      if (paramList.isSublist("TransferFactories") == false) {
        factory = Build2<T>(paramList, factoryMapIn, factoryManagersIn);

      } else {
        RCP<Teuchos::ParameterList>       paramListNonConst = rcp(new Teuchos::ParameterList(paramList));
        RCP<const Teuchos::ParameterList> transferFactories = rcp(new Teuchos::ParameterList(*sublist(paramListNonConst, "TransferFactories")));

        paramListNonConst->remove("TransferFactories");

        factory = Build2<T>(*paramListNonConst, factoryMapIn, factoryManagersIn);

        for (Teuchos::ParameterList::ConstIterator param = transferFactories->begin(); param != transferFactories->end(); ++param) {
          RCP<const FactoryBase> p = BuildFactory(transferFactories->entry(param), factoryMapIn, factoryManagersIn);
          factory->AddTransferFactory(p);
        }
      }

      return factory;
    }

    //! CoupledAggregationFactory
    RCP<FactoryBase> BuildCoupledAggregationFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<CoupledAggregationFactory> factory = Build<CoupledAggregationFactory>(paramList, factoryMapIn, factoryManagersIn);

      if (paramList.isParameter("aggregation: ordering"))
        factory->SetOrdering(paramList.get<std::string>("aggregation: ordering"));

      if (paramList.isParameter("aggregation: max selected neighbors"))
        factory->SetMaxNeighAlreadySelected(paramList.get<int>("aggregation: max selected neighbors"));

      if (paramList.isParameter("Phase3AggCreation"))
        factory->SetPhase3AggCreation(paramList.get<double>("Phase3AggCreation"));

      if(paramList.isParameter("aggregation: min agg size"))
        factory->SetMinNodesPerAggregate(paramList.get<int>("aggregation: min agg size"));

      return factory;
    }

    //! UncoupledAggregationFactory
    RCP<FactoryBase> BuildUncoupledAggregationFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<UncoupledAggregationFactory> factory = Build<UncoupledAggregationFactory>(paramList, factoryMapIn, factoryManagersIn);

      ParameterList paramListWithFactories(paramList); // copy  (*might* also avoid indicating that parameter entry is used)
      paramListWithFactories.remove("factory", false);

      // Read the RCP<Factory> parameters of the class T
      RCP<const ParameterList> validParamList = factory->GetValidParameterList();
      for (ParameterList::ConstIterator param = validParamList->begin(); param != validParamList->end(); ++param) {
        const std::string & pName = validParamList->name(param);

        if (validParamList->isType< RCP<const FactoryBase> >(pName) && paramList.isParameter(pName)) {
          // Generate or get factory described by param
          RCP<const FactoryBase> generatingFact = BuildFactory(paramList.getEntry(pName), factoryMapIn, factoryManagersIn);

          // Replace <std::string> or sub-list entry by an RCP<Factory> in paramListWithFactories
          paramListWithFactories.remove(pName);
          paramListWithFactories.set(pName, generatingFact);
        }

        if (pName == "ParameterList" && validParamList->isType<RCP<const ParameterList> >(pName) && paramList.isParameter(pName)) {
          // NOTE: we cannot use
          //     subList = sublist(rcpFromRef(paramList), pName)
          // here as that would result in sublist also being a reference to a temporary object.
          // The resulting dereferencing in the corresponding factory would then segfault
          RCP<const ParameterList> subList = Teuchos::sublist(rcp(new ParameterList(paramList)), pName);
          paramListWithFactories.set(pName, subList);
        }
      }

      // Configure the factory
      factory->SetParameterList(paramListWithFactories);

      return factory;
    }

    //! TrilinosSmoother
    // Parameter List Parsing:
    //     <ParameterList name="smootherFact1">
    //       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
    //       <Parameter name="verbose" type="string" value="Warnings"/>
    //       <Parameter name="type" type="string" value="RELAXATION"/>
    //       <ParameterList name="ParameterList">
    //       ...
    //       </ParameterList>
    //     </ParameterList>
    RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new TrilinosSmoother())));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "TrilinosSmoother", Exceptions::RuntimeError, "");

      // Is it true? TEUCHOS_TEST_FOR_EXCEPTION(!paramList.isParameter("type"), Exceptions::RuntimeError, "TrilinosSmoother: parameter 'type' is mandatory");
      // type="" is default in TrilinosSmoother, but what happen then?

      std::string type="";            if(paramList.isParameter("type"))          type    = paramList.get<std::string>("type");
      int         overlap=0;          if(paramList.isParameter("overlap"))       overlap = paramList.get<int>        ("overlap");
      // std::string verbose;         if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params;  if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new TrilinosSmoother(type, params, overlap))));
    }

    RCP<FactoryBase> BuildDirectSolver(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new DirectSolver()), Teuchos::null));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "DirectSolver", Exceptions::RuntimeError, "");

      std::string type;              if(paramList.isParameter("type"))          type = paramList.get<std::string>("type");
      // std::string verbose;        if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new DirectSolver(type, params)), Teuchos::null));
    }

#ifdef HAVE_MUELU_EXPERIMENTAL
    template <class T> // T must implement the Factory interface
    RCP<FactoryBase> BuildBlockedSmoother(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      // read in sub lists
      RCP<ParameterList> paramListNonConst = rcp(new ParameterList(paramList));

      // internal vector of factory managers
      std::vector<RCP<FactoryManager> > facManagers;

      // loop over all "block%i" sublists in parameter list
      int blockid = 1;
      bool blockExists = true;
      while (blockExists == true) {
        std::stringstream ss;
        ss << "block" << blockid;

        if(paramList.isSublist(ss.str()) == true) {
          // we either have a parameter group or we have a list of factories in here
          RCP<const ParameterList> b = rcp(new ParameterList(*sublist(paramListNonConst, ss.str())));

          RCP<FactoryManager> M = Teuchos::null;

          if (b->isParameter("group")) {
            // use a factory manager
            std::string facManagerName = b->get< std::string >("group");
            TEUCHOS_TEST_FOR_EXCEPTION(factoryManagersIn.count(facManagerName) != 1, Exceptions::RuntimeError, "Factory manager has not been found. Please check the spelling of the factory managers in your xml file.");
            RCP<FactoryManagerBase> Mb = factoryManagersIn.find(facManagerName)->second;
            M = Teuchos::rcp_dynamic_cast<FactoryManager>(Mb);
            TEUCHOS_TEST_FOR_EXCEPTION(M==Teuchos::null, Exceptions::RuntimeError, "Failed to cast FactoryManagerBase object to FactoryManager.");
          } else {
            // read in the list of factories
            M = rcp(new FactoryManager());
            for (ParameterList::ConstIterator param = b->begin(); param != b->end(); ++param) {
              RCP<const FactoryBase> p = BuildFactory(b->entry(param), factoryMapIn, factoryManagersIn);
              M->SetFactory(b->name(param),p);
            }
          }

          // add factory manager to internal vector of factory managers
          M->SetIgnoreUserData(true);
          facManagers.push_back(M);
          paramListNonConst->remove(ss.str());
          blockid++;
        } else {
          blockExists = false;
          break;
        }

      }

      // create a new blocked smoother
      RCP<T> bs = Build2<T>(*paramListNonConst, factoryMapIn, factoryManagersIn);

      // important: set block factory for A here! TODO think about this in more detail
      bs->SetFactory("A", MueLu::NoFactory::getRCP());

      for (int i = 0; i<Teuchos::as<int>(facManagers.size()); i++) {
        bs->AddFactoryManager(facManagers[i],i);
      }

      return rcp(new SmootherFactory(bs));
    }

    RCP<FactoryBase> BuildBlockedDirectSolver(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      //if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new BlockedDirectSolver())));

      /*TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "DirectSolver", Exceptions::RuntimeError, "");

      std::string type;              if(paramList.isParameter("type"))          type = paramList.get<std::string>("type");
      // std::string verbose;        if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new DirectSolver(type, params))));*/
    }

    RCP<FactoryBase> BuildBlockedPFactory(const Teuchos::ParameterList& paramList, const FactoryMap& factoryMapIn, const FactoryManagerMap& factoryManagersIn) const {
      RCP<BlockedPFactory> pfac = rcp(new BlockedPFactory());

      // read in sub lists
      RCP<ParameterList> paramListNonConst = rcp(new ParameterList(paramList));

      // internal vector of factory managers
      std::vector<RCP<FactoryManager> > facManagers;

      // loop over all "block%i" sublists in parameter list
      int blockid = 1;
      bool blockExists = true;
      while (blockExists == true) {
        std::stringstream ss;
        ss << "block" << blockid;

        if(paramList.isSublist(ss.str()) == true) {
          // we either have a parameter group or we have a list of factories in here
          RCP<const ParameterList> b = rcp(new ParameterList(*sublist(paramListNonConst, ss.str())));

          RCP<FactoryManager> M = Teuchos::null;

          if (b->isParameter("group")) {
            // use a factory manager
            std::string facManagerName = b->get< std::string >("group");
            TEUCHOS_TEST_FOR_EXCEPTION(factoryManagersIn.count(facManagerName) != 1, Exceptions::RuntimeError, "Factory manager has not been found. Please check the spelling of the factory managers in your xml file.");
            RCP<FactoryManagerBase> Mb = factoryManagersIn.find(facManagerName)->second;
            M = Teuchos::rcp_dynamic_cast<FactoryManager>(Mb);
            TEUCHOS_TEST_FOR_EXCEPTION(M==Teuchos::null, Exceptions::RuntimeError, "Failed to cast FactoryManagerBase object to FactoryManager.");
          } else {
            // read in the list of factories
            M = rcp(new FactoryManager());
            for (ParameterList::ConstIterator param = b->begin(); param != b->end(); ++param) {
              RCP<const FactoryBase> p = BuildFactory(b->entry(param), factoryMapIn, factoryManagersIn);
              M->SetFactory(b->name(param),p);
            }
          }

          // add factory manager to internal vector of factory managers
          M->SetIgnoreUserData(true);
          facManagers.push_back(M);
          paramListNonConst->remove(ss.str());
          blockid++;
        } else {
          blockExists = false;
          break;
        }

      }

      // build BlockedPFactory (without sub block information)
      pfac = Build2<BlockedPFactory>(*paramListNonConst, factoryMapIn, factoryManagersIn);

      // add FactoryManager objects
      for(size_t i = 0; i<facManagers.size(); i++) {
        pfac->AddFactoryManager(facManagers[i]); // add factory manager
      }

      return pfac;
    }
#endif /* #ifdef HAVE_MUELU_EXPERIMENTAL */
  }; // class
} // namespace MueLu

#define MUELU_FACTORYFACTORY_SHORT
#endif // MUELU_FACTORYFACTORY_DECL_HPP

  // TODO: handle factory parameters
  // TODO: parameter validator
  // TODO: static
  // TODO: default parameters should not be duplicated here and on the Factory (ex: default for overlap (=0) is defined both here and on TrilinosSmoother constructors)

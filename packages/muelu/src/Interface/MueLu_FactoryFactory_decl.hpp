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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
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
#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

#include "MueLu_AmalgamationFactory.hpp" //TMP
#include "MueLu_CoalesceDropFactory.hpp" //TMP
#include "MueLu_RAPFactory.hpp" //TMP
#include "MueLu_RebalanceAcFactory.hpp" //TMP
#include "MueLu_TransPFactory.hpp" //TMP
#include "MueLu_GenericRFactory.hpp" //TMP
#include "MueLu_SaPFactory.hpp" //TMP
#include "MueLu_PgPFactory.hpp" //TMP
#include "MueLu_TrilinosSmoother.hpp" //TMP
#include "MueLu_SmootherFactory.hpp" //TMP
#include "MueLu_TentativePFactory.hpp" //TMP
#include "MueLu_UCAggregationFactory.hpp" //TMP
#include "MueLu_UncoupledAggregationFactory.hpp" //TMP
#include "MueLu_DirectSolver.hpp" //TMP
#include "MueLu_Exceptions.hpp" //TMP
#include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_RepartitionFactory.hpp"

namespace MueLu {

/*! class FactoryFactory

    @brief Factory that can generate other factories from 


*/
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class FactoryFactory {
#undef MUELU_FACTORYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; // TODO: remove

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
    RCP<const FactoryBase> BuildFactory(const Teuchos::ParameterEntry & param, const FactoryMap & factoryMapIn) const {

      //TODO: add test restricted keyword

      // Find factory
      std::string factoryName;
      Teuchos::ParameterList paramList;
      if (!param.isList()) {
        factoryName = Teuchos::getValue<std::string>(param);
      } else {
        paramList = Teuchos::getValue<Teuchos::ParameterList>(param);
        factoryName = paramList.get<std::string>("factory");
      }

      // TODO: see how Teko handles this (=> register factories).
      if (factoryName == "AmalgamationFactory") {
        return Build<AmalgamationFactory>(paramList, factoryMapIn);
      }
      if (factoryName == "CoalesceDropFactory") {
        return Build<CoalesceDropFactory>(paramList, factoryMapIn);
      }
      if (factoryName == "TentativePFactory") {
        return Build<TentativePFactory>(paramList, factoryMapIn);
      }
      if (factoryName == "SaPFactory") {
        return BuildSaPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "TransPFactory") {
        return Build<TransPFactory>(paramList, factoryMapIn);
      }
      if (factoryName == "RAPFactory") {
        return BuildRAPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "RebalanceAcFactory") {
        return Build<RebalanceAcFactory>(paramList, factoryMapIn);
      }
      if (factoryName == "UCAggregationFactory") {
        return BuildUCAggregationFactory(paramList, factoryMapIn);
      }
      if (factoryName == "UncoupledAggregationFactory") {
        return BuildUncoupledAggregationFactory(paramList, factoryMapIn);
      }
      if (factoryName == "TrilinosSmoother") {
        return BuildTrilinosSmoother(paramList, factoryMapIn);
      }
      if (factoryName == "DirectSolver") {
        return BuildDirectSolver(paramList, factoryMapIn);
      }
      if (factoryName == "MultiVectorTransferFactory") {
        return BuildMultiVectorTransferFactory(paramList, factoryMapIn);
      }
      if (factoryName == "ZoltanInterface") {
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
        return Build<ZoltanInterface>(paramList, factoryMapIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a ZoltanInterface object: Zoltan is disabled: HAVE_MUELU_ZOLTAN && HAVE_MPI == false.");
#endif // HAVE_MUELU_ZOLTAN && HAVE_MPI
      }

      if (factoryName == "RepartitionFactory") {
#ifdef HAVE_MPI
        return BuildRepartitionFactory(paramList, factoryMapIn);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory:BuildFactory(): Cannot create a RepartitionFactory object: HAVE_MPI == false.");
#endif // HAVE_MPI
      }
      if (factoryName == "RebalanceTransferFactory") {
        return BuildRebalanceTransferFactory(paramList, factoryMapIn);
      }

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
    RCP<T> Build(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      // TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "T", Exceptions::RuntimeError, "");

      RCP<T> factory = rcp(new T());

      const char* strarray[] = {"A", "P", "R", "Graph", "UnAmalgamationInfo", "Aggregates", "Nullspace", "TransferFactory"};
      std::vector<std::string> v(strarray, strarray + arraysize(strarray));
      for (std::vector<std::string>::iterator it = v.begin() ; it != v.end(); ++it) {
        if (paramList.isParameter(*it)) { factory->SetFactory(*it, BuildFactory(paramList.getEntry(*it), factoryMapIn)); }
      }

      return factory;
    }


#define MUELU_FACTORY_PARAM2(name)                                      \
    if (paramList.isParameter(name)) { factory->SetFactory(name, BuildFactory(paramList.getEntry(name), factoryMapIn)); }

    //! SaPFactory
    RCP<FactoryBase> BuildSaPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      RCP<SaPFactory> factory = Build<SaPFactory>(paramList, factoryMapIn);
      if (paramList.isParameter("DampingFactor")) factory->SetDampingFactor(paramList.get<Scalar>("DampingFactor"));

      return factory;
    }

    //! RAPFactory
    RCP<FactoryBase> BuildRAPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      RCP<RAPFactory> factory = Build<RAPFactory>(paramList, factoryMapIn);

      if (paramList.isSublist("TransferFactories")) {
        Teuchos::ParameterList transferList = paramList.sublist("TransferFactories");
        for (Teuchos::ParameterList::ConstIterator param = transferList.begin(); param != transferList.end(); ++param) {
          RCP<const FactoryBase> p = BuildFactory(transferList.entry(param), factoryMapIn);
          factory->AddTransferFactory(p);
        }
      }

      return factory;
    }

    //! UCAggregationFactory
    RCP<FactoryBase> BuildUCAggregationFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      RCP<UCAggregationFactory> factory = Build<UCAggregationFactory>(paramList, factoryMapIn);

      if(paramList.isParameter("Ordering")) {
        std::string orderingStr = paramList.get<std::string>("Ordering");
        Ordering ordering;
        if (orderingStr == "Natural")
          ordering = NATURAL;
        else if (orderingStr == "Random")
          ordering = RANDOM;
        else if (orderingStr == "Graph")
          ordering = GRAPH;
        else TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory::BuildUCAggregationFactory()::Unknown Ordering type");

        factory->SetOrdering(ordering);
      }

      if(paramList.isParameter("MaxNeighAlreadySelected")) {
        factory->SetMaxNeighAlreadySelected(paramList.get<int>("MaxNeighAlreadySelected"));
      }

      if(paramList.isParameter("Phase3AggCreation")) {
        factory->SetPhase3AggCreation(paramList.get<double>("Phase3AggCreation"));
      }

      if(paramList.isParameter("MinNodesPerAggregate")) {
        factory->SetMinNodesPerAggregate(paramList.get<int>("MinNodesPerAggregate"));
      }

      return factory;
    }

    //! UncoupledAggregationFactory
    RCP<FactoryBase> BuildUncoupledAggregationFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      RCP<UncoupledAggregationFactory> factory = Build<UncoupledAggregationFactory>(paramList, factoryMapIn);

      if(paramList.isParameter("Ordering")) {
        std::string orderingStr = paramList.get<std::string>("Ordering");
        Ordering ordering;
        if (orderingStr == "Natural")
          ordering = NATURAL;
        else if (orderingStr == "Random")
          ordering = RANDOM;
        else if (orderingStr == "Graph")
          ordering = GRAPH;
        else TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory::BuildUncoupledAggregationFactory()::Unknown Ordering type");

        factory->SetOrdering(ordering);
      }

      if(paramList.isParameter("MaxNeighAlreadySelected")) {
        factory->SetMaxNeighAlreadySelected(paramList.get<int>("MaxNeighAlreadySelected"));
      }

      if(paramList.isParameter("MinNodesPerAggregate")) {
        factory->SetMinNodesPerAggregate(paramList.get<int>("MinNodesPerAggregate"));
      }

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
    RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
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

    RCP<FactoryBase> BuildDirectSolver(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new DirectSolver())));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "DirectSolver", Exceptions::RuntimeError, "");

      std::string type;              if(paramList.isParameter("type"))          type = paramList.get<std::string>("type");
      // std::string verbose;        if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new DirectSolver(type, params))));
    }

    RCP<FactoryBase> BuildMultiVectorTransferFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "MultiVectorTransferFactory", Exceptions::RuntimeError, "");

      std::string vectorName="";      vectorName = paramList.get<std::string>("vectorName");
      std::string restrictionName=""; restrictionName = paramList.get<std::string>("restrictionName");
      RCP<Factory> factory = rcp(new MultiVectorTransferFactory(vectorName, restrictionName));
      MUELU_FACTORY_PARAM2(restrictionName);

      return factory;
    }

#ifdef HAVE_MPI
    RCP<FactoryBase> BuildRepartitionFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "RepartitionFactory", Exceptions::RuntimeError, "");

      int minRowsPerProc=2000;     if(paramList.isParameter("minRowsPerProc"))   minRowsPerProc = paramList.get<int>("minRowsPerProc");
      double nonzeroImbalance=1.2; if(paramList.isParameter("nonzeroImbalance")) nonzeroImbalance = paramList.get<double>("nonzeroImbalance");
      int startLevel=1;            if(paramList.isParameter("startLevel"))       startLevel = paramList.get<int>("startLevel");
      int diffusive=0;             if(paramList.isParameter("diffusive"))        diffusive = paramList.get<int>("diffusive");

      RCP<Factory> factory = rcp(new RepartitionFactory(minRowsPerProc, nonzeroImbalance, startLevel, diffusive));

      MUELU_FACTORY_PARAM2("A");
      MUELU_FACTORY_PARAM2("Partition");

      return factory;
    }
#endif

    RCP<FactoryBase> BuildRebalanceTransferFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "RebalanceTransferFactory", Exceptions::RuntimeError, "");

      std::string type; type = paramList.get<std::string>("type");
      if (type == "Interpolation") {
        RCP<Factory> factory = rcp(new RebalanceTransferFactory(MueLu::INTERPOLATION));
        MUELU_FACTORY_PARAM2("Importer");
        MUELU_FACTORY_PARAM2("A");
        MUELU_FACTORY_PARAM2("P");
        return factory;
      } else if (type == "Restriction") {
        RCP<Factory> factory = rcp(new RebalanceTransferFactory(MueLu::RESTRICTION));
        MUELU_FACTORY_PARAM2("Importer");
        MUELU_FACTORY_PARAM2("A");
        MUELU_FACTORY_PARAM2("R");
        MUELU_FACTORY_PARAM2("Nullspace");
        return factory;
      } else {
        TEUCHOS_TEST_FOR_EXCEPT(1);
      }
    }
  }; // class

} // namespace MueLu

#define MUELU_FACTORYFACTORY_SHORT
#endif // MUELU_FACTORYFACTORY_DECL_HPP

// TODO: handle factory parameters
// TODO: parameter validator
// TODO: static
// TODO: default parameters should not be duplicated here and on the Factory (ex: default for overlap (=0) is defined both here and on TrilinosSmoother constructors)

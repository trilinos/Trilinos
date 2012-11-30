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
#include "MueLu_RepartitionAcFactory.hpp" //TMP
#include "MueLu_TransPFactory.hpp" //TMP
#include "MueLu_GenericRFactory.hpp" //TMP
#include "MueLu_SaPFactory.hpp" //TMP
#include "MueLu_PgPFactory.hpp" //TMP
#include "MueLu_TrilinosSmoother.hpp" //TMP
#include "MueLu_SmootherFactory.hpp" //TMP
#include "MueLu_TentativePFactory.hpp" //TMP
#include "MueLu_UCAggregationFactory.hpp" //TMP
#include "MueLu_DirectSolver.hpp" //TMP
#include "MueLu_Exceptions.hpp" //TMP
#include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_PermutedTransferFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_RepartitionFactory.hpp"

namespace MueLu {

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

      // TODO: see how Teko handles this.
      if (factoryName == "AmalgamationFactory") {
        return BuildAmalgamationFactory(paramList, factoryMapIn);
      }
      if (factoryName == "CoalesceDropFactory") {
        return BuildCoalesceDropFactory(paramList, factoryMapIn);
      }
      if (factoryName == "TentativePFactory") {
        return BuildTentativePFactory(paramList, factoryMapIn);
      }
      if (factoryName == "SaPFactory") {
        return BuildSaPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "TransPFactory") {
        return BuildTransPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "RAPFactory") {
        return BuildRAPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "RepartitionAcFactory") {
        return BuildRepartitionAcFactory(paramList, factoryMapIn);
      }
      if (factoryName == "UCAggregationFactory") {
        return BuildUCAggregationFactory(paramList, factoryMapIn);
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
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
      if (factoryName == "ZoltanInterface") {
        return BuildZoltanInterface(paramList, factoryMapIn);
      }
#endif
#ifdef HAVE_MPI
      if (factoryName == "RepartitionFactory") {
        return BuildRepartitionFactory(paramList, factoryMapIn);
      }
#endif
      if (factoryName == "PermutedTransferFactory") {
        return BuildPermutedTransferFactory(paramList, factoryMapIn);
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

#define MUELU_FACTORY_PARAM2(name)                                      \
    if (paramList.isParameter(name)) { factory->SetFactory(name, BuildFactory(paramList.getEntry(name), factoryMapIn)); }

    //! AmalgamationFactory
    RCP<FactoryBase> BuildAmalgamationFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      RCP<Factory> factory = rcp(new AmalgamationFactory());
      MUELU_FACTORY_PARAM2("A");
      return factory;
    }

    //! CoalesceDropFactory
    RCP<FactoryBase> BuildCoalesceDropFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      RCP<Factory> factory = rcp(new CoalesceDropFactory());
      MUELU_FACTORY_PARAM2("A");
      MUELU_FACTORY_PARAM2("UnAmalgamationInfo");
      return factory;
    }

    //! TentativePFactory
    RCP<FactoryBase> BuildTentativePFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      RCP<Factory> factory = rcp(new TentativePFactory());
      MUELU_FACTORY_PARAM2("Aggregates");
      MUELU_FACTORY_PARAM2("Nullspace");
      MUELU_FACTORY_PARAM2("A");
      MUELU_FACTORY_PARAM2("UnAmalgamationInfo");
      return factory;
    }

    //! SaPFactory
    RCP<FactoryBase> BuildSaPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new SaPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "SaPFactory", Exceptions::RuntimeError, "");
      RCP<SaPFactory> factory = rcp(new SaPFactory());
      MUELU_FACTORY_PARAM2("A");
      MUELU_FACTORY_PARAM2("P");

      if (paramList.isParameter("DampingFactor")) factory->SetDampingFactor(paramList.get<Scalar>("DampingFactor"));

      return factory;
    }

    //! PgPFactory
    RCP<FactoryBase> BuildPgPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new PgPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "PgPFactory", Exceptions::RuntimeError, "");
      RCP<PgPFactory> factory = rcp(new PgPFactory());
      MUELU_FACTORY_PARAM2("A");
      MUELU_FACTORY_PARAM2("P");

      return factory;
    }

    //! TransPFactory
    RCP<FactoryBase> BuildTransPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "TransPFactory", Exceptions::RuntimeError, "");
      RCP<Factory> factory = rcp(new TransPFactory());
      MUELU_FACTORY_PARAM2("P");
      return factory;
    }

    //! GenericRFactory
    RCP<FactoryBase> BuildGenericRFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "GenericRFactory", Exceptions::RuntimeError, "");
      RCP<Factory> factory = rcp(new GenericRFactory());
      MUELU_FACTORY_PARAM2("P");

      return factory;
    }

    //! RAPFactory
    RCP<FactoryBase> BuildRAPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new RAPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "RAPFactory", Exceptions::RuntimeError, "");

      RCP<RAPFactory> factory = rcp(new RAPFactory());
      MUELU_FACTORY_PARAM2("P");
      MUELU_FACTORY_PARAM2("R");
      MUELU_FACTORY_PARAM2("A");

      if (paramList.isSublist("TransferFactories")) {
        Teuchos::ParameterList transferList = paramList.sublist("TransferFactories");
        for (Teuchos::ParameterList::ConstIterator param = transferList.begin(); param != transferList.end(); ++param) {
          RCP<const FactoryBase> p = BuildFactory(transferList.entry(param), factoryMapIn);
          factory->AddTransferFactory(p);
        }
      }

      return factory;
    }

    //! RepartitionAcFactory
    RCP<FactoryBase> BuildRepartitionAcFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new RepartitionAcFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "RepartitionAcFactory", Exceptions::RuntimeError, "");

      RCP<RepartitionAcFactory> factory = rcp(new RepartitionAcFactory());
      MUELU_FACTORY_PARAM2("P"); //TODO
      MUELU_FACTORY_PARAM2("R");
      MUELU_FACTORY_PARAM2("A");

      return factory;
    }

    //! UCAggregationFactory
    RCP<FactoryBase> BuildUCAggregationFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new UCAggregationFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "UCAggregationFactory", Exceptions::RuntimeError, "");

      RCP<UCAggregationFactory> factory = rcp(new UCAggregationFactory());
      MUELU_FACTORY_PARAM2("Graph");

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

#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
    RCP<FactoryBase> BuildZoltanInterface(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "ZoltanInterface", Exceptions::RuntimeError, "");

      RCP<Factory> factory = rcp(new ZoltanInterface());
      MUELU_FACTORY_PARAM2("A");
      MUELU_FACTORY_PARAM2("TransferFactory");
      return factory;
    }
#endif

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

    RCP<FactoryBase> BuildPermutedTransferFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "PermutedTransferFactory", Exceptions::RuntimeError, "");

      std::string type; type = paramList.get<std::string>("type");
      if (type == "Interpolation") {
        RCP<Factory> factory = rcp(new PermutedTransferFactory(MueLu::INTERPOLATION));
        MUELU_FACTORY_PARAM2("Importer");
        MUELU_FACTORY_PARAM2("A");
        MUELU_FACTORY_PARAM2("P");
        return factory;
      } else if (type == "Restriction") {
        RCP<Factory> factory = rcp(new PermutedTransferFactory(MueLu::RESTRICTION));
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

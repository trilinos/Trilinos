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
#ifndef MUELU_FACTORYMANAGER_DEF_HPP
#define MUELU_FACTORYMANAGER_DEF_HPP

#include <Teuchos_ParameterList.hpp>

// Headers for factories used by default:
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_AggregateQualityEstimateFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_InitialBlockNumberFactory.hpp"
#include "MueLu_LineDetectionFactory.hpp"
// #include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_NoFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_RepartitionHeuristicFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_ScaledNullspaceFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_StructuredAggregationFactory.hpp"
#include "MueLu_HybridAggregationFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_InterfaceMappingTransferFactory.hpp"
#include "MueLu_InterfaceAggregationFactory.hpp"
#include "MueLu_InverseApproximationFactory.hpp"

#include "MueLu_AmalgamationFactory_kokkos.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_CoarseMapFactory_kokkos.hpp"
#include "MueLu_CoordinatesTransferFactory_kokkos.hpp"
#include "MueLu_NullspaceFactory_kokkos.hpp"
#include "MueLu_SaPFactory_kokkos.hpp"
#include "MueLu_TentativePFactory_kokkos.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"

#include "MueLu_FactoryManager_decl.hpp"


namespace MueLu {

#define MUELU_KOKKOS_FACTORY(varName, oldFactory, newFactory)   \
  (!useKokkos_) ? SetAndReturnDefaultFactory(varName, rcp(new oldFactory())) : \
                  SetAndReturnDefaultFactory(varName, rcp(new newFactory()));

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory) {
    factoryTable_[varName] = factory;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const FactoryBase> FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetFactory(const std::string& varName) const {
    if (factoryTable_.count(varName)) {
      // Search user provided factories
      return factoryTable_.find(varName)->second;
    }

    // Search/create default factory for this name
    return GetDefaultFactory(varName);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<FactoryBase> FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetFactoryNonConst(const std::string& varName) {
    return Teuchos::rcp_const_cast<FactoryBase>(GetFactory(varName));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasFactory(const std::string& varName) const {
    if (factoryTable_.count(varName)) return true;
    return false;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const FactoryBase> FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetDefaultFactory(const std::string& varName) const {
    if (defaultFactoryTable_.count(varName)) {
      // The factory for this name was already created (possibly, for previous level, if we reuse factory manager)
      return defaultFactoryTable_.find(varName)->second;

    } else {
      // No factory was created for this name, but we may know which one to create
      if (varName == "A")                               return SetAndReturnDefaultFactory(varName, rcp(new RAPFactory()));
      if (varName == "Ainv")                            return SetAndReturnDefaultFactory(varName, rcp(new InverseApproximationFactory()));
      if (varName == "RAP Pattern")                     return GetFactory("A");
      if (varName == "AP Pattern")                      return GetFactory("A");
      if (varName == "Ptent")                           return MUELU_KOKKOS_FACTORY(varName, TentativePFactory, TentativePFactory_kokkos);
      if (varName == "P") {
        // GetFactory("Ptent"): we need to use the same factory instance for both "P" and "Nullspace"
        RCP<Factory> factory;
        if (useKokkos_)
          factory = rcp(new SaPFactory_kokkos());
        else
          factory = rcp(new SaPFactory());
        factory->SetFactory("P", GetFactory("Ptent"));
        return SetAndReturnDefaultFactory(varName, factory);
      }
      if (varName == "Nullspace") {
        // GetFactory("Ptent"): we need to use the same factory instance for both "P" and "Nullspace"
        RCP<Factory> factory;
        if (useKokkos_)
          factory = rcp(new NullspaceFactory_kokkos());
        else
          factory = rcp(new NullspaceFactory());
        factory->SetFactory("Nullspace", GetFactory("Ptent"));
        return SetAndReturnDefaultFactory(varName, factory);
      }
      if (varName == "Scaled Nullspace")                return SetAndReturnDefaultFactory(varName, rcp(new ScaledNullspaceFactory()));

      if (varName == "Coordinates")                     return GetFactory("Ptent");
      if (varName == "Node Comm")                       return GetFactory("Ptent");

      if (varName == "R")                               return SetAndReturnDefaultFactory(varName, rcp(new TransPFactory()));
      if (varName == "RfromPfactory")                   return  GetFactory("P");
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
      if (varName == "Partition")                       return SetAndReturnDefaultFactory(varName, rcp(new ZoltanInterface()));
#endif //ifdef HAVE_MPI

      if (varName == "Importer") {
#ifdef HAVE_MPI
                                                        return SetAndReturnDefaultFactory(varName, rcp(new RepartitionFactory()));
#else
                                                        return SetAndReturnDefaultFactory(varName, NoFactory::getRCP());
#endif
      }
      if (varName == "number of partitions") {
#ifdef HAVE_MPI
                                                        return SetAndReturnDefaultFactory(varName, rcp(new RepartitionHeuristicFactory()));
#else
                                                        return SetAndReturnDefaultFactory(varName, NoFactory::getRCP());
#endif
      }
      if (varName == "repartition: heuristic target rows per process") return GetFactory("number of partitions");

      if (varName == "Graph")                           return MUELU_KOKKOS_FACTORY(varName, CoalesceDropFactory, CoalesceDropFactory_kokkos);
      if (varName == "UnAmalgamationInfo")              return MUELU_KOKKOS_FACTORY(varName, AmalgamationFactory, AmalgamationFactory_kokkos);
      if (varName == "Aggregates")                      return MUELU_KOKKOS_FACTORY(varName, UncoupledAggregationFactory, UncoupledAggregationFactory_kokkos);
      if (varName == "AggregateQualities")              return SetAndReturnDefaultFactory(varName, rcp(new AggregateQualityEstimateFactory()));
      if (varName == "CoarseMap")                       return MUELU_KOKKOS_FACTORY(varName, CoarseMapFactory, CoarseMapFactory_kokkos);
      if (varName == "DofsPerNode")                     return GetFactory("Graph");
      if (varName == "Filtering")                       return GetFactory("Graph");
      if (varName == "BlockNumber")                     return SetAndReturnDefaultFactory(varName, rcp(new InitialBlockNumberFactory()));
      if (varName == "LineDetection_VertLineIds")       return SetAndReturnDefaultFactory(varName, rcp(new LineDetectionFactory()));
      if (varName == "LineDetection_Layers")            return GetFactory("LineDetection_VertLineIds");
      if (varName == "CoarseNumZLayers")                return GetFactory("LineDetection_VertLineIds");
      
      // Structured
      if (varName == "structuredInterpolationOrder")    return SetAndReturnDefaultFactory(varName, rcp(new StructuredAggregationFactory()));

      // Non-Galerkin
      if (varName == "K")                               return GetFactory("A");
      if (varName == "M")                               return GetFactory("A");
      if (varName == "Mdiag")                           return GetFactory("A");
      if (varName == "cfl-based shift array")           return GetFactory("A");

      // Same factory for both Pre and Post Smoother. Factory for key "Smoother" can be set by users.
      if (varName == "PreSmoother")                     return GetFactory("Smoother");
      if (varName == "PostSmoother")                    return GetFactory("Smoother");

      if (varName == "Ppattern") {
        RCP<PatternFactory> PpFact = rcp(new PatternFactory);
        PpFact->SetFactory("P", GetFactory("Ptent"));
        return SetAndReturnDefaultFactory(varName, PpFact);
      }
      if (varName == "Constraint")                      return SetAndReturnDefaultFactory(varName, rcp(new ConstraintFactory()));

      if (varName == "Smoother") {
        Teuchos::ParameterList smootherParamList;
        smootherParamList.set("relaxation: type",           "Symmetric Gauss-Seidel");
        smootherParamList.set("relaxation: sweeps",         Teuchos::OrdinalTraits<LO>::one());
        smootherParamList.set("relaxation: damping factor", Teuchos::ScalarTraits<Scalar>::one());
        return SetAndReturnDefaultFactory(varName, rcp(new SmootherFactory(rcp(new TrilinosSmoother("RELAXATION", smootherParamList)))));
      }
      if (varName == "CoarseSolver")                    return SetAndReturnDefaultFactory(varName, rcp(new SmootherFactory(rcp(new DirectSolver()), Teuchos::null)));

      if (varName == "DualNodeID2PrimalNodeID")         return SetAndReturnDefaultFactory(varName, rcp(new InterfaceMappingTransferFactory()));
      if (varName == "CoarseDualNodeID2PrimalNodeID")   return SetAndReturnDefaultFactory(varName, rcp(new InterfaceAggregationFactory()));
#ifdef HAVE_MUELU_INTREPID2
      // If we're asking for it, find who made P
      if (varName == "pcoarsen: element to node map")   return GetFactory("P");
#endif

      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::FactoryManager::GetDefaultFactory(): No default factory available for building '" + varName + "'.");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const FactoryBase> FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetAndReturnDefaultFactory(const std::string& varName, const RCP<const FactoryBase>& factory) const {
    TEUCHOS_TEST_FOR_EXCEPTION(factory.is_null(), Exceptions::RuntimeError, "The default factory for building '" << varName << "' is null");

    GetOStream(Runtime1) << "Using default factory (" << factory->ShortClassName() <<"["<<factory->GetID()<<"]) for building '" << varName << "'." << std::endl;

    defaultFactoryTable_[varName] = factory;

    return defaultFactoryTable_[varName];
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Print() const {
    std::map<std::string, RCP<const FactoryBase> >::const_iterator it;
    Teuchos::FancyOStream& fancy = GetOStream(Debug);
    //auto & fancy = std::cout;// For debugging


    fancy << "Users factory table (factoryTable_):" << std::endl;
    for (it = factoryTable_.begin(); it != factoryTable_.end(); it++) {
      fancy << "  " << it->first << " -> ";
      if (it->second.get() == NoFactory::get()) fancy << "NoFactory";
      else if (!it->second.get()) fancy<< "NULL";
      else {
        fancy << it->second.get()->ShortClassName()<<"["<<it->second.get()->GetID()<<"]";
#ifdef HAVE_MUELU_DEBUG
        fancy<<"("<<Teuchos::toString(it->second.get()) <<")";
#endif
      }
      fancy<< std::endl;
    }

    fancy << "Default factory table (defaultFactoryTable_):" << std::endl;
    for (it = defaultFactoryTable_.begin(); it != defaultFactoryTable_.end(); it++) {
      fancy << "  " << it->first << " -> ";
      if (it->second.get() == NoFactory::get()) fancy << "NoFactory";
      else if (!it->second.get()) fancy<< "NULL";
      else {
        fancy << it->second.get()->ShortClassName()<<"["<<it->second.get()->GetID()<<"]";
#ifdef HAVE_MUELU_DEBUG
        fancy<<"("<<Teuchos::toString(it->second.get()) <<")";
#endif
      }
      fancy<< std::endl;
    }

  }

#ifdef HAVE_MUELU_DEBUG
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ResetDebugData() const {
    std::map<std::string, RCP<const FactoryBase> >::const_iterator it;

    for (it = factoryTable_.begin(); it != factoryTable_.end(); it++)
      if (!it->second.is_null())
        it->second->ResetDebugData();

    for (it = defaultFactoryTable_.begin(); it != defaultFactoryTable_.end(); it++)
      if (!it->second.is_null())
        it->second->ResetDebugData();
  }
#endif


#undef MUELU_KOKKOS_FACTORY

} // namespace MueLu

//TODO: add operator[]
//TODO: should we use a parameterList instead of a std::map? It might be useful to tag which factory have been used and report unused factory.
//TODO: add an option 'NoDefault' to check if we are using any default factory.
//TODO: use Teuchos::ConstNonConstObjectContainer to allow user to modify factories after a GetFactory()

#endif // MUELU_FACTORYMANAGER_DEF_HPP

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
// Type definitions for templated classes (generally graph-related) that do not require a scalar.

#include <Xpetra_UseShortNamesOrdinal.hpp>

#ifdef MUELU_GRAPH_SHORT
typedef MueLu::Graph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Graph;
#endif

#ifdef MUELU_AGGREGATES_SHORT
typedef MueLu::Aggregates<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> Aggregates;
#endif

#ifdef MUELU_LOCALAGGREGATIONALGORITHM_SHORT
typedef MueLu::LocalAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LocalAggregationAlgorithm;
#endif

#ifdef MUELU_LEFTOVERAGGREGATIONALGORITHM_SHORT
typedef MueLu::LeftoverAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> LeftoverAggregationAlgorithm;
#endif

#ifdef MUELU_UCAGGREGATIONFACTORY_SHORT
typedef MueLu::UCAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UCAggregationFactory;
#endif

#ifdef MUELU_UCAGGREGATIONCOMMHELPER_SHORT
typedef MueLu::UCAggregationCommHelper<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> UCAggregationCommHelper;
#endif

#ifdef MUELU_PRFACTORY_SHORT
typedef MueLu::PRFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> PRFactory;
#endif

#ifdef MUELU_ZOLTANINTERFACE_SHORT
typedef MueLu::ZoltanInterface<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ZoltanInterface;
#endif

#ifdef MUELU_AMALGAMATIONINFO_SHORT
typedef MueLu::AmalgamationInfo<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> AmalgamationInfo;
#endif

#ifdef MUELU_CHEAPAGGREGATIONALGORITHM_SHORT
typedef MueLu::CheapAggregationAlgorithm<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> CheapAggregationAlgorithm;
#endif

#ifdef MUELU_EXPERIMENTALAGGREGATIONFACTORY_SHORT
typedef MueLu::ExperimentalAggregationFactory<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> ExperimentalAggregationFactory;
#endif

#ifdef MUELU_SINGLELEVELFACTORYBASE_SHORT
typedef MueLu::SingleLevelFactoryBase SingleLevelFactoryBase;
#endif

#ifdef MUELU_TWOLEVELFACTORYBASE_SHORT
typedef MueLu::TwoLevelFactoryBase TwoLevelFactoryBase;
#endif

#ifdef MUELU_RFACTORY_SHORT
typedef MueLu::RFactory RFactory;
#endif

#ifdef MUELU_AMESOSSMOOTHER_SHORT
typedef MueLu::AmesosSmoother AmesosSmoother;
#endif

#ifdef MUELU_IFPACKSMOOTHER_SHORT
typedef MueLu::IfpackSmoother IfpackSmoother;
#endif

#ifdef MUELU_FACTORYBASE_SHORT
typedef MueLu::FactoryBase FactoryBase;
#endif

#ifdef MUELU_FACTORYMANAGERBASE_SHORT
typedef MueLu::FactoryManagerBase FactoryManagerBase;
#endif

#ifdef MUELU_LEVEL_SHORT
typedef MueLu::Level Level;
#endif

#ifdef MUELU_PFACTORY_SHORT
typedef MueLu::PFactory PFactory;
#endif

#ifdef MUELU_TWOKEYMAP_SHORT
typedef MueLu::TwoKeyMap TwoKeyMap;
#endif

#ifdef MUELU_VARIABLECONTAINER_SHORT
typedef MueLu::VariableContainer VariableContainer;
#endif


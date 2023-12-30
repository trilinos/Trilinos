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
#ifndef MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_DECL_HPP
#define MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "MueLu_AggregationPhase2bAlgorithm_kokkos_fwd.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_AggregationAlgorithmBase_kokkos.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

namespace MueLu {
/*!
  @class AggregationPhase2bAlgorithm class.
  @brief Add leftovers to existing aggregates
  @ingroup Aggregation

  ### Idea ###
  In phase 2b non-aggregated nodes are added to existing aggregates.
  All neighbors of the unaggregated node are checked and the corresponding
  aggregate weight is increased. The unaggregated node is added to the aggregate
  with the best weight. A simple penalty strategy makes sure that the non-aggregated
  nodes are added to different aggregates.
  The routine runs twice to cover non-aggregate nodes which have a node distance
  of two to existing aggregates. Assuming that the node distance is not greater
  than 3 (the aggregate diameter size), running the algorithm only twice should
  be sufficient.

  ### Comments ###
  Only nodes with state READY are changed to AGGREGATED. There are no aggregation criteria considered. Especially the aggregation: max agg size criterion is ignored.
  This is not a problem, since after the previous aggregation phases one should not be able to build too large aggregates.
*/

template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationPhase2bAlgorithm_kokkos : public MueLu::AggregationAlgorithmBase_kokkos<LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

 public:
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;
  using memory_space    = typename LWGraph_kokkos::memory_space;

  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  AggregationPhase2bAlgorithm_kokkos(const RCP<const FactoryBase>& /* graphFact */ = Teuchos::null) {}

  //! Destructor.
  virtual ~AggregationPhase2bAlgorithm_kokkos() {}

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */

  void BuildAggregates(const ParameterList& params,
                       const LWGraph_kokkos& graph,
                       Aggregates& aggregates,
                       Kokkos::View<unsigned*, device_type>& aggStat,
                       LO& numNonAggregatedNodes) const;

  void BuildAggregatesRandom(const ParameterList& params,
                             const LWGraph_kokkos& graph,
                             Aggregates& aggregates,
                             Kokkos::View<unsigned*, device_type>& aggStat,
                             LO& numNonAggregatedNodes) const;

  void BuildAggregatesDeterministic(const ParameterList& params,
                                    const LWGraph_kokkos& graph,
                                    Aggregates& aggregates,
                                    Kokkos::View<unsigned*, device_type>& aggStat,
                                    LO& numNonAggregatedNodes) const;
  //@}

  std::string description() const { return "Phase 2b (expansion)"; }
};

}  // namespace MueLu

#define MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_SHORT
#endif  // MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_DECL_HPP

// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONALGORITHMBASE_HPP_
#define MUELU_AGGREGATIONALGORITHMBASE_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Types.hpp"

namespace MueLu {

/*!
     @class AggregationAlgorithmBase
     @brief Pure virtual base class for all MueLu aggregation algorithms

     @ingroup MueLuBaseClasses
 */
template <class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class AggregationAlgorithmBase : public BaseClass {
#undef MUELU_AGGREGATIONALGORITHMBASE_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"
 public:
  using LWGraphHostType = LWGraph;
  using AggStatHostType = Kokkos::View<unsigned*, typename LWGraphHostType::device_type>;

  using LWGraphType = LWGraph_kokkos;
  using AggStatType = Kokkos::View<unsigned*, typename LWGraphType::device_type>;

  //! @name Constructors/Destructors
  //@{

  //! Destructor.
  virtual ~AggregationAlgorithmBase() {}

  //@}

  //! @name Build routines
  //@{

  //! BuildAggregatesNonKokkos routine.
  virtual void BuildAggregatesNonKokkos(const Teuchos::ParameterList& params,
                                        const LWGraphHostType& graph,
                                        Aggregates& aggregates,
                                        AggStatHostType& aggStat,
                                        LO& numNonAggregatedNodes) const = 0;

  //! BuildAggregates routine.
  virtual void BuildAggregates(const Teuchos::ParameterList& params,
                               const LWGraphType& graph,
                               Aggregates& aggregates,
                               AggStatType& aggStat,
                               LO& numNonAggregatedNodes) const = 0;
  //@}
};

}  // namespace MueLu

#define MUELU_AGGREGATIONALGORITHMBASE_SHORT
#endif /* MUELU_AGGREGATIONALGORITHMBASE_HPP_ */

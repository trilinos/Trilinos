// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_TpetraCrsGraphAdapter.hpp
    \brief Defines TpetraCrsGraphAdapter class.
*/

#ifndef _ZOLTAN2_TPETRACRSGRAPHADAPTER_HPP_
#define _ZOLTAN2_TPETRACRSGRAPHADAPTER_HPP_

#include <Zoltan2_PartitioningHelpers.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_TpetraRowGraphAdapter.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <string>

namespace Zoltan2 {

/*!  \brief Provides access for Zoltan2 to Tpetra::CrsGraph data.

    \todo test for memory alloc failure when we resize a vector
    \todo we assume FillComplete has been called.  Should we support
                objects that are not FillCompleted.
*/

template <typename User, typename UserCoord = User>
class TpetraCrsGraphAdapter : public TpetraRowGraphAdapter<User, UserCoord> {

public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using scalar_t = typename InputTraits<User>::scalar_t;
  using offset_t = typename InputTraits<User>::offset_t;
  using lno_t = typename InputTraits<User>::lno_t;
  using gno_t = typename InputTraits<User>::gno_t;
  using part_t = typename InputTraits<User>::part_t;
  using node_t = typename InputTraits<User>::node_t;
  using user_t = User;
  using userCoord_t = UserCoord;

  using Base = GraphAdapter<User, UserCoord>;
  using RowGraph = TpetraRowGraphAdapter<User, UserCoord>;
#endif

  /*! \brief Constructor for graph with no weights or coordinates.
   *  \param ingraph the Tpetra::CrsGraph
   *  \param numVtxWeights  the number of weights per vertex (default = 0)
   *  \param numEdgeWeights the number of weights per edge  (default = 0)
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  TpetraCrsGraphAdapter(const RCP<const User> &graph, int nVtxWeights = 0,
                        int nEdgeWeights = 0);

  /*! \brief Access to user's graph
   */
  RCP<const User> getUserGraph() const { return this->graph_; }

  template <typename Adapter>
  void applyPartitioningSolution(
      const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const;

  template <typename Adapter>
  void applyPartitioningSolution(
      const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User, typename UserCoord>
TpetraCrsGraphAdapter<User, UserCoord>::TpetraCrsGraphAdapter(
    const RCP<const User> &graph, int nVtxWgts, int nEdgeWgts)
    : TpetraRowGraphAdapter<User>(nVtxWgts, nEdgeWgts, graph) {
  auto adjIdsHost = graph->getLocalIndicesHost();

  auto adjIdsGlobalHost =
      typename Base::IdsHostView("adjIdsGlobalHost", adjIdsHost.extent(0));
  auto colMap = graph->getColMap();

  // Convert to global IDs using Tpetra::Map
  Kokkos::parallel_for("adjIdsGlobalHost",
                       Kokkos::RangePolicy<Kokkos::HostSpace::execution_space>(
                           0, adjIdsGlobalHost.extent(0)),
                       [=](const int i) {
                         adjIdsGlobalHost(i) =
                             colMap->getGlobalElement(adjIdsHost(i));
                       });

  auto adjIdsDevice = Kokkos::create_mirror_view_and_copy(
      typename Base::device_t(), adjIdsGlobalHost);

  this->adjIdsDevice_ = adjIdsDevice;
  this->offsDevice_ = graph->getLocalRowPtrsDevice();

  if (this->nWeightsPerVertex_ > 0) {

    this->vertexWeightsDevice_ = typename Base::WeightsDeviceView(
        "vertexWeightsDevice_", graph->getLocalNumRows(),
        this->nWeightsPerVertex_);

    this->vertexDegreeWeightsHost_ = typename Base::VtxDegreeHostView(
        "vertexDegreeWeightsHost_", this->nWeightsPerVertex_);

    for (int i = 0; i < this->nWeightsPerVertex_; ++i) {
      this->vertexDegreeWeightsHost_(i) = false;
    }
  }

  if (this->nWeightsPerEdge_) {
    this->edgeWeightsDevice_ = typename Base::WeightsDeviceView(
        "nWeightsPerEdge_", graph->getLocalNumRows(), this->nWeightsPerEdge_);
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
template <typename Adapter>
void TpetraCrsGraphAdapter<User, UserCoord>::applyPartitioningSolution(
    const User &in, User *&out,
    const PartitioningSolution<Adapter> &solution) const {
  TpetraRowGraphAdapter<User, UserCoord>::applyPartitioningSolution(in, out,
                                                                    solution);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
template <typename Adapter>
void TpetraCrsGraphAdapter<User, UserCoord>::applyPartitioningSolution(
    const User &in, RCP<User> &out,
    const PartitioningSolution<Adapter> &solution) const {
  TpetraRowGraphAdapter<User, UserCoord>::applyPartitioningSolution(in, out,
                                                                    solution);
}

} // namespace Zoltan2

#endif

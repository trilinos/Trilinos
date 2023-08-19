// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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

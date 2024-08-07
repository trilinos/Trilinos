// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_TpetraRowGraphAdapter.hpp
    \brief Defines TpetraRowGraphAdapter class.
*/

#ifndef _ZOLTAN2_TPETRAROWGRAPHADAPTER_HPP_
#define _ZOLTAN2_TPETRAROWGRAPHADAPTER_HPP_

#include "Kokkos_DualView.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include <Tpetra_RowGraph.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>
#include <Zoltan2_StridedData.hpp>
#include <string>

namespace Zoltan2 {

/*!  \brief Provides access for Zoltan2 to Tpetra::RowGraph data.

    \todo test for memory alloc failure when we resize a vector
    \todo we assume FillComplete has been called.  Should we support
                objects that are not FillCompleted.

    The template parameter is the user's input object:
     \li Tpetra::CrsGraph
     \li Tpetra::RowGraph
     \li Epetra_CrsGraph

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.
*/

template <typename User, typename UserCoord = User>
class TpetraRowGraphAdapter : public GraphAdapter<User, UserCoord> {

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
#endif

  /*! \brief Constructor for graph with no weights or coordinates.
   *  \param ingraph the  Tpetra::RowGraph
   *  \param numVtxWeights  the number of weights per vertex (default = 0)
   *  \param numEdgeWeights the number of weights per edge  (default = 0)
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  TpetraRowGraphAdapter(const RCP<const User> &ingraph, int nVtxWeights = 0,
                        int nEdgeWeights = 0);

  /*! \brief Provide a pointer to weights for the primary entity type.
   *    \param val A pointer to the weights for index \c idx.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th entity for index \idx.
   *    \param idx A number from 0 to one less than
   *          weight idx specified in the constructor.
   *
   *  The order of the weights should match the order that
   *  entities appear in the input data structure.
   */

  void setWeights(const scalar_t *val, int stride, int idx);

  /*! \brief Provide a device view of weights for the primary entity type.
   *    \param val A view to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          weight idx specified in the constructor.
   *
   *  The order of the weights should match the order that
   *  entities appear in the input data structure.
   */
  void setWeightsDevice(typename Base::ConstWeightsDeviceView1D val, int idx);

  /*! \brief Provide a host view of weights for the primary entity type.
   *    \param val A view to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          weight idx specified in the constructor.
   *
   *  The order of the weights should match the order that
   *  entities appear in the input data structure.
   */
  void setWeightsHost(typename Base::ConstWeightsHostView1D val, int idx);

  /*! \brief Provide a pointer to vertex weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th vertex for index \idx.
   *    \param idx A number from 0 to one less than
   *          number of vertex weights specified in the constructor.
   *
   *  The order of the vertex weights should match the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getLocalElementList()
   *     \endcode
   */

  void setVertexWeights(const scalar_t *val, int stride, int idx);

  /*! \brief Provide a device view to vertex weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          number of vertex weights specified in the constructor.
   *
   *  The order of the vertex weights should match the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getLocalElementList()
   *     \endcode
   */
  void setVertexWeightsDevice(typename Base::ConstWeightsDeviceView1D val,
                              int idx);

  /*! \brief Provide a host view to vertex weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param idx A number from 0 to one less than
   *          number of vertex weights specified in the constructor.
   *
   *  The order of the vertex weights should match the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getLocalElementList()
   *     \endcode
   */
  void setVertexWeightsHost(typename Base::ConstWeightsHostView1D val, int idx);

  /*! \brief Specify an index for which the weight should be
              the degree of the entity
   *    \param idx Zoltan2 will use the entity's
   *         degree as the entity weight for index \c idx.
   */
  void setWeightIsDegree(int idx);

  /*! \brief Specify an index for which the vertex weight should be
              the degree of the vertex
   *    \param idx Zoltan2 will use the vertex's
   *         degree as the vertex weight for index \c idx.
   */
  void setVertexWeightIsDegree(int idx);

  /*! \brief Provide a pointer to edge weights.
   *    \param val A pointer to the weights for index \c idx.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th edge for index \idx.
   *    \param dim A number from 0 to one less than the number
   *          of edge weights specified in the constructor.
   *
   *  The order of the edge weights should follow the order that the
   *  the vertices and edges appear in the input data structure.
   *
   *  By vertex:
   *     \code
   *       TheGraph->getRowMap()->getLocalElementList()
   *     \endcode
   *
   *  Then by vertex neighbor:
   *     \code
   *       TheGraph->getLocalRowView(vertexNum, neighborList);
   *     \endcode
   */

  void setEdgeWeights(const scalar_t *val, int stride, int idx);

  /*! \brief Provide a device view to edge weights.
   *  \param val A pointer to the weights for index \c idx.
   *  \param idx A number from 0 to one less than the number
   *             of edge weights specified in the constructor.
   */
  void setEdgeWeightsDevice(typename Base::ConstWeightsDeviceView1D val,
                            int idx);

  /*! \brief Provide a host view to edge weights.
   *  \param val A pointer to the weights for index \c idx.
   *  \param idx A number from 0 to one less than the
   *             number of edge weights specified in the constructor.
   */
  void setEdgeWeightsHost(typename Base::ConstWeightsHostView1D val, int idx);

  ////////////////////////////////////////////////////
  // The GraphAdapter interface.
  ////////////////////////////////////////////////////

  // TODO:  Assuming rows == objects;
  // TODO:  Need to add option for columns or nonzeros?
  size_t getLocalNumVertices() const override;

  void getVertexIDsView(const gno_t *&ids) const override;

  void
  getVertexIDsDeviceView(typename Base::ConstIdsDeviceView &ids) const override;

  void
  getVertexIDsHostView(typename Base::ConstIdsHostView &ids) const override;

  size_t getLocalNumEdges() const override;

  void getEdgesView(const offset_t *&offsets,
                    const gno_t *&adjIds) const override;

  void
  getEdgesDeviceView(typename Base::ConstOffsetsDeviceView &offsets,
                     typename Base::ConstIdsDeviceView &adjIds) const override;

  void getEdgesHostView(typename Base::ConstOffsetsHostView &offsets,
                        typename Base::ConstIdsHostView &adjIds) const override;

  int getNumWeightsPerVertex() const override;

  void getVertexWeightsView(const scalar_t *&weights, int &stride,
                            int idx) const override;

  void getVertexWeightsDeviceView(typename Base::WeightsDeviceView1D &weights,
                                  int idx = 0) const override;

  void getVertexWeightsDeviceView(
      typename Base::WeightsDeviceView &weights) const override;

  void getVertexWeightsHostView(typename Base::WeightsHostView1D &weights,
                                int idx = 0) const override;

  void getVertexWeightsHostView(
      typename Base::WeightsHostView &weights) const override;

  bool useDegreeAsVertexWeight(int idx) const override;

  int getNumWeightsPerEdge() const override;

  void getEdgeWeightsView(const scalar_t *&weights, int &stride,
                          int idx) const override;

  void getEdgeWeightsDeviceView(typename Base::WeightsDeviceView1D &weights,
                                int idx = 0) const override;

  void getEdgeWeightsDeviceView(
      typename Base::WeightsDeviceView &weights) const override;

  void getEdgeWeightsHostView(typename Base::WeightsHostView1D &weights,
                              int idx = 0) const override;

  void getEdgeWeightsHostView(
      typename Base::WeightsHostView &weights) const override;

  template <typename Adapter>
  void applyPartitioningSolution(
      const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const;

  template <typename Adapter>
  void applyPartitioningSolution(
      const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const;

protected:
  // Useb by TpetraCrsGraphAdapter
  TpetraRowGraphAdapter(int nVtxWgts, int nEdgeWgts,
                        const RCP<const User> &graph)
      : graph_(graph), nWeightsPerVertex_(nVtxWgts),
        nWeightsPerEdge_(nEdgeWgts) {}

  RCP<const User> graph_;

  typename Base::ConstOffsetsHostView offsHost_;
  typename Base::ConstIdsHostView adjIdsHost_;

  typename Base::ConstIdsDeviceView adjIdsDevice_;
  typename Base::ConstOffsetsDeviceView offsDevice_;

  int nWeightsPerVertex_;
  ArrayRCP<StridedData<lno_t, scalar_t>> vertexWeights_;
  typename Base::WeightsDeviceView vertexWeightsDevice_;
  typename Base::VtxDegreeHostView vertexDegreeWeightsHost_;

  int nWeightsPerEdge_;
  ArrayRCP<StridedData<lno_t, scalar_t>> edgeWeights_;
  typename Base::WeightsDeviceView edgeWeightsDevice_;

  virtual RCP<User> doMigration(const User &from, size_t numLocalRows,
                                const gno_t *myNewRows) const;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User, typename UserCoord>
TpetraRowGraphAdapter<User, UserCoord>::TpetraRowGraphAdapter(
    const RCP<const User> &ingraph, int nVtxWgts, int nEdgeWgts)
    : graph_(ingraph), nWeightsPerVertex_(nVtxWgts),
      nWeightsPerEdge_(nEdgeWgts), edgeWeights_() {
  using strided_t = StridedData<lno_t, scalar_t>;
  using localInds_t = typename User::nonconst_local_inds_host_view_type;

  const auto nvtx = graph_->getLocalNumRows();
  const auto nedges = graph_->getLocalNumEntries();
  // Diff from CrsMatrix
  const auto maxNumEntries = graph_->getLocalMaxNumRowEntries();

  // Unfortunately we have to copy the offsets and edge Ids
  // because edge Ids are not usually stored in vertex id order.

  adjIdsHost_ = typename Base::ConstIdsHostView("adjIdsHost_", nedges);
  offsHost_ = typename Base::ConstOffsetsHostView("offsHost_", nvtx + 1);

  localInds_t nbors("nbors", maxNumEntries);

  for (size_t v = 0; v < nvtx; v++) {
    size_t numColInds = 0;
    graph_->getLocalRowCopy(v, nbors, numColInds); // Diff from CrsGraph

    offsHost_(v + 1) = offsHost_(v) + numColInds;
    for (offset_t e = offsHost_(v), i = 0; e < offsHost_(v + 1); e++) {
      adjIdsHost_(e) = graph_->getColMap()->getGlobalElement(nbors(i++));
    }
  }

  // Since there's no direct getter of offsets and edges in device view,
  // we have to deep copy here
  offsDevice_ =
      Kokkos::create_mirror_view_and_copy(typename Base::device_t(), offsHost_);
  adjIdsDevice_ = Kokkos::create_mirror_view_and_copy(typename Base::device_t(),
                                                      adjIdsHost_);

  if (nWeightsPerVertex_ > 0) {
    vertexWeights_ =
        arcp(new strided_t[nWeightsPerVertex_], 0, nWeightsPerVertex_, true);

    vertexWeightsDevice_ = typename Base::WeightsDeviceView(
        "vertexWeightsDevice_", nvtx, nWeightsPerVertex_);

    vertexDegreeWeightsHost_ = typename Base::VtxDegreeHostView(
        "vertexDegreeWeightsHost_", nWeightsPerVertex_);

    for (int i = 0; i < nWeightsPerVertex_; ++i) {
      vertexDegreeWeightsHost_(i) = false;
    }
  }

  if (nWeightsPerEdge_ > 0) {
    edgeWeights_ =
        arcp(new strided_t[nWeightsPerEdge_], 0, nWeightsPerEdge_, true);

    edgeWeightsDevice_ = typename Base::WeightsDeviceView(
        "nWeightsPerEdge_", graph_->getLocalNumRows(), nWeightsPerEdge_);
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setWeights(
    const scalar_t *weightVal, int stride, int idx) {
  if (this->getPrimaryEntityType() == GRAPH_VERTEX)
    setVertexWeights(weightVal, stride, idx);
  else
    setEdgeWeights(weightVal, stride, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setWeightsDevice(
    typename Base::ConstWeightsDeviceView1D val, int idx) {
  if (this->getPrimaryEntityType() == GRAPH_VERTEX)
    setVertexWeightsDevice(val, idx);
  else
    setEdgeWeightsDevice(val, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setWeightsHost(
    typename Base::ConstWeightsHostView1D val, int idx) {
  if (this->getPrimaryEntityType() == GRAPH_VERTEX)
    setVertexWeightsHost(val, idx);
  else
    setEdgeWeightsHost(val, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setVertexWeights(
    const scalar_t *weightVal, int stride, int idx) {
  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid vertex weight index: " + std::to_string(idx));

  size_t nvtx = getLocalNumVertices();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nvtx * stride, false);
  vertexWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setVertexWeightsDevice(
    typename Base::ConstWeightsDeviceView1D weights, int idx) {

  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid vertex weight index: " + std::to_string(idx));

  AssertCondition(vertexWeightsDevice_.extent(0) == weights.extent(0),
                  "Invalid sizes!");

  Kokkos::parallel_for(
      vertexWeightsDevice_.extent(0), KOKKOS_CLASS_LAMBDA(const int vertexID) {
        vertexWeightsDevice_(vertexID, idx) = weights(vertexID);
      });

  Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setVertexWeightsHost(
    typename Base::ConstWeightsHostView1D weightsHost, int idx) {
  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid vertex weight index: " + std::to_string(idx));

  auto weightsDevice = Kokkos::create_mirror_view_and_copy(
      typename Base::device_t(), weightsHost);

  setVertexWeightsDevice(weightsDevice, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setWeightIsDegree(int idx) {
  AssertCondition(this->getPrimaryEntityType() == GRAPH_VERTEX,
                  "setWeightIsNumberOfNonZeros is supported only for vertices");

  setVertexWeightIsDegree(idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setVertexWeightIsDegree(int idx) {
  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid vertex weight index.");

  vertexDegreeWeightsHost_(idx) = true;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setEdgeWeights(
    const scalar_t *weightVal, int stride, int idx) {
  typedef StridedData<lno_t, scalar_t> input_t;

  AssertCondition((idx >= 0) and (idx < nWeightsPerEdge_),
                  "Invalid edge weight index" + std::to_string(idx));

  size_t nedges = getLocalNumEdges();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nedges * stride, false);
  edgeWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setEdgeWeightsDevice(
    typename Base::ConstWeightsDeviceView1D weights, int idx) {
  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid edge weight index.");

  AssertCondition(edgeWeightsDevice_.extent(0) == weights.extent(0),
                  "Invalid sizes!");

  Kokkos::parallel_for(
      edgeWeightsDevice_.extent(0), KOKKOS_CLASS_LAMBDA(const int vertexID) {
        edgeWeightsDevice_(vertexID, idx) = weights(vertexID);
      });

  Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::setEdgeWeightsHost(
    typename Base::ConstWeightsHostView1D weightsHost, int idx) {
  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid edge weight index.");

  auto weightsDevice = Kokkos::create_mirror_view_and_copy(
      typename Base::device_t(), weightsHost);

  setEdgeWeightsDevice(weightsDevice);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
size_t TpetraRowGraphAdapter<User, UserCoord>::getLocalNumVertices() const {
  return graph_->getLocalNumRows();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexIDsView(
    const gno_t *&ids) const {
  ids = NULL;
  if (getLocalNumVertices())
    ids = graph_->getRowMap()->getLocalElementList().getRawPtr();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexIDsDeviceView(
    typename Base::ConstIdsDeviceView &ids) const {

  // TODO: Making a  ConstIdsDeviceView LayoutLeft would proably remove the
  //       need of creating tmpIds
  auto idsDevice = graph_->getRowMap()->getMyGlobalIndices();
  auto tmpIds = typename Base::IdsDeviceView("", idsDevice.extent(0));

  Kokkos::deep_copy(tmpIds, idsDevice);

  ids = tmpIds;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexIDsHostView(
    typename Base::ConstIdsHostView &ids) const {
  // TODO: Making a  ConstIdsDeviceView LayoutLeft would proably remove the
  //       need of creating tmpIds
  auto idsDevice = graph_->getRowMap()->getMyGlobalIndices();
  auto tmpIds = typename Base::IdsHostView("", idsDevice.extent(0));

  Kokkos::deep_copy(tmpIds, idsDevice);

  ids = tmpIds;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
size_t TpetraRowGraphAdapter<User, UserCoord>::getLocalNumEdges() const {
  return graph_->getLocalNumEntries();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgesView(
    const offset_t *&offsets, const gno_t *&adjIds) const {
  offsets = offsHost_.data();
  adjIds = (getLocalNumEdges() ? adjIdsHost_.data() : NULL);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgesDeviceView(
    typename Base::ConstOffsetsDeviceView &offsets,
    typename Base::ConstIdsDeviceView &adjIds) const {

  offsets = offsDevice_;
  adjIds = adjIdsDevice_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgesHostView(
    typename Base::ConstOffsetsHostView &offsets,
    typename Base::ConstIdsHostView &adjIds) const {

  auto hostIDs = Kokkos::create_mirror_view(adjIdsDevice_);
  Kokkos::deep_copy(hostIDs, adjIdsDevice_);
  adjIds = hostIDs;

  auto hostOffsets = Kokkos::create_mirror_view(offsDevice_);
  Kokkos::deep_copy(hostOffsets, offsDevice_);
  offsets = hostOffsets;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
int TpetraRowGraphAdapter<User, UserCoord>::getNumWeightsPerVertex() const {
  return nWeightsPerVertex_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexWeightsView(
    const scalar_t *&weights, int &stride, int idx) const {

  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid vertex weight index.");

  size_t length;
  vertexWeights_[idx].getStridedList(length, weights, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexWeightsDeviceView(
    typename Base::WeightsDeviceView1D &weights, int idx) const {
  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid vertex weight index.");

  const auto size = vertexWeightsDevice_.extent(0);
  weights = typename Base::WeightsDeviceView1D("weights", size);

  Kokkos::parallel_for(
      size, KOKKOS_CLASS_LAMBDA(const int id) {
        weights(id) = vertexWeightsDevice_(id, idx);
      });

  Kokkos::fence();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexWeightsDeviceView(
    typename Base::WeightsDeviceView &weights) const {

  weights = vertexWeightsDevice_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexWeightsHostView(
    typename Base::WeightsHostView1D &weights, int idx) const {
  AssertCondition((idx >= 0) and (idx < nWeightsPerVertex_),
                  "Invalid vertex weight index.");

  auto weightsDevice = typename Base::WeightsDeviceView1D(
      "weights", vertexWeightsDevice_.extent(0));
  getVertexWeightsDeviceView(weightsDevice, idx);

  weights = Kokkos::create_mirror_view(weightsDevice);
  Kokkos::deep_copy(weights, weightsDevice);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getVertexWeightsHostView(
    typename Base::WeightsHostView &weights) const {

  weights = Kokkos::create_mirror_view(vertexWeightsDevice_);
  Kokkos::deep_copy(weights, vertexWeightsDevice_);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
bool TpetraRowGraphAdapter<User, UserCoord>::useDegreeAsVertexWeight(
    int idx) const {
  return vertexDegreeWeightsHost_(idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
int TpetraRowGraphAdapter<User, UserCoord>::getNumWeightsPerEdge() const {
  return nWeightsPerEdge_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgeWeightsView(
    const scalar_t *&weights, int &stride, int idx) const {
  AssertCondition((idx >= 0) and (idx < nWeightsPerEdge_),
                  "Invalid edge weight index.");

  size_t length;
  edgeWeights_[idx].getStridedList(length, weights, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgeWeightsDeviceView(
    typename Base::WeightsDeviceView1D &weights, int idx) const {

  weights = Kokkos::subview(edgeWeightsDevice_, Kokkos::ALL, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgeWeightsDeviceView(
    typename Base::WeightsDeviceView &weights) const {

  weights = edgeWeightsDevice_;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgeWeightsHostView(
    typename Base::WeightsHostView1D &weights, int idx) const {

  auto weightsDevice = Kokkos::subview(edgeWeightsDevice_, Kokkos::ALL, idx);
  weights = Kokkos::create_mirror_view(weightsDevice);
  Kokkos::deep_copy(weights, weightsDevice);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void TpetraRowGraphAdapter<User, UserCoord>::getEdgeWeightsHostView(
    typename Base::WeightsHostView &weights) const {

  weights = Kokkos::create_mirror_view(edgeWeightsDevice_);
  Kokkos::deep_copy(weights, edgeWeightsDevice_);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
template <typename Adapter>
void TpetraRowGraphAdapter<User, UserCoord>::applyPartitioningSolution(
    const User &in, User *&out,
    const PartitioningSolution<Adapter> &solution) const {
  // Get an import list (rows to be received)
  size_t numNewVtx;
  ArrayRCP<gno_t> importList;
  try {
    numNewVtx =
        Zoltan2::getImportList<Adapter, TpetraRowGraphAdapter<User, UserCoord>>(
            solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new graph.
  RCP<User> outPtr = doMigration(in, numNewVtx, importList.getRawPtr());
  out = outPtr.get();
  outPtr.release();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
template <typename Adapter>
void TpetraRowGraphAdapter<User, UserCoord>::applyPartitioningSolution(
    const User &in, RCP<User> &out,
    const PartitioningSolution<Adapter> &solution) const {
  // Get an import list (rows to be received)
  size_t numNewVtx;
  ArrayRCP<gno_t> importList;
  try {
    numNewVtx =
        Zoltan2::getImportList<Adapter, TpetraRowGraphAdapter<User, UserCoord>>(
            solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new graph.
  out = doMigration(in, numNewVtx, importList.getRawPtr());
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
RCP<User> TpetraRowGraphAdapter<User, UserCoord>::doMigration(
    const User &from, size_t numLocalRows, const gno_t *myNewRows) const {
  typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
  typedef Tpetra::CrsGraph<lno_t, gno_t, node_t> tcrsgraph_t;

  // We cannot create a Tpetra::RowGraph, unless the underlying type is
  // something we know (like Tpetra::CrsGraph).
  // If the underlying type is something different, the user probably doesn't
  // want a Tpetra::CrsGraph back, so we throw an error.

  // Try to cast "from" graph to a TPetra::CrsGraph
  // If that fails we throw an error.
  // We could cast as a ref which will throw std::bad_cast but with ptr
  // approach it might be clearer what's going on here
  const tcrsgraph_t *pCrsGraphSrc = dynamic_cast<const tcrsgraph_t *>(&from);

  if (!pCrsGraphSrc) {
    throw std::logic_error("TpetraRowGraphAdapter cannot migrate data for "
                           "your RowGraph; it can migrate data only for "
                           "Tpetra::CrsGraph.  "
                           "You can inherit from TpetraRowGraphAdapter and "
                           "implement migration for your RowGraph.");
  }

  // source map
  const RCP<const map_t> &smap = from.getRowMap();
  int oldNumElts = smap->getLocalNumElements();
  gno_t numGlobalRows = smap->getGlobalNumElements();
  gno_t base = smap->getMinAllGlobalIndex();

  // target map
  ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
  const RCP<const Teuchos::Comm<int>> &comm = from.getComm();
  RCP<const map_t> tmap = rcp(new map_t(numGlobalRows, rowList, base, comm));

  // importer
  Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

  // number of entries in my new rows
  typedef Tpetra::Vector<gno_t, lno_t, gno_t, node_t> vector_t;
  vector_t numOld(smap);
  vector_t numNew(tmap);
  for (int lid = 0; lid < oldNumElts; lid++) {
    numOld.replaceGlobalValue(smap->getGlobalElement(lid),
                              from.getNumEntriesInLocalRow(lid));
  }
  numNew.doImport(numOld, importer, Tpetra::INSERT);

  size_t numElts = tmap->getLocalNumElements();
  ArrayRCP<const gno_t> nnz;
  if (numElts > 0)
    nnz = numNew.getData(0); // hangs if vector len == 0

  ArrayRCP<const size_t> nnz_size_t;

  if (numElts && sizeof(gno_t) != sizeof(size_t)) {
    size_t *vals = new size_t[numElts];
    nnz_size_t = arcp(vals, 0, numElts, true);
    for (size_t i = 0; i < numElts; i++) {
      vals[i] = static_cast<size_t>(nnz[i]);
    }
  } else {
    nnz_size_t = arcp_reinterpret_cast<const size_t>(nnz);
  }

  // target graph
  RCP<tcrsgraph_t> G = rcp(new tcrsgraph_t(tmap, nnz_size_t()));

  G->doImport(*pCrsGraphSrc, importer, Tpetra::INSERT);
  G->fillComplete();
  return Teuchos::rcp_dynamic_cast<User>(G);
}

} // namespace Zoltan2

#endif

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

#ifndef MUELU_ZOLTAN2GRAPHADAPTER_HPP_
#define MUELU_ZOLTAN2GRAPHADAPTER_HPP_

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_ZOLTAN2)

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Xpetra_Map.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include "MueLu_GraphBase.hpp"

// Zoltab2 InputTraits for MueLu Graph objects
namespace Zoltan2 {

template <typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<MueLu::GraphBase<LocalOrdinal, GlobalOrdinal, Node> > {
  typedef Zoltan2::default_scalar_t scalar_t;
  typedef LocalOrdinal lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t part_t;
  typedef Node node_t;
  static inline std::string name() { return "MueLu::Graph"; }

  Z2_STATIC_ASSERT_TYPES  // validate the types
};
}  // end namespace Zoltan2

namespace MueLu {

template <typename User, typename UserCoord = User>
class MueLuGraphBaseAdapter : public Zoltan2::GraphAdapter<User, UserCoord> {
 public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Zoltan2::InputTraits<User>::scalar_t scalar_t;
  typedef typename Zoltan2::InputTraits<User>::offset_t offset_t;
  typedef typename Zoltan2::InputTraits<User>::lno_t lno_t;
  typedef typename Zoltan2::InputTraits<User>::gno_t gno_t;
  typedef typename Zoltan2::InputTraits<User>::part_t part_t;
  typedef typename Zoltan2::InputTraits<User>::node_t node_t;
  typedef User xgraph_t;
  typedef User user_t;
  typedef UserCoord userCoord_t;
#endif

  //! MueLu::GraphBase Compatibility Layer
  const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const { return graph_->GetComm(); }
  const Teuchos::RCP<const Xpetra::Map<lno_t, gno_t, node_t> > getRowMap() const { return graph_->GetDomainMap(); }
  const RCP<const Xpetra::Map<lno_t, gno_t, node_t> > getColMap() const {
    // For some GraphBases' this is a ColMap, in others it is a seperate map that is
    // only non-null in parallel.
    Teuchos::RCP<const Xpetra::Map<lno_t, gno_t, node_t> > map = graph_->GetImportMap();
    if (map.is_null()) map = graph_->GetDomainMap();
    return map;
  }
  size_t getLocalNumEntries() const { return graph_->GetNodeNumEdges(); }
  size_t getLocalNumRows() const { return getRowMap()->getLocalNumElements(); }
  size_t getLocalNumCols() const { return getColMap()->getLocalNumElements(); }

  void getLocalRowView(lno_t LocalRow, Teuchos::ArrayView<const lno_t> &indices) const {
    indices = graph_->getNeighborVertices(LocalRow);
  }

  /*! \brief Destructor
   */
  ~MueLuGraphBaseAdapter() {}

  /*! \brief Constructor for graph with no weights or coordinates.
   *  \param ingraph the Epetra_CrsGraph, Tpetra::CrsGraph or Xpetra::CrsGraph
   *  \param numVtxWeights  the number of weights per vertex (default = 0)
   *  \param numEdgeWeights the number of weights per edge  (default = 0)
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  MueLuGraphBaseAdapter(const RCP<const User> &ingraph,
                        int nVtxWeights = 0, int nEdgeWeights = 0);

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

  /*! \brief Access to Xpetra-wrapped user's graph.
   */
  RCP<const xgraph_t> getXpetraGraph() const { return graph_; }

  /*! \brief Access to user's graph
   */
  RCP<const User> getUserGraph() const { return ingraph_; }

  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  ////////////////////////////////////////////////////
  // The GraphAdapter interface.
  ////////////////////////////////////////////////////

  // TODO:  Assuming rows == objects;
  // TODO:  Need to add option for columns or nonzeros?
  size_t getLocalNumVertices() const { return getLocalNumRows(); }

  void getVertexIDsView(const gno_t *&ids) const {
    ids = NULL;
    if (getLocalNumVertices())
      ids = getRowMap()->getLocalElementList().getRawPtr();
  }

  size_t getLocalNumEdges() const { return getLocalNumEntries(); }

  void getEdgesView(const offset_t *&offsets, const gno_t *&adjIds) const {
    offsets = offs_.getRawPtr();
    adjIds  = (getLocalNumEdges() ? adjids_.getRawPtr() : NULL);
  }

  int getNumWeightsPerVertex() const { return nWeightsPerVertex_; }

  void getVertexWeightsView(const scalar_t *&weights, int &stride,
                            int idx) const {
    if (idx < 0 || idx >= nWeightsPerVertex_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vertex weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }

    size_t length;
    vertexWeights_[idx].getStridedList(length, weights, stride);
  }

  bool useDegreeAsVertexWeight(int idx) const { return vertexDegreeWeight_[idx]; }

  int getNumWeightsPerEdge() const { return nWeightsPerEdge_; }

  void getEdgeWeightsView(const scalar_t *&weights, int &stride, int idx) const {
    if (idx < 0 || idx >= nWeightsPerEdge_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid edge weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }

    size_t length;
    edgeWeights_[idx].getStridedList(length, weights, stride);
  }

  template <typename Adapter>
  void applyPartitioningSolution(const User &in, User *&out,
                                 const Zoltan2::PartitioningSolution<Adapter> &solution) const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, std::invalid_argument, "applyPartitionlingSolution not implemeneted");
  }

  template <typename Adapter>
  void applyPartitioningSolution(const User &in, RCP<User> &out,
                                 const Zoltan2::PartitioningSolution<Adapter> &solution) const {
    TEUCHOS_TEST_FOR_EXCEPTION(1, std::invalid_argument, "applyPartitionlingSolution not implemeneted");
  }

 private:
  RCP<const User> ingraph_;
  RCP<const xgraph_t> graph_;
  RCP<const Teuchos::Comm<int> > comm_;

  ArrayRCP<const offset_t> offs_;
  ArrayRCP<const gno_t> adjids_;

  int nWeightsPerVertex_;
  ArrayRCP<Zoltan2::StridedData<lno_t, scalar_t> > vertexWeights_;
  ArrayRCP<bool> vertexDegreeWeight_;

  int nWeightsPerEdge_;
  ArrayRCP<Zoltan2::StridedData<lno_t, scalar_t> > edgeWeights_;

  int coordinateDim_;
  ArrayRCP<Zoltan2::StridedData<lno_t, scalar_t> > coords_;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User, typename UserCoord>
MueLuGraphBaseAdapter<User, UserCoord>::MueLuGraphBaseAdapter(
    const RCP<const User> &ingraph, int nVtxWgts, int nEdgeWgts)
  : ingraph_(ingraph)
  , graph_()
  , comm_()
  , offs_()
  , adjids_()
  , nWeightsPerVertex_(nVtxWgts)
  , vertexWeights_()
  , vertexDegreeWeight_()
  , nWeightsPerEdge_(nEdgeWgts)
  , edgeWeights_()
  , coordinateDim_(0)
  , coords_() {
  typedef Zoltan2::StridedData<lno_t, scalar_t> input_t;
  graph_ = ingraph;

  comm_         = getRowMap()->getComm();
  size_t nvtx   = getLocalNumRows();
  size_t nedges = getLocalNumEntries();

  // Unfortunately we have to copy the offsets and edge Ids
  // because edge Ids are not usually stored in vertex id order.
  size_t n = nvtx + 1;
  offs_.resize(n);
  offset_t *offs = const_cast<offset_t *>(offs_.getRawPtr());
  gno_t *adjids  = 0;
  if (nedges > 0) {
    adjids_.resize(nedges);
    adjids = const_cast<gno_t *>(adjids_.getRawPtr());
  }

  offs[0] = 0;
  for (size_t v = 0; v < nvtx; v++) {
    ArrayView<const lno_t> nbors;
    getLocalRowView(v, nbors);
    offs[v + 1] = offs[v] + nbors.size();
    for (offset_t e = offs[v], i = 0; e < offs[v + 1]; e++) {
      adjids[e] = getColMap()->getGlobalElement(nbors[i++]);
    }
  }

  if (nWeightsPerVertex_ > 0) {
    vertexWeights_ =
        arcp(new input_t[nWeightsPerVertex_], 0, nWeightsPerVertex_, true);
    vertexDegreeWeight_ =
        arcp(new bool[nWeightsPerVertex_], 0, nWeightsPerVertex_, true);
    for (int i = 0; i < nWeightsPerVertex_; i++)
      vertexDegreeWeight_[i] = false;
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void MueLuGraphBaseAdapter<User, UserCoord>::setWeights(
    const scalar_t *weightVal, int stride, int idx) {
  if (this->getPrimaryEntityType() == Zoltan2::GRAPH_VERTEX)
    setVertexWeights(weightVal, stride, idx);
  else
    setEdgeWeights(weightVal, stride, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void MueLuGraphBaseAdapter<User, UserCoord>::setVertexWeights(
    const scalar_t *weightVal, int stride, int idx) {
  typedef Zoltan2::StridedData<lno_t, scalar_t> input_t;

  if (idx < 0 || idx >= nWeightsPerVertex_) {
    std::ostringstream emsg;
    emsg << __FILE__ << ":" << __LINE__
         << "  Invalid vertex weight index " << idx << std::endl;
    throw std::runtime_error(emsg.str());
  }

  size_t nvtx = getLocalNumVertices();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nvtx * stride, false);
  vertexWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void MueLuGraphBaseAdapter<User, UserCoord>::setWeightIsDegree(
    int idx) {
  if (this->getPrimaryEntityType() == Zoltan2::GRAPH_VERTEX)
    setVertexWeightIsDegree(idx);
  else {
    std::ostringstream emsg;
    emsg << __FILE__ << "," << __LINE__
         << " error:  setWeightIsNumberOfNonZeros is supported only for"
         << " vertices" << std::endl;
    throw std::runtime_error(emsg.str());
  }
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void MueLuGraphBaseAdapter<User, UserCoord>::setVertexWeightIsDegree(
    int idx) {
  if (idx < 0 || idx >= nWeightsPerVertex_) {
    std::ostringstream emsg;
    emsg << __FILE__ << ":" << __LINE__
         << "  Invalid vertex weight index " << idx << std::endl;
    throw std::runtime_error(emsg.str());
  }

  vertexDegreeWeight_[idx] = true;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
void MueLuGraphBaseAdapter<User, UserCoord>::setEdgeWeights(
    const scalar_t *weightVal, int stride, int idx) {
  typedef Zoltan2::StridedData<lno_t, scalar_t> input_t;

  if (idx < 0 || idx >= nWeightsPerEdge_) {
    std::ostringstream emsg;
    emsg << __FILE__ << ":" << __LINE__
         << "  Invalid edge weight index " << idx << std::endl;
    throw std::runtime_error(emsg.str());
  }

  size_t nedges = getLocalNumEdges();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nedges * stride, false);
  edgeWeights_[idx] = input_t(weightV, stride);
}

}  // namespace MueLu

#endif  // HAVE_MUELU_ZOLTAN2

#endif

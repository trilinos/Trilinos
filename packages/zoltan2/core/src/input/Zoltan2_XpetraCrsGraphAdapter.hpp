// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsGraphAdapter.hpp
    \brief Defines XpetraCrsGraphAdapter class.
*/

#ifndef _ZOLTAN2_XPETRACRSGRAPHADAPTER_HPP_
#define _ZOLTAN2_XPETRACRSGRAPHADAPTER_HPP_

#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_PartitioningHelpers.hpp>

#include <Xpetra_CrsGraph.hpp>

namespace Zoltan2 {

/*!  \brief Provides access for Zoltan2 to Xpetra::CrsGraph data.

    \todo test for memory alloc failure when we resize a vector
    \todo we assume FillComplete has been called.  Should we support
                objects that are not FillCompleted.

    The template parameter is the user's input object:
     \li Tpetra::CrsGraph
     \li Xpetra::CrsGraph
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

template <typename User, typename UserCoord=User>
  class XpetraCrsGraphAdapter : public GraphAdapter<User,UserCoord> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::offset_t    offset_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
  typedef User user_t;
  typedef UserCoord userCoord_t;

  using Base = GraphAdapter<User,UserCoord>;
#endif

  /*! \brief Constructor for graph with no weights or coordinates.
   *  \param ingraph the Epetra_CrsGraph, Tpetra::CrsGraph or Xpetra::CrsGraph
   *  \param numVtxWeights  the number of weights per vertex (default = 0)
   *  \param numEdgeWeights the number of weights per edge  (default = 0)
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphAdapter(const RCP<const User> &ingraph,
                        int nVtxWeights=0, int nEdgeWeights=0);

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
  void setWeightsDevice(typename Base::ConstWeightsDeviceView& val, int idx) {}
  void setWeightsHost(typename Base::ConstWeightsHostView& val, int idx) {}

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
  void setVertexWeightsDevice(typename Base::ConstWeightsDeviceView& val, int idx);
  void setVertexWeightsHost(typename Base::ConstWeightsHostView& val, int idx);

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
  void setEdgeWeightsDevice(typename Base::ConstWeightsDeviceView& val, int idx);
  void setEdgeWeightsHost(typename Base::ConstWeightsHostView& val, int idx);

  /*! \brief Access to Xpetra-wrapped user's graph.
   */
  RCP<const xgraph_t> getXpetraGraph() const { return graph_; }

  /*! \brief Access to user's graph
   */
  RCP<const User> getUserGraph() const { return ingraph_; }

  ////////////////////////////////////////////////////
  // The GraphAdapter interface.
  ////////////////////////////////////////////////////

  // TODO:  Assuming rows == objects;
  // TODO:  Need to add option for columns or nonzeros?
  size_t getLocalNumVertices() const override { return graph_->getLocalNumRows(); }

  void getVertexIDsView(const gno_t *&ids) const override
  {
    ids = NULL;
    if (getLocalNumVertices())
      ids = graph_->getRowMap()->getLocalElementList().getRawPtr();
  }

  size_t getLocalNumEdges() const override { return graph_->getLocalNumEntries(); }

  void getEdgesView(const offset_t *&offsets, const gno_t *&adjIds) const override
  {
    offsets = offs_.getRawPtr();
    adjIds = (getLocalNumEdges() ? adjids_.getRawPtr() : NULL);
  }

  int getNumWeightsPerVertex() const override { return nWeightsPerVertex_;}

  void getVertexWeightsView(const scalar_t *&weights, int &stride,
                            int idx) const override
  {
    if(idx<0 || idx >= nWeightsPerVertex_)
    {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vertex weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }


    size_t length;
    vertexWeights_[idx].getStridedList(length, weights, stride);
  }

  bool useDegreeAsVertexWeight(int idx) const override {return vertexDegreeWeight_[idx];}

  int getNumWeightsPerEdge() const override { return nWeightsPerEdge_;}

  void getEdgeWeightsView(const scalar_t *&weights, int &stride, int idx) const override
  {
    if(idx<0 || idx >= nWeightsPerEdge_)
    {
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
      const PartitioningSolution<Adapter> &solution) const;

  template <typename Adapter>
    void applyPartitioningSolution(const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const;

private:

  RCP<const User > ingraph_;
  RCP<const xgraph_t > graph_;
  RCP<const Comm<int> > comm_;

  ArrayRCP<const offset_t> offs_;
  ArrayRCP<const gno_t> adjids_;

  int nWeightsPerVertex_;
  ArrayRCP<StridedData<lno_t, scalar_t> > vertexWeights_;
  ArrayRCP<bool> vertexDegreeWeight_;

  int nWeightsPerEdge_;
  ArrayRCP<StridedData<lno_t, scalar_t> > edgeWeights_;

  int coordinateDim_;
  ArrayRCP<StridedData<lno_t, scalar_t> > coords_;

};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User, typename UserCoord>
  XpetraCrsGraphAdapter<User,UserCoord>::XpetraCrsGraphAdapter(
    const RCP<const User> &ingraph, int nVtxWgts, int nEdgeWgts):
      ingraph_(ingraph), graph_(), comm_() , offs_(), adjids_(),
      nWeightsPerVertex_(nVtxWgts), vertexWeights_(), vertexDegreeWeight_(),
      nWeightsPerEdge_(nEdgeWgts), edgeWeights_(),
      coordinateDim_(0), coords_()
{
  typedef StridedData<lno_t,scalar_t> input_t;

  try {
    graph_ = rcp_const_cast<const xgraph_t>(
           XpetraTraits<User>::convertToXpetra(rcp_const_cast<User>(ingraph)));
  }
  Z2_FORWARD_EXCEPTIONS

  comm_ = graph_->getComm();
  size_t nvtx = graph_->getLocalNumRows();
  size_t nedges = graph_->getLocalNumEntries();

  // Unfortunately we have to copy the offsets and edge Ids
  // because edge Ids are not usually stored in vertex id order.

  size_t n = nvtx + 1;
  offset_t *offs = new offset_t [n];

  if (!offs)
  {
    std::cerr << "Error: " << __FILE__ << ", " << __LINE__<< std::endl;
    std::cerr << n << " objects" << std::endl;
    throw std::bad_alloc();
  }

  gno_t *adjids = NULL;
  if (nedges)
  {
    adjids = new gno_t [nedges];

    if (!adjids)
    {
      std::cerr << "Error: " << __FILE__ << ", " << __LINE__<< std::endl;
      std::cerr << nedges << " objects" << std::endl;
      throw std::bad_alloc();
    }
  }

  offs[0] = 0;
  for (size_t v=0; v < nvtx; v++){
    ArrayView<const lno_t> nbors;
    graph_->getLocalRowView(v, nbors);
    offs[v+1] = offs[v] + nbors.size();
    for (offset_t e=offs[v], i=0; e < offs[v+1]; e++)
      adjids[e] = graph_->getColMap()->getGlobalElement(nbors[i++]);
  }

  offs_ = arcp(offs, 0, n, true);
  adjids_ = arcp(adjids, 0, nedges, true);

  if (nWeightsPerVertex_ > 0) {
    vertexWeights_ =
          arcp(new input_t[nWeightsPerVertex_], 0, nWeightsPerVertex_, true);
    vertexDegreeWeight_ =
          arcp(new bool[nWeightsPerVertex_], 0, nWeightsPerVertex_, true);
    for (int i=0; i < nWeightsPerVertex_; i++)
      vertexDegreeWeight_[i] = false;
  }

  if (nWeightsPerEdge_ > 0)
    edgeWeights_ = arcp(new input_t[nWeightsPerEdge_], 0, nWeightsPerEdge_, true);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsGraphAdapter<User,UserCoord>::setWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  if (this->getPrimaryEntityType() == GRAPH_VERTEX)
    setVertexWeights(weightVal, stride, idx);
  else
    setEdgeWeights(weightVal, stride, idx);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsGraphAdapter<User,UserCoord>::setVertexWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  typedef StridedData<lno_t,scalar_t> input_t;

  if(idx<0 || idx >= nWeightsPerVertex_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vertex weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }

  size_t nvtx = getLocalNumVertices();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nvtx*stride, false);
  vertexWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsGraphAdapter<User,UserCoord>::setWeightIsDegree(
    int idx)
{
  if (this->getPrimaryEntityType() == GRAPH_VERTEX)
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
  void XpetraCrsGraphAdapter<User,UserCoord>::setVertexWeightIsDegree(
    int idx)
{
  if(idx<0 || idx >= nWeightsPerVertex_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vertex weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }

  vertexDegreeWeight_[idx] = true;
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  void XpetraCrsGraphAdapter<User,UserCoord>::setEdgeWeights(
    const scalar_t *weightVal, int stride, int idx)
{
  typedef StridedData<lno_t,scalar_t> input_t;

  if(idx<0 || idx >= nWeightsPerEdge_)
  {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid edge weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
  }

  size_t nedges = getLocalNumEdges();
  ArrayRCP<const scalar_t> weightV(weightVal, 0, nedges*stride, false);
  edgeWeights_[idx] = input_t(weightV, stride);
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  template<typename Adapter>
    void XpetraCrsGraphAdapter<User,UserCoord>::applyPartitioningSolution(
      const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewVtx;
  ArrayRCP<gno_t> importList;
  try{
    numNewVtx = Zoltan2::getImportList<Adapter,
                                       XpetraCrsGraphAdapter<User,UserCoord> >
                                      (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new graph.
  RCP<User> outPtr = XpetraTraits<User>::doMigration(in, numNewVtx,
                                                     importList.getRawPtr());
  out = outPtr.get();
  outPtr.release();
}

////////////////////////////////////////////////////////////////////////////
template <typename User, typename UserCoord>
  template<typename Adapter>
    void XpetraCrsGraphAdapter<User,UserCoord>::applyPartitioningSolution(
      const User &in, RCP<User> &out,
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list (rows to be received)
  size_t numNewVtx;
  ArrayRCP<gno_t> importList;
  try{
    numNewVtx = Zoltan2::getImportList<Adapter,
                                       XpetraCrsGraphAdapter<User,UserCoord> >
                                      (solution, this, importList);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Move the rows, creating a new graph.
  out = XpetraTraits<User>::doMigration(in, numNewVtx,
                                        importList.getRawPtr());
}

}  //namespace Zoltan2

#endif

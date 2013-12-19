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

/*! \file Zoltan2_XpetraCrsGraphAdapter.hpp
    \brief Defines XpetraCrsGraphAdapter class.
*/

#ifndef _ZOLTAN2_XPETRACRSGRAPHADAPTER_HPP_
#define _ZOLTAN2_XPETRACRSGRAPHADAPTER_HPP_

#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Util.hpp>

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

template <typename User>
  class XpetraCrsGraphAdapter : public GraphAdapter<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
  typedef GraphAdapter<User> base_adapter_t;
  typedef User user_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraCrsGraphAdapter() { }

  /*! \brief Constructor for graph with no weights or coordinates.
   *  \param ingraph the Epetra_CrsGraph, Tpetra::CrsGraph or Xpetra::CrsGraph
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphAdapter(const RCP<const User> &ingraph);

  /*! \brief Constructor for graph with weights but no coordinates.
   *  \param ingraph the Epetra_CrsGraph, Tpetra::CrsGraph or Xpetra::CrsGraph
   *  \param vWeights  a list of pointers to vertex weights.
   *      The number of weights per graph vertex is assumed to be
   *      \c vWeights.size().
   *  \param vWeightStrides  a list of strides for the \c vWeights.
   *     The weight for weight dimension \c n for vertex \c k should be
   *     found at <tt>vWeights[n][vWeightStrides[n] * k]</tt>.
   *     If \c vWeightStrides.size() is zero, it is assumed all strides are one.
   *  \param eWeights  a list of pointers to edge weights.
   *      The number of weights per edge is assumed to be
   *      \c eWeights.size().
   *  \param eWeightStrides  a list of strides for the \c eWeights.
   *     The weight for weight dimension \c n for edge \c k should be
   *     found at <tt>eWeights[n][eWeightStrides[n] * k]</tt>.
   *     If \c eWeightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The order of the vertex weights should match the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   *
   *  The order of the edge weights should follow the order that the
   *  the vertices and edges appear in the input data structure.
   *
   *  By vertex:
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   *
   *  Then by vertex neighbor:
   *     \code
   *       TheGraph->getLocalRowView(vertexNum, neighborList);
   *     \endcode
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphAdapter(const RCP<const User> &ingraph,
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides);

  /*! \brief Constructor for graph with weights and vertex coordinates.
   *  \param ingraph the Epetra_CrsGraph, Tpetra::CrsGraph or Xpetra::CrsGraph
   *  \param vWeights  a list of pointers to vertex weights.
   *      The number of weights per graph vertex is assumed to be
   *      \c vWeights.size().
   *  \param vWeightStrides  a list of strides for the \c vWeights.
   *     The weight for weight dimension \c n for vertex \c k should be
   *     found at <tt>vWeights[n][vWeightStrides[n] * k]</tt>.
   *     If \c vWeightStrides.size() is zero, it is assumed all strides are one.
   *  \param eWeights  a list of pointers to edge weights.
   *      The number of weights per edge is assumed to be
   *      \c eWeights.size().
   *  \param eWeightStrides  a list of strides for the \c eWeights.
   *     The weight for weight dimension \c n for edge \c k should be
   *     found at <tt>eWeights[n][eWeightStrides[n] * k]</tt>.
   *     If \c eWeightStrides.size() is zero, it is assumed all strides are one.
   *  \param coords  a list of pointers to vertex coordinates.
   *      The coordinate dimension is assumed to be \c coords.size().
   *  \param coordStrides  a list of strides for the \c coords.
   *     The coordinate for dimension \c n for vertex \c k should be
   *     found at <tt>coords[n][coordStrides[n] * k]</tt>.
   *     If \c coordStrides.size() is zero, it is assumed all strides are one.
   *
   *  The order of the vertex weights and coordinates should coorespond to
   *  the order that vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   *
   *  The order of the edge weights should follow the order that the
   *  the vertices and edges appear in the input data structure.
   *
   *  By vertex:
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   *
   *  Then by vertex neighbor:
   *     \code
   *       TheGraph->getLocalRowView(vertexNum, neighborList);
   *     \endcode
   *
   * Most adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphAdapter(const RCP<const User> &ingraph,
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides,
    vector<const scalar_t *> &coords,  vector<int> &coordStrides);

  /*! \brief Provide a pointer to one dimension of the vertex weights.
   *    \param val A pointer to the weights for dimension \c dim.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th vertex for dimension \dim.
   *    \param idx A number from 0 to one less than 
   *          vertex weight dimension specified in the constructor.
   *
   *  The order of the vertex weights should match the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   */

  void setVertexWeights(const scalar_t *val, int stride, int idx);

  /*! \brief Provide a pointer to one dimension of the edge weights.
   *    \param val A pointer to the weights for dimension \c dim.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th edge for dimension \dim.
   *    \param dim A number from 0 to one less than 
   *          edge weight dimension specified in the constructor.
   *
   *  The order of the edge weights should follow the order that the
   *  the vertices and edges appear in the input data structure.
   *
   *  By vertex:
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   *
   *  Then by vertex neighbor:
   *     \code
   *       TheGraph->getLocalRowView(vertexNum, neighborList);
   *     \endcode
   */

  void setEdgeWeights(const scalar_t *val, int stride, int idx);

  /*! \brief Provide a pointer to one dimension of the vertex coordinates.
   *    \param dim A number from 0 to one less than 
   *          vertex coordinate dimension specified in the constructor.
   *    \param val A pointer to the coordinates for dimension \c dim.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the coordinate for the
   *             \c n th vertex.
   *
   *  The order of the vertex coordinates should coorespond to the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   */

  void setVertexCoordinates(int dim, const scalar_t *val, int stride);

  /*! \brief Access to Xpetra-wrapped user's graph.
   */ 
  RCP<const xgraph_t> getXpetraGraph() const
  {
    return graph_;
  }

  /*! \brief Access to user's graph 
   */ 
  RCP<const User> getUserGraph() const
  {
    return ingraph_;
  }

  ////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////

  ////////////////////////////////////////////////////
  // The GraphAdapter interface.
  ////////////////////////////////////////////////////

  // TODO:  Assuming rows == objects; 
  // TODO:  Need to add option for columns or nonzeros?
  size_t getLocalNumVertices() const { return graph_->getNodeNumRows(); }

  void getVertexIDsView(const gid_t *&ids) const 
  {
    ids = NULL;
    if (getLocalNumVertices())
      ids = graph_->getRowMap()->getNodeElementList().getRawPtr();
  }

  size_t getLocalNumEdges() const { return graph_->getNodeNumEntries(); }

  void getEdgeView(const lno_t *&offsets, const gid_t *&adjIds) const
  {
    adjIds = NULL;
    offsets = NULL;
    if (getLocalNumVertices()) {
      offsets = offs_.getRawPtr();
      adjIds = adjids_.getRawPtr();
    }
  }

  int getNumWeightsPerVertex() const { return vertexWeightDim_;}

  void getVertexWeightsView(const scalar_t *&weights, int &stride,
                            int idx) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid weight index",
      idx >= 0 && idx < vertexWeightDim_, BASIC_ASSERTION);
    size_t length;
    vertexWeights_[idx].getStridedList(length, weights, stride);
  }


  int getNumWeightsPerEdge() const { return edgeWeightDim_;}

  void getEdgeWeightsView(const scalar_t *&weights, int &stride, int idx) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid weight index",
      idx >= 0 && idx < edgeWeightDim_, BASIC_ASSERTION);
    size_t length;
    edgeWeights_[idx].getStridedList(length, weights, stride);
  }


  int getCoordinateDimension() const { return coordinateDim_; }

  void getVertexCoordinatesView(const scalar_t *&coords, int &stride,
                                int idx) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, 
      "invalid coordinate dimension",
      idx >= 0 && idx < coordinateDim_, BASIC_ASSERTION);
    size_t length;
    coords_[idx].getStridedList(length, coords, stride);
  }

  template <typename Adapter>
    void applyPartitioningSolution(const User &in, User *&out,
      const PartitioningSolution<Adapter> &solution) const;

private:

  void initializeData(
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides,
    vector<const scalar_t *> &coords,  vector<int> &coordStrides);

  RCP<const User > ingraph_;
  RCP<const xgraph_t > graph_;
  RCP<const Comm<int> > comm_;

  ArrayRCP<const lno_t> offs_;
  ArrayRCP<const gid_t> adjids_;

  int vertexWeightDim_;
  ArrayRCP<StridedData<lno_t, scalar_t> > vertexWeights_;

  int edgeWeightDim_;
  ArrayRCP<StridedData<lno_t, scalar_t> > edgeWeights_;

  int coordinateDim_;
  ArrayRCP<StridedData<lno_t, scalar_t> > coords_;

  // A default Environment for error messages.  User-written
  // InputAdapter classes can use some other error return convention
  // if desired.
  RCP<const Environment> env_;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User>
  XpetraCrsGraphAdapter<User>::XpetraCrsGraphAdapter(
    const RCP<const User> &ingraph):
      ingraph_(ingraph), graph_(), comm_() , offs_(), adjids_(),
      vertexWeightDim_(0), vertexWeights_(),
      edgeWeightDim_(0), edgeWeights_(),
      coordinateDim_(0), coords_(),
      env_(rcp(new Environment))
{
  vector<const scalar_t *> emptyValues;
  vector<int> emptyStrides;

  initializeData(emptyValues, emptyStrides, emptyValues, emptyStrides,
    emptyValues, emptyStrides);
}

template <typename User>
  XpetraCrsGraphAdapter<User>::XpetraCrsGraphAdapter(
  const RCP<const User> &ingraph,
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides):
      ingraph_(ingraph), graph_(), comm_() , offs_(), adjids_(),
      vertexWeightDim_(vWeights.size()), vertexWeights_(),
      edgeWeightDim_(eWeights.size()), edgeWeights_(),
      coordinateDim_(0), coords_(),
      env_(rcp(new Environment))
{
  vector<const scalar_t *> emptyValues;
  vector<int> emptyStrides;

  initializeData(vWeights, vWeightStrides, eWeights, eWeightStrides,
    emptyValues, emptyStrides);
}

template <typename User>
  XpetraCrsGraphAdapter<User>::XpetraCrsGraphAdapter(
    const RCP<const User> &ingraph,
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides,
    vector<const scalar_t *> &coords,  vector<int> &coordStrides):
      ingraph_(ingraph), graph_(), comm_() , offs_(), adjids_(),
      vertexWeightDim_(vWeights.size()), vertexWeights_(),
      edgeWeightDim_(eWeights.size()), edgeWeights_(),
      coordinateDim_(coords.size()), coords_(),
      env_(rcp(new Environment))
{
  initializeData(vWeights, vWeightStrides, eWeights, eWeightStrides,
    coords, coordStrides);
}

template <typename User>
  void XpetraCrsGraphAdapter<User>::initializeData(
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides,
    vector<const scalar_t *> &coords,  vector<int> &coordStrides)
{
  env_->localInputAssertion(__FILE__, __LINE__, 
    "invalid number of dimensions", 
    vertexWeightDim_ >= 0 && edgeWeightDim_ >= 0 && coordinateDim_ >= 0, 
    BASIC_ASSERTION);

  typedef StridedData<lno_t,scalar_t> input_t;

  if (vertexWeightDim_)
    vertexWeights_ = 
      arcp(new input_t [vertexWeightDim_], 0, vertexWeightDim_, true);

  if (edgeWeightDim_)
    edgeWeights_ = 
      arcp(new input_t [edgeWeightDim_], 0, edgeWeightDim_, true);

  if (coordinateDim_)
    coords_ = 
      arcp(new input_t [coordinateDim_], 0, coordinateDim_, true);

  graph_ = XpetraTraits<User>::convertToXpetra(ingraph_);
  comm_ = graph_->getComm();
  size_t nvtx = graph_->getNodeNumRows();
  size_t nedges = graph_->getNodeNumEntries();

  // Unfortunately we have to copy the offsets and edge Ids
  // because edge Ids are not usually stored in vertex id order.

  size_t n = nvtx + 1;
  lno_t *offs = new lno_t [n];
  env_->localMemoryAssertion(__FILE__, __LINE__, n, offs);

  gid_t *adjids = NULL;
  if (nedges){
    adjids = new gid_t [nedges];
    env_->localMemoryAssertion(__FILE__, __LINE__, nedges, adjids);
  }

  offs[0] = 0;
  for (size_t v=0; v < nvtx; v++){
    ArrayView<const lno_t> nbors;
    graph_->getLocalRowView(v, nbors);
    offs[v+1] = offs[v] + nbors.size();
    for (lno_t e=offs[v], i=0; e < offs[v+1]; e++)
      adjids[e] = graph_->getColMap()->getGlobalElement(nbors[i++]);
  }

  offs_ = arcp(offs, 0, n, true);
  adjids_ = arcp(adjids, 0, nedges, true);

  int stride = 1;
  for (int dim=0; dim < coordinateDim_; dim++){
    if (coordStrides.size())
      stride = coordStrides[dim];
    ArrayRCP<const scalar_t> coordV(coords[dim], 0, nvtx, false);
    coords_[dim] = input_t(coordV, stride);
  }

  stride = 1;
  for (int dim=0; dim < vertexWeightDim_; dim++){
    if (vWeightStrides.size())
      stride = vWeightStrides[dim];
    ArrayRCP<const scalar_t> wgtV(vWeights[dim], 0, nvtx, false);
    vertexWeights_[dim] = input_t(wgtV, stride);
  }

  stride = 1;
  for (int dim=0; dim < edgeWeightDim_; dim++){
    if (eWeightStrides.size())
      stride = eWeightStrides[dim];
    ArrayRCP<const scalar_t> ewgtV(eWeights[dim], 0, nedges, false);
    edgeWeights_[dim] = input_t(ewgtV, stride);
  }
}

template <typename User>
  template<typename Adapter>
    void XpetraCrsGraphAdapter<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list

  size_t len = solution.getLocalNumberOfIds();
  const gid_t *gids = solution.getIdList();
  const partId_t *parts = solution.getPartList();
  ArrayRCP<gid_t> gidList = arcp(const_cast<gid_t *>(gids), 0, len, false);
  ArrayRCP<partId_t> partList = arcp(const_cast<partId_t *>(parts), 0, len, 
    false);

  ArrayRCP<lno_t> dummyIn;
  ArrayRCP<gid_t> importList;
  ArrayRCP<lno_t> dummyOut;
  size_t numNewVtx;
  const RCP<const Comm<int> > comm = graph_->getComm();

  try{
    numNewVtx = solution.convertSolutionToImportList(
      0, dummyIn, importList, dummyOut);
  }
  Z2_FORWARD_EXCEPTIONS;

  RCP<const User> inPtr = rcp(&in, false);

  RCP<const User> outPtr = XpetraTraits<User>::doMigration(
   inPtr, numNewVtx, importList.getRawPtr());

  out = const_cast<User *>(outPtr.get());
  outPtr.release();
}
  
}  //namespace Zoltan2
  
#endif

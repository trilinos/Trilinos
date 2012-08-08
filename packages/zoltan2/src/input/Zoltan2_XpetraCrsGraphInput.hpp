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

/*! \file Zoltan2_XpetraCrsGraphInput.hpp
    \brief Defines XpetraCrsGraphInput class.
*/

#ifndef _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_

#include <Zoltan2_GraphInput.hpp>
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
  class XpetraCrsGraphInput : public GraphInput<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
  typedef GraphInput<User>       base_adapter_t;
  typedef User user_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraCrsGraphInput() { }

  /*! \brief Constructor for graph with no weights or coordinates.
   *  \param ingraph the Epetra_CrsGraph, Tpetra::CrsGraph or Xpetra::CrsGraph
   *
   * Most input adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphInput(const RCP<const User> &ingraph);

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
   * Most input adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphInput(const RCP<const User> &ingraph,
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
   * Most input adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphInput(const RCP<const User> &ingraph,
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides,
    vector<const scalar_t *> &coords,  vector<int> &coordStrides);

  /*! \brief Provide a pointer to one dimension of the vertex weights.
   *    \param dim A number from 0 to one less than 
   *          vertex weight dimension specified in the constructor.
   *    \param val A pointer to the weights for dimension \c dim.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th vertex for dimension \dim.
   *
   *  The order of the vertex weights should match the order that
   *  vertices appear in the input data structure.
   *     \code
   *       TheGraph->getRowMap()->getNodeElementList()
   *     \endcode
   */

  void setVertexWeights(int dim, const scalar_t *val, int stride);

  /*! \brief Provide a pointer to one dimension of the edge weights.
   *    \param dim A number from 0 to one less than 
   *          edge weight dimension specified in the constructor.
   *    \param val A pointer to the weights for dimension \c dim.
   *    \param stride    A stride for the \c val array.  If \stride is
   *             \c k, then val[n * k] is the weight for the
   *             \c n th edge for dimension \dim.
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

  void setEdgeWeights(int dim, const scalar_t *val, int stride);

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
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  string inputAdapterName()const { return string("XpetraCrsGraph");}

  size_t getLocalNumberOfObjects() const { return getLocalNumberOfVertices();}

  int getNumberOfWeightsPerObject() const { return 0;}

  size_t getObjectWeights(int dim, const scalar_t *&wgt, int &stride) const
  {
    return getVertexWeights(dim, wgt, stride);
  }

  ////////////////////////////////////////////////////
  // The GraphInput interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumberOfVertices() const { 
    return graph_->getNodeNumRows(); 
  }

  size_t getLocalNumberOfEdges() const { 
    return graph_->getNodeNumEntries();
  }

  int getVertexWeightDimension() const { 
    return vertexWeightDim_;
  }

  int getEdgeWeightDimension() const { 
    return edgeWeightDim_;
  }

  int getCoordinateDimension() const { 
    return coordinateDim_;
  }

  size_t getVertexListView(const gid_t *&ids,
    const lno_t *&offsets, const gid_t *& edgeId) const
  {
    size_t nvtx = getLocalNumberOfVertices();
    ids = edgeId = NULL;
    offsets = NULL;

    if (nvtx){
      ids = graph_->getRowMap()->getNodeElementList().getRawPtr();
      offsets = offs_.getRawPtr();
      edgeId = eids_.getRawPtr();
    }
    
    return nvtx;
  }

  size_t getVertexWeights(int dim,
    const scalar_t *&weights, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid weight dimension",
      dim >= 0 && dim < vertexWeightDim_, BASIC_ASSERTION);

    size_t length;
    vertexWeights_[dim].getStridedList(length, weights, stride);
    return length;
  }

  size_t getEdgeWeights(int dim,
    const scalar_t *&weights, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid weight dimension",
      dim >= 0 && dim < edgeWeightDim_, BASIC_ASSERTION);

    size_t length;
    edgeWeights_[dim].getStridedList(length, weights, stride);
    return length;
  }

  size_t getVertexCoordinates(int dim,
    const scalar_t *&coords, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, 
      "invalid coordinate dimension",
      dim >= 0 && dim < coordinateDim_, BASIC_ASSERTION);

    size_t length;
    coords_[dim].getStridedList(length, coords, stride);
    return length;
  }

 /*! \brief Repartition a graph that has the same structure as
   *   the graph that instantiated this input adapter.
   */
  template<typename Adapter>
    size_t applyPartitioningSolution(const User &in, User *&out,
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
  ArrayRCP<const gid_t> eids_;

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
  XpetraCrsGraphInput<User>::XpetraCrsGraphInput(
    const RCP<const User> &ingraph):
      ingraph_(ingraph), graph_(), comm_() , offs_(), eids_(),
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
  XpetraCrsGraphInput<User>::XpetraCrsGraphInput(
  const RCP<const User> &ingraph,
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides):
      ingraph_(ingraph), graph_(), comm_() , offs_(), eids_(),
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
  XpetraCrsGraphInput<User>::XpetraCrsGraphInput(
    const RCP<const User> &ingraph,
    vector<const scalar_t *> &vWeights,  vector<int> &vWeightStrides,
    vector<const scalar_t *> &eWeights,  vector<int> &eWeightStrides,
    vector<const scalar_t *> &coords,  vector<int> &coordStrides):
      ingraph_(ingraph), graph_(), comm_() , offs_(), eids_(),
      vertexWeightDim_(vWeights.size()), vertexWeights_(),
      edgeWeightDim_(eWeights.size()), edgeWeights_(),
      coordinateDim_(coords.size()), coords_(),
      env_(rcp(new Environment))
{
  initializeData(vWeights, vWeightStrides, eWeights, eWeightStrides,
    coords, coordStrides);
}

template <typename User>
  void XpetraCrsGraphInput<User>::initializeData(
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

  gid_t *eids = NULL;
  if (nedges){
    eids = new gid_t [nedges];
    env_->localMemoryAssertion(__FILE__, __LINE__, nedges, eids);
  }

  offs[0] = 0;
  for (size_t v=0; v < nvtx; v++){
    ArrayView<const lno_t> nbors;
    graph_->getLocalRowView(v, nbors);
    offs[v+1] = offs[v] + nbors.size();
    for (lno_t e=offs[v], i=0; e < offs[v+1]; e++)
      eids[e] = graph_->getColMap()->getGlobalElement(nbors[i++]);
  }

  offs_ = arcp(offs, 0, n, true);
  eids_ = arcp(eids, 0, nedges, true);

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
    size_t XpetraCrsGraphInput<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<Adapter> &solution) const
{
  // Get an import list

  Zoltan2::Environment env;
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

  gno_t lsum = numNewVtx;
  gno_t gsum = 0;
  reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum, &gsum);

  RCP<const User> inPtr = rcp(&in, false);

  RCP<const User> outPtr = XpetraTraits<User>::doMigration(
   inPtr, lsum, importList.getRawPtr());

  out = const_cast<User *>(outPtr.get());
  outPtr.release();
  return numNewVtx;
}
  
}  //namespace Zoltan2
  
#endif

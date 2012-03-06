// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsGraphInput.hpp
    \brief Defines XpetraCrsGraphInput class.
*/

#ifndef _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_XPETRACRSGRAPHINPUT_HPP_

#include <Zoltan2_GraphInput.hpp>
#include <Zoltan2_StridedInput.hpp>
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
*/

template <typename User>
class XpetraCrsGraphInput : public GraphInput<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
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

  /*! \brief Constructor
   *  \param ingraph the Epetra_CrsGraph, Tpetra::CrsGraph or Xpetra::CrsGraph
   *  \param vertexWeightDim  the number of weights per vertex 
            that will be supplied in setVertexWeights()
   *  \param edgeWeightDim  the number of weights per edge 
            that will be supplied in setEdgeWeights()
   *  \param vertexCoordinateDim the number of coordinates per vertex 
            that will be supplied in setCoordinates()
   *
   * Most input adapters do not have RCPs in their interface.  This
   * one does because the user is obviously a Trilinos user.
   */

  XpetraCrsGraphInput(const RCP<const User> &ingraph, 
    int vertexWeightDim=0, int edgeWeightDim=0, int vertexCoordinateDim=0);

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

  std::string inputAdapterName()const { return std::string("XpetraCrsGraph");}

  size_t getLocalNumberOfObjects() const { return getLocalNumberOfVertices();}

  int getNumberOfWeightsPerObject() const { return 0;}

  ////////////////////////////////////////////////////
  // The GraphInput interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumberOfVertices() const { 
    return graph_->getNodeNumRows(); 
  }

  global_size_t getGlobalNumberOfVertices() const { 
    return graph_->getGlobalNumRows(); 
  }

  size_t getLocalNumberOfEdges() const { 
    return graph_->getNodeNumEntries();
  }

  global_size_t getGlobalNumberOfEdges() const { 
    return graph_->getGlobalNumEntries();
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
    vertexWeights_[dim]->getStridedList(length, weights, stride);
    return length;
  }

  size_t getEdgeWeights(int dim,
    const scalar_t *&weights, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid weight dimension",
      dim >= 0 && dim < edgeWeightDim_, BASIC_ASSERTION);

    size_t length;
    edgeWeights_[dim]->getStridedList(length, weights, stride);
    return length;
  }

  size_t getVertexCoordinates(int dim,
    const scalar_t *&coords, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, 
      "invalid coordinate dimension",
      dim >= 0 && dim < coordinateDim_, BASIC_ASSERTION);

    size_t length;
    coords_[dim]->getStridedList(length, coords, stride);
    return length;
  }

 /*! \brief Repartition a graph that has the same structure as
   *   the graph that instantiated this input adapter.
   */
  template<typename User2>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<User2> &solution) const;

private:

  RCP<const User > ingraph_;
  RCP<const xgraph_t > graph_;
  RCP<const Comm<int> > comm_;

  ArrayRCP<const lno_t> offs_;
  ArrayRCP<const gid_t> eids_;

  int vertexWeightDim_;
  Array<RCP<StridedInput<lno_t, scalar_t> > > vertexWeights_;

  int edgeWeightDim_;
  Array<RCP<StridedInput<lno_t, scalar_t> > > edgeWeights_;

  int coordinateDim_;
  Array<RCP<StridedInput<lno_t, scalar_t> > > coords_;

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
    const RCP<const User> &ingraph, 
    int vertexWeightDim, 
    int edgeWeightDim, 
    int vertexCoordinateDim):
      ingraph_(ingraph), graph_(), comm_() , offs_(), eids_(),
      vertexWeightDim_(vertexWeightDim), vertexWeights_(vertexWeightDim),
      edgeWeightDim_(edgeWeightDim), edgeWeights_(edgeWeightDim),
      coordinateDim_(vertexCoordinateDim), coords_(vertexCoordinateDim),
      env_(rcp(new Environment))
{
  env_->localInputAssertion(__FILE__, __LINE__, 
    "invalid number of dimensions", 
    vertexWeightDim >= 0 && edgeWeightDim >= 0 && vertexCoordinateDim >= 0, 
    BASIC_ASSERTION);

  graph_ = XpetraTraits<User>::convertToXpetra(ingraph);
  comm_ = graph_->getComm();
  size_t nvtx = graph_->getNodeNumRows();
  size_t nedges = graph_->getNodeNumEntries();
  Environment env;

  // Unfortunately we have to copy the offsets and edge Ids
  // because edge Ids are not usually stored in vertex id order.

  size_t n = nvtx + 1;
  lno_t *offs = new lno_t [n];
  env.localMemoryAssertion(__FILE__, __LINE__, n, offs);

  gid_t *eids = NULL;
  if (nedges){
    eids = new gid_t [nedges];
    env.localMemoryAssertion(__FILE__, __LINE__, nedges, eids);
  }

  offs[0] = 0;
  for (lno_t v=0; v < nvtx; v++){
    ArrayView<const lno_t> nbors;
    graph_->getLocalRowView(v, nbors);
    offs[v+1] = offs[v] + nbors.size();
    for (lno_t e=offs[v], i=0; e < offs[v+1]; e++)
      eids[e] = graph_->getColMap()->getGlobalElement(nbors[i++]);
  }

  offs_ = arcp(offs, 0, n, true);
  eids_ = arcp(eids, 0, nedges, true);
}

template <typename User>
  void XpetraCrsGraphInput<User>::setVertexCoordinates(int dim, 
    const scalar_t *coordVal, int stride)
{
  typedef StridedInput<lno_t,scalar_t> input_t;

  env_->localInputAssertion(__FILE__, __LINE__, 
    "invalid vertex coordinate dimension", 
    dim >= 0 && dim < coordinateDim_, BASIC_ASSERTION);

  size_t nvtx = getLocalNumberOfVertices();

  coords_[dim] = 
    rcp<input_t>(
      new input_t(ArrayView<const scalar_t>(coordVal, nvtx), stride));
}

template <typename User>
  void XpetraCrsGraphInput<User>::setEdgeWeights(int dim, 
    const scalar_t *val, int stride)
{
  typedef StridedInput<lno_t,scalar_t> input_t;

  env_->localInputAssertion(__FILE__, __LINE__, 
    "invalid edge weight dimension", 
    dim >= 0 && dim < edgeWeightDim_, BASIC_ASSERTION);

  size_t nedges = getLocalNumberOfEdges();

  edgeWeights_[dim] = 
    rcp<input_t>(
      new input_t(ArrayView<const scalar_t>(val, nedges), stride));
}

template <typename User>
  void XpetraCrsGraphInput<User>::setVertexWeights(int dim, 
    const scalar_t *val, int stride)
{
  typedef StridedInput<lno_t,scalar_t> input_t;

  env_->localInputAssertion(__FILE__, __LINE__, 
    "invalid vertex weight dimension", 
    dim >= 0 && dim < vertexWeightDim_, BASIC_ASSERTION);

  size_t nvtx = getLocalNumberOfVertices();

  vertexWeights_[dim] = 
    rcp<input_t>(
      new input_t(ArrayView<const scalar_t>(val, nvtx), stride));
}

template <typename User>
  template<typename User2>
    size_t XpetraCrsGraphInput<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<User2> &solution) const
{
  // Get an import list

  Zoltan2::Environment env;
  size_t len = solution.getNumberOfIds();
  const gid_t *gids = solution.getGlobalIdList();
  const size_t *parts = solution.getPartList();
  ArrayRCP<gid_t> gidList = arcp(const_cast<gid_t *>(gids), 0, len, false);
  ArrayRCP<size_t> partList = arcp(const_cast<size_t *>(parts), 0, len, false);

  ArrayRCP<lno_t> dummyIn;
  ArrayRCP<gid_t> importList;
  ArrayRCP<lno_t> dummyOut;
  size_t numNewVtx;
  const RCP<const Comm<int> > comm = graph_->getComm();

  try{
    numNewVtx = convertPartListToImportList<gid_t, lno_t, lno_t>(
      *comm, partList, gidList, dummyIn, importList, dummyOut);
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

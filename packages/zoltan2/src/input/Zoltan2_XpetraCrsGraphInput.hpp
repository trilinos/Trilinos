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
    std::vector<const scalar_t *> &vWeights,  std::vector<int> &vWeightStrides,
    std::vector<const scalar_t *> &eWeights,  std::vector<int> &eWeightStrides);

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
    std::vector<const scalar_t *> &vWeights,  std::vector<int> &vWeightStrides,
    std::vector<const scalar_t *> &eWeights,  std::vector<int> &eWeightStrides,
    std::vector<const scalar_t *> &coords,  std::vector<int> &coordStrides);

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

  void initializeData(
    std::vector<const scalar_t *> &vWeights,  std::vector<int> &vWeightStrides,
    std::vector<const scalar_t *> &eWeights,  std::vector<int> &eWeightStrides,
    std::vector<const scalar_t *> &coords,  std::vector<int> &coordStrides);

  RCP<const User > ingraph_;
  RCP<const xgraph_t > graph_;
  RCP<const Comm<int> > comm_;

  ArrayRCP<const lno_t> offs_;
  ArrayRCP<const gid_t> eids_;

  int vertexWeightDim_;
  Array<RCP<StridedData<lno_t, scalar_t> > > vertexWeights_;

  int edgeWeightDim_;
  Array<RCP<StridedData<lno_t, scalar_t> > > edgeWeights_;

  int coordinateDim_;
  Array<RCP<StridedData<lno_t, scalar_t> > > coords_;

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
      vertexWeightDim_(0), vertexWeights_(0),
      edgeWeightDim_(0), edgeWeights_(0),
      coordinateDim_(0), coords_(0),
      env_(rcp(new Environment))
{
  std::vector<const scalar_t *> emptyValues;
  std::vector<int> emptyStrides;

  initializeData(emptyValues, emptyStrides, emptyValues, emptyStrides,
    emptyValues, emptyStrides);
}

template <typename User>
  XpetraCrsGraphInput<User>::XpetraCrsGraphInput(const RCP<const User> &ingraph,
    std::vector<const scalar_t *> &vWeights,  std::vector<int> &vWeightStrides,
    std::vector<const scalar_t *> &eWeights,  std::vector<int> &eWeightStrides):
      ingraph_(ingraph), graph_(), comm_() , offs_(), eids_(),
      vertexWeightDim_(vWeights.size()), vertexWeights_(vWeights.size()),
      edgeWeightDim_(eWeights.size()), edgeWeights_(eWeights.size()),
      coordinateDim_(0), coords_(0),
      env_(rcp(new Environment))
{
  std::vector<const scalar_t *> emptyValues;
  std::vector<int> emptyStrides;

  initializeData(vWeights, vWeightStrides, eWeights, eWeightStrides,
    emptyValues, emptyStrides);
}

template <typename User>
  XpetraCrsGraphInput<User>::XpetraCrsGraphInput(const RCP<const User> &ingraph,
    std::vector<const scalar_t *> &vWeights,  std::vector<int> &vWeightStrides,
    std::vector<const scalar_t *> &eWeights,  std::vector<int> &eWeightStrides,
    std::vector<const scalar_t *> &coords,  std::vector<int> &coordStrides):
      ingraph_(ingraph), graph_(), comm_() , offs_(), eids_(),
      vertexWeightDim_(vWeights.size()), vertexWeights_(vWeights.size()),
      edgeWeightDim_(eWeights.size()), edgeWeights_(eWeights.size()),
      coordinateDim_(coords.size()), coords_(coords.size()),
      env_(rcp(new Environment))
{
  initializeData(vWeights, vWeightStrides, eWeights, eWeightStrides,
    coords, coordStrides);
}

template <typename User>
  void XpetraCrsGraphInput<User>::initializeData(
    std::vector<const scalar_t *> &vWeights,  std::vector<int> &vWeightStrides,
    std::vector<const scalar_t *> &eWeights,  std::vector<int> &eWeightStrides,
    std::vector<const scalar_t *> &coords,  std::vector<int> &coordStrides)
{
  typedef StridedData<lno_t,scalar_t> input_t;
  env_->localInputAssertion(__FILE__, __LINE__, 
    "invalid number of dimensions", 
    vertexWeightDim_ >= 0 && edgeWeightDim_ >= 0 && coordinateDim_ >= 0, 
    BASIC_ASSERTION);

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
  for (lno_t v=0; v < nvtx; v++){
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
    coords_[dim] = 
      rcp<input_t>(
        new input_t(ArrayView<const scalar_t>(coords[dim], nvtx), stride));
  }

  stride = 1;
  for (int dim=0; dim < vertexWeightDim_; dim++){
    if (vWeightStrides.size())
      stride = vWeightStrides[dim];
    vertexWeights_[dim] = 
      rcp<input_t>(
        new input_t(ArrayView<const scalar_t>(vWeights[dim], nvtx), stride));
  }

  stride = 1;
  for (int dim=0; dim < edgeWeightDim_; dim++){
    if (eWeightStrides.size())
      stride = eWeightStrides[dim];
    edgeWeights_[dim] = 
      rcp<input_t>(
        new input_t(ArrayView<const scalar_t>(eWeights[dim], nedges), stride));
  }
}

template <typename User>
  template<typename User2>
    size_t XpetraCrsGraphInput<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<User2> &solution) const
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
    numNewVtx = convertSolutionToImportList<User2, lno_t>(
      solution, dummyIn, importList, dummyOut);
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

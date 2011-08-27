// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GraphInput.hpp

    \brief The abstract interface for a graph input adapter.
*/


#ifndef _ZOLTAN2_GRAPHINPUT_HPP_
#define _ZOLTAN2_GRAPHINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphInput
    \brief The GraphInput is the abstract base class for weighted graph input adapters.

    The Graph accessor methods defined here mimic those of Tpetra::CrsGraph and Tpetra::Map.

    These public methods define the graph adapter interface to Zoltan2 models.  They
    uses the internal local and global IDs used by Zoltan2.  However the concrete subclasses 
    of GraphInput will define in addition an interface to the application providing
    the graph.  That interface will use the application's local and global ID types.

    TODO: It may be necessary to put a migration interface at this level.
*/

template<typename Scalar, typename LNO, typename GNO>
  class GraphInput : public InputAdapter{
private:

public:

  /*! Returns the global number of vertices in the graph.
   */
  virtual GNO getGlobalNumVertices() const = 0;

  /*! Returns the global number edges in the graph.
   */
  virtual GNO getGlobalNumEdges() const = 0;

  /*! Returns the number vertices on this process.
   */
  virtual LNO getLocalNumVertices() const = 0;

  /*! Returns the number edges on this process.
   */
  virtual LNO getLocalNumEdges() const = 0;

  virtual LNO getVertexWeightDim() const = 0;
  virtual LNO getEdgeWeightDim() const = 0;

  /*! Returns the minimum global vertex ID on this process.
   */
  virtual GNO getMinGlobalVertex() const = 0;

  /*! Returns the maximum global vertex ID on this process.
   */
  virtual GNO getMaxGlobalVertex() const = 0;

  /*! Returns the minimum local global vertex ID everywhere.
   */
  virtual GNO getMinAllGlobalVertex() const = 0;

  /*! Returns the maximum local global vertex ID everywhere.
   */
  virtual GNO getMaxAllGlobalVertex() const = 0;

  /*! Translate vertex ID from global to local.
      \param vtxID The global ID of a vertex on this process.
      \result the local ID for the specified vertex.
   */
  virtual LNO getLocalVertex(GNO vtxID) const = 0;

  /*! Translate vertex ID from local to global.
      \param vtxID The local ID of a vertex on this process.
      \result the global ID for the specified vertex.
   */
  virtual GNO getGlobalVertex(LNO vtxID) const = 0;

  /*! Find the owner of and local ID of a list vertices.
      \param vtxID   A list of global vertex IDs
      \param nodeID  A view of memory allocated by the caller
      \param nodeVtxLID A view of memory allocated by the caller
      \result nodeID  A list of the process rank for the owner corresponding to each vertex in vtxID
      \result nodeVtxLID  A list of the local ID at the owning process for each vertex vtxID

       This method must be called by all processes in the communicator.
   */
  virtual void getRemoteVertexList const Teuchos::ArrayView<const GNO> &vtxID,
    const Teuchos::ArrayView<const int> &nodeID,
    const Teuchos::ArrayView<const LNO> &nodeVtxLID) const = 0;

  /*! Find the process rank for the owner of each vertex in a list.
      \param vtxID   A list of global vertex IDs.
      \param nodeID  A view of memory allocated by the caller.
      \result nodeID  A list of the process rank for the owner corresponding to each vertex in vtxID.

       This method must be called by all processes in the communicator.
   */
  virtual void getRemoteVertexList const Teuchos::ArrayView<const GNO> &vtxID,
    const Teuchos::ArrayView<const int> &nodeID) const = 0;

  /*! Get the list of vertex global IDs and their weights.
      \param ids a view of an array allocated by the caller.
      \param wgts a view of an array allocated by the caller.
      \result ids a list of the global IDs of each vertex on this process.
      \result wgts a list of the weight or weights associated with each vertex in the ids list.  Weights are
                     listed by vertex by weight component.

      This method throws std::runtime_error if there are no weights associated with the IDs.
   */
  virtual void getGlobalVertexList(Teuchos::ArrayView<const GN0> &ids, Teuchos::ArrayView<const Scalar> &wgt) const = 0;

  /*! Get the list of vertex global IDs.
      \param ids a view of an array allocated by the caller.
      \result ids a list of the global IDs of each vertex on this process.
   */
  virtual void getGlobalVertexList(Teuchos::ArrayView<const GN0> &ids) const = 0;

  /*! Does the supplied local vertex ID exist on this process.
      \result true if the ID is a vertex local ID on this process.
      \result false if the ID is not a vertex local ID on this process.
   */
  virtual bool isNodeLocalVertex(LNO vtxID) const = 0;

  /*! Does the supplied global vertex ID exist on this process.
      \result true if the ID is a vertex global ID on this process.
      \result false if the ID is not a vertex global ID on this process.
   */
  virtual bool isNodeGlobalVertex(GNO vtxID) const = 0;

  /*! Are the vertex global IDs contiguous across process?
      \result true if global IDs process contiguously from one process to the next
      \result false if global IDs are not contiguous
   */
  virtual bool isContiguous() const = 0;

  /*! Are the graph distributed across more than this one process.
      \result true if the graph is distributed.
      \result false if the graph is not distributed.
   */
  virtual bool isDistributed() const = 0;

  /*! Returns the sum of the number of neighbors of all vertices
   */
  virtual GNO getGlobalNumberOfNbors() const = 0;

  /*! Returns the sum of the number of neighbors of vertices on this process
   */
  virtual LNO getLocalNumberOfNbors() const = 0;

  /*! Returns the number of neighbors for a given vertex.
      \param vtxID The global ID of a vertex owned by this process.
   */
  virtual LNO getNumNborsOfGlobalVertex(GNO vtxID) const = 0;

  /*! Returns the number of neighbors for a given vertex.
      \param vtxID The local ID of a vertex.
   */
  virtual LNO getNumNborsOfLocalVertex(LNO vtxID) const = 0;

  /*! Returns the global number of self-edges.
   */
  virtual GNO getGlobalNumDiags() const = 0;

  /*! Returns the local number of self-edges.
   */
  virtual LNO getLocalNumDiags() const = 0;

  /*! Returns the global maximum vertex degree.
   */
  virtual LNO getGlobalMaxNumVertexNbors() const = 0;

  /*! Returns the local maximum vertex degree.
   */
  virtual LNO getLocalMaxNumVertexNbors() const = 0;

  /*! Vertex weights are optional in a user supplied graph.
      \return true if vertex weights are available
      \return false if vertex weights are not available
   */
  virtual bool hasVertexWeights() const = 0;

  /*! Edge weights are optional in a user supplied graph.
      \return true if edge weights are available
      \return false if edge weights are not available
   */
  virtual bool hasEdgeWeights() const = 0;

  /*! Obtain a copy of the edge IDs of the input vertex global ID
      \param vtxID  global ID for a vertex on this process
      \param edgeID user allocated array for edge IDs
      \result edgeID global edge IDs are copied to user's array
      \result numEdges  the number of edgeIDs written to the array
   */
  virtual void getGlobalVertexCopy(GNO vtxID, const Teuchos::ArrayView<GNO> &edgeID, LNO &numEdges) const = 0;

  /*! Obtain a copy of the edge IDs of the input vertex local ID
      \param vtxID  a vertex local ID 
      \param edgeID user allocated array for edge IDs
      \result edgeID global edge IDs are copied to user's array
      \result numEdges  the number of edgeIDs written to the array
   */
  virtual void getLocalVertexCopy(LNO vtxID, const Teuchos::ArrayView<GNO> &edgeID, LNO &numEdges) const = 0;

  /*! Obtain a read-only view of the edge IDs of the input vertex
      \param vtxID  global ID for a vertex on this process
      \param edgeID user supplied ArrayView object
      \result edgeID the ArrayView is set to a read-only view of the edge IDs
   */
  virtual void getGlobalVertexView(GNO vtxID, Teuchos::ArrayView<const GNO> &edgeID) const = 0;

  /*! Obtain a read-only view of the edge IDs of the input vertex
      \param vtxID  a local vertex ID
      \param edgeID user supplied ArrayView object
      \result edgeID the ArrayView is set to a read-only view of the edge IDs
   */
  virtual void getLocalVertexView(LNO vtxID, Teuchos::ArrayView<const GNO> &edgeID) const = 0;


};
  
  
}  //namespace Zoltan2
  
#endif

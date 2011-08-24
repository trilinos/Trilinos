// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GraphAdapter_decl.hpp

    \brief The abstract interface for a graph input adapter.
*/


#ifndef _ZOLTAN2_GRAPHADAPTER_DECL_HPP_
#define _ZOLTAN2_GRAPHADAPTER_DECL_HPP_

#include <Teuchos_ArrayView.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphAdapter
    \brief The GraphAdapter is the abstract base class for weighted graph input adapters.

    The data type that the caller uses for local IDs and global IDs may be different
    than the type that Zoltan2 uses for these values.  The first template variable is
    the data type used by the caller for vertex and edge weights.  The following two
    are the data type used by the caller for its local and global IDs.  The
    final two are the data type used by Zoltan2 for local and global IDs.  For these
    the caller must choose an integral type that is large enough to represent the space
    of the caller's local and global identifiers respectively.  For example, if the
    caller's global IDs are std::pair<int, int>, this caller should set the GNO
    to int or long (or long long if needed and if supported).

    Zoltan2 is more efficient if the callers global IDs are an integral type.

    The NODE template parameter should be set if the caller is using Kokkos for
    multicore/manycore optimizations.

    The Graph accessor methods defined here mimic those of Tpetra::CrsGraph and Tpetra::Map.
*/

template<typename Scalar, typename LNO, typename GNO>
  class GraphAdapter : public InputAdapter {
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

  /*! Edge IDs may or may not be supplied in increasing order in queries.
      \return true if edge IDs are supplied in increasing order.
      \return false if edge IDs may not be supplied in increasing order.
   */
  virtual bool isSorted() const = 0;

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

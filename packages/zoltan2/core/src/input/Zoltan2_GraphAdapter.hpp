// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_GraphAdapter.hpp
    \brief Defines the GraphAdapter interface.
*/

#ifndef _ZOLTAN2_GRAPHADAPTER_HPP_
#define _ZOLTAN2_GRAPHADAPTER_HPP_

#include <Zoltan2_Adapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>

namespace Zoltan2 {

/*!  \brief Enumerated entity type for graphs:  Vertices or Edges
 */
enum GraphEntityType { GRAPH_VERTEX, GRAPH_EDGE };

/*!  \brief GraphAdapter defines the interface for graph-based user data.

    Adapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t vertex and edge weights
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c node_t   is a Kokkos Node type

    The Kokkos node type can be safely ignored.

    The template parameter \c User is a user-defined data type
    which, through a traits mechanism, provides the actual data types
    with which the Zoltan2 library will be compiled.
    \c User may be the actual class or structure used by application to
    represent a vector, or it may be the helper class BasicUserTypes.
    See InputTraits for more information.

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
class GraphAdapter : public AdapterWithCoordsWrapper<User, UserCoord> {
private:
  /// Enum to represent the primary entity type in the graph (vertex or edge).
  /// This is the entity that will be partitioned, ordered, colored, matched,
  /// etc.
  enum GraphEntityType primaryEntityType_ = GRAPH_VERTEX;

  /// Enum to represent the adjacency entity type in the graph (edge or vertex).
  /// This typically refers to an entity type opposite to the
  /// primaryEntityType_.
  enum GraphEntityType adjacencyEntityType_ = GRAPH_EDGE;

  /// Pointer to a VectorAdapter containing coordinates of objects of type
  /// primaryEntityType_. This is an optional attribute. \sa
  /// haveCoordinateInput_
  VectorAdapter<UserCoord> *coordinateInput_ = nullptr;

  /// Flag indicating whether the coordinate input is provided.
  /// When this flag is set to true, coordinateInput_ should not be null.
  bool haveCoordinateInput_ = false;

public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using scalar_t = typename InputTraits<User>::scalar_t;
  using lno_t = typename InputTraits<User>::lno_t;
  using gno_t = typename InputTraits<User>::gno_t;
  using node_t = typename InputTraits<User>::node_t;
  using offset_t = typename InputTraits<User>::offset_t;
  using user_t = User;
  using userCoord_t = UserCoord;
  using base_adapter_t = GraphAdapter<User, UserCoord>;
  using Base = AdapterWithCoordsWrapper<User, UserCoord>;
  using VtxDegreeHostView = Kokkos::View<bool *, Kokkos::HostSpace>;
  using device_t = typename node_t::device_type;
#endif

  enum BaseAdapterType adapterType() const override { return GraphAdapterType; }

  ////////////////////////////////////////////////////////////////////////////
  // Methods to be defined in derived classes.

  /*! \brief Returns the number of vertices on this process.
   */
  virtual size_t getLocalNumVertices() const = 0;

  /*! \brief Returns the number of edges on this process.
   */
  virtual size_t getLocalNumEdges() const = 0;

  /*! \brief Sets pointers to this process' graph entries.
      \param vertexIds will on return a pointer to vertex global Ids
   */
  virtual void getVertexIDsView(const gno_t *&vertexIds) const = 0;

  /*! \brief Sets pointers to this process' graph entries.
      \param vertexIds will on return a device Kokkos::View with vertex global
     Ids
   */
  virtual void
  getVertexIDsDeviceView(typename Base::ConstIdsDeviceView &vertexIds) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Sets pointers to this process' graph entries.
      \param vertexIds will on return a host Kokkos::View with vertex global Ids
   */
  virtual void
  getVertexIDsHostView(typename Base::ConstIdsHostView &vertexIds) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Gets adjacency lists for all vertices in a compressed
             sparse row (CSR) format.
      \param offsets is an array of size getLocalNumVertices() + 1.
         The neighboring vertices for vertexId[i]
         begin at adjIds[offsets[i]].
          The last element of offsets is the size of the adjIds array.
      \param adjIds on return will point to the array of adjacent vertices for
         for each vertex.
   */
  virtual void getEdgesView(const offset_t *&offsets,
                            const gno_t *&adjIds) const = 0;

  /*! \brief Gets adjacency lists for all vertices in a compressed
             sparse row (CSR) format.
      \param offsets is device Kokkos::View of size getLocalNumVertices() + 1.
         The neighboring vertices for vertexId[i]
         begin at adjIds[offsets[i]].
         The last element of offsets is the size of the adjIds array.
      \param adjIds Device Kokkos::View of adjacent vertices for for each
     vertex.
   */
  virtual void
  getEdgesDeviceView(typename Base::ConstOffsetsDeviceView &offsets,
                     typename Base::ConstIdsDeviceView &adjIds) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Gets adjacency lists for all vertices in a compressed
             sparse row (CSR) format.
      \param offsets is host Kokkos::View of size getLocalNumVertices() + 1.
         The neighboring vertices for vertexId[i]
         begin at adjIds[offsets[i]].
         The last element of offsets is the size of the adjIds array.
      \param adjIds Host Kokkos::View of adjacent vertices for for each vertex.
   */
  virtual void getEdgesHostView(typename Base::ConstOffsetsHostView &offsets,
                                typename Base::ConstIdsHostView &adjIds) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Returns the number (0 or greater) of weights per vertex
   */
  virtual int getNumWeightsPerVertex() const { return 0; }

  /*! \brief  Provide a pointer to the vertex weights, if any.
      \param weights is the list of weights of the given index for
           the vertices returned in getVertexIDsView().
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than getNumWeightsPerVertex().
   */
  virtual void getVertexWeightsView(const scalar_t *&weights, int &stride,
                                    int /* idx */ = 0) const {
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a device view of the vertex weights, if any.
      \param weights is the list of weights of the given index for
           the vertices returned in getVertexIDsView().
      \param idx ranges from zero to one less than getNumWeightsPerVertex().
   */
  virtual void
  getVertexWeightsDeviceView(typename Base::WeightsDeviceView1D &weights,
                             int /* idx */ = 0) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a device view of the vertex weights, if any.
      \param weights is the view of all the weights for the vertices returned in getVertexIDsView().
   */
  virtual void
  getVertexWeightsDeviceView(typename Base::WeightsDeviceView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a host view of the vertex weights, if any.
      \param weights is the list of weights of the given index for
           the vertices returned in getVertexIDsView().
      \param idx ranges from zero to one less than getNumWeightsPerVertex().
   */
  virtual void
  getVertexWeightsHostView(typename Base::WeightsHostView1D &weights,
                           int /* idx */ = 0) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a host view of the vertex weights, if any.
      \param weights is the list of all the weights for the vertices returned in getVertexIDsView()
   */
  virtual void
  getVertexWeightsHostView(typename Base::WeightsHostView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Indicate whether vertex weight with index idx should be the
   *         global degree of the vertex
   */
  virtual bool useDegreeAsVertexWeight(int /* idx */) const { return false; }

  /*! \brief Returns the number (0 or greater) of edge weights.
   */
  virtual int getNumWeightsPerEdge() const { return 0; }

  /*! \brief  Provide a pointer to the edge weights, if any.
      \param weights is the list of weights of the given index for
           the edges returned in getEdgeView().
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than getNumWeightsPerEdge().
   */
  virtual void getEdgeWeightsView(const scalar_t *&weights, int &stride,
                                  int /* idx */ = 0) const {
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a device view of the edge weights, if any.
      \param weights is the list of weights of the given index for
           the edges returned in getEdgeView().
      \param idx ranges from zero to one less than getNumWeightsPerEdge().
   */
  virtual void
  getEdgeWeightsDeviceView(typename Base::WeightsDeviceView1D &weights,
                           int /* idx */ = 0) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a device view of the edge weights, if any.
      \param weights is the list of weights for the edges returned in getEdgeView().
   */
  virtual void
  getEdgeWeightsDeviceView(typename Base::WeightsDeviceView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a host view of the edge weights, if any.
      \param weights is the list of weights of the given index for
           the edges returned in getEdgeView().
      \param idx ranges from zero to one less than getNumWeightsPerEdge().
   */
  virtual void
  getEdgeWeightsHostView(typename Base::WeightsHostView1D &weights,
                         int /* idx */ = 0) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief  Provide a host view of the edge weights, if any.
      \param weights is the list of weights for the edges returned in getEdgeView().
   */
  virtual void
  getEdgeWeightsHostView(typename Base::WeightsHostView &weights) const {
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Allow user to provide additional data that contains coordinate
   *         info associated with the MatrixAdapter's primaryEntityType_.
   *         Associated data must have the same parallel distribution and
   *         ordering of entries as the primaryEntityType_.
   *
   *  \param coordData is a pointer to a VectorAdapter with the user's
   *         coordinate data.
   */
  void setCoordinateInput(VectorAdapter<UserCoord> *coordData) override {
    coordinateInput_ = coordData;
    haveCoordinateInput_ = true;
  }

  /*! \brief Indicate whether coordinate information has been set for this
   *         MatrixAdapter
   */
  bool coordinatesAvailable() const { return haveCoordinateInput_; }

  /*! \brief Obtain the coordinate data registered by the user.
   *  \return pointer a VectorAdapter with the user's coordinate data.
   */
  VectorAdapter<UserCoord> *getCoordinateInput() const override {
    return coordinateInput_;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Implementations of base-class methods

  /*! \brief Returns the entity to be partitioned, ordered, colored, etc.
   *  Valid values are GRAPH_VERTEX or GRAPH_EDGE.
   */
  inline enum GraphEntityType getPrimaryEntityType() const {
    return this->primaryEntityType_;
  }

  /*! \brief Sets the primary entity type.  Called by algorithm based on
   *  parameter value in parameter list from application.
   *  Also sets to adjacencyEntityType_ to something reasonable:  opposite of
   *  primaryEntityType_.
   */
  void setPrimaryEntityType(const std::string &typestr) {
    if (typestr == "vertex") {
      this->primaryEntityType_ = GRAPH_VERTEX;
      this->adjacencyEntityType_ = GRAPH_EDGE;
    } else if (typestr == "edge") {
      this->primaryEntityType_ = GRAPH_EDGE;
      this->adjacencyEntityType_ = GRAPH_VERTEX;
    } else {
      AssertCondition(true, "Invalid GraphEntityType (" + typestr +
                                "). Valid values are 'vertex' and 'edge'");
    }
  }

  /*! \brief Returns the entity that describes adjacencies between the
   *  entities to be partitioned, ordered, colored, etc.
   *  Valid values are GRAPH_VERTEX or GRAPH_EDGE.
   */
  inline enum GraphEntityType getAdjacencyEntityType() const {
    return this->adjacencyEntityType_;
  }

  /*! \brief Sets the adjacency entity type.  Called by algorithm based on
   *  parameter value in parameter list from application.
   *  Also sets to primaryEntityType_ to something reasonable:  opposite of
   *  adjacencyEntityType_.
   */
  void setAdjacencyEntityType(const std::string &typestr) {
    if (typestr == "vertex") {
      this->adjacencyEntityType_ = GRAPH_VERTEX;
      this->primaryEntityType_ = GRAPH_EDGE;
    } else if (typestr == "edge") {
      this->adjacencyEntityType_ = GRAPH_EDGE;
      this->primaryEntityType_ = GRAPH_VERTEX;
    } else {
      AssertCondition(true, "Invalid GraphEntityType (" + typestr +
                                "). Valid values are 'vertex' and 'edge'");
    }
  }

  // Functions from the BaseAdapter interface
  size_t getLocalNumIDs() const override {
    if (getPrimaryEntityType() == GRAPH_VERTEX)
      return getLocalNumVertices();
    else
      return getLocalNumEdges();
  }

  void getIDsView(const gno_t *&Ids) const override {
    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getIDsView not yet supported for graph edges.");

    getVertexIDsView(Ids);
  }

  void getIDsDeviceView(typename Base::ConstIdsDeviceView &Ids) const override {
    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getIDsDeviceView not yet supported for graph edges.");

    getVertexIDsDeviceView(Ids);
  }

  void getIDsHostView(typename Base::ConstIdsHostView &Ids) const override {
    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getIDsHostView not yet supported for graph edges.");

    getVertexIDsHostView(Ids);
  }

  int getNumWeightsPerID() const override {
    if (getPrimaryEntityType() == GRAPH_VERTEX)
      return getNumWeightsPerVertex();
    else
      return getNumWeightsPerEdge();
  }

  void getWeightsView(const scalar_t *&wgt, int &stride,
                      int idx = 0) const override {

    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getWeightsView not yet supported for graph edges.");

    getVertexWeightsView(wgt, stride, idx);
  }

  void getWeightsHostView(typename Base::WeightsHostView1D &hostWgts,
                          int idx = 0) const override {
    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getWeightsHostView not yet supported for graph edges.");

    getVertexWeightsHostView(hostWgts, idx);
  }

  void getWeightsHostView(typename Base::WeightsHostView &hostWgts) const override {
    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getWeightsHostView not yet supported for graph edges.");

    getVertexWeightsHostView(hostWgts);
  }

  void getWeightsDeviceView(typename Base::WeightsDeviceView1D &deviceWgts,
                            int idx = 0) const override {
    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getWeightsDeviceView not yet supported for graph edges.");

    getVertexWeightsDeviceView(deviceWgts, idx);
  }

  void getWeightsDeviceView(typename Base::WeightsDeviceView &deviceWgts) const override {
    AssertCondition(getPrimaryEntityType() == GRAPH_VERTEX,
                    "getWeightsDeviceView not yet supported for graph edges.");

    getVertexWeightsDeviceView(deviceWgts);
  }

  bool useDegreeAsWeight(int idx) const {
    AssertCondition(this->getPrimaryEntityType() == GRAPH_VERTEX,
                    "useDegreeAsWeight not yet supported for graph edges.");

    return useDegreeAsVertexWeight(idx);
  }
};

} // namespace Zoltan2

#endif

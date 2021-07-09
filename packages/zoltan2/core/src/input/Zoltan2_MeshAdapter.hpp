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

/*! \file Zoltan2_MeshAdapter.hpp
    \brief Defines the MeshAdapter interface.
*/


#ifndef _ZOLTAN2_MESHADAPTER_HPP_
#define _ZOLTAN2_MESHADAPTER_HPP_

#include <Zoltan2_Adapter.hpp>
#include "TpetraExt_MatrixMatrix.hpp"

namespace Zoltan2 {

  /*!  \brief Enumerate entity types for meshes:  Regions, Faces, Edges, or
   *                                              Vertices
   */

enum MeshEntityType {
  MESH_VERTEX,
  MESH_EDGE,
  MESH_FACE,
  MESH_REGION
};

  /*!  \brief Enumerate entity topology types for meshes:
   *          points,lines,polygons,triangles,quadrilaterals,
   *          polyhedrons, tetrahedrons, hexhedrons, prisms, or pyramids
   */

enum EntityTopologyType {
  POINT,         // a 0D entity (e.g. a vertex)
  LINE_SEGMENT,  // a 1D entity (e.g. an edge)
  POLYGON,       // a general 2D entity
  TRIANGLE,      // a specific 2D entity bounded by 3 edge entities
  QUADRILATERAL, // a specific 2D entity bounded by 4 edge entities
  POLYHEDRON,    // a general 3D entity
  TETRAHEDRON,   // a specific 3D entity bounded by 4 triangle entities
  HEXAHEDRON,    // a specific 3D entity bounded by 6 quadrilateral
                 // entities
  PRISM,         // a specific 3D entity bounded by a combination of 3
                 //quadrilateral entities and 2 triangle entities
  PYRAMID        // a specific 3D entity bounded by a combination of 1
                 // quadrilateral entity and 4 triangle entities
};

/*!  \brief MeshAdapter defines the interface for mesh input.

    Adapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t entity and adjacency weights
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c node_t   is a Kokkos CPU node

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

template <typename User>
class MeshAdapter : public BaseAdapter<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::offset_t offset_t;
  typedef typename InputTraits<User>::lno_t lno_t;
  typedef typename InputTraits<User>::gno_t gno_t;
  typedef typename InputTraits<User>::part_t part_t;
  typedef typename InputTraits<User>::node_t node_t;
  typedef User user_t;
  typedef User userCoord_t;
  typedef MeshAdapter<User> base_adapter_t;
#endif

  enum BaseAdapterType adapterType() const {return MeshAdapterType;}

  /*! \brief Destructor
   */
  virtual ~MeshAdapter() {};

  // Default MeshEntityType is MESH_REGION with MESH_FACE-based adjacencies and
  // second adjacencies and coordinates
  MeshAdapter() : primaryEntityType(MESH_REGION),
                  adjacencyEntityType(MESH_FACE),
                  secondAdjacencyEntityType(MESH_FACE) {};

  ////////////////////////////////////////////////////////////////////////////
  // Methods to be defined in derived classes.

  /*! \brief Provide a pointer to the entity topology types
      \param Types will on return point to the list of entity topology types
      for this process.
  */
  virtual bool areEntityIDsUnique(MeshEntityType etype) const
  {
    return etype==this->getPrimaryEntityType();
  }

  /*! \brief Returns the global number of mesh entities of MeshEntityType
   */
  //virtual size_t getGlobalNumOf(MeshEntityType etype) const = 0;

  /*! \brief Returns the number of mesh entities on this process.
   */
  virtual size_t getLocalNumOf(MeshEntityType etype) const = 0;


  /*! \brief Provide a pointer to this process' identifiers.
      \param Ids will on return point to the list of the global Ids for this
       process.
  */
  virtual void getIDsViewOf(MeshEntityType etype,
                            gno_t const *&Ids) const = 0;


  /*! \brief Provide a pointer to the entity topology types
      \param Types will on return point to the list of entity topology types
      for this process.
  */
  virtual void getTopologyViewOf(MeshEntityType etype,
                                 enum EntityTopologyType const *&Types) const
  {
    Types = NULL;
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Return the number of weights per entity.
   *  \return the count of weights, zero or more per entity.
   *   If the number of weights is zero, then we assume that the entities
   *   are equally weighted.
   */
  virtual int getNumWeightsPerOf(MeshEntityType etype) const { return 0; }

  /*! \brief Provide a pointer to one of the number of this process'
                optional entity weights.

      \param weights on return will contain a list of the weights for the
               number specified.

      \param stride on return will indicate the stride of the weights list.

       The k'th weight is located at weights[stride*k].

      \param idx is a value ranging from zero to one less than
                   getNumWeightsPerEntityID()
  */
  virtual void getWeightsViewOf(MeshEntityType etype,
     const scalar_t *&weights, int &stride, int idx = 0) const
  {
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED
  }


  /*! \brief Return dimension of the entity coordinates, if any.
   *
   *  Some algorithms can partition mesh entities using geometric coordinate
   *    information
   *
   *  Some algorithms can use geometric entity coordinate
   *    information if it is present.
   */
  virtual int getDimension() const { return 0; }

  /*! \brief Provide a pointer to one dimension of entity coordinates.
      \param coords  points to a list of coordinate values for the dimension.
      \param stride  describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].
      \param coordDim  is a value from 0 to one less than
         getEntityCoordinateDimension() specifying which dimension is
         being provided in the coords list.
  */
  virtual void getCoordinatesViewOf(MeshEntityType etype,
    const scalar_t *&coords, int &stride, int coordDim) const
  {
    coords = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED
  }


  /*! \brief Returns whether a first adjacency combination is available.
   */
  virtual bool availAdjs(MeshEntityType source, MeshEntityType target) const {
    return false;
  }


  /*! \brief Returns the number of first adjacencies on this process.
   */
  virtual size_t getLocalNumAdjs(MeshEntityType source,
                                 MeshEntityType target) const { return 0;}


  /*! \brief Sets pointers to this process' mesh first adjacencies.
      \param source
      \param offsets is an array of size getLocalNumOf() + 1.
         The first adjacency Ids for Ids[i] (returned in
         getIDsViewOf()) begin at adjacencyIds[offsets[i]].
          The last element of offsets
          is the size of the adjacencyIds array.
      \param adjacencyIds on return will point to the global first adjacency
         Ids for each entity.
  */
  virtual void getAdjsView(MeshEntityType source, MeshEntityType target,
     const offset_t *&offsets, const gno_t *& adjacencyIds) const
  {
    offsets = NULL;
    adjacencyIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED
  }


  /*! \brief Returns whether a second adjacency combination is available.
   *   If combination is not available in the MeshAdapter, Zoltan2 will
   *   compute them, using A^T A, where A is matrix of first adjacencies.
   */
  virtual bool avail2ndAdjs(MeshEntityType sourcetarget,
                            MeshEntityType through) const
  {
    return false;
  }

  /*! \brief if avail2ndAdjs(), returns the number of second adjacencies
   *   on this process.
   */
  virtual size_t getLocalNum2ndAdjs(MeshEntityType sourcetarget,
                                    MeshEntityType through) const
  {
    return 0;
  }

  /*! \brief if avail2ndAdjs(), set pointers to this process' second adjacencies
      \param sourcetarget
      \param offsets is an array of size getLocalNumOf() + 1.
         The second adjacency Ids for Ids[i] (returned in
         getIDsViewOf()) begin at adjacencyIds[offsets[i]].
          The last element of offsets
          is the size of the adjacencyIds array.
      \param adjacencyIds on return will point to the global second adjacency
         Ids for each entity.
   */
  virtual void get2ndAdjsView(MeshEntityType sourcetarget,
                              MeshEntityType through,
                              const offset_t *&offsets,
                              const gno_t *&adjacencyIds) const
  {
    offsets = NULL;
    adjacencyIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED
  }

  /*! \brief Returns the number (0 or greater) of weights per second adjacency.
   *  Note:  second-adjacency weights may be used only if avail2ndAdjs().
   */
  virtual int getNumWeightsPer2ndAdj(MeshEntityType sourcetarget,
                                     MeshEntityType through) const { return 0;}


  /*! \brief  Provide a pointer to the second adjacency weights, if any.
   *  Note:  second-adjacency weights may be used only if avail2ndAdjs().

      \param weights is the list of weights of the given number for
           the second adjacencies returned in get2ndAdjsView().
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than
                   getNumWeightsPer2ndAdj().
   */
  virtual void get2ndAdjWeightsView(MeshEntityType sourcetarget,
                                    MeshEntityType through,
                                    const scalar_t *&weights,
                                    int &stride,
                                    int idx) const
  {
    weights = NULL;
    stride = 0;
    Z2_THROW_NOT_IMPLEMENTED
  }

  ////////////////////////////////////////////////////////////////////////////
  // Implementations of base-class methods

  /*! \brief Returns the entity to be partitioned, ordered, colored, etc.
   */
  inline enum MeshEntityType getPrimaryEntityType() const {
    return this->primaryEntityType;
  }

  /*! \brief Returns the entity that describes adjacencies between the
   *  entities to be partitioned, ordered, colored, etc.
   *  That is, a primaryEntityType that contains an adjacencyEntityType are
   *  adjacent.
   */
  inline enum MeshEntityType getAdjacencyEntityType() const {
    return this->adjacencyEntityType;
  }

  /*! \brief Returns the entity that describes second adjacencies between the
   *  entities to be partitioned, ordered, colored, etc.
   *  That is, two primaryEntityType that share a secondAdjacencyEntityType
   *  are adjacent.
   */
  inline enum MeshEntityType getSecondAdjacencyEntityType() const {
    return this->secondAdjacencyEntityType;
  }

  /*! \brief Sets the primary, adjacency, and second adjacency entity types.
   *  Called by algorithm based on parameter values in parameter list from
   *  application.  Also sets primaryEntityType, adjacencyEntityType, and
   *  secondAdjacencyEntityType to something reasonable:  primaryEntityType not
   *  adjacencyEntityType or secondAdjacencyEntityType.
   */
  void setEntityTypes(std::string ptypestr, std::string atypestr,
                      std::string satypestr) {

    if (ptypestr != atypestr && ptypestr != satypestr) {
      if (ptypestr == "region")
        this->primaryEntityType = MESH_REGION;
      else if (ptypestr == "face")
        this->primaryEntityType = MESH_FACE;
      else if (ptypestr == "edge")
        this->primaryEntityType = MESH_EDGE;
      else if (ptypestr == "vertex")
        this->primaryEntityType = MESH_VERTEX;
      else {
        std::ostringstream emsg;
        emsg << __FILE__ << "," << __LINE__
             << " error:  Invalid MeshEntityType " << ptypestr << std::endl;
        emsg << "Valid values: region  face  edge  vertex" << std::endl;
        throw std::runtime_error(emsg.str());
      }

      if (atypestr == "region")
        this->adjacencyEntityType = MESH_REGION;
      else if (atypestr == "face")
        this->adjacencyEntityType = MESH_FACE;
      else if (atypestr == "edge")
        this->adjacencyEntityType = MESH_EDGE;
      else if (atypestr == "vertex")
        this->adjacencyEntityType = MESH_VERTEX;
      else {
        std::ostringstream emsg;
        emsg << __FILE__ << "," << __LINE__
             << " error:  Invalid MeshEntityType " << atypestr << std::endl;
        emsg << "Valid values: region  face  edge  vertex" << std::endl;
        throw std::runtime_error(emsg.str());
      }

      if (satypestr == "region")
        this->secondAdjacencyEntityType = MESH_REGION;
      else if (satypestr == "face")
        this->secondAdjacencyEntityType = MESH_FACE;
      else if (satypestr == "edge")
        this->secondAdjacencyEntityType = MESH_EDGE;
      else if (satypestr == "vertex")
        this->secondAdjacencyEntityType = MESH_VERTEX;
      else {
        std::ostringstream emsg;
        emsg << __FILE__ << "," << __LINE__
             << " error:  Invalid MeshEntityType " << satypestr << std::endl;
        emsg << "Valid values: region  face  edge  vertex" << std::endl;
        throw std::runtime_error(emsg.str());
      }
    }
    else {
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  PrimaryEntityType " << ptypestr
           << " matches AdjacencyEntityType " << atypestr
           << " or SecondAdjacencyEntityType " << satypestr << std::endl;
      throw std::runtime_error(emsg.str());
    }
  }

  /*! \brief Optional method allowing the idx-th weight of entity type etype
   *  to be set as the number of neighbors (the degree) of the entity
   *  Default is false; user can change in his MeshAdapter implementation.
   */
  virtual bool useDegreeAsWeightOf(MeshEntityType etype, int idx) const
  {
    return false;
  }

  ///////////////////////////////////////////
  // Functions from the BaseAdapter interface
  size_t getLocalNumIDs() const {
    return getLocalNumOf(getPrimaryEntityType());
  }

  void getIDsView(const gno_t *&Ids) const {
    getIDsViewOf(getPrimaryEntityType(), Ids);
  }

  void getIDsKokkosView(Kokkos::View<const gno_t *,
    typename node_t::device_type> &ids) const
  {
    Kokkos::View<gno_t *, typename node_t::device_type>
      kokkos_ids("gids", getLocalNumIDs());
    auto host_kokkos_ids = Kokkos::create_mirror_view(kokkos_ids);

    const gno_t * gnos;
    getIDsView(gnos);
    for(size_t n = 0; n < getLocalNumIDs(); ++n) {
      host_kokkos_ids(n) = gnos[n];
    }
    Kokkos::deep_copy(kokkos_ids, host_kokkos_ids);
    ids = kokkos_ids;
  }

  int getNumWeightsPerID() const {
    return getNumWeightsPerOf(getPrimaryEntityType());
  }

  void getWeightsView(const scalar_t *&wgt, int &stride, int idx = 0) const {
    getWeightsViewOf(getPrimaryEntityType(), wgt, stride, idx);
  }

  void getCoordinatesView(const scalar_t *&coords, int &stride,
                          int coordDim) const
  {
    getCoordinatesViewOf(getPrimaryEntityType(), coords, stride, coordDim);
  }

  inline void getCoordinatesKokkosView(
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type> & elements) const
  {
    // coordinates in MJ are LayoutLeft since Tpetra Multivector gives LayoutLeft
    Kokkos::View<scalar_t **, Kokkos::LayoutLeft, typename node_t::device_type>
      kokkos_coordinates("pamgen coords", getLocalNumIDs(), getDimension());
    auto host_temp_values = Kokkos::create_mirror_view(kokkos_coordinates);
    const scalar_t * coords;
    for(int dim = 0; dim < getDimension(); ++dim) {
      int stride = -1;
      getCoordinatesView(coords, stride, dim);
      for(size_t n = 0; n < getLocalNumIDs(); ++n) {
        host_temp_values(n, dim) = coords[n*stride];
      }
    }
    Kokkos::deep_copy(kokkos_coordinates, host_temp_values);
    elements = kokkos_coordinates;
  }

  bool useDegreeAsWeight(int idx) const
  {
    return useDegreeAsWeightOf(getPrimaryEntityType(), idx);
  }

private:
  enum MeshEntityType primaryEntityType; // Entity type
                                         // to be partitioned, ordered,
                                         // colored, matched, etc.
  enum MeshEntityType adjacencyEntityType; // Entity type defining first-order
                                           // adjacencies; adjacencies are of
                                           // this type.
  enum MeshEntityType secondAdjacencyEntityType; // Bridge entity type
                                                 // defining second-order
                                                 // adjacencies.
};

}  //namespace Zoltan2

#endif

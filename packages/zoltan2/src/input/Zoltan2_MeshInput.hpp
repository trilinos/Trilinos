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
// Questions? Contact Vitus Leung (vjleung@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_MeshAdapter.hpp
    \brief Defines the MeshAdapter interface.
*/


#ifndef _ZOLTAN2_MESHADAPTER_HPP_
#define _ZOLTAN2_MESHADAPTER_HPP_

#include <Zoltan2_InputAdapter.hpp>

#include <string>

namespace Zoltan2 {

enum MeshEntityType {
  MESH_REGION,
  MESH_FACE,
  MESH_EDGE,
  MESH_VERTEX
};
/*!  \brief MeshAdapter defines the interface for mesh input.

    Adapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t entity and adjacency weights
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c gid_t    application global Ids
    \li \c node_t is a sub class of KokkosClassic::StandardNodeMemoryModel

    See IdentifierTraits to understand why the user's global ID type (\c gid_t)
    may differ from that used by Zoltan2 (\c gno_t).

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
  class MeshAdapter : public BaseAdapter {
private:
  enum MeshEntityType primaryEntityType; // Entity to be partitioned, ordered,
                                         // colored, matched, etc.
  enum MeshEntityType adjacencyEntityType; // Entity describing adjacencies
  //KDD Do we need a 2nd-adjacency Entity Type to use as "through"?

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  enum BaseAdapterType adapterType() const {return MeshAdapterType;}

  // Default MeshEntityType is MESH_REGION with MESH_FACE-based adjacencies
  MeshAdapter() : primaryEntityType(MESH_REGION),
                  adjacencyEntityType(MESH_FACE) {};

  /*! \brief Destructor
   */
  virtual ~MeshAdapter() {};

  /*! \brief Returns the number of mesh entities on this process.
   */
  virtual size_t getLocalNumOf(MeshEntityType etype) const = 0;


  /*! \brief Return the number of weights per entity.
   *  \return the count of weights, zero or more per entity.
   *   If the number of weights is zero, then we assume that the entities
   *   are equally weighted.
   */
  virtual int getNumWeightsPerOf(MeshEntityType etype) const = 0;


  /*! \brief Provide a pointer to this process' identifiers.
      \param Ids will on return point to the list of the global Ids for this
       process.
      \return The number of ids in the Ids list.

      Some algorithms can partition a simple list of weighted identifiers
        with no geometry or topology provided.
  */
  virtual size_t getIDsViewOf(MeshEntityType etype,
                              gid_t const *&Ids) const = 0;


  /*! \brief Provide a pointer to one of the number of this process'
                optional entity weights.

      \param weights on return will contain a list of the weights for the
               number specified.  If weights for
           this number are to be uniform for all entities in the
           global problem, the \c weights should be a NULL pointer.

      \param stride on return will indicate the stride of the weights list.

       The k'th weight is located at weights[stride*k].

      \param idx is a value ranging from zero to one less than
                   getNumWeightsPerEntityID()

       \return The number of values in the weights list.  This may be greater
          than the number of entities, because the stride may be greater
          than one.

      Zoltan2 does not copy your data.  The data pointed to by weights
      must remain valid for the lifetime of this Adapter.
  */
  virtual void getWeightsViewOf(MeshEntityType etype,
     const scalar_t *&weights, int &stride, int idx = 0) const = 0;


  /*! \brief Copy the weights for a specific entity into array provided by
                Zoltan2.

      \param id specific entity.  Works only if id is local to this process.

      \param wgts array provided by Zoltan2.
                wgts has size getNumWeightsPerID().

      \return not sure what is needed, if anything.
  */
//KDD  For now, haven't been implementing functions for individual IDs.
//KDD  virtual size_t getEntityIDWeights(gid_t id, scalar_t *wgts) const = 0;


  /*! \brief Return dimension of the entity coordinates, if any.
   *
   *  Some algorithms can partition mesh entities using geometric coordinate
   *    information
   *
   *  Some algorithms can use geometric entity coordinate
   *    information if it is present.
   */
//KDD Assumed uniform dimension of all entity types in other adapters.
//KDD Should we do the same here?
  virtual int getDimensionOf(MeshEntityType etype) const = 0;

  /*! \brief Provide a pointer to one dimension of entity coordinates.
      \param coords  points to a list of coordinate values for the dimension.
      \param stride  describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].
      \param coordDim  is a value from 0 to one less than
         getEntityCoordinateDimension() specifying which dimension is
         being provided in the coords list.

       \return The length of the \c coords list.  This may be more than
              getLocalNumEntityIDs() because the \c stride
              may be more than one.

      Zoltan2 does not copy your data.  The data pointed to coords
      must remain valid for the lifetime of this Adapter.
  */
  virtual void getCoordinatesViewOf(MeshEntityType etype,
    const scalar_t *&coords, int &stride, int coordDim) const = 0;


  /*! \brief Copy the coordinates for a specific ID into array
                provided by Zoltan2.

      \param id specific ID.  Works only if id is local to this        process.

      \param coords array provided by Zoltan2.

      \return not sure what is needed, if anything.
  */
//KDD  For now, haven't been implementing functions for individual IDs.
//KDD  virtual size_t getEntityCoordinates(gid_t id, scalar_t *coords) const = 0;


  /*! \brief Returns whether a first adjacency combination is available.
   */
  virtual bool availAdjs(MeshEntityType source, MeshEntityType target);


  /*! \brief Returns the number of first adjacencies on this process.
   */
  virtual size_t getLocalNumAdjs(MeshEntityType source,
                                 MeshEntityType target) const = 0;


  /*! \brief Returns the number of first adjacencies for this entity
      \param source
   */
//KDD  For now, haven't been implementing functions for individual IDs.
//KDD  virtual size_t getEntityNumAdjs(MeshEntityType source, MeshEntityType target,
//KDD                                  gid_t id) const = 0;


  /*! \brief Sets pointers to this process' mesh first adjacencies.
      \param source
      \param offsets is an array of size getLocalNumOf() + 1.
         The first adjacency Ids for Ids[i] (returned in
         getIDsViewOf()) begin at adjacencyIds[offsets[i]].
          The last element of offsets
          is the size of the adjacencyIds array.
      \param adjacencyIds on return will point to the global first adjacency
         Ids for each entity.
       \return The number of ids in the adjacencyIds list.

      Zoltan2 does not copy your data.  The data pointed to by
      offsets and adjacencyIds
      must remain valid for the lifetime of this Adapter.
   */
//KDD Since the source objects are assumed to be gotten from getIDsViewOf(),
//KDD is the source MeshEntityType understood here?
//KDD What about the target?
  virtual size_t getAdjsView(MeshEntityType source, MeshEntityType target,
     const lno_t *&offsets, const gid_t *& adjacencyIds) const = 0;


  /*! \brief Copy all the first adjacency IDs for a specific entity into arrays
         provided by Zoltan2.
      \param source
      \param id specific entity.  Works only if id is local to this process.
      \param adjacencyIds array provided by Zoltan2.
      \return not sure what is needed, if anything.
  */
//KDD  For now, haven't been implementing functions for individual IDs.
//KDD  virtual size_t getEntityAdjs(MeshEntityType source, MeshEntityType target,
//KDD                               gid_t id, gid_t *adjacencyIds) const = 0;


  /*! \brief Returns whether a second adjacency combination is available.
   */
  virtual bool avail2ndAdjs(MeshEntityType sourcetarget, MeshEntityType through);


  /*! \brief Returns the number of second adjacencies on this process.
   *
   *  Some algorithms can partition a graph of mesh entities
   *
   *  Parameters will specify algorithm options:
   *   balance_entity_type==MeshEntityType, adjacency_through==MeshEntityType
   */
  virtual size_t getLocalNum2ndAdjs(MeshEntityType sourcetarget,
                                    MeshEntityType through) const = 0;


  /*! \brief Return number of second adjacencies to specific entity.
   *  \param sourcetarget
   *  \param id specific entity.  Works only if id is local to this process.
   *  \return number of second adjacendies to entity.
   */
//KDD  For now, haven't been implementing functions for individual IDs.
//KDD  virtual lno_t getEntityNum2ndAdjs(MeshEntityType sourcetarget,
//KDD     MeshEntityType through, gid_t id) const = 0;


  /*! \brief Sets pointers to this process' mesh second adjacencies.
      \param sourcetarget
      \param offsets is an array of size getLocalNumEntityIDs() + 1.
         The second adjacency Ids for Ids[i] (returned in
         getLocalEntityIDsView()) begin at adjacencyIds[offsets[i]].
          The last element of offsets
          is the size of the adjacencyIds array.
      \param adjacencyIds on return will point to the global second adjacency
         Ids for each entity.
       \return The number of ids in the adjacencyIds list.

      Zoltan2 does not copy your data.  The data pointed to by
      offsets and adjacencyIds
      must remain valid for the lifetime of this Adapter.
   */
// TODO:  Later may allow user to not implement second adjacencies and, if we want them,
// TODO:  we compute A^T A, where A is matrix of first adjacencies.
//KDD Since the source objects are assumed to be gotten from getIDsViewOf(),
//KDD is the sourcetarget MeshEntityType understood here?
//KDD What about the through MeshEntityType?
  virtual size_t get2ndAdjsView(MeshEntityType sourcetarget,
     MeshEntityType through, const lno_t *&offsets,
     const gid_t *& adjacencyIds) const = 0;


  /*! \brief Copy the second adjacency IDs for all adjacencies incident to a
         specific entity into array provided by Zoltan2.
      \param sourcetarget
      \param id specific entity.  Works only if id is local to this process.
      \param nborIds array provided by Zoltan2.
      \return not sure what is needed, if anything.
  */
//KDD  For now, haven't been implementing functions for individual IDs.
//KDD  virtual size_t getEntity2ndAdjs(MeshEntityType sourcetarget,
//KDD     MeshEntityType through, gid_t id, gid_t *nborIds) const = 0;


  /*! \brief Returns the number (0 or greater) of weights per second adjacency.
   */
  virtual int getNumWeightsPer2ndAdj(MeshEntityType sourcetarget,
                                     MeshEntityType through) const = 0;


  /*! \brief  Provide a pointer to the second adjacency weights, if any.

      \param weights is the list of weights of the given number for
           the second adjacencies returned in getLocal2ndAdjsView().
      \param stride The k'th weight is located at weights[stride*k]
      \param idx ranges from zero to one less than
                   getNumWeightsPer2ndAdj().
       \return The number of weights listed, which should be the same
               as the number of adjacencies in getLocal2ndAdjsView().

      Zoltan2 does not copy your data.  The data pointed to by weights
      must remain valid for the lifetime of this Adapter.
   */
//KDD Since the source objects are assumed to be gotten from getIDsViewOf(),
//KDD is the sourcetarget MeshEntityType understood here?
//KDD What about the through MeshEntityType?
  virtual void get2ndAdjWeightsView(MeshEntityType sourcetarget,
     MeshEntityType through, const scalar_t *&weights, int &stride,
     int idx) const = 0;


  /*! \brief  Copy the second adjacency weights incident to a specific entity
          id into array provided by Zoltan2.
      \param sourcetarget
      \param id specific entity id.  Works only if id is local to this process.
      \param ewgts array provided by Zoltan2.
      \return not sure what is needed, if anything.
  */
//KDD  For now, haven't been implementing functions for individual IDs.
//KDD  virtual size_t getEntity2ndAdjWeights(MeshEntityType sourcetarget,
//KDD     MeshEntityType through, gid_t id, scalar_t *ewgts);


  /*! \brief Apply a partitioning problem solution to an input.
   *
   *  This is not a required part of the MeshAdapter interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *           of global Ids to processes that matches that of the input
   *           data that instantiated this Adapter.
   *  \param out On return this should point to a newly created object
   *            with the specified partitioning.
   *  \param solution  The Solution object created by a Problem should
   *      be supplied as the third argument.  It must have been templated
   *      on user data that has the same global ID distribution as this
   *      user data.
   *  \return   Returns the number of local Ids in the new partitioning.
   */

  template <typename Adapter>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<Adapter> &solution) const
  {
    return 0;
  }


  ////////////////////////////////////////////////////////////////////////////
  // Implementations of base-class methods

  /*! \brief Returns the entity to be partitioned, ordered, colored, etc.
   */
  inline enum MeshEntityType getPrimaryEntityType() const {
    return this->primaryEntityType;
  }

  /*! \brief Sets the primary entity type.  Called by algorithm based on
   *  parameter value in parameter list from application.
   *  Also sets to adjacencyEntityType to something reasonable:  opposite of
   *  primaryEntityType.
   */
  void setPrimaryEntityType(string typestr) {
    if (typestr == "region")
      this->primaryEntityType = MESH_REGION;
    else if (typestr == "face")
      this->primaryEntityType = MESH_FACE;
    else if (typestr == "edge")
      this->primaryEntityType = MESH_EDGE;
    else if (typestr == "vertex")
      this->primaryEntityType = MESH_VERTEX;
    else {
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  Invalid MeshEntityType " << typestr << std::endl;
      emsg << "Valid values: region  face  edge  vertex" << std::endl;
      throw std::runtime_error(emsg.str());
    }
  }

  /*! \brief Returns the entity that describes adjacencies between the
   *  entities to be partitioned, ordered, colored, etc.
   *  That is, two primaryEntityType that share an adjacencyEntityType are
   *  adjacent.
   *  KDD:  Is Adjacency a poorly chosen name here?  Is it overloaded?
   */
  inline enum MeshEntityType getAdjacencyEntityType() const {
    return this->adjacencyEntityType;
  }

  /*! \brief Sets the adjacency entity type.  Called by algorithm based on
   *  parameter value in parameter list from application.
   *  Also sets to primaryEntityType to something reasonable:  opposite of
   *  adjacencyEntityType.
   *  KDD:  Is Adjacency a poorly chosen name here?  Is it overloaded?
   */
  void setAdjacencyEntityType(string typestr) {
    if (typestr == "region")
      this->adjacencyEntityType = MESH_REGION;
    else if (typestr == "face")
      this->adjacencyEntityType = MESH_FACE;
    else if (typestr == "edge")
      this->adjacencyEntityType = MESH_EDGE;
    else if (typestr == "vertex")
      this->adjacencyEntityType = MESH_VERTEX;
    else {
      std::ostringstream emsg;
      emsg << __FILE__ << "," << __LINE__
           << " error:  Invalid MeshEntityType " << typestr << std::endl;
      emsg << "Valid values: region  face  edge  vertex" << std::endl;
      throw std::runtime_error(emsg.str());
    }
  }

  // Functions from the BaseAdapter interface
  size_t getLocalNum() const {
    return getLocalNumOf(getPrimaryEntityType());
  }

  size_t getIDsView(const gid_t *&Ids) const {
    return getIDsViewOf(getPrimaryEntityType());
  }

  int getNumWeightsPerID() const {
    return getNumWeightsPerOf(getPrimaryEntityType());
  }

  void getWeightsView(const scalar_t *&wgt, int &stride, int idx = 0) const {
    getWeightsViewOf(getPrimaryEntityType());
  }

};

}  //namespace Zoltan2

#endif

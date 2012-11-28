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

/*! \file Zoltan2_MeshInput.hpp
    \brief Defines the MeshInput interface.
*/


#ifndef _ZOLTAN2_MESHINPUT_HPP_
#define _ZOLTAN2_MESHINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <string>

namespace Zoltan2 {

/*!  \brief MeshInput defines the interface for mesh input adapters.

    InputAdapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t entity and adjacency weights 
    \li \c lno_t    local indices and local counts
    \li \c gno_t    global indices and global counts
    \li \c gid_t    application global Ids
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel

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
  class MeshInput : public InputAdapter {
private:

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  enum InputAdapterType inputAdapterType() const {return MeshAdapterType;}

  /*! \brief Destructor
   */
  virtual ~MeshInput() {};

  /*! \brief Returns the number of mesh entities on this process.
   */

  virtual size_t getLocalNumEntityIDs(entityType etype) const = 0;


  /*! \brief Return the number of weights per entity.
   *   \return the count of weights, zero or more per entity.
   *   If the number of weights is zero, then we assume that the entities
   *   are equally weighted.
   */
  virtual int getNumWeightsPerEntityID(entityType etype) const = 0;


  /*! \brief Provide a pointer to this process' identifiers.

      \param Ids will on return point to the list of the global Ids for this
        process.

       \return The number of ids in the Ids list.
    
      Some algorithms can partition a simple list of weighted identifiers
        with no geometry or topology provided.
  */

  virtual size_t getLocalEntityIDsView(entityType etype,
				       gid_t const *&Ids) const = 0;


  /*! \brief Provide a pointer to one of the number of this process'
                optional entity weights.

      \param weights on return will contain a list of the weights for the
               number specified.  If weights for
	   this number are to be uniform for all entities in the
	   global problem, the \c weights should be a NULL pointer.

      \param stride on return will indicate the stride of the weights list.

       The k'th weight is located at weights[stride*k].

      \param number is a value ranging from zero to one less than
                   getNumWeightsPerEntityID()

       \return The number of values in the weights list.  This may be greater
          than the number of entities, because the stride may be greater
	  than one.

      Zoltan2 does not copy your data.  The data pointed to by weights
      must remain valid for the lifetime of this InputAdapter.
  */

  virtual size_t getLocalEntityIDsWeightsView(entityType etype,
     const scalar_t *&weights, int &stride, int number) const = 0;

  /*! \brief Copy the weights for a specific entity into array provided by
                Zoltan2.

      \param id specific entity.  Works only if id is local to this process.

      \param wgts array provided by Zoltan2.
		wgts has size getNumWeightsPerID().

      \return not sure what is needed, if anything.
  */

  virtual size_t getWeights(gid_t id, scalar_t *wgts) const = 0;


  /*! \brief Return dimension of the entity coordinates, if any.
   *
   *  Some algorithms can partition mesh entities using geometric coordinate
   *    information
   *
   *  Some algorithms can use geometric entity coordinate
   *    information if it is present.
   */
  virtual int getEntityCoordinateDimension(entityType etype) const = 0;


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
      must remain valid for the lifetime of this InputAdapter.
  */

  virtual size_t getLocalEntityCoordinatesView(entityType etype,
    const scalar_t *&coords, int &stride, int coordDim) const = 0;


  /*! \brief Copy the coordinates for a specific ID into array
                provided by Zoltan2.

      \param id specific ID.  Works only if id is local to this	process.

      \param coords array provided by Zoltan2.

      \return not sure what is needed, if anything.
  */
  virtual size_t getCoordinates(gid_t id, scalar_t *coords) const = 0;


  /*! \brief Returns whether a first adjacency combination is available.
   */
  virtual bool availAdjacency(entityType source, entityType target);


  /*! \brief Returns the number of first adjacencies on this process.
   */
  virtual size_t getLocalNumAdjacencyIDs(entityType source,
					 entityType target) const = 0;


  /*! \brief Sets pointers to this process' mesh first adjacencies.
      \param source
      \param offsets is an array of size getLocalNumEntityIDs() + 1.
         The first adjacency Ids for Ids[i] (returned in
	 getLocalEntityIDsView()) begin at adjacencyIds[offsets[i]].
          The last element of offsets
          is the size of the adjacencyIds array.
      \param adjacencyIds on return will point to the global first adjacency
         Ids for each entity.
       \return The number of ids in the adjacencyIds list.

      Zoltan2 does not copy your data.  The data pointed to by
      offsets and adjacencyIds
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getLocalAdjacencyIDsView(entityType source, entityType target,
     const lno_t *&offsets, const gid_t *& adjacencyIds) const = 0;


  /*! \brief Returns whether a second adjacency combination is available.
   */
  virtual bool avail2ndAdjacency(entityType sourcetarget, entityType through);


  /*! \brief Returns the number of second adjacencies on this process.
   *
   *  Some algorithms can partition a graph of mesh entities
   *
   *  Parameters will specify algorithm options:
   *   balance_entity_type==entityType, adjacency_through==entityType
   */
  virtual size_t getLocalNum2ndAdjacencyIDs(entityType sourcetarget,
					    entityType through) const = 0;

  /*! \brief Return number of second adjacencies to specific entity.
   *  \param sourcetarget
   *  \param id specific entity.  Works only if id is local to this process.
   *  \return number of second adjacendies to entity.
   */
  virtual lno_t getEntityNum2ndAdjacencyIDs(entityType sourcetarget,
     entityType through, gid_t id) const = 0;

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
      must remain valid for the lifetime of this InputAdapter.
   */

// TODO:  Later may allow user to not implement second adjacencies and, if we want them,
// TODO:  we compute A^T A, where A is matrix of first adjacencies.

  virtual size_t getLocal2ndAdjacencyIDsView(entityType sourcetarget,
     entityType through, const lno_t *&offsets,
     const gid_t *& adjacencyIds) const = 0;


  /*! \brief Returns the number (0 or greater) of weights per second adjacency.
   */
  virtual int getNumWeightsPer2ndAdjacencyID(entityType sourcetarger,
					     entityType through) const = 0;


  /*! \brief  Provide a pointer to the second adjacency weights, if any.

      \param weights is the list of weights of the given number for
           the second adjacencies returned in getLocal2ndAdjacencyIDsView().
      \param stride The k'th weight is located at weights[stride*k]
      \param number ranges from zero to one less than 
                   getNumWeightsPer2ndAdjacencyID().
       \return The number of weights listed, which should be the same
               as the number of adjacencies in getLocal2ndAdjacencyIDsView().

      Zoltan2 does not copy your data.  The data pointed to by weights
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getLocal2ndAdjacencyIDWeightsView(entityType sourcetarget,
     entityType through, const scalar_t *&weights, int &stride,
     int number) const = 0;

  /*! \brief Apply a partitioning problem solution to an input.  
   *
   *  This is not a required part of the MeshInput interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *           of global Ids to processes that matches that of the input
   *           data that instantiated this InputAdapter.
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
  
};
  
}  //namespace Zoltan2
  
#endif

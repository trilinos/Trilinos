// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GraphInput.hpp
    \brief Defines the GraphInput interface.
*/


#ifndef _ZOLTAN2_GRAPHINPUT_HPP_
#define _ZOLTAN2_GRAPHINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

/*!  \brief GraphInput defines the interface for graph input adapters.

    InputAdapter objects provide access for Zoltan2 to the user's data.
    Many built-in adapters are already defined for common data structures,
    such as Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t vertex and edge weights 
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

\todo add weights
\todo how do we want to handle graphs that have coordinates for the
        vertices? 
*/

template <typename User>
  class GraphInput : public InputAdapter {
private:

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
#endif

  enum InputAdapterType inputAdapterType() const {return GraphAdapterType;}

  /*! \brief Destructor
   */
  virtual ~GraphInput() {};

  /*! \brief Returns the number vertices on this process.
   */
  virtual size_t getLocalNumVertices() const = 0;

  /*! \brief Returns the global number vertices.
   */
  virtual global_size_t getGlobalNumVertices() const = 0;

  /*! \brief Returns the number edges on this process.
   */
  virtual size_t getLocalNumEdges() const = 0;

  /*! \brief Returns the global number edges.
   */
  virtual global_size_t getGlobalNumEdges() const = 0;

#if 0
  /*! Returns the dimension (0 or greater) of vertex weights.
   */
  virtual int getVertexWeightDim() const = 0;

  /*! Returns the dimension (0 or greater) of edge weights.
   */
  virtual int getEdgeWeightDim() const = 0;

  /*! Returns the dimension (0 to 3) of vertex coordinates.
   */
  virtual int getCoordinateDim() const = 0;
#endif

  /*! \brief Sets pointers to this process' graph entries.
      \param vertexIds will on return a pointer to vertex global Ids
      \param offsets is an array of size numVertices + 1.  
         The edge Ids for vertexId[i] begin at edgeIds[offsets[i]].  
          The last element of offsets
          is the size of the edgeIds array.
      \param edgeIds on return will point to the global edge Ids for
         for each vertex.
       \return The number of ids in the vertexIds list.
   */

  virtual size_t getVertexListView(const gid_t *&vertexIds, 
    const lno_t *&offsets, const gid_t *& edgeIds) const = 0; 

  /*! \brief Apply a partitioning problem solution to an input.  
   *
   *  This is not a required part of the GraphInput interface. However
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

  template <typename User2>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<User2> &solution)
  {
    return 0;
  }
  
};
  
}  //namespace Zoltan2
  
#endif

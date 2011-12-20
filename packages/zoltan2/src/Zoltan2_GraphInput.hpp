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
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

/*!  \brief GraphInput defines the interface for input adapters for
            graphs that may have vertex and edge weight. 
*/

template <typename User>
  class GraphInput : public InputAdapter {
private:

public:

  /*! \brief Data type for weights. */
  typedef typename InputTraits<User>::scalar_t scalar_t;

  /*! \brief Data type for Zoltan2's indices and local counts. */
  typedef typename InputTraits<User>::lno_t    lno_t;

  /*! \brief Data type for Zoltan2's global identifiers and counts. */
  typedef typename InputTraits<User>::gno_t    gno_t;

  /*! \brief The application's global identifier data type. */
  typedef typename InputTraits<User>::gid_t    gid_t;

  /*! \brief The Kokkos node type. */
  typedef typename InputTraits<User>::node_t   node_t;

  enum InputAdapterType inputAdapterType() const {return GraphAdapterType;}

  typedef User user_t;

  /*! \brief Pure virtual destructor
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

  /*! \brief Apply the solution to a partitioning problem to an input.  
   *
   *  This is not a required part of the GraphInput interface.  However
   *  if the PartitioningProblem::redistribute() method is called, it 
   *  will use this method to redistribute the data.  If the user has 
   *  no intention of calling redistribute(), then it is not necessary to 
   *  define applyPartitioningSolution in the InputAdapter.
   *
   *  \param in  An input object with a structure and assignment of
   *           of global Ids to processes that matches that of the input
   *           data that instantiated this InputAdapter.
   *  \param out On return this should point to a newly created object 
   *            with the specified partitioning.
   *  \return   Returns the number of local Ids in the new partitioning.
   *
   * TODO - A solution needs to be more than a list of partitions, but
   *   also how those partitions map to processes.  For now it's
   *   process "p" gets part "p".
   */

  size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<User> &solution)
  {
    return 0;
  }
  
};
  
}  //namespace Zoltan2
  
#endif

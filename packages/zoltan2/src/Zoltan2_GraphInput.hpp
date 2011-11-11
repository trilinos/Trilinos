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

#include <string>
#include <Zoltan2_InputAdapter.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphInput
    \brief GraphInput defines the interface for input adapters for
            graphs that may have vertex and edge weight. 

    The Graph accessor methods defined here mimic those of 
    Tpetra::CrsGraph and Tpetra::Map.

    LID: the type for the application's local Ids
    GID: the type for the application's global Ids
    LNO: the integral type that Zoltan2 will use for local counters.
    GNO: the integral type that Zoltan2 will use for the global 
      counts and identifiers.  It needs to be large enough for the
      problem's number of objects.
*/

template <typename User>
  class GraphInput : public InputAdapter {
private:

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;

  // adapterType == GraphAdapterType
  // Function must return one of Zoltan2's enumerated types in InputAdapter
  // User should not rewrite this function.
  enum InputAdapterType inputAdapterType() {return GraphAdapterType;}

  /*! Pure virtual destructor
   */
  virtual ~GraphInput() {};

  /*! Returns the number vertices on this process.
   */
  virtual size_t getLocalNumVertices() const = 0;

  /*! Returns the global number vertices.
   */
  virtual global_size_t getGlobalNumVertices() const = 0;

  /*! Returns the number edges on this process.
   */
  virtual size_t getLocalNumEdges() const = 0;

  /*! Returns the global number edges.
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

  /*! Sets pointers to this process' graph entries.
      \param vertexIds will on return a pointer to vertex global Ids
      \param localIds can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  If localIds are omitted and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local Id order.
      \param offsets is an array of size numVertices + 1.  
         The edge Ids for vertexId[i] begin at edgeIds[offsets[i]].  
          The last element of offsets
          is the size of the edgeIds array.
      \param edgeIds on return will point to the global edge Ids for
         for each vertex.
       \return The number of ids in the vertexIds list.
   */

  virtual size_t getVertexListView(const gid_t *&vertexIds, 
    const lid_t *&localIds, 
    const lno_t *&offsets, const gid_t *& edgeIds) const = 0; 

  /*! Apply the solution to a partitioning problem to an input.  
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
   *  \param numIds  The number of ids in the gid and partition lists.
   *  \param numParts  The global number of partitions.  Partitions are 
   *     numbered from 0 through numParts-1. 
   *  \param gid     A list of object global Ids.
   *  \param lid     A corresponding list of object local Ids, if the
   *      InputAdapter had supplied local Ids.
   *  \param partition  A corresponding list of partitions.  gid[i]
   *            has been assigned to partition[i].
   *  \return   Returns the number of local Ids in the new partitioning.
   *
   * TODO - A solution needs to be more than a list of partitions, but
   *   also how those partitions map to processes.  For now it's
   *   process "p" gets part "p".
   */

  size_t applyPartitioningSolution(const User &in, User *&out,
    lno_t numIds, lno_t numParts, const gid_t *gid, 
    const lid_t *lid, const lno_t *partition)
  {
    return 0;
  }
  
};
  
}  //namespace Zoltan2
  
#endif

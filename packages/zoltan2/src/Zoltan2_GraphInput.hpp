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
  class GraphInput : public InputAdapter<User> {
private:

public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

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

  /*! Returns the dimension (0 or greater) of vertex weights.
   */
  virtual int getVertexWeightDim() const = 0;

  /*! Returns the dimension (0 or greater) of edge weights.
   */
  virtual int getEdgeWeightDim() const = 0;

  /*! Returns the dimension (0 to 3) of vertex coordinates.
   */
  virtual int getCoordinateDim() const = 0;

  /*! Returns list of this process' vertex Ids and their weights.
      \param Ids will on return hold the list of the global Ids for 
        each vertex on this process.
      \param localIds can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  If localIds are omitted and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local Id order.
      \param xyz can, optionally, on return hold coordinates for the
        vertices in order by vertex by coordinate.
      \param wgts can, optionallly, on return hold list of the 
         weight or weights associated with each vertex.  Weights 
         are listed by vertex by weight component.  
   */

  virtual void getVertexListCopy(std::vector<gid_t> &Ids, 
    std::vector<lid_t> &localIds,
    std::vector<scalar_t> &xyz,
    std::vector<scalar_t> &wgts) const = 0;

  /*! Sets pointers to this process' vertex Ids and their weights.
      If this optional call is defined in the adapter, it can save a memory
      copy of application data.
      \param Ids will on return point to the list of the global Ids for 
        each vertex on this process.
      \param localIds can, optionally, on return point to a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list. If localIds is NULL and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local ID order.
      \param xyz will on return point to a list coordinates for
         each vertex in the Ids list.  Coordinates are listed by 
         vertex by component.  
      \param wgts will on return point to a list of the weight or weights 
         associated with each vertex in the Ids list.  Weights are listed by 
         vertex by weight component.  
       \return The number of ids in the Ids list.
   */

  lno_t getVertexListView(gid_t *&Ids, lid_t *&localIds,
     scalar_t *&xyz, scalar_t *&wgts)
  {
    Ids = NULL;
    localIds = NULL;
    xyz = NULL;
    wgts = NULL;
    return 0;
  }

  /*! Return a list of the edge Ids of the input vertex global ID
      \param ID  global ID for a vertex on this process
      \param localID  app's local ID, if any, associated with this vertex
      \param edgeID on return will contain the list of edge Ids
      \param wgts on return contains the weights, if any associated with the
         edges. Weights are listed by edge by weight component.
   */
  virtual void getVertexEdgeCopy(gid_t ID, lid_t localID,
    std::vector<gid_t> &edgeID, std::vector<scalar_t> &wgts) const = 0;

  /*! Obtain a read-only view, if possible, of the edge Ids of the 
      input vertex.
      \param ID  global ID for a vertex on this process
      \param localID  if input adapter supplied local Ids, this
         is that localID
      \param edgeID on return will point to the list of edge global Ids.
      \param wgts on return points to the weights, if any associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeID list.
   */
  lno_t getVertexEdgeView(gid_t ID, lid_t localID, gid_t *&edgeID,
    scalar_t * &wgts) const
  {
    edgeID = NULL;
    return 0;
  }
};
  
  
}  //namespace Zoltan2
  
#endif

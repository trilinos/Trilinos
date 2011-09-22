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

    TODO: It may be necessary to put a migration interface at this level.

    Scalar: This data type is used for weights and coordinates.
    AppLID: the type for the application's local IDs
    AppGID: the type for the application's global IDs
    LNO: the integral type that Zoltan2 will use for local counters.
    GNO: the integral type that Zoltan2 will use for the global 
      counts and identifiers.  It needs to be large enough for the
      problem's number of objects.
*/

template<typename Scalar, typename AppLID, typename AppGID,
  typename LNO=int, typename GNO=AppGID>
  class GraphInput : public InputAdapter {
private:

public:

  // TODO - what to do with destructor
  //virtual ~GraphInput();

  /*! return a name that identifies the concrete adapter
   */

  virtual std::string inputAdapterName() const = 0;

  /*! Returns the number vertices on this process.
   */
  virtual LNO getLocalNumVertices() const = 0;

  /*! Returns the number edges on this process.
   */
  virtual LNO getLocalNumEdges() const = 0;

  /*! Returns the dimension (0 or greater) of vertex weights.
   */
  virtual int getVertexWeightDim() const = 0;

  /*! Returns the dimension (0 or greater) of edge weights.
   */
  virtual int getEdgeWeightDim() const = 0;

  /*! Returns the dimension (0 to 3) of vertex coordinates.
   */
  virtual int getCoordinateDim() const = 0;

  /*! Returns list of this process' vertex IDs and their weights.
      \param IDs will on return hold the list of the global IDs for 
        each vertex on this process.
      \param localIDs can, optionally, on return hold a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.
      \param xyz can, optionally, on return hold coordinates for the
        vertices in order by vertex by coordinate.
      \param wgts can, optionallly, on return hold list of the 
         weight or weights associated with each vertex.  Weights 
         are listed by vertex by weight component.  
   */

  virtual void getVertexListCopy(std::vector<AppGID> &IDs, 
    std::vector<AppLID> &localIDs,
    std::vector<Scalar> &xyz,
    std::vector<Scalar> &wgts) const = 0;

  /*! Sets pointers to this process' vertex IDs and their weights.
      If this optional call is defined in the adapter, it can save a memory
      copy of application data.
      \param IDs will on return point to the list of the global IDs for 
        each vertex on this process.
      \param localIDs can, optionally, on return point to a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list.  It's a copy since the common case is
        that they are not stored by the application but rather taken from
        indices or pointers.
      \param xyz will on return point to a list coordinates for
         each vertex in the IDs list.  Coordinates are listed by 
         vertex by component.  
      \param wgts will on return point to a list of the weight or weights 
         associated with each vertex in the IDs list.  Weights are listed by 
         vertex by weight component.  
       \return The number of ids in the IDs list.
   */

  LNO getVertexListView(AppGID const *&IDs, 
    std::vector<AppLID> &localIDs, Scalar const *&xyz, Scalar const *&wgts)
  {
    IDs = NULL;
    localIDs = NULL;
    xyz = NULL;
    wgts = NULL;
    return 0;
  }

  /*! Return a list of the edge IDs of the input vertex global ID
      \param ID  global ID for a vertex on this process
      \param localID  app's local ID, if any, associated with this vertex
      \param edgeID on return will contain the list of edge IDs
      \param wgts on return contains the weights, if any associated with the
         edges. Weights are listed by edge by weight component.
   */
  virtual void getVertexEdgeCopy(AppGID ID, AppLID localID,
    std::vector<AppGID> &edgeID, std::vector<Scalar> &wgts) const = 0;

  /*! Obtain a read-only view, if possible, of the edge IDs of the 
      input vertex.
      \param ID  global ID for a vertex on this process
      \param localID  app's local ID, if any, associated with this vertex
      \param edgeID on return will point to the list of edge global IDs.
      \param wgts on return points to the weights, if any associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeID list.
   */
  LNO getVertexEdgeView(AppGID ID, AppLID localID, AppGID const *&edgeID,
    Scalar const * &wgts) const
  {
    edgeID = NULL;
    return 0;
  }

  // TODO - do we want the application to specify properties of the
  //  graph?  isDirected(), isBipartite(), ...
};
  
  
}  //namespace Zoltan2
  
#endif

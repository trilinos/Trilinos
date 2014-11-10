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

#include "Tpetra_DefaultPlatform.hpp"

namespace Zoltan2 {
  
  /*!  \brief Enumerate entity types for meshes:  Regions, Faces, Edges, or 
   *                                              Vertices
   */
  
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
    \li \c zgid_t    application global Ids
    \li \c node_t is a sub class of KokkosClassic::StandardNodeMemoryModel

    See IdentifierTraits to understand why the user's global ID type
    (\c zgid_t) may differ from that used by Zoltan2 (\c gno_t).

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
  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::zgid_t    zgid_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;
  typedef User userCoord_t;
  typedef int LO;
  typedef int GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  typedef Tpetra::Map<LO, GO, Node>      map_type;
  typedef Tpetra::CrsGraph<LO, GO, Node> sparse_graph_type;
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
  
  /*! \brief Returns the number of mesh entities on this process.
   */
  virtual size_t getLocalNumOf(MeshEntityType etype) const = 0;
  
  
  /*! \brief Provide a pointer to this process' identifiers.
      \param Ids will on return point to the list of the global Ids for this
       process.
  */
  virtual void getIDsViewOf(MeshEntityType etype,
                            zgid_t const *&Ids) const = 0;


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
    Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
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
    Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
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
     const lno_t *&offsets, const zgid_t *& adjacencyIds) const 
  {
    offsets = NULL;
    adjacencyIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
  }


  /*! \brief Returns whether a second adjacency combination is available.
   */
  virtual bool avail2ndAdjs(MeshEntityType sourcetarget, 
			    MeshEntityType through) const {
    if (!availAdjs(sourcetarget, through))
      return false;
    else {
      return false;

      using Tpetra::DefaultPlatform;
      using Tpetra::global_size_t;
      using Teuchos::Array;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::rcp;

      // Get the default communicator and Kokkos Node instance
      RCP<const Comm<int> > comm =
	DefaultPlatform::getDefaultPlatform ().getComm ();
      RCP<Node> node = DefaultPlatform::getDefaultPlatform ().getNode ();

      lno_t const *offsets=NULL;
      zgid_t const *adjacencyIds=NULL;
      getAdjsView(sourcetarget, through, offsets, adjacencyIds);

      zgid_t const *Ids=NULL;
      getIDsViewOf(MESH_VERTEX, Ids);

      int LocalNumIDs = getLocalNumIDs();
      int LocalNumAdjs = getLocalNumAdjs(sourcetarget, through);

      Array<GO> adjsGIDs;
      Array<GO> GIDs;
      RCP<const map_type> adjsMapG;
      RCP<const map_type> MapG;

      // Count owned nodes
      int adjsNodes = getLocalNumOf(MESH_VERTEX);

      // Build a list of the ADJS global ids...
      adjsGIDs.resize (adjsNodes);
      for (int i = 0; i < adjsNodes; ++i) {
	adjsGIDs[i] = as<int> (Ids[i]);
      }

      getIDsView(Ids);

      // Build a list of the global ids...
      GIDs.resize (LocalNumIDs);
      for (int i = 0; i < LocalNumIDs; ++i) {
	GIDs[i] = as<int> (Ids[i]);
      }

      //Generate adjs Map for nodes.
      adjsMapG = rcp (new map_type (-1, adjsGIDs (), 0, comm, node));

      //Generate Map for elements.
      MapG = rcp (new map_type (-1, GIDs (), 0, comm, node));

      RCP<sparse_graph_type> adjsGraph;
      RCP<sparse_graph_type> adjsGraphTranspose;

      // Construct Tpetra::CrsGraph objects.
      adjsGraph = rcp (new sparse_graph_type (adjsMapG, 0));
      adjsGraphTranspose = rcp (new sparse_graph_type (MapG, 0));

      for (int i = 0; i < LocalNumIDs; ++i) {

	int Row = Ids[i];
	//globalRow for Tpetra Graph
	global_size_t globalRowT = as<global_size_t> (Row);
	int globalRow = as<int> (Row);
	//create ArrayView globalRow object for Tpetra
	ArrayView<int> globalRowAV = Teuchos::arrayView (&globalRow, 1);

	int NumAdjs;
	if (i + 1 < LocalNumIDs) {
	  NumAdjs = offsets[i+1];
	} else {
	  NumAdjs = LocalNumAdjs;
	}

	for (int j = offsets[i]; j < NumAdjs; ++j) {
	  int Col = adjacencyIds[j];
	  //globalCol for Tpetra Graph
	  global_size_t globalColT = as<global_size_t> (Col);
	  int globalCol = as<int> (Col);
	  //create ArrayView globalCol object for Tpetra
	  ArrayView<int> globalColAV = Teuchos::arrayView (&globalCol, 1);

	  //Update Tpetra adjs Graph
	  adjsGraph->insertGlobalIndices (globalRowT, globalColAV);
	  adjsGraphTranspose->insertGlobalIndices (globalColT, globalColAV);
	}// *** col loop ***
      }// *** row loop ***

      // Fill-complete adjs Graph.
      adjsGraph->fillComplete ();

      //Fill-complete adjs Graph Transpose.
      adjsGraphTranspose->fillComplete ();
    }
  }


  /*! \brief Returns the number of second adjacencies on this process.
   *
   *  Parameters will specify algorithm options:
   *   balance_entity_type==MeshEntityType, adjacency_through==MeshEntityType
   */
  virtual size_t getLocalNum2ndAdjs(MeshEntityType sourcetarget,
                                    MeshEntityType through) const {
    if (!avail2ndAdjs(sourcetarget, through))
      return 0;
    else {
      return nadj_;
    }
  }

  /*! \brief Sets pointers to this process' mesh second adjacencies.
      \param sourcetarget
      \param offsets is an array of size getLocalNumOf() + 1.
         The second adjacency Ids for Ids[i] (returned in
         getIDsViewOf()) begin at adjacencyIds[offsets[i]].
          The last element of offsets
          is the size of the adjacencyIds array.
      \param adjacencyIds on return will point to the global second adjacency
         Ids for each entity.
   */
// TODO:  Later may allow user to not implement second adjacencies and,
// TODO:  if we want them, we compute A^T A, where A is matrix of first adjacencies.
  virtual void get2ndAdjsView(MeshEntityType sourcetarget,
                              MeshEntityType through,
                              const lno_t *&offsets,
                              const zgid_t *&adjacencyIds) const
  {
    if (!avail2ndAdjs(sourcetarget, through)) {
      offsets = NULL;
      adjacencyIds = NULL;
      Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
    } else {
      offsets = start_;
      adjacencyIds =  adj_;
    }
  }


  /*! \brief Returns the number (0 or greater) of weights per second adjacency.
   */
  virtual int getNumWeightsPer2ndAdj(MeshEntityType sourcetarget,
                                     MeshEntityType through) const { return 0;}


  /*! \brief  Provide a pointer to the second adjacency weights, if any.

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
    Z2_THROW_NOT_IMPLEMENTED_IN_ADAPTER
  }

//KDD What if we wanted to provide weights with respect to first adjacencies?
//KDD Should we add functions for that?
//VJL Yes.

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
   *  KDD:  Is Adjacency a poorly chosen name here?  Is it overloaded?
   *  VJL:  Maybe
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
   *  KDD:  Is Adjacency a poorly chosen name here?  Is it overloaded?
   *  VJL:  Maybe
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
  
  void getIDsView(const zgid_t *&Ids) const {
    getIDsViewOf(getPrimaryEntityType(), Ids);
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
  lno_t *start_;
  zgid_t *adj_;
  size_t nadj_;
};
  
}  //namespace Zoltan2

#endif

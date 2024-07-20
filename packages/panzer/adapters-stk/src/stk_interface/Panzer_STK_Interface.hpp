// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_Interface_hpp__
#define __Panzer_STK_Interface_hpp__

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>

#include "Kokkos_Core.hpp"

#include <Shards_CellTopology.hpp>
#include <Shards_CellTopologyData.h>

#include <PanzerAdaptersSTK_config.hpp>
#include <Kokkos_ViewFactory.hpp>

#include <unordered_map>

#ifdef PANZER_HAVE_IOSS
#include <stk_io/StkMeshIoBroker.hpp>
#endif

#ifdef PANZER_HAVE_PERCEPT
namespace percept {
  class PerceptMesh;
  class URP_Heterogeneous_3D;
}
#endif

namespace panzer_stk {

class PeriodicBC_MatcherBase;

/** Pure virtial base class that builds a basic element. To be
  * overidden with several types of elements.
  */
class ElementDescriptor {
public:
   ElementDescriptor(stk::mesh::EntityId gid,const std::vector<stk::mesh::EntityId> & nodes);
   virtual ~ElementDescriptor();

   stk::mesh::EntityId getGID() const { return gid_; }
   const std::vector<stk::mesh::EntityId> & getNodes() const { return nodes_; }
protected:
   stk::mesh::EntityId gid_;
   std::vector<stk::mesh::EntityId> nodes_;

   ElementDescriptor();
};

/** Constructor function for building the element descriptors.
  */
Teuchos::RCP<ElementDescriptor>
buildElementDescriptor(stk::mesh::EntityId elmtId,std::vector<stk::mesh::EntityId> & nodes);

class STK_Interface {
public:
   typedef double ProcIdData; // ECC: Not sure why?
   typedef stk::mesh::Field<double> SolutionFieldType;
   typedef stk::mesh::Field<double> VectorFieldType;
   typedef stk::mesh::Field<ProcIdData> ProcIdFieldType;

   // some simple exception classes
   struct ElementBlockException : public std::logic_error
   { ElementBlockException(const std::string & what) : std::logic_error(what) {} };

   struct SidesetException : public std::logic_error
   { SidesetException(const std::string & what) : std::logic_error(what) {} };

   struct EdgeBlockException : public std::logic_error
   { EdgeBlockException(const std::string & what) : std::logic_error(what) {} };

   struct FaceBlockException : public std::logic_error
   { FaceBlockException(const std::string & what) : std::logic_error(what) {} };

   STK_Interface();

   /** Default constructor
     */
   STK_Interface(unsigned dim);

   STK_Interface(Teuchos::RCP<stk::mesh::MetaData> metaData);

   // functions called before initialize
   //////////////////////////////////////////

   /** Add an element block with a string name
     */
   void addElementBlock(const std::string & name,const CellTopologyData * ctData);

   /** Add an edge block with a string name
     */
//   void addEdgeBlock(const std::string & name,const CellTopologyData * ctData);
   void addEdgeBlock(const std::string & elemBlockName,
                     const std::string & edgeBlockName,
                     const stk::topology & topology);
   void addEdgeBlock(const std::string & elemBlockName,
                     const std::string & edgeBlockName,
                     const CellTopologyData * ctData);

   /** Add a face block with a string name
     */
   void addFaceBlock(const std::string & elemBlockName,
                     const std::string & faceBlockName,
                     const stk::topology & topology);
   void addFaceBlock(const std::string & elemBlockName,
                     const std::string & faceBlockName,
                     const CellTopologyData * ctData);

   /** Add a side set with a string name
     */
   void addSideset(const std::string & name,const CellTopologyData * ctData);

   /** Add a node set with a string name
     */
   void addNodeset(const std::string & name);

   /** Add a solution field
     */
   void addSolutionField(const std::string & fieldName,const std::string & blockId);

   /** Add a solution field
     */
   void addCellField(const std::string & fieldName,const std::string & blockId);

   /** Add an edge field
     */
   void addEdgeField(const std::string & fieldName,const std::string & blockId);

   /** Add a face field
     */
   void addFaceField(const std::string & fieldName,const std::string & blockId);

   /** Add a solution field for coordinates with a particular prefix, force it
     * to be outputed as a to be mesh displacement field. This
     * is really only relevant for I/O and how the field is stored internally in the mesh.
     *
     * \param[in] blockId Element block to use
     * \param[in] coordPrefix Prefix for solution coordinate field (assumes using "X","Y","Z" as postfix)
     * \param[in] dispPrefix Prefix for displacment (output) field
     */
   void addMeshCoordFields(const std::string & blockId,
                           const std::vector<std::string> & coordField,
                           const std::string & dispPrefix);

   /** Add a vector of strings to the Information Records block.  Each string
     * will be it's own record.  The info records will be deduped before they
     * are added to IOSS.
     */
   void addInformationRecords(const std::vector<std::string> & info_records);

   //////////////////////////////////////////

   /** Initialize the mesh with the current dimension This also calls
     * commit on the meta data causing it to be frozen. Information
     * about elements blocks has to be commited before this. If
     * parallel machine has already been specified through <code>instantiateBulkData</code>
     * that communicator is used. Otherwise a new copy is constructed and
     * will be used through out this mesh object's lifetime.
     *
     * \param[in] parallelMach Communicator
     * \param[in] setupIO If set to true and IOSS is enabled, the output mesh will be initialized.
     * \param[in] buildRefinementSupport If true, build percept uniform refinement objects.
     */
   void initialize(stk::ParallelMachine parallelMach,bool setupIO=true,
                  const bool buildRefinementSupport = false);

   /** Build a bulk data object but don't do anything with it.
     * If parallel machine has already been specified through <code>initialize</code>
     * that communicator is used. Otherwise a new copy is constructed and will
     * be used by default throughout the lifetime of this object.
     */
   void instantiateBulkData(stk::ParallelMachine parallelMach);

   // functions to manage and manipulate bulk data
   //////////////////////////////////////////

   /** Put the bulk data manager in modification mode.
     */
   void beginModification();

   /** Take the bulk data manager out of modification mode.
     */
   void endModification();

   /** Add a node to the mesh with a specific set of coordinates to the mesh.
     *
     * \pre <code>coord.size()==getDimension()</code>
     * \pre <code>isModifiable()==true</code>
     */
   void addNode(stk::mesh::EntityId gid, const std::vector<double> & coord);

   void addElement(const Teuchos::RCP<ElementDescriptor> & ed,stk::mesh::Part * block);

   void addEdges();

   void addFaces();

   /** Adds an entity to a specified side set.
     */
   void addEntityToSideset(stk::mesh::Entity entity,stk::mesh::Part * sideset);

   /** Adds an entity to a specified node set.
     */
   void addEntityToNodeset(stk::mesh::Entity entity,stk::mesh::Part * nodeset);

   /** Adds an entity to a specified edge block.
     */
   void addEntityToEdgeBlock(stk::mesh::Entity entity,stk::mesh::Part * edgeblock);
   /** Adds a vector of entities to a specified edge block.
     */
   void addEntitiesToEdgeBlock(std::vector<stk::mesh::Entity> entities,stk::mesh::Part * edgeblock);

   /** Adds an entity to a specified face block.
     */
   void addEntityToFaceBlock(stk::mesh::Entity entity,stk::mesh::Part * faceblock);
   /** Adds a vector of entities to a specified face block.
     */
   void addEntitiesToFaceBlock(std::vector<stk::mesh::Entity> entities,stk::mesh::Part * faceblock);

   // Methods to interrogate the mesh topology and structure
   //////////////////////////////////////////

   /** Grab the coordinates field
     */
   const VectorFieldType & getCoordinatesField() const
   { return *coordinatesField_; }

   /** Grab the edges field
     */
   const VectorFieldType & getEdgesField() const
   { return *edgesField_; }

   const VectorFieldType & getFacesField() const
   { return *facesField_; }

   /** Look up a global node and get the coordinate.
     */
   const double * getNodeCoordinates(stk::mesh::EntityId nodeId) const;

   /** Look up a global node and get the coordinate.
     */
   const double * getNodeCoordinates(stk::mesh::Entity node) const;

   /** Get subcell global IDs
     */
   void getSubcellIndices(unsigned entityRank,stk::mesh::EntityId elementId,
                          std::vector<stk::mesh::EntityId> & subcellIds) const;

   /** Get a vector of elements owned by this processor
     */
   void getMyElements(std::vector<stk::mesh::Entity> & elements) const;

   /** Get a vector of elements owned by this processor on a particular block ID
     */
   void getMyElements(const std::string & blockID,std::vector<stk::mesh::Entity> & elements) const;

   /** Get a vector of elements that share an edge/face with an owned element. Note that these elements
     * are not owned.
     */
   void getNeighborElements(std::vector<stk::mesh::Entity> & elements) const;

   /** Get a vector of elements not owned by this processor but in a particular block
     */
   void getNeighborElements(const std::string & blockID,std::vector<stk::mesh::Entity> & elements) const;

   /** Get a vector of edges owned by this processor
     */
   void getMyEdges(std::vector<stk::mesh::Entity> & edges) const;

   /** Get Entities corresponding to the edge block requested.
     * The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] edgeBlockName Name of edge block
     * \param[in,out] edges Vector of entities containing the requested edges.
     */
   void getMyEdges(const std::string & edgeBlockName,std::vector<stk::mesh::Entity> & edges) const;

   /** Get Entities corresponding to the locally owned part of the edge block requested. This also limits
     * the entities to be in a particular element block. The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] edgeBlockName Name of edge block
     * \param[in] blockName Name of block
     * \param[in,out] edges Vector of entities containing the requested edges.
     */
   void getMyEdges(const std::string & edgeBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & edges) const;

   /** Get Entities corresponding to the locally owned part of the edge block requested.
     * The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] edgeBlockName Name of edge block
     * \param[in,out] edges Vector of entities containing the requested edges.
     */
   void getAllEdges(const std::string & edgeBlockName,std::vector<stk::mesh::Entity> & edges) const;

   /** Get Entities corresponding to the edge block requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] edgeBlockName Name of edge block
     * \param[in] blockName Name of block
     * \param[in,out] edges Vector of entities containing the requested edges.
     */
   void getAllEdges(const std::string & edgeBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & edges) const;

   /** Get a vector of faces owned by this processor
     */
   void getMyFaces(std::vector<stk::mesh::Entity> & faces) const;

   /** Get Entities corresponding to the face block requested.
     * The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] faceBlockName Name of face block
     * \param[in,out] faces Vector of entities containing the requested faces.
     */
   void getMyFaces(const std::string & faceBlockName,std::vector<stk::mesh::Entity> & faces) const;

   /** Get Entities corresponding to the locally owned part of the face block requested. This also limits
     * the entities to be in a particular element block. The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] faceBlockName Name of face block
     * \param[in] blockName Name of block
     * \param[in,out] faces Vector of entities containing the requested faces.
     */
   void getMyFaces(const std::string & faceBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & faces) const;

   /** Get Entities corresponding to the locally owned part of the face block requested.
     * The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] faceBlockName Name of face block
     * \param[in,out] faces Vector of entities containing the requested faces.
     */
   void getAllFaces(const std::string & faceBlockName,std::vector<stk::mesh::Entity> & faces) const;

   /** Get Entities corresponding to the face block requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be a dimension
     * lower than <code>getDimension()</code>.
     *
     * \param[in] faceBlockName Name of face block
     * \param[in] blockName Name of block
     * \param[in,out] faces Vector of entities containing the requested faces.
     */
   void getAllFaces(const std::string & faceBlockName,const std::string & blockName,std::vector<stk::mesh::Entity> & faces) const;

   /** Get Entities corresponding to the side set requested.
     * The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMySides(const std::string & sideName,std::vector<stk::mesh::Entity> & sides) const;

   /** Get Entities corresponding to the locally owned part of the side set requested. This also limits
     * the entities to be in a particular element block. The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in] blockName Name of block
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMySides(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity> & sides) const;

   /** Get Entities corresponding to the locally owned part of the side set requested.
     * The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getAllSides(const std::string & sideName,std::vector<stk::mesh::Entity> & sides) const;

   /** Get Entities corresponding to the side set requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in] blockName Name of block
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getAllSides(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity> & sides) const;

   /** Get Entities corresponding to the node set requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be of dimension
     * <code>0</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in] blockName Name of block
     * \param[in,out] nodes Vector of entities containing the requested nodes.
     */

   void getMyNodes(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity> & nodes) const;

   /**
    * Searches for connected entity by rank and relation id. Returns
    * invalid entity on failure.
    *
    * \param[in] src The handle to the source entity (the 'from' part of the relation)
    * \param[in] tgt_rank The entity rank of relations to search
    * \param[in] rel_id The id of the relation to search for
    */
   stk::mesh::Entity findConnectivityById(stk::mesh::Entity src, stk::mesh::EntityRank tgt_rank, unsigned rel_id) const;

   // Utility functions
   //////////////////////////////////////////

  /**
   *  \brief Write this mesh and associated fields to the given output file.
   *
   *  \note This will only write a single timestep at time = 0.
   *
   *  \param[in] filename The name of the output Exodus file.
   *  \param[in] append If set to true, the output will be appended to the output Exodus file. If set to false, output file will be overwritten. Default is false.
   */
  void
  writeToExodus(const std::string& filename,
                const bool append = false);

  /**
   *  \brief Set up an output Exodus file for writing results.
   *
   *  Create an output mesh associated with the given `filename` and add the
   *  fields to it.
   *
   *  \note No writing is actually done at this stage.  That happens with a
   *        call to `writeToExodus(double timestep)`.
   *
   *  \param[in] filename The name of the output Exodus file.
   *  \param[in] append If set to true, the output will be appended to the output Exodus file. If set to false, output file will be overwritten. Default is false.
   *  \param[in] append_after_restart_time If set to true, instead of appending to the end of the Exodus file, this option will append new writes after a specified time and overwrite subsequent time steps that were in the file. Allows users to restart anywhere in the exodus file and write consistently. If set to false, the next write will append to the end of the file. Default is false.
   *  \param[in] restart_time If append_after_restart_time is true, this is the time value to append after. Otherwise this value is ignored.
   *
   *  \throws `std::logic_error` If the `STK_Interface` does not yet have a MPI
   *                             communicator.
   */
  void
  setupExodusFile(const std::string& filename,
                  const bool append = false,
                  const bool append_after_restart_time = false,
                  const double restart_time = 0.0);

  /**
   *  \brief Write this mesh and associated fields at the given `timestep`.
   *
   *  Write this mesh and the current state of the associated fields for time =
   *  `timestep` to the output Exodus file created via `setupExodusFile()`.
   *
   *  \note `setupExodusFile(const std::string& filename)` must be called
   *        before any invocations of this routine.
   *
   *  \param[in] timestep The simulation time at which you're writing the
   *                      results.
   */
  void
  writeToExodus(
    double timestep);

  /**
   *  \brief Add an `int` global variable to the information to be written to
   *         the Exodus output file.
   *
   *  This allows you to add any arbitrary global data (that is, data not
   *  associated with nodes, elements, etc.) to the Exodus output file tagged
   *  with a `string` key.
   *
   *  \note Data is not actually written to the output Exodus file until
   *        `writeToExodus()` is called.
   *
   *  \param[in] key   The name of the global variable you'll be adding.  Note
   *                   that keys should be unique if you're adding multiple
   *                   global variables with repeated calls to
   *                   `addGlobalToExodus()`.  Repeated calls with identical
   *                   keys will result in the value associated with the key
   *                   simply being overwritten.
   *  \param[in] value The value of the global variable you'll be adding.
   */
  void
  addGlobalToExodus(
    const std::string& key,
    const int&         value);

  /**
   *  \brief Add a `double` global variable to the information to be written to
   *         the Exodus output file.
   *
   *  This allows you to add any arbitrary global data (that is, data not
   *  associated with nodes, elements, etc.) to the Exodus output file tagged
   *  with a `string` key.
   *
   *  \note Data is not actually written to the output Exodus file until
   *        `writeToExodus()` is called.
   *
   *  \param[in] key   The name of the global variable you'll be adding.  Note
   *                   that keys should be unique if you're adding multiple
   *                   global variables with repeated calls to
   *                   `addGlobalToExodus()`.  Repeated calls with identical
   *                   keys will result in the value associated with the key
   *                   simply being overwritten.
   *  \param[in] value The value of the global variable you'll be adding.
   */
  void
  addGlobalToExodus(
    const std::string& key,
    const double&      value);

  /**
   *  \brief Add a `std::vector<int>` global variable to the information to be
   *         written to the Exodus output file.
   *
   *  This allows you to add any arbitrary global data (that is, data not
   *  associated with nodes, elements, etc.) to the Exodus output file tagged
   *  with a `string` key.
   *
   *  \note Data is not actually written to the output Exodus file until
   *        `writeToExodus()` is called.
   *
   *  \param[in] key   The name of the global variable you'll be adding.  Note
   *                   that keys should be unique if you're adding multiple
   *                   global variables with repeated calls to
   *                   `addGlobalToExodus()`.  Repeated calls with identical
   *                   keys will result in the value associated with the key
   *                   simply being overwritten.
   *  \param[in] value The value of the global variable you'll be adding.
   */
  void
  addGlobalToExodus(
    const std::string&      key,
    const std::vector<int>& value);

  /**
   *  \brief Add a `std::vector<double>` global variable to the information to
   *         be written to the Exodus output file.
   *
   *  This allows you to add any arbitrary global data (that is, data not
   *  associated with nodes, elements, etc.) to the Exodus output file tagged
   *  with a `string` key.
   *
   *  \note Data is not actually written to the output Exodus file until
   *        `writeToExodus()` is called.
   *
   *  \param[in] key   The name of the global variable you'll be adding.  Note
   *                   that keys should be unique if you're adding multiple
   *                   global variables with repeated calls to
   *                   `addGlobalToExodus()`.  Repeated calls with identical
   *                   keys will result in the value associated with the key
   *                   simply being overwritten.
   *  \param[in] value The value of the global variable you'll be adding.
   */
  void
  addGlobalToExodus(
    const std::string&         key,
    const std::vector<double>& value);

   // Accessor functions
   //////////////////////////////////////////

   //! get the comm associated with this mesh
   Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

   Teuchos::RCP<stk::mesh::BulkData> getBulkData() const { return bulkData_; }
   Teuchos::RCP<stk::mesh::MetaData> getMetaData() const { return metaData_; }

#ifdef PANZER_HAVE_PERCEPT
  //! Get the uniformly refined PerceptMesh object.
  Teuchos::RCP<percept::PerceptMesh> getRefinedMesh() const
  { TEUCHOS_ASSERT(Teuchos::nonnull(refinedMesh_)); return refinedMesh_; }
#endif

   bool isWritable() const;

   bool isModifiable() const
   {  if(bulkData_==Teuchos::null) return false;
      return bulkData_->in_modifiable_state(); }

   //! get the dimension
   unsigned getDimension() const
   { return dimension_; }

   //! get the block count
   std::size_t getNumElementBlocks() const
   { return elementBlocks_.size(); }

   /** Get a vector containing the names of the element blocks.
     * This function always returns the current set of element blocks
     * in lexiographic order (uses the sorting built into the std::map).
     * This method can only be called after <code>initialize</code>.
     *
     * \param[in,out] names Vector of names of the element blocks.
     */
   void getElementBlockNames(std::vector<std::string> & names) const;

   /** Get a vector containing the names of the side sets.
     * This function always returns the current set of side sets
     * in lexiographic order (uses the sorting built into the std::map).
     * This method can only be called after <code>initialize</code>.
     *
     * \param[in,out] names Vector of names of the side sets
     */
   void getSidesetNames(std::vector<std::string> & name) const;

   /** Get a vector containing the names of the node sets.
     * This function always returns the current set of node sets
     * in lexiographic order (uses the sorting built into the std::map).
     * This method can only be called after <code>initialize</code>.
     *
     * \param[in,out] names Vector of names of the node sets.
     */
   void getNodesetNames(std::vector<std::string> & name) const;

   /** Get a vector containing the names of the edge blocks.
     * This function always returns the current set of edge blocks
     * in lexiographic order (uses the sorting built into the std::map).
     * This method can only be called after <code>initialize</code>.
     *
     * \param[in,out] names Vector of names of the edge blocks.
     */
   void getEdgeBlockNames(std::vector<std::string> & names) const;

   /** Get a vector containing the names of the face blocks.
     * This function always returns the current set of face blocks
     * in lexiographic order (uses the sorting built into the std::map).
     * This method can only be called after <code>initialize</code>.
     *
     * \param[in,out] names Vector of names of the face blocks.
     */
   void getFaceBlockNames(std::vector<std::string> & names) const;

   //! Get a pointer to the locally owned part
   stk::mesh::Part * getOwnedPart() const
   { return &getMetaData()->locally_owned_part(); } // I don't like the pointer access here, but it will do for now!

   //! get the block part
   stk::mesh::Part * getElementBlockPart(const std::string & name) const
   {
      std::map<std::string, stk::mesh::Part*>::const_iterator itr = elementBlocks_.find(name);   // Element blocks
      if(itr==elementBlocks_.end()) return 0;
      return itr->second;
   }

   //! get the block part
   stk::mesh::Part * getEdgeBlock(const std::string & name) const
   {
      std::map<std::string, stk::mesh::Part*>::const_iterator itr = edgeBlocks_.find(name);   // edge blocks
      if(itr==edgeBlocks_.end()) return 0;
      return itr->second;
   }

   //! get the block part
   stk::mesh::Part * getFaceBlock(const std::string & name) const
   {
      std::map<std::string, stk::mesh::Part*>::const_iterator itr = faceBlocks_.find(name);   // face blocks
      if(itr==faceBlocks_.end()) return 0;
      return itr->second;
   }

   //! get the side set count
   std::size_t getNumSidesets() const
   { return sidesets_.size(); }

   stk::mesh::Part * getSideset(const std::string & name) const
   {
     auto itr = sidesets_.find(name);
     return (itr != sidesets_.end()) ? itr->second : nullptr;
   }

   //! get the side set count
   std::size_t getNumNodesets() const
   { return nodesets_.size(); }

   stk::mesh::Part * getNodeset(const std::string & name) const
   {
     auto itr = nodesets_.find(name);
     return (itr != nodesets_.end()) ? itr->second : nullptr;
   }

   //! get the global counts for the entity of specified rank
   std::size_t getEntityCounts(unsigned entityRank) const;

   //! get max entity ID of type entityRank
   stk::mesh::EntityId getMaxEntityId(unsigned entityRank) const;

   // Utilities
   //////////////////////////////////////////

   //! get a set of elements sharing a single node
   void getElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity> & elements) const;

   //! get a list of node ids for nodes connected to an element
   void getNodeIdsForElement(stk::mesh::Entity element,std::vector<stk::mesh::EntityId> & nodeIds) const;

   /** Get set of element sharing a single node and its local node id.
     */
   void getOwnedElementsSharingNode(stk::mesh::Entity node,std::vector<stk::mesh::Entity> & elements,
                                    std::vector<int> & relIds) const;

   /** Get set of element sharing a single node and its local node id.
     */
   void getOwnedElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity> & elements,
                                    std::vector<int> & relIds, unsigned int matchType) const;


   //! get a set of elements sharing multiple nodes
   void getElementsSharingNodes(const std::vector<stk::mesh::EntityId> nodeId,std::vector<stk::mesh::Entity> & elements) const;

   //! force the mesh to build subcells: edges and faces
   void buildSubcells();

   /** Get an elements local index
     */
   std::size_t elementLocalId(stk::mesh::Entity elmt) const;

   /** Get an elements local index
     */
   std::size_t elementLocalId(stk::mesh::EntityId gid) const;

   /** Get an elements global index
     */
   inline stk::mesh::EntityId elementGlobalId(std::size_t lid) const
   { return bulkData_->identifier((*orderedElementVector_)[lid]); }

   /** Get an elements global index
     */
   inline stk::mesh::EntityId elementGlobalId(stk::mesh::Entity elmt) const
   { return bulkData_->identifier(elmt); }

   /** Is an edge local to this processor?
     */
   bool isEdgeLocal(stk::mesh::Entity edge) const;

   /** Is an edge local to this processor?
     */
   bool isEdgeLocal(stk::mesh::EntityId gid) const;

   /** Get an edge's local index
     */
   std::size_t edgeLocalId(stk::mesh::Entity elmt) const;

   /** Get an edge's local index
     */
   std::size_t edgeLocalId(stk::mesh::EntityId gid) const;

   /** Get an edge's global index
     */
   inline stk::mesh::EntityId edgeGlobalId(std::size_t lid) const
   { return bulkData_->identifier((*orderedEdgeVector_)[lid]); }

   /** Get an edge's global index
     */
   inline stk::mesh::EntityId edgeGlobalId(stk::mesh::Entity edge) const
   { return bulkData_->identifier(edge); }

   /** Is a face local to this processor?
     */
   bool isFaceLocal(stk::mesh::Entity face) const;

   /** Is a face local to this processor?
     */
   bool isFaceLocal(stk::mesh::EntityId gid) const;

   /** Get a face's local index
     */
   std::size_t faceLocalId(stk::mesh::Entity elmt) const;

   /** Get a face's local index
     */
   std::size_t faceLocalId(stk::mesh::EntityId gid) const;

   /** Get a face's global index
     */
   inline stk::mesh::EntityId faceGlobalId(std::size_t lid) const
   { return bulkData_->identifier((*orderedFaceVector_)[lid]); }

   /** Get a face's global index
     */
   inline stk::mesh::EntityId faceGlobalId(stk::mesh::Entity face) const
   { return bulkData_->identifier(face); }

  /** Get an Entity's parallel owner (process rank)
   */
  inline unsigned entityOwnerRank(stk::mesh::Entity entity) const
  { return bulkData_->parallel_owner_rank(entity); }

  /** Check if entity handle is valid
   */
  inline bool isValid(stk::mesh::Entity entity) const
  { return bulkData_->is_valid(entity); }

   /**  Get the containing block ID of this element.
     */
   std::string containingBlockId(stk::mesh::Entity elmt) const;

   /** Get the stk mesh field pointer associated with a particular solution value
     * Assumes there is a field associated with "fieldName,blockId" pair. If none
     * is found an exception (std::runtime_error) is raised.
     */
   stk::mesh::Field<double> * getSolutionField(const std::string & fieldName,
                                               const std::string & blockId) const;

   /** Get the stk mesh field pointer associated with a particular value
     * Assumes there is a field associated with "fieldName,blockId" pair. If none
     * is found an exception (std::runtime_error) is raised.
     */
   stk::mesh::Field<double> * getCellField(const std::string & fieldName,
                                           const std::string & blockId) const;

   /** Get the stk mesh field pointer associated with a particular value
     * Assumes there is a field associated with "fieldName,blockId" pair. If none
     * is found an exception (std::runtime_error) is raised.
     */
   stk::mesh::Field<double> * getEdgeField(const std::string & fieldName,
                                           const std::string & blockId) const;

   /** Get the stk mesh field pointer associated with a particular value
     * Assumes there is a field associated with "fieldName,blockId" pair. If none
     * is found an exception (std::runtime_error) is raised.
     */
   stk::mesh::Field<double> * getFaceField(const std::string & fieldName,
                                           const std::string & blockId) const;

   ProcIdFieldType * getProcessorIdField() { return processorIdField_; }

   //! Has <code>initialize</code> been called on this mesh object?
   bool isInitialized() const { return initialized_; }

   /** Get Vector of element entities ordered by their LID, returns an RCP so that
     * it is easily stored by the caller.
     */
   Teuchos::RCP<const std::vector<stk::mesh::Entity> > getElementsOrderedByLID() const;

   /** Writes a particular field to an array. Notice this is setup to work with
     * the worksets associated with Panzer.
     *
     * \param[in] fieldName Name of field to be filled
     * \param[in] blockId Name of block this set of elements belongs to
     * \param[in] localElementIds Local element IDs for this set of solution values
     * \param[in] solutionValues An two dimensional array object sized by (Cells,Basis Count)
     *
     * \note The block ID is not strictly needed in this context. However forcing the
     *       user to provide it does permit an additional level of safety. The implicit
     *       assumption is that the elements being "set" are part of the specified block.
     *       This prevents the need to perform a null pointer check on the field data, because
     *       the STK_Interface construction of the fields should force it to be nonnull...
     */
   template <typename ArrayT>
   void setSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                             const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues,double scaleValue=1.0);

   /** Reads a particular field into an array. Notice this is setup to work with
     * the worksets associated with Panzer.
     *
     * \param[in] fieldName Name of field to be filled
     * \param[in] blockId Name of block this set of elements belongs to
     * \param[in] localElementIds Local element IDs for this set of solution values
     * \param[in] solutionValues An two dimensional array object sized by (Cells,Basis Count)
     *
     * \note The block ID is not strictly needed in this context. However forcing the
     *       user to provide it does permit an additional level of safety. The implicit
     *       assumption is that the elements being retrieved are part of the specified block.
     *       This prevents the need to perform a null pointer check on the field data, because
     *       the STK_Interface construction of the fields should force it to be nonnull...
     */
   template <typename ArrayT>
   void getSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                             const std::vector<std::size_t> & localElementIds,ArrayT & solutionValues) const;

   /** Writes a particular field to a cell array. Notice this is setup to work with
     * the worksets associated with Panzer.
     *
     * \param[in] fieldName Name of field to be filled
     * \param[in] blockId Name of block this set of elements belongs to
     * \param[in] localElementIds Local element IDs for this set of solution values
     * \param[in] solutionValues A one dimensional array object sized by (Cells)
     *
     * \note The block ID is not strictly needed in this context. However forcing the
     *       user to provide it does permit an additional level of safety. The implicit
     *       assumption is that the elements being "set" are part of the specified block.
     *       This prevents the need to perform a null pointer check on the field data, because
     *       the STK_Interface construction of the fields should force it to be nonnull...
     */
   template <typename ArrayT>
   void setCellFieldData(const std::string & fieldName,const std::string & blockId,
                         const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues,double scaleValue=1.0);

   /** Get Vector of edge entities ordered by their LID, returns an RCP so that
     * it is easily stored by the caller.
     */
   Teuchos::RCP<const std::vector<stk::mesh::Entity> > getEdgesOrderedByLID() const;

   /** Get Vector of face entities ordered by their LID, returns an RCP so that
     * it is easily stored by the caller.
     */
   Teuchos::RCP<const std::vector<stk::mesh::Entity> > getFacesOrderedByLID() const;

   /** Writes a particular field to an edge array. Notice this is setup to work with
     * the worksets associated with Panzer.
     *
     * \param[in] fieldName Name of field to be filled
     * \param[in] blockId Name of block this set of elements belongs to
     * \param[in] localEdgeIds Local edge IDs for this set of solution values
     * \param[in] edgeValues A one dimensional array object sized by (edges)
     *
     * \note The block ID is not strictly needed in this context. However forcing the
     *       user to provide it does permit an additional level of safety. The implicit
     *       assumption is that the elements being "set" are part of the specified block.
     *       This prevents the need to perform a null pointer check on the field data, because
     *       the STK_Interface construction of the fields should force it to be nonnull...
     */
   template <typename ArrayT>
   void setEdgeFieldData(const std::string & fieldName,const std::string & blockId,
                         const std::vector<std::size_t> & localEdgeIds,const ArrayT & edgeValues,double scaleValue=1.0);

   /** Writes a particular field to a face array. Notice this is setup to work with
     * the worksets associated with Panzer.
     *
     * \param[in] fieldName Name of field to be filled
     * \param[in] blockId Name of block this set of elements belongs to
     * \param[in] localFaceIds Local face IDs for this set of solution values
     * \param[in] faceValues A one dimensional array object sized by (faces)
     *
     * \note The block ID is not strictly needed in this context. However forcing the
     *       user to provide it does permit an additional level of safety. The implicit
     *       assumption is that the elements being "set" are part of the specified block.
     *       This prevents the need to perform a null pointer check on the field data, because
     *       the STK_Interface construction of the fields should force it to be nonnull...
     */
   template <typename ArrayT>
   void setFaceFieldData(const std::string & fieldName,const std::string & blockId,
                         const std::vector<std::size_t> & localFaceIds,const ArrayT & faceValues,double scaleValue=1.0);

   //////////// TO BE DEPRECATED

   /** Get vertices associated with a number of elements of the same geometry.
     *
     * \param[in] localIds Element local IDs to construct vertices
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVertices(const std::vector<std::size_t> & localIds, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry.
     *
     * \param[in] elements Element entities to construct vertices
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVertices(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry.
     *
     * \param[in] localIds Element local IDs to construct vertices
     * \param[in] eBlock Element block the elements are in
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVertices(const std::vector<std::size_t> & localIds,const std::string & eBlock, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry.
     *
     * \param[in] elements Element entities to construct vertices
     * \param[in] eBlock Element block the elements are in
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVertices(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry. The vertex array will not be resized.
     *
     * \param[in] localIds Element local IDs to construct vertices
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVerticesNoResize(const std::vector<std::size_t> & localIds, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry. The vertex array will not be resized.
     *
     * \param[in] elements Element entities to construct vertices
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVerticesNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry. The vertex array will not be resized.
     *
     * \param[in] localIds Element local IDs to construct vertices
     * \param[in] eBlock Element block the elements are in
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVerticesNoResize(const std::vector<std::size_t> & localIds,const std::string & eBlock, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry. The vertex array will not be resized.
     *
     * \param[in] elements Element entities to construct vertices
     * \param[in] eBlock Element block the elements are in
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVerticesNoResize(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & vertices) const;

   ///////////// END TO BE DEPRECATED

   /** Get nodes associated with a number of elements of the same geometry.
     *
     * \param[in] localIds Element local IDs to construct nodes
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */

   template <typename ArrayT>
   void getElementNodes(const std::vector<std::size_t> & localIds, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry.
     *
     * \param[in] elements Element entities to construct nodes
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodes(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry.
     *
     * \param[in] localIds Element local IDs to construct nodes
     * \param[in] eBlock Element block the elements are in
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodes(const std::vector<std::size_t> & localIds,const std::string & eBlock, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry.
     *
     * \param[in] elements Element entities to construct nodes
     * \param[in] eBlock Element block the elements are in
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodes(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry. The node array will not be resized.
     *
     * \param[in] localIds Element local IDs to construct nodes
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodesNoResize(const std::vector<std::size_t> & localIds, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry. The node array will not be resized.
     *
     * \param[in] elements Element entities to construct nodes
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodesNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry. The node array will not be resized.
     *
     * \param[in] localIds Element local IDs to construct nodes
     * \param[in] eBlock Element block the elements are in
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodesNoResize(const std::vector<std::size_t> & localIds,const std::string & eBlock, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry. The node array will not be resized.
     *
     * \param[in] elements Element entities to construct nodes
     * \param[in] eBlock Element block the elements are in
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodesNoResize(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & nodes) const;

   // const stk::mesh::FEMInterface & getFEMInterface() const
   // { return *femPtr_; }

   stk::mesh::EntityRank getElementRank() const { return stk::topology::ELEMENT_RANK; }
   stk::mesh::EntityRank getSideRank() const { return metaData_->side_rank(); }
   stk::mesh::EntityRank getFaceRank() const { return stk::topology::FACE_RANK; }
   stk::mesh::EntityRank getEdgeRank() const { return stk::topology::EDGE_RANK; }
   stk::mesh::EntityRank getNodeRank() const { return stk::topology::NODE_RANK; }

   /** Build fields and parts from the meta data
     */
   void initializeFromMetaData();

   /** Setup local element IDs
     */
   void buildLocalElementIDs();

   /** Setup local edge IDs
     */
   void buildLocalEdgeIDs();

   /** Setup local face IDs
     */
   void buildLocalFaceIDs();

   /** Return a vector containing all the periodic boundary conditions.
     */
   const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > &
   getPeriodicBCVector() const
   { return periodicBCs_; }

   /** Return a vector containing all the periodic boundary conditions.
     */
   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > &
   getPeriodicBCVector()
   { return periodicBCs_; }

   /** Return a flag indicating if the bounding box search is used
    * when matching periodic Ids. */ 
   const bool & useBoundingBoxSearch() const
   { return useBBoxSearch_; }

   /** Set the periodic search flag. Indicates if the bounding box search is used */
   void setBoundingBoxSearchFlag(const bool & searchFlag)
   { useBBoxSearch_ = searchFlag; return; }

   /** Add a periodic boundary condition.
     *
     * \note This does not actually change the underlying mesh.
     *       The object itself only communciates the matched IDs (currently nodes)
     */
   void addPeriodicBC(const Teuchos::RCP<const PeriodicBC_MatcherBase> & bc)
   { periodicBCs_.push_back(bc); }

   /** Add a set of periodic boundary conditions.
     *
     * \note This does not actually change the underlying mesh.
     *       The object itself only communciates the matched IDs (currently nodes)
     */
   void addPeriodicBCs(const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & bc_vec)
   { periodicBCs_.insert(periodicBCs_.end(),bc_vec.begin(),bc_vec.end()); }

   /** Pairs DOFs on periodic entities 
    */
   std::pair<Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >, Teuchos::RCP<std::vector<unsigned int> > >
   getPeriodicNodePairing() const;

   /** check for a valid block id
     */
   bool validBlockId(const std::string & blockId) const;

   /** Print a brief description of the STK mesh to a stream.
     */
   void print(std::ostream & os) const;

   /** Print a brief description of the STK mesh to a stream.
     */
   void printMetaData(std::ostream & os) const;

   /** Get the cell topology from the element block.
     */
   Teuchos::RCP<const shards::CellTopology> getCellTopology(const std::string & eBlock) const;

   /** Get the value of the time-stamp the last time this object was written to Exodus (default 0.0)
     *
     * \note The initial state time is completely disconnected from the current state time
     */
   double getCurrentStateTime() const { return currentStateTime_; }

   /** Get the value of the time-stamp when this object was created (default 0.0)
     * or when setInitialStateTime was called.
     *
     * \note The initial state time is completely disconnected from the current state time
     */
   double getInitialStateTime() const { return initialStateTime_; }

   /** Set the value of the initial time-stamp
     *
     * \note The initial state time is completely disconnected from the current state time
     */
   void setInitialStateTime(double value) { initialStateTime_ = value; }

   /** Rebalance using zoltan. Note that this will void the local element ids.
     */
   void rebalance(const Teuchos::ParameterList & params);

   /** Set the weight for a particular element block. Larger means more costly
     * to assemble and evaluate.
     */
   void setBlockWeight(const std::string & blockId,double weight)
   { blockWeights_[blockId] = weight; }

   /** When coordinates are returned in the getElementNodes
     * method, extract coordinates using a specified field (not the intrinsic coordinates)
     * where available (where unavailable the intrinsic coordinates are used.
     * Note that this does not change the behavior of getNodeCoordinates.
     * This is set to false by default.
     */
   void setUseFieldCoordinates(bool useFieldCoordinates)
   { useFieldCoordinates_ = useFieldCoordinates; }

   /** Return the use field coordinates flag */
   bool getUseFieldCoordinates() const
   { return useFieldCoordinates_; }

   /** Use lower case (or not) for I/O */
   void setUseLowerCaseForIO(bool useLowerCase)
   { useLowerCase_ = useLowerCase; }

   /** Use lower case (or not) for I/O */
   bool getUseLowerCaseForIO() const
   { return useLowerCase_; }

   ///////////////// TO BE DEPRECATED....

   /** Get vertices associated with a number of elements of the same geometry, note that a coordinate field
     * will be used (if not is available an exception will be thrown).
     *
     * \param[in] elements Element entities to construct vertices
     * \param[in] eBlock Element block the elements are in
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVertices_FromField(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & vertices) const;

   template <typename ArrayT>
   void getElementVertices_FromFieldNoResize(const std::vector<stk::mesh::Entity> & elements,
                                             const std::string & eBlock, ArrayT & vertices) const;

   /** Get vertices associated with a number of elements of the same geometry. This access the true mesh coordinates
     * array.
     *
     * \param[in] elements Element entities to construct vertices
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVertices_FromCoords(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const;

   template <typename ArrayT>
   void getElementVertices_FromCoordsNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const;

   /////////////////// END TO BE DEPRECATED

   /** Get nodes associated with a number of elements of the same geometry, note that a coordinate field
     * will be used (if not is available an exception will be thrown).
     *
     * \param[in] elements Element entities to construct nodes
     * \param[in] eBlock Element block the elements are in
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodes_FromField(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & nodes) const;

   template <typename ArrayT>
   void getElementNodes_FromFieldNoResize(const std::vector<stk::mesh::Entity> & elements,
                                             const std::string & eBlock, ArrayT & nodes) const;

   /** Get nodes associated with a number of elements of the same geometry. This access the true mesh coordinates
     * array.
     *
     * \param[in] elements Element entities to construct nodes
     * \param[out] nodes Output array that will be sized (<code>localIds.size()</code>,#Nodes,#Dim)
     *
     * \note If not all elements have the same number of nodes an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementNodes_FromCoords(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const;

   template <typename ArrayT>
   void getElementNodes_FromCoordsNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const;

  /** Uniformly refine the mesh using Percept
   *
   * \param[in] numberOfLevels Number of uniform refinement levels to apply. Must be >=1.
   * \param[in] deleteParentElements If true, deletes the parent elements from the mesh to save memory.
   */
  void refineMesh(const int numberOfLevels, const bool deleteParentElements);

public: // static operations
   static const std::string coordsString;
   static const std::string nodesString;
   static const std::string edgesString;
   static const std::string edgeBlockString;
   static const std::string faceBlockString;
   static const std::string facesString;

protected:

   /** Compute global entity counts.
     */
   void buildEntityCounts();

   /** Compute global entity counts.
     */
   void buildMaxEntityIds();

   /** Initialize STK fields using a map (allocate space for storage and writing)
     * to a specific entity rank.
     */
   void initializeFieldsInSTK(const std::map<std::pair<std::string,std::string>,SolutionFieldType*> & nameToField,
                              bool setupIO);

   /** Build a safely handled Teuchos MPI communicator from a parallel machine.
     * This object asserts ownership of the communicator so that we can gurantee
     * existence so the outer user can do whatever they'd like with the original.
     */
   Teuchos::RCP<Teuchos::MpiComm<int> > getSafeCommunicator(stk::ParallelMachine parallelMach) const;

   /** In a pure local operation apply the user specified block weights for each
     * element block to the field that defines the load balance weighting. This
     * uses the blockWeights_ member to determine the user value that has been
     * set for a particular element block.
     */
   void applyElementLoadBalanceWeights();

   /** Determine if a particular field in an element block is a displacement field. Return
     * the displacement field name in the requested argument if so.
     */
   bool isMeshCoordField(const std::string & eBlock,const std::string & fieldName,int & axis) const;

   /** Writes a particular field to an array as a coordinate diplacement. Notice this is setup to work with
     * the worksets associated with Panzer.
     *
     * \param[in] fieldName Name of field to be filled
     * \param[in] blockId Name of block this set of elements belongs to
     * \param[in] dimension What coordinate axis to write to
     * \param[in] localElementIds Local element IDs for this set of solution values
     * \param[in] solutionValues An two dimensional array object sized by (Cells,Basis Count)
     *
     * \note The block ID is not strictly needed in this context. However forcing the
     *       user to provide it does permit an additional level of safety. The implicit
     *       assumption is that the elements being "set" are part of the specified block.
     *       This prevents the need to perform a null pointer check on the field data, because
     *       the STK_Interface construction of the fields should force it to be nonnull...
     */
   template <typename ArrayT>
   void setDispFieldData(const std::string & fieldName,const std::string & blockId,int axis,
                         const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues);

   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > periodicBCs_;
   bool useBBoxSearch_ = false;  // TODO swap this to change default periodic BC search (see also PeriodicBC_Parser.cpp)

   Teuchos::RCP<stk::mesh::MetaData> metaData_;
   Teuchos::RCP<stk::mesh::BulkData> bulkData_;
#ifdef PANZER_HAVE_PERCEPT
  Teuchos::RCP<percept::PerceptMesh> refinedMesh_;
  Teuchos::RCP<percept::URP_Heterogeneous_3D> breakPattern_;
#endif

   std::map<std::string, stk::mesh::Part*> elementBlocks_;  // Element blocks
   std::map<std::string, stk::mesh::Part*> sidesets_;       // Side sets
   std::map<std::string, stk::mesh::Part*> nodesets_;       // Node sets
   std::map<std::string, stk::mesh::Part*> edgeBlocks_;     // Edge blocks
   std::map<std::string, stk::mesh::Part*> faceBlocks_;     // Face blocks

   std::map<std::string, Teuchos::RCP<shards::CellTopology> > elementBlockCT_;

   // for storing/accessing nodes
   stk::mesh::Part * nodesPart_;
   std::vector<stk::mesh::Part*> nodesPartVec_;
   stk::mesh::Part * edgesPart_;
   std::vector<stk::mesh::Part*> edgesPartVec_;
   stk::mesh::Part * facesPart_;
   std::vector<stk::mesh::Part*> facesPartVec_;

   VectorFieldType * coordinatesField_;
   VectorFieldType * edgesField_;
   VectorFieldType * facesField_;
   ProcIdFieldType * processorIdField_;
   SolutionFieldType * loadBalField_;

   // maps field names to solution field stk mesh handles
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToSolution_;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToCellField_;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToEdgeField_;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToFaceField_;

   // use a set to maintain a list of unique information records
   std::set<std::string> informationRecords_;

   unsigned dimension_;

   bool initialized_;

   // how many elements, faces, edges, and nodes are there globally
   std::vector<std::size_t> entityCounts_;

   // what is maximum entity ID
   std::vector<stk::mesh::EntityId> maxEntityId_;

   unsigned procRank_;
   std::size_t currentLocalId_;

   Teuchos::RCP<Teuchos::MpiComm<int> > mpiComm_;

   double initialStateTime_; // the time stamp at the time this object was constructed (default 0.0)
   double currentStateTime_; // the time stamp set by the user when writeToExodus is called (default 0.0)

#ifdef PANZER_HAVE_IOSS
   // I/O support
   Teuchos::RCP<stk::io::StkMeshIoBroker> meshData_;
   int meshIndex_;

  /**
   *  \brief An enumeration to indicate to `globalToExodus()` whether it should
   *         be adding or writing the global variables to the mesh database.
   */
  enum class GlobalVariable
  {
    ADD,
    WRITE
  }; // end of enum class GlobalVariable

  /**
   *  \brief Add or write global variables to the mesh database.
   *
   *  This routine either adds or writes, depending on the value of `flag`, any
   *  global variables added via `addGlobalToExodus()` to the mesh database.
   *
   *  \param[in] flag Either `GlobalVariable::ADD` or `GlobalVariable::WRITE`,
   *                  indicating that the global variables should be added or
   *                  written to the mesh database, respectively.
   *
   *  \throws `std::invalid_argument` If a global variable is not an `int`,
   *                                  `double`, `std::string`,
   *                                  `std::vector<int>`,
   *                                  `std::vector<double>`, or
   *                                  `std::vector<std::string>`.
   */
  void
  globalToExodus(
    const GlobalVariable& flag);

  /**
   *  \brief The global variable(s) to be added to the Exodus output.
   */
  Teuchos::ParameterList globalData_;
#endif

   // uses lazy evaluation
   mutable Teuchos::RCP<std::vector<stk::mesh::Entity> > orderedElementVector_;

   // uses lazy evaluation
   mutable Teuchos::RCP<std::vector<stk::mesh::Entity> > orderedEdgeVector_;

   // uses lazy evaluation
   mutable Teuchos::RCP<std::vector<stk::mesh::Entity> > orderedFaceVector_;

   // for element block weights
   std::map<std::string,double> blockWeights_;

   std::unordered_map<stk::mesh::EntityId,std::size_t> localIDHash_;
   std::unordered_map<stk::mesh::EntityId,std::size_t> localEdgeIDHash_;
   std::unordered_map<stk::mesh::EntityId,std::size_t> localFaceIDHash_;

   // Store mesh displacement fields by element block. This map
   // goes like this meshCoordFields_[eBlock][axis_index] => coordinate FieldName
   // goes like this meshDispFields_[eBlock][axis_index] => displacement FieldName
   std::map<std::string,std::vector<std::string> > meshCoordFields_; // coordinate  fields written by user
   std::map<std::string,std::vector<std::string> > meshDispFields_;  // displacement fields, output to exodus

   bool useFieldCoordinates_;

   bool useLowerCase_;

   // Object describing how to sort a vector of elements using
   // local ID as the key, very short lived object
   class LocalIdCompare {
   public:
     LocalIdCompare(const STK_Interface * mesh) : mesh_(mesh) {}

     // Compares two stk mesh entities based on local ID
     bool operator() (stk::mesh::Entity a,stk::mesh::Entity b)
     { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b);}

   private:
     const STK_Interface * mesh_;
   };
};

template <typename ArrayT>
void STK_Interface::setSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                                         const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues,double scaleValue)
{
   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());
   auto solutionValues_h = Kokkos::create_mirror_view(solutionValues);
   Kokkos::deep_copy(solutionValues_h, solutionValues);

   int field_axis = -1;
   if(isMeshCoordField(blockId,fieldName,field_axis)) {
     setDispFieldData(fieldName,blockId,field_axis,localElementIds,solutionValues_h);
     return;
   }

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity element = elements[localId];

      // loop over nodes set solution values
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* nodes = bulkData_->begin_nodes(element);
      for(std::size_t i=0; i<num_nodes; ++i) {
        stk::mesh::Entity node = nodes[i];

        double * solnData = stk::mesh::field_data(*field,node);
        // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
        solnData[0] = scaleValue*solutionValues_h(cell,i);
      }
   }
}

template <typename ArrayT>
void STK_Interface::setDispFieldData(const std::string & fieldName,const std::string & blockId,int axis,
                                     const std::vector<std::size_t> & localElementIds,const ArrayT & dispValues)
{
   TEUCHOS_ASSERT(axis>=0); // sanity check

   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);
   const VectorFieldType & coord_field = this->getCoordinatesField();

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity element = elements[localId];

      // loop over nodes set solution values
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* nodes = bulkData_->begin_nodes(element);
      for(std::size_t i=0; i<num_nodes; ++i) {
        stk::mesh::Entity node = nodes[i];

        double * solnData = stk::mesh::field_data(*field,node);
        double * coordData = stk::mesh::field_data(coord_field,node);
        // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
        solnData[0] = dispValues(cell,i)-coordData[axis];
      }
   }
}

template <typename ArrayT>
void STK_Interface::getSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                                         const std::vector<std::size_t> & localElementIds,ArrayT & solutionValues) const
{
   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

   solutionValues = Kokkos::createDynRankView(solutionValues,
					      "solutionValues",
					      localElementIds.size(),
					      bulkData_->num_nodes(elements[localElementIds[0]]));

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity element = elements[localId];

      // loop over nodes set solution values
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* nodes = bulkData_->begin_nodes(element);
      for(std::size_t i=0; i<num_nodes; ++i) {
        stk::mesh::Entity node = nodes[i];

        double * solnData = stk::mesh::field_data(*field,node);
        // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
        solutionValues(cell,i) = solnData[0];
      }
   }
}

template <typename ArrayT>
void STK_Interface::setCellFieldData(const std::string & fieldName,const std::string & blockId,
                                     const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues,double scaleValue)
{
   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

   SolutionFieldType * field = this->getCellField(fieldName,blockId);

   auto solutionValues_h = Kokkos::create_mirror_view(solutionValues);
   Kokkos::deep_copy(solutionValues_h, solutionValues);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity element = elements[localId];

      double * solnData = stk::mesh::field_data(*field,element);
      TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
      solnData[0] = scaleValue*solutionValues_h.access(cell,0);
   }
}

template <typename ArrayT>
void STK_Interface::setEdgeFieldData(const std::string & fieldName,const std::string & blockId,
                                     const std::vector<std::size_t> & localEdgeIds,const ArrayT & edgeValues,double scaleValue)
{
   const std::vector<stk::mesh::Entity> & edges = *(this->getEdgesOrderedByLID());

   SolutionFieldType * field = this->getEdgeField(fieldName,blockId);

   auto edgeValues_h = Kokkos::create_mirror_view(edgeValues);
   Kokkos::deep_copy(edgeValues_h, edgeValues);

   for(std::size_t idx=0;idx<localEdgeIds.size();idx++) {
      std::size_t localId = localEdgeIds[idx];
      stk::mesh::Entity edge = edges[localId];

      double * solnData = stk::mesh::field_data(*field,edge);
      TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
      solnData[0] = scaleValue*edgeValues_h.access(idx,0);
   }
}

template <typename ArrayT>
void STK_Interface::setFaceFieldData(const std::string & fieldName,const std::string & blockId,
                                     const std::vector<std::size_t> & localFaceIds,const ArrayT & faceValues,double scaleValue)
{
   const std::vector<stk::mesh::Entity> & faces = *(this->getFacesOrderedByLID());

   SolutionFieldType * field = this->getFaceField(fieldName,blockId);

   auto faceValues_h = Kokkos::create_mirror_view(faceValues);
   Kokkos::deep_copy(faceValues_h, faceValues);

   for(std::size_t idx=0;idx<localFaceIds.size();idx++) {
      std::size_t localId = localFaceIds[idx];
      stk::mesh::Entity face = faces[localId];

      double * solnData = stk::mesh::field_data(*field,face);
      TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
      solnData[0] = scaleValue*faceValues_h.access(idx,0);
   }
}

//////// TO BE DEPRECATED...

template <typename ArrayT>
void STK_Interface::getElementVertices(const std::vector<std::size_t> & localElementIds, ArrayT & vertices) const
{
   if(!useFieldCoordinates_) {
     //
     // gather from the intrinsic mesh coordinates (non-lagrangian)
     //

     const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

     // convert to a vector of entity objects
     std::vector<stk::mesh::Entity> selected_elements;
     for(std::size_t cell=0;cell<localElementIds.size();cell++)
       selected_elements.push_back(elements[localElementIds[cell]]);

     getElementVertices_FromCoords(selected_elements,vertices);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementVertices: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementVertices(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const
{
   if(!useFieldCoordinates_) {
     getElementVertices_FromCoords(elements,vertices);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementVertices: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementVertices(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & vertices) const
{
   if(!useFieldCoordinates_) {
     getElementVertices_FromCoords(elements,vertices);
   }
   else {
     getElementVertices_FromField(elements,eBlock,vertices);
   }
}

template <typename ArrayT>
void STK_Interface::getElementVertices(const std::vector<std::size_t> & localElementIds,const std::string & eBlock, ArrayT & vertices) const
{
   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

   // convert to a vector of entity objects
   std::vector<stk::mesh::Entity> selected_elements;
   for(std::size_t cell=0;cell<localElementIds.size();cell++)
     selected_elements.push_back(elements[localElementIds[cell]]);

   if(!useFieldCoordinates_) {
     getElementVertices_FromCoords(selected_elements,vertices);
   }
   else {
     getElementVertices_FromField(selected_elements,eBlock,vertices);
   }
}

template <typename ArrayT>
void STK_Interface::getElementVerticesNoResize(const std::vector<std::size_t> & localElementIds, ArrayT & vertices) const
{
   if(!useFieldCoordinates_) {
     //
     // gather from the intrinsic mesh coordinates (non-lagrangian)
     //

     const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

     // convert to a vector of entity objects
     std::vector<stk::mesh::Entity> selected_elements;
     for(std::size_t cell=0;cell<localElementIds.size();cell++)
       selected_elements.push_back(elements[localElementIds[cell]]);

     getElementVertices_FromCoordsNoResize(selected_elements,vertices);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementVerticesNoResize: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementVerticesNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const
{
   if(!useFieldCoordinates_) {
     getElementVertices_FromCoordsNoResize(elements,vertices);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementVerticesNoResize: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementVerticesNoResize(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & vertices) const
{
   if(!useFieldCoordinates_) {
     getElementVertices_FromCoordsNoResize(elements,vertices);
   }
   else {
     getElementVertices_FromFieldNoResize(elements,eBlock,vertices);
   }
}

template <typename ArrayT>
void STK_Interface::getElementVerticesNoResize(const std::vector<std::size_t> & localElementIds,const std::string & eBlock, ArrayT & vertices) const
{
   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

   // convert to a vector of entity objects
   std::vector<stk::mesh::Entity> selected_elements;
   for(std::size_t cell=0;cell<localElementIds.size();cell++)
     selected_elements.push_back(elements[localElementIds[cell]]);

   if(!useFieldCoordinates_) {
     getElementVertices_FromCoordsNoResize(selected_elements,vertices);
   }
   else {
     getElementVertices_FromFieldNoResize(selected_elements,eBlock,vertices);
   }
}

template <typename ArrayT>
void STK_Interface::getElementVertices_FromCoords(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const
{
   // nothing to do! silently return
   if(elements.size() == 0) {
     vertices = Kokkos::createDynRankView(vertices, "vertices", 0, 0, 0);
     return;
   }

   //
   // gather from the intrinsic mesh coordinates (non-lagrangian)
   //

   // get *master* cell toplogy...(belongs to first element)
   const auto masterVertexCount
     = stk::mesh::get_cell_topology(bulkData_->bucket(elements[0]).topology()).getCellTopologyData()->vertex_count;

   // allocate space
   vertices = Kokkos::createDynRankView(vertices, "vertices", elements.size(), masterVertexCount,getDimension());
   auto vertices_h = Kokkos::create_mirror_view(vertices);
   Kokkos::deep_copy(vertices_h, vertices);

   // loop over each requested element
   const auto dim = getDimension();
   for(std::size_t cell = 0; cell < elements.size(); cell++) {
      const auto element = elements[cell];
      TEUCHOS_ASSERT(element != 0);

      const auto vertexCount
        = stk::mesh::get_cell_topology(bulkData_->bucket(element).topology()).getCellTopologyData()->vertex_count;
      TEUCHOS_TEST_FOR_EXCEPTION(vertexCount != masterVertexCount, std::runtime_error,
                         "In call to STK_Interface::getElementVertices all elements "
                         "must have the same vertex count!");

      // loop over all element nodes
      const size_t num_nodes = bulkData_->num_nodes(element);
      auto const* nodes = bulkData_->begin_nodes(element);
      TEUCHOS_TEST_FOR_EXCEPTION(num_nodes!=masterVertexCount,std::runtime_error,
                         "In call to STK_Interface::getElementVertices cardinality of "
                                 "element node relations must be the vertex count!");
      for(std::size_t node = 0; node < num_nodes; ++node) {
        const double * coord = getNodeCoordinates(nodes[node]);

        // set each dimension of the coordinate
        for(unsigned d=0;d<dim;d++)
          vertices_h(cell,node,d) = coord[d];
      }
   }
   Kokkos::deep_copy(vertices, vertices_h);
}

template <typename ArrayT>
void STK_Interface::getElementVertices_FromCoordsNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & vertices) const
{
   // nothing to do! silently return
   if(elements.size()==0) {
      return;
   }

   //
   // gather from the intrinsic mesh coordinates (non-lagrangian)
   //

   // get *master* cell toplogy...(belongs to first element)
   unsigned masterVertexCount
     = stk::mesh::get_cell_topology(bulkData_->bucket(elements[0]).topology()).getCellTopologyData()->vertex_count;

   // loop over each requested element
   unsigned dim = getDimension();
   auto vertices_h = Kokkos::create_mirror_view(vertices);
   for(std::size_t cell=0;cell<elements.size();cell++) {
      stk::mesh::Entity element = elements[cell];
      TEUCHOS_ASSERT(element!=0);

      unsigned vertexCount
        = stk::mesh::get_cell_topology(bulkData_->bucket(element).topology()).getCellTopologyData()->vertex_count;
      TEUCHOS_TEST_FOR_EXCEPTION(vertexCount!=masterVertexCount,std::runtime_error,
                         "In call to STK_Interface::getElementVertices all elements "
                         "must have the same vertex count!");

      // loop over all element nodes
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* nodes = bulkData_->begin_nodes(element);
      TEUCHOS_TEST_FOR_EXCEPTION(num_nodes!=masterVertexCount,std::runtime_error,
                         "In call to STK_Interface::getElementVertices cardinality of "
                                 "element node relations must be the vertex count!");
      for(std::size_t node=0; node<num_nodes; ++node) {
        const double * coord = getNodeCoordinates(nodes[node]);

        // set each dimension of the coordinate
        for(unsigned d=0;d<dim;d++)
          vertices_h(cell,node,d) = coord[d];
      }
   }
   Kokkos::deep_copy(vertices, vertices_h);
}

template <typename ArrayT>
void STK_Interface::getElementVertices_FromField(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & vertices) const
{
   TEUCHOS_ASSERT(useFieldCoordinates_);

   // nothing to do! silently return
   if(elements.size()==0) {
     vertices = Kokkos::createDynRankView(vertices,"vertices",0,0,0);
      return;
   }

   // get *master* cell toplogy...(belongs to first element)
   unsigned masterVertexCount
     = stk::mesh::get_cell_topology(bulkData_->bucket(elements[0]).topology()).getCellTopologyData()->vertex_count;

   // allocate space
   vertices = Kokkos::createDynRankView(vertices,"vertices",elements.size(),masterVertexCount,getDimension());
   auto vertices_h = Kokkos::create_mirror_view(vertices);
   std::map<std::string,std::vector<std::string> >::const_iterator itr = meshCoordFields_.find(eBlock);
   if(itr==meshCoordFields_.end()) {
     // no coordinate field set for this element block
     TEUCHOS_ASSERT(false);
   }

   const std::vector<std::string> & coordField = itr->second;
   std::vector<SolutionFieldType*> fields(getDimension());
   for(std::size_t d=0;d<fields.size();d++) {
     fields[d] = this->getSolutionField(coordField[d],eBlock);
   }

   for(std::size_t cell=0;cell<elements.size();cell++) {
      stk::mesh::Entity element = elements[cell];

      // loop over nodes set solution values
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* nodes = bulkData_->begin_nodes(element);
      for(std::size_t i=0; i<num_nodes; ++i) {
        stk::mesh::Entity node = nodes[i];

        const double * coord = getNodeCoordinates(node);

        for(unsigned d=0;d<getDimension();d++) {
          double * solnData = stk::mesh::field_data(*fields[d],node);

          // recall mesh field coordinates are stored as displacements
          // from the mesh coordinates, make sure to add them together
	vertices_h(cell,i,d) = solnData[0]+coord[d];
        }
      }
   }
   Kokkos::deep_copy(vertices, vertices_h);
}

template <typename ArrayT>
void STK_Interface::getElementVertices_FromFieldNoResize(const std::vector<stk::mesh::Entity> & elements,
                                                         const std::string & eBlock, ArrayT & vertices) const
{
   TEUCHOS_ASSERT(useFieldCoordinates_);

   // nothing to do! silently return
   if(elements.size()==0) {
      return;
   }

   std::map<std::string,std::vector<std::string> >::const_iterator itr = meshCoordFields_.find(eBlock);
   if(itr==meshCoordFields_.end()) {
     // no coordinate field set for this element block
     TEUCHOS_ASSERT(false);
   }

   const std::vector<std::string> & coordField = itr->second;
   std::vector<SolutionFieldType*> fields(getDimension());
   for(std::size_t d=0;d<fields.size();d++) {
     fields[d] = this->getSolutionField(coordField[d],eBlock);
   }

   for(std::size_t cell=0;cell<elements.size();cell++) {
      stk::mesh::Entity element = elements[cell];

      // loop over nodes set solution values
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* nodes = bulkData_->begin_nodes(element);
      for(std::size_t i=0; i<num_nodes; ++i) {
        stk::mesh::Entity node = nodes[i];

        const double * coord = getNodeCoordinates(node);

        for(unsigned d=0;d<getDimension();d++) {
          double * solnData = stk::mesh::field_data(*fields[d],node);

          // recall mesh field coordinates are stored as displacements
          // from the mesh coordinates, make sure to add them together
          vertices(cell,i,d) = solnData[0]+coord[d];
        }
      }
   }
}

//////////// END TO BE DEPRECATED

template <typename ArrayT>
void STK_Interface::getElementNodes(const std::vector<std::size_t> & localElementIds, ArrayT & nodes) const
{
   if(!useFieldCoordinates_) {
     //
     // gather from the intrinsic mesh coordinates (non-lagrangian)
     //

     const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

     // convert to a vector of entity objects
     std::vector<stk::mesh::Entity> selected_elements;
     for(std::size_t cell=0;cell<localElementIds.size();cell++)
       selected_elements.push_back(elements[localElementIds[cell]]);

     getElementNodes_FromCoords(selected_elements,nodes);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementNodes: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodes(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const
{
   if(!useFieldCoordinates_) {
     getElementNodes_FromCoords(elements,nodes);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementNodes: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodes(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & nodes) const
{
   if(!useFieldCoordinates_) {
     getElementNodes_FromCoords(elements,nodes);
   }
   else {
     getElementNodes_FromField(elements,eBlock,nodes);
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodes(const std::vector<std::size_t> & localElementIds,const std::string & eBlock, ArrayT & nodes) const
{
   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

   // convert to a vector of entity objects
   std::vector<stk::mesh::Entity> selected_elements;
   for(std::size_t cell=0;cell<localElementIds.size();cell++)
     selected_elements.push_back(elements[localElementIds[cell]]);

   if(!useFieldCoordinates_) {
     getElementNodes_FromCoords(selected_elements,nodes);
   }
   else {
     getElementNodes_FromField(selected_elements,eBlock,nodes);
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodesNoResize(const std::vector<std::size_t> & localElementIds, ArrayT & nodes) const
{
   if(!useFieldCoordinates_) {
     //
     // gather from the intrinsic mesh coordinates (non-lagrangian)
     //

     const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

     // convert to a vector of entity objects
     std::vector<stk::mesh::Entity> selected_elements;
     for(std::size_t cell=0;cell<localElementIds.size();cell++)
       selected_elements.push_back(elements[localElementIds[cell]]);

     getElementNodes_FromCoordsNoResize(selected_elements,nodes);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementNodesNoResize: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodesNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const
{
   if(!useFieldCoordinates_) {
     getElementNodes_FromCoordsNoResize(elements,nodes);
   }
   else {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                "STK_Interface::getElementNodesNoResize: Cannot call this method when field coordinates are used "
                                "without specifying an element block.");
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodesNoResize(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & nodes) const
{
   if(!useFieldCoordinates_) {
     getElementNodes_FromCoordsNoResize(elements,nodes);
   }
   else {
     getElementNodes_FromFieldNoResize(elements,eBlock,nodes);
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodesNoResize(const std::vector<std::size_t> & localElementIds,const std::string & eBlock, ArrayT & nodes) const
{
   const std::vector<stk::mesh::Entity> & elements = *(this->getElementsOrderedByLID());

   // convert to a vector of entity objects
   std::vector<stk::mesh::Entity> selected_elements;
   for(std::size_t cell=0;cell<localElementIds.size();cell++)
     selected_elements.push_back(elements[localElementIds[cell]]);

   if(!useFieldCoordinates_) {
     getElementNodes_FromCoordsNoResize(selected_elements,nodes);
   }
   else {
     getElementNodes_FromFieldNoResize(selected_elements,eBlock,nodes);
   }
}

template <typename ArrayT>
void STK_Interface::getElementNodes_FromCoords(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const
{
   // nothing to do! silently return
   if(elements.size() == 0) {
     nodes = Kokkos::createDynRankView(nodes, "nodes", 0, 0, 0);
     return;
   }

   //
   // gather from the intrinsic mesh coordinates (non-lagrangian)
   //

   // get *master* cell toplogy...(belongs to first element)
   const auto masterNodeCount
     = stk::mesh::get_cell_topology(bulkData_->bucket(elements[0]).topology()).getCellTopologyData()->node_count;

   // allocate space
   nodes = Kokkos::createDynRankView(nodes, "nodes", elements.size(), masterNodeCount,getDimension());
   auto nodes_h = Kokkos::create_mirror_view(nodes);
   Kokkos::deep_copy(nodes_h, nodes);

   // loop over each requested element
   const auto dim = getDimension();
   for(std::size_t cell = 0; cell < elements.size(); cell++) {
      const auto element = elements[cell];
      TEUCHOS_ASSERT(element != 0);

      const auto nodeCount
        = stk::mesh::get_cell_topology(bulkData_->bucket(element).topology()).getCellTopologyData()->node_count;
      TEUCHOS_TEST_FOR_EXCEPTION(nodeCount != masterNodeCount, std::runtime_error,
                         "In call to STK_Interface::getElementNodes all elements "
                         "must have the same node count!");

      // loop over all element nodes
      const size_t num_nodes = bulkData_->num_nodes(element);
      auto const* elem_nodes = bulkData_->begin_nodes(element);
      TEUCHOS_TEST_FOR_EXCEPTION(num_nodes!=masterNodeCount,std::runtime_error,
                         "In call to STK_Interface::getElementNodes cardinality of "
                                 "element node relations must be the node count!");
      for(std::size_t node = 0; node < num_nodes; ++node) {
        const double * coord = getNodeCoordinates(elem_nodes[node]);

        // set each dimension of the coordinate
        for(unsigned d=0;d<dim;d++)
          nodes_h(cell,node,d) = coord[d];
      }
   }
   Kokkos::deep_copy(nodes, nodes_h);
}

template <typename ArrayT>
void STK_Interface::getElementNodes_FromCoordsNoResize(const std::vector<stk::mesh::Entity> & elements, ArrayT & nodes) const
{
   // nothing to do! silently return
   if(elements.size()==0) {
      return;
   }

   //
   // gather from the intrinsic mesh coordinates (non-lagrangian)
   //

   // get *master* cell toplogy...(belongs to first element)
   unsigned masterNodeCount
     = stk::mesh::get_cell_topology(bulkData_->bucket(elements[0]).topology()).getCellTopologyData()->node_count;

   // loop over each requested element
   unsigned dim = getDimension();
   auto nodes_h = Kokkos::create_mirror_view(nodes);
   for(std::size_t cell=0;cell<elements.size();cell++) {
      stk::mesh::Entity element = elements[cell];
      TEUCHOS_ASSERT(element!=0);

      unsigned nodeCount
        = stk::mesh::get_cell_topology(bulkData_->bucket(element).topology()).getCellTopologyData()->node_count;
      TEUCHOS_TEST_FOR_EXCEPTION(nodeCount!=masterNodeCount,std::runtime_error,
                         "In call to STK_Interface::getElementNodes all elements "
                         "must have the same node count!");

      // loop over all element nodes
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* elem_nodes = bulkData_->begin_nodes(element);
      TEUCHOS_TEST_FOR_EXCEPTION(num_nodes!=masterNodeCount,std::runtime_error,
                         "In call to STK_Interface::getElementNodes cardinality of "
                                 "element node relations must be the node count!");
      for(std::size_t node=0; node<num_nodes; ++node) {
        const double * coord = getNodeCoordinates(elem_nodes[node]);

        // set each dimension of the coordinate
        for(unsigned d=0;d<dim;d++)
          nodes_h(cell,node,d) = coord[d];
      }
   }
   Kokkos::deep_copy(nodes, nodes_h);
}

template <typename ArrayT>
void STK_Interface::getElementNodes_FromField(const std::vector<stk::mesh::Entity> & elements,const std::string & eBlock, ArrayT & nodes) const
{
   TEUCHOS_ASSERT(useFieldCoordinates_);

   // nothing to do! silently return
   if(elements.size()==0) {
     nodes = Kokkos::createDynRankView(nodes,"nodes",0,0,0);
      return;
   }

   // get *master* cell toplogy...(belongs to first element)
   unsigned masterNodeCount
     = stk::mesh::get_cell_topology(bulkData_->bucket(elements[0]).topology()).getCellTopologyData()->node_count;

   // allocate space
   nodes = Kokkos::createDynRankView(nodes,"nodes",elements.size(),masterNodeCount,getDimension());
   auto nodes_h = Kokkos::create_mirror_view(nodes);
   std::map<std::string,std::vector<std::string> >::const_iterator itr = meshCoordFields_.find(eBlock);
   if(itr==meshCoordFields_.end()) {
     // no coordinate field set for this element block
     TEUCHOS_ASSERT(false);
   }

   const std::vector<std::string> & coordField = itr->second;
   std::vector<SolutionFieldType*> fields(getDimension());
   for(std::size_t d=0;d<fields.size();d++) {
     fields[d] = this->getSolutionField(coordField[d],eBlock);
   }

   // loop over elements
   for(std::size_t cell=0;cell<elements.size();cell++) {
      stk::mesh::Entity element = elements[cell];

      // loop over nodes set solution values
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* elem_nodes = bulkData_->begin_nodes(element);
      for(std::size_t i=0; i<num_nodes; ++i) {
        stk::mesh::Entity node = elem_nodes[i];

        const double * coord = getNodeCoordinates(node);

        for(unsigned d=0;d<getDimension();d++) {
          double * solnData = stk::mesh::field_data(*fields[d],node);

          // recall mesh field coordinates are stored as displacements
          // from the mesh coordinates, make sure to add them together
	        nodes_h(cell,i,d) = solnData[0]+coord[d];
        }
      }
   }
   Kokkos::deep_copy(nodes, nodes_h);
}

template <typename ArrayT>
void STK_Interface::getElementNodes_FromFieldNoResize(const std::vector<stk::mesh::Entity> & elements,
                                                         const std::string & eBlock, ArrayT & nodes) const
{
   TEUCHOS_ASSERT(useFieldCoordinates_);

   // nothing to do! silently return
   if(elements.size()==0) {
      return;
   }

   std::map<std::string,std::vector<std::string> >::const_iterator itr = meshCoordFields_.find(eBlock);
   if(itr==meshCoordFields_.end()) {
     // no coordinate field set for this element block
     TEUCHOS_ASSERT(false);
   }

   const std::vector<std::string> & coordField = itr->second;
   std::vector<SolutionFieldType*> fields(getDimension());
   for(std::size_t d=0;d<fields.size();d++) {
     fields[d] = this->getSolutionField(coordField[d],eBlock);
   }

   // loop over elements
   for(std::size_t cell=0;cell<elements.size();cell++) {
      stk::mesh::Entity element = elements[cell];

      // loop over nodes set solution values
      const size_t num_nodes = bulkData_->num_nodes(element);
      stk::mesh::Entity const* elem_nodes = bulkData_->begin_nodes(element);
      for(std::size_t i=0; i<num_nodes; ++i) {
        stk::mesh::Entity node = elem_nodes[i];

        const double * coord = getNodeCoordinates(node);

        for(unsigned d=0;d<getDimension();d++) {
          double * solnData = stk::mesh::field_data(*fields[d],node);

          // recall mesh field coordinates are stored as displacements
          // from the mesh coordinates, make sure to add them together
          nodes(cell,i,d) = solnData[0]+coord[d];
        }
      }
   }
}

}

#endif

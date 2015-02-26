// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_STK_Interface_hpp__
#define __Panzer_STK_Interface_hpp__

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <Shards_CellTopology.hpp>
#include <Shards_CellTopologyData.h>

#include <Panzer_STK_config.hpp>

#include <boost/unordered_map.hpp>

#ifdef HAVE_IOSS
#include <stk_io/MeshReadWriteUtils.hpp>
#endif

namespace panzer_stk_classic {

class PeriodicBC_MatcherBase;

/** Pure virtial base class that builds a basic element. To be
  * overidden with several types of elements.
  */ 
class ElementDescriptor {
public:
   ElementDescriptor(stk_classic::mesh::EntityId gid,const std::vector<stk_classic::mesh::EntityId> & nodes);
   virtual ~ElementDescriptor();

   stk_classic::mesh::EntityId getGID() const { return gid_; }
   const std::vector<stk_classic::mesh::EntityId> & getNodes() const { return nodes_; }
protected:
   stk_classic::mesh::EntityId gid_;
   std::vector<stk_classic::mesh::EntityId> nodes_;

   ElementDescriptor();
};

/** Constructor function for building the element descriptors.
  */ 
Teuchos::RCP<ElementDescriptor> 
buildElementDescriptor(stk_classic::mesh::EntityId elmtId,std::vector<stk_classic::mesh::EntityId> & nodes);

class STK_Interface {
public:
   typedef double ProcIdData; // ECC: Not sure why?
   typedef stk_classic::mesh::Field<double> SolutionFieldType;
   typedef stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> VectorFieldType;
   typedef stk_classic::mesh::Field<ProcIdData> ProcIdFieldType;

   // some simple exception classes
   struct ElementBlockException : public std::logic_error 
   { ElementBlockException(const std::string & what) : std::logic_error(what) {} };

   struct SidesetException : public std::logic_error 
   { SidesetException(const std::string & what) : std::logic_error(what) {} };

   STK_Interface();

   /** Default constructor
     */
   STK_Interface(unsigned dim);

   // functions called before initialize
   //////////////////////////////////////////

   /** Add an element block with a string name
     */
   void addElementBlock(const std::string & name,const CellTopologyData * ctData);

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

   //////////////////////////////////////////

   /** Initialize the mesh with the current dimension This also calls
     * commit on the meta data causing it to be frozen. Information
     * about elements blocks has to be commited before this. If 
     * parallel machine has already been specified through <code>instantiateBulkData</code>
     * that communicator is used. Otherwise a new copy is constructed and
     * will be used through out this mesh object's lifetime.
     */
   void initialize(stk_classic::ParallelMachine parallelMach,bool setupIO=true);

   /** Build a bulk data object but don't do anything with it.
     * If parallel machine has already been specified through <code>initialize</code>
     * that communicator is used. Otherwise a new copy is constructed and will
     * be used by default throughout the lifetime of this object.
     */
   void instantiateBulkData(stk_classic::ParallelMachine parallelMach);

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
   void addNode(stk_classic::mesh::EntityId gid, const std::vector<double> & coord);

   void addElement(const Teuchos::RCP<ElementDescriptor> & ed,stk_classic::mesh::Part * block);

   void addEdges();

   /** Addes an entity to a specified side set.
     */
   void addEntityToSideset(stk_classic::mesh::Entity & entity,stk_classic::mesh::Part * sideset);

   /** Addes an entity to a specified node set.
     */
   void addEntityToNodeset(stk_classic::mesh::Entity & entity,stk_classic::mesh::Part * nodeset);

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

   /** Look up a global node and get the coordinate.
     */
   const double * getNodeCoordinates(stk_classic::mesh::EntityId nodeId) const;

   /** Look up a global node and get the coordinate.
     */
   const double * getNodeCoordinates(stk_classic::mesh::Entity * node) const;

   /** Get subcell global IDs
     */
   void getSubcellIndices(unsigned entityRank,stk_classic::mesh::EntityId elementId,
                          std::vector<stk_classic::mesh::EntityId> & subcellIds) const;

   /** Get a vector of elements owned by this processor
     */
   void getMyElements(std::vector<stk_classic::mesh::Entity*> & elements) const;

   /** Get a vector of elements owned by this processor on a particular block ID
     */
   void getMyElements(const std::string & blockID,std::vector<stk_classic::mesh::Entity*> & elements) const;

   /** Get a vector of elements that share an edge/face with an owned element. Note that these elements
     * are not owned.
     */
   void getNeighborElements(std::vector<stk_classic::mesh::Entity*> & elements) const;

   /** Get a vector of elements not owned by this processor but in a particular block
     */
   void getNeighborElements(const std::string & blockID,std::vector<stk_classic::mesh::Entity*> & elements) const;

   /** Get Entities corresponding to the side set requested. 
     * The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMySides(const std::string & sideName,std::vector<stk_classic::mesh::Entity*> & sides) const;

   /** Get Entities corresponding to the side set requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in] blockName Name of block
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMySides(const std::string & sideName,const std::string & blockName,std::vector<stk_classic::mesh::Entity*> & sides) const;

   /** Get Entities corresponding to the node set requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be ofdimension
     * <code>0</code>.
     *
     * \param[in] nodesetName Name of side set
     * \param[in] blockName Name of block
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMyNodes(const std::string & sideName,const std::string & blockName,std::vector<stk_classic::mesh::Entity*> & nodes) const;

   // Utility functions
   //////////////////////////////////////////

   /** Write this mesh to exodus. This is a one shot function
     * that will write to a particular exodus output file.
     */
   void writeToExodus(const std::string & filename);

   /** This simply sets up a transient exodus file for writing.
     * No work is performed at this stage. This is used
     * in combination with <code>writeToExodus(double timestep)</code>.
     */
   void setupTransientExodusFile(const std::string & filename);

   /** Write this timestep to the exodus file specified in the
     * <code>setupTransientExodusFile</code>. This uses the
     * current state of the STK fields as the time step.
     */
   void writeToExodus(double timestep);

   // Accessor functions
   //////////////////////////////////////////
 
   Teuchos::RCP<stk_classic::mesh::BulkData> getBulkData() const { return bulkData_; }
   Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData> getMetaData() const { return metaData_; }

   bool isWritable() const;

   bool isModifiable() const
   {  if(bulkData_==Teuchos::null) return false; 
      return bulkData_->synchronized_state()==stk_classic::mesh::BulkData::MODIFIABLE; }

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
     * \param[in,out] names Vector of names of the element blocks.
     */
   void getSidesetNames(std::vector<std::string> & name) const;

   /** Get a vector containing the names of the node sets.
     * This function always returns the current set of node sets
     * in lexiographic order (uses the sorting built into the std::map).
     * This method can only be called after <code>initialize</code>.
     *
     * \param[in,out] names Vector of names of the element blocks.
     */
   void getNodesetNames(std::vector<std::string> & name) const;

   //! Get a pointer to the locally owned part
   stk_classic::mesh::Part * getOwnedPart() const
   { return &getMetaData()->locally_owned_part(); } // I don't like the pointer access here, but it will do for now!

   //! get the block count
   stk_classic::mesh::Part * getElementBlockPart(const std::string & name) const
   { 
      std::map<std::string, stk_classic::mesh::Part*>::const_iterator itr = elementBlocks_.find(name);   // Element blocks
      if(itr==elementBlocks_.end()) return 0;
      return itr->second; 
   }

   //! get the side set count
   std::size_t getNumSidesets() const
   { return sidesets_.size(); }

   stk_classic::mesh::Part * getSideset(const std::string & name) const
   { return sidesets_.find(name)->second; }

   //! get the side set count
   std::size_t getNumNodesets() const
   { return nodesets_.size(); }

   stk_classic::mesh::Part * getNodeset(const std::string & name) const
   { return nodesets_.find(name)->second; }

   //! get the global counts for the entity of specified rank
   std::size_t getEntityCounts(unsigned entityRank) const;

   //! get max entity ID of type entityRank
   stk_classic::mesh::EntityId getMaxEntityId(unsigned entityRank) const;

   // Utilities
   //////////////////////////////////////////

   //! get a set of elements sharing a single node
   void getElementsSharingNode(stk_classic::mesh::EntityId nodeId,std::vector<stk_classic::mesh::Entity *> & elements) const;

   /** Get set of element sharing a single node and its local node id.
     */
   void getOwnedElementsSharingNode(stk_classic::mesh::Entity * node,std::vector<stk_classic::mesh::Entity *> & elements,
                                                        std::vector<int> & localNodeId) const;

   /** Get set of element sharing a single node and its local node id.
     */
   void getOwnedElementsSharingNode(stk_classic::mesh::EntityId nodeId,std::vector<stk_classic::mesh::Entity *> & elements,
                                                          std::vector<int> & localNodeId, unsigned int matchType) const;


   //! get a set of elements sharing multiple nodes
   void getElementsSharingNodes(const std::vector<stk_classic::mesh::EntityId> nodeId,std::vector<stk_classic::mesh::Entity *> & elements) const;

   //! force the mesh to build subcells: edges and faces
   void buildSubcells();

   /** Get an elements local index
     */
   std::size_t elementLocalId(stk_classic::mesh::Entity * elmt) const;

   /** Get an elements local index
     */
   std::size_t elementLocalId(stk_classic::mesh::EntityId gid) const;

   /** Get an elements global index
     */
   inline stk_classic::mesh::EntityId elementGlobalId(std::size_t lid) const
   { return (*orderedElementVector_)[lid]->identifier(); }

   /** Get an elements global index
     */
   inline stk_classic::mesh::EntityId elementGlobalId(stk_classic::mesh::Entity * elmt) const
   { return elmt->identifier(); }

   /**  Get the containing block ID of this element.
     */ 
   std::string containingBlockId(stk_classic::mesh::Entity * elmt);

   /** Get the stk mesh field pointer associated with a particular solution value
     * Assumes there is a field associated with "fieldName,blockId" pair. If none
     * is found an exception (std::runtime_error) is raised.
     */
   stk_classic::mesh::Field<double> * getSolutionField(const std::string & fieldName,
                                               const std::string & blockId) const;

   /** Get the stk mesh field pointer associated with a particular value
     * Assumes there is a field associated with "fieldName,blockId" pair. If none
     * is found an exception (std::runtime_error) is raised.
     */
   stk_classic::mesh::Field<double> * getCellField(const std::string & fieldName,
                                           const std::string & blockId) const;

   ProcIdFieldType * getProcessorIdField() { return processorIdField_; }

   //! Has <code>initialize</code> been called on this mesh object?
   bool isInitialized() const { return initialized_; }

   /** Get Vector of element entities ordered by their LID, returns an RCP so that
     * it is easily stored by the caller.
     */
   Teuchos::RCP<const std::vector<stk_classic::mesh::Entity*> > getElementsOrderedByLID() const;

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
   void getElementVertices(const std::vector<stk_classic::mesh::Entity*> & elements, ArrayT & vertices) const;

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
   void getElementVertices(const std::vector<stk_classic::mesh::Entity*> & elements,const std::string & eBlock, ArrayT & vertices) const;

   // const stk_classic::mesh::fem::FEMInterface & getFEMInterface() const 
   // { return *femPtr_; }

   stk_classic::mesh::EntityRank getElementRank() const { return metaData_->element_rank(); }
   stk_classic::mesh::EntityRank getSideRank() const { return metaData_->side_rank(); }
   stk_classic::mesh::EntityRank getFaceRank() const { return metaData_->face_rank(); }
   stk_classic::mesh::EntityRank getEdgeRank() const { return metaData_->edge_rank(); }
   stk_classic::mesh::EntityRank getNodeRank() const { return metaData_->node_rank(); }

   /** Build fields and parts from the meta data
     */
   void initializeFromMetaData();

   /** Setup local element IDs
     */
   void buildLocalElementIDs();

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

   /** When coordinates are returned in the getElementVertices
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

public: // static operations
   static const std::string coordsString;
   static const std::string nodesString;
   static const std::string edgesString;

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
                             stk_classic::mesh::EntityRank rank,bool setupIO);

   /** Build a safely handled Teuchos MPI communicator from a parallel machine.
     * This object asserts ownership of the communicator so that we can gurantee
     * existence so the outer user can do whatever they'd like with the original.
     */
   Teuchos::RCP<Teuchos::MpiComm<int> > getSafeCommunicator(stk_classic::ParallelMachine parallelMach) const;

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
   void getElementVertices_FromField(const std::vector<stk_classic::mesh::Entity*> & elements,const std::string & eBlock, ArrayT & vertices) const;

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
   void getElementVertices_FromCoords(const std::vector<stk_classic::mesh::Entity*> & elements, ArrayT & vertices) const;

   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > periodicBCs_;

   Teuchos::RCP<stk_classic::mesh::fem::FEMMetaData> metaData_;
   Teuchos::RCP<stk_classic::mesh::BulkData> bulkData_;

   std::map<std::string, stk_classic::mesh::Part*> elementBlocks_;   // Element blocks
   std::map<std::string, stk_classic::mesh::Part*> sidesets_; // Side sets 
   std::map<std::string, stk_classic::mesh::Part*> nodesets_; // Node sets 
   std::map<std::string, Teuchos::RCP<shards::CellTopology> > elementBlockCT_;

   // for storing/accessing nodes
   stk_classic::mesh::Part * nodesPart_;
   std::vector<stk_classic::mesh::Part*> nodesPartVec_;
   stk_classic::mesh::Part * edgesPart_;
   std::vector<stk_classic::mesh::Part*> edgesPartVec_;

   VectorFieldType * coordinatesField_;
   VectorFieldType * edgesField_;
   ProcIdFieldType * processorIdField_;
   SolutionFieldType * loadBalField_;
   
   // maps field names to solution field stk mesh handles
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToSolution_;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToCellField_;

   unsigned dimension_;

   bool initialized_;

   // how many elements, faces, edges, and nodes are there globally
   std::vector<std::size_t> entityCounts_;

   // what is maximum entity ID
   std::vector<stk_classic::mesh::EntityId> maxEntityId_;

   unsigned procRank_;
   std::size_t currentLocalId_;

   Teuchos::RCP<Teuchos::MpiComm<int> > mpiComm_;

   double initialStateTime_; // the time stamp at the time this object was constructed (default 0.0)
   double currentStateTime_; // the time stamp set by the user when writeToExodus is called (default 0.0)

#ifdef HAVE_IOSS
   // I/O support
   Teuchos::RCP<stk_classic::io::MeshData> meshData_;
#endif

   // uses lazy evaluation
   mutable Teuchos::RCP<std::vector<stk_classic::mesh::Entity*> > orderedElementVector_;

   // for element block weights
   std::map<std::string,double> blockWeights_;

   boost::unordered_map<stk_classic::mesh::EntityId,std::size_t> localIDHash_;

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
     bool operator() (stk_classic::mesh::Entity * a,stk_classic::mesh::Entity * b) 
     { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b);}
   
   private:
     const STK_Interface * mesh_;
   };
};

template <typename ArrayT>
void STK_Interface::setSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                                         const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues,double scaleValue)
{
   const std::vector<stk_classic::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   int field_axis = -1;
   if(isMeshCoordField(blockId,fieldName,field_axis)) {
     setDispFieldData(fieldName,blockId,field_axis,localElementIds,solutionValues); 
     return;
   }

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk_classic::mesh::Entity * element = elements[localId];

      // loop over nodes set solution values
      stk_classic::mesh::PairIterRelation relations = element->relations(getNodeRank());
      for(std::size_t i=0;i<relations.size();++i) {
         stk_classic::mesh::Entity * node = relations[i].entity();

         double * solnData = stk_classic::mesh::field_data(*field,*node);
         // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
         solnData[0] = scaleValue*solutionValues(cell,i);
      }
   }
}

template <typename ArrayT>
void STK_Interface::setDispFieldData(const std::string & fieldName,const std::string & blockId,int axis,
                                     const std::vector<std::size_t> & localElementIds,const ArrayT & dispValues)
{
   TEUCHOS_ASSERT(axis>=0); // sanity check

   const std::vector<stk_classic::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);
   const VectorFieldType & coord_field = this->getCoordinatesField();

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk_classic::mesh::Entity * element = elements[localId];

      // loop over nodes set solution values
      stk_classic::mesh::PairIterRelation relations = element->relations(getNodeRank());
      for(std::size_t i=0;i<relations.size();++i) {
         stk_classic::mesh::Entity * node = relations[i].entity();

         double * solnData = stk_classic::mesh::field_data(*field,*node);
         double * coordData = stk_classic::mesh::field_data(coord_field,*node);
         // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
         solnData[0] = dispValues(cell,i)-coordData[axis];
      }
   }
}

template <typename ArrayT>
void STK_Interface::getSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                                         const std::vector<std::size_t> & localElementIds,ArrayT & solutionValues) const
{
   const std::vector<stk_classic::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   solutionValues.resize(localElementIds.size(),elements[localElementIds[0]]->relations(getNodeRank()).size());

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk_classic::mesh::Entity * element = elements[localId];

      // loop over nodes set solution values
      stk_classic::mesh::PairIterRelation relations = element->relations(getNodeRank());
      for(std::size_t i=0;i<relations.size();++i) {
         stk_classic::mesh::Entity * node = relations[i].entity();

         double * solnData = stk_classic::mesh::field_data(*field,*node);
         // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
         solutionValues(cell,i) = solnData[0]; 
      }
   }
}

template <typename ArrayT>
void STK_Interface::setCellFieldData(const std::string & fieldName,const std::string & blockId,
                                     const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues,double scaleValue)
{
   const std::vector<stk_classic::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   SolutionFieldType * field = this->getCellField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk_classic::mesh::Entity * element = elements[localId];

      double * solnData = stk_classic::mesh::field_data(*field,*element);
      TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
      solnData[0] = scaleValue*solutionValues(cell,0);
   }
}

template <typename ArrayT>
void STK_Interface::getElementVertices(const std::vector<std::size_t> & localElementIds, ArrayT & vertices) const
{
   if(!useFieldCoordinates_) {
     //
     // gather from the intrinsic mesh coordinates (non-lagrangian)
     //

     const std::vector<stk_classic::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());
    
     // convert to a vector of entity objects
     std::vector<stk_classic::mesh::Entity*> selected_elements;
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
void STK_Interface::getElementVertices(const std::vector<stk_classic::mesh::Entity*> & elements, ArrayT & vertices) const
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
void STK_Interface::getElementVertices(const std::vector<stk_classic::mesh::Entity*> & elements,const std::string & eBlock, ArrayT & vertices) const
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
   const std::vector<stk_classic::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   // convert to a vector of entity objects
   std::vector<stk_classic::mesh::Entity*> selected_elements;
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
void STK_Interface::getElementVertices_FromCoords(const std::vector<stk_classic::mesh::Entity*> & elements, ArrayT & vertices) const
{
   // nothing to do! silently return
   if(elements.size()==0) {
      vertices.resize(0,0,0);
      return;
   }

   //
   // gather from the intrinsic mesh coordinates (non-lagrangian)
   //

   // get *master* cell toplogy...(belongs to first element)
   unsigned masterVertexCount 
      = stk_classic::mesh::fem::get_cell_topology(elements[0]->bucket()).getCellTopologyData()->vertex_count;

   // allocate space
   vertices.resize(elements.size(),masterVertexCount,getDimension());

   // loop over each requested element
   unsigned dim = getDimension();
   for(std::size_t cell=0;cell<elements.size();cell++) {
      stk_classic::mesh::Entity * element = elements[cell];
      TEUCHOS_ASSERT(element!=0);
 
      unsigned vertexCount 
         = stk_classic::mesh::fem::get_cell_topology(element->bucket()).getCellTopologyData()->vertex_count;
      TEUCHOS_TEST_FOR_EXCEPTION(vertexCount!=masterVertexCount,std::runtime_error,
                         "In call to STK_Interface::getElementVertices all elements "
                         "must have the same vertex count!");

      // loop over all element nodes
      stk_classic::mesh::PairIterRelation nodes = element->relations(getNodeRank());
      TEUCHOS_TEST_FOR_EXCEPTION(nodes.size()!=masterVertexCount,std::runtime_error,
                         "In call to STK_Interface::getElementVertices cardinality of "
                         "element node relations must be the vertex count!");
      for(std::size_t node=0;node<nodes.size();++node) {
         const double * coord = getNodeCoordinates(nodes[node].entity()->identifier());

         // set each dimension of the coordinate
         for(unsigned d=0;d<dim;d++)
            vertices(cell,node,d) = coord[d];
      }
   }
}

template <typename ArrayT>
void STK_Interface::getElementVertices_FromField(const std::vector<stk_classic::mesh::Entity*> & elements,const std::string & eBlock, ArrayT & vertices) const
{
   TEUCHOS_ASSERT(useFieldCoordinates_);

   // nothing to do! silently return
   if(elements.size()==0) {
      vertices.resize(0,0,0);
      return;
   }

   // get *master* cell toplogy...(belongs to first element)
   unsigned masterVertexCount 
      = stk_classic::mesh::fem::get_cell_topology(elements[0]->bucket()).getCellTopologyData()->vertex_count;

   // allocate space
   vertices.resize(elements.size(),masterVertexCount,getDimension());

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
      stk_classic::mesh::Entity * element = elements[cell];

      // loop over nodes set solution values
      stk_classic::mesh::PairIterRelation relations = element->relations(getNodeRank());
      for(std::size_t i=0;i<relations.size();++i) {
         stk_classic::mesh::Entity * node = relations[i].entity();

         const double * coord = getNodeCoordinates(node->identifier());

         for(int d=0;d<getDimension();d++) {
           double * solnData = stk_classic::mesh::field_data(*fields[d],*node);
 
           // recall mesh field coordinates are stored as displacements
           // from the mesh coordinates, make sure to add them together
           vertices(cell,i,d) = solnData[0]+coord[d]; 
         }
      }
   }
}
 
}

#endif

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

#ifdef HAVE_IOSS
#include <stk_io/MeshReadWriteUtils.hpp>
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
   typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;
   typedef stk::mesh::Field<ProcIdData> ProcIdFieldType;
   typedef stk::mesh::Field<std::size_t> LocalIdFieldType;

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

   //////////////////////////////////////////

   /** Initialize the mesh with the current dimension This also calls
     * commit on the meta data causing it to be frozen. Information
     * about elements blocks has to be commited before this. If 
     * parallel machine has already been specified through <code>instantiateBulkData</code>
     * that communicator is used. Otherwise a new copy is constructed and
     * will be used through out this mesh object's lifetime.
     */
   void initialize(stk::ParallelMachine parallelMach,bool setupIO=true);

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

   /** Addes an entity to a specified side set.
     */
   void addEntityToSideset(stk::mesh::Entity & entity,stk::mesh::Part * sideset);

   /** Addes an entity to a specified node set.
     */
   void addEntityToNodeset(stk::mesh::Entity & entity,stk::mesh::Part * nodeset);

   // Methods to interrogate the mesh topology and structure
   //////////////////////////////////////////

   /** Grab the coordinates field 
     */
   const VectorFieldType & getCoordinatesField() const
   { return *coordinatesField_; }

   /** Look up a global node and get the coordinate.
     */
   const double * getNodeCoordinates(stk::mesh::EntityId nodeId) const;

   /** Look up a global node and get the coordinate.
     */
   const double * getNodeCoordinates(stk::mesh::Entity * node) const;

   /** Get subcell global IDs
     */
   void getSubcellIndices(unsigned entityRank,stk::mesh::EntityId elementId,
                          std::vector<stk::mesh::EntityId> & subcellIds) const;

   /** Get a vector of elements owned by this processor
     */
   void getMyElements(std::vector<stk::mesh::Entity*> & elements) const;

   /** Get a vector of elements owned by this processor on a particular block ID
     */
   void getMyElements(const std::string & blockID,std::vector<stk::mesh::Entity*> & elements) const;

   /** Get Entities corresponding to the side set requested. 
     * The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMySides(const std::string & sideName,std::vector<stk::mesh::Entity*> & sides) const;

   /** Get Entities corresponding to the side set requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be a dimension
     * lower then <code>getDimension()</code>.
     *
     * \param[in] sideName Name of side set
     * \param[in] blockName Name of block
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMySides(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity*> & sides) const;

   /** Get Entities corresponding to the node set requested. This also limits the entities
     * to be in a particular element block. The Entites in the vector should be ofdimension
     * <code>0</code>.
     *
     * \param[in] nodesetName Name of side set
     * \param[in] blockName Name of block
     * \param[in,out] sides Vector of entities containing the requested sides.
     */
   void getMyNodes(const std::string & sideName,const std::string & blockName,std::vector<stk::mesh::Entity*> & nodes) const;

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
 
   Teuchos::RCP<stk::mesh::BulkData> getBulkData() const { return bulkData_; }
   Teuchos::RCP<stk::mesh::fem::FEMMetaData> getMetaData() const { return metaData_; }

   bool isWritable() const;

   bool isModifiable() const
   {  if(bulkData_==Teuchos::null) return false; 
      return bulkData_->synchronized_state()==stk::mesh::BulkData::MODIFIABLE; }

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

   //! get the block count
   stk::mesh::Part * getElementBlockPart(const std::string & name) const
   { 
      std::map<std::string, stk::mesh::Part*>::const_iterator itr = elementBlocks_.find(name);   // Element blocks
      if(itr==elementBlocks_.end()) return 0;
      return itr->second; 
   }

   //! get the side set count
   std::size_t getNumSidesets() const
   { return sidesets_.size(); }

   stk::mesh::Part * getSideset(const std::string & name) const
   { return sidesets_.find(name)->second; }

   //! get the side set count
   std::size_t getNumNodesets() const
   { return nodesets_.size(); }

   stk::mesh::Part * getNodeset(const std::string & name) const
   { return nodesets_.find(name)->second; }

   //! get the global counts for the entity of specified rank
   std::size_t getEntityCounts(unsigned entityRank) const;

   //! get max entity ID of type entityRank
   stk::mesh::EntityId getMaxEntityId(unsigned entityRank) const;

   // Utilities
   //////////////////////////////////////////

   //! get a set of elements sharing a single node
   void getElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity *> & elements) const;

   /** Get set of element sharing a single node and its local node id.
     */
   void getOwnedElementsSharingNode(stk::mesh::Entity * node,std::vector<stk::mesh::Entity *> & elements,
                                                        std::vector<int> & localNodeId) const;

   /** Get set of element sharing a single node and its local node id.
     */
   void getOwnedElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity *> & elements,
                                                          std::vector<int> & localNodeId) const;


   //! get a set of elements sharing multiple nodes
   void getElementsSharingNodes(const std::vector<stk::mesh::EntityId> nodeId,std::vector<stk::mesh::Entity *> & elements) const;

   //! force the mesh to build subcells: edges and faces
   void buildSubcells();

   /** Get an elements local index
     */
   std::size_t elementLocalId(stk::mesh::Entity * elmt) const;

   /** Get an elements local index
     */
   std::size_t elementLocalId(stk::mesh::EntityId gid) const;

   /** Get an elements global index
     */
   inline stk::mesh::EntityId elementGlobalId(std::size_t lid) const
   { return (*orderedElementVector_)[lid]->identifier(); }

   /** Get an elements global index
     */
   inline stk::mesh::EntityId elementGlobalId(stk::mesh::Entity * elmt) const
   { return elmt->identifier(); }

   /**  Get the containing block ID of this element.
     */ 
   std::string containingBlockId(stk::mesh::Entity * elmt);

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

   ProcIdFieldType * getProcessorIdField() { return processorIdField_; }

   //! Has <code>initialize</code> been called on this mesh object?
   bool isInitialized() const { return initialized_; }

   /** Get Vector of element entities ordered by their LID, returns an RCP so that
     * it is easily stored by the caller.
     */
   Teuchos::RCP<const std::vector<stk::mesh::Entity*> > getElementsOrderedByLID() const;

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
                             const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues);

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
                         const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues);

   /** Get vertices associated with a number of elements of the same geometry.
     *
     * \param[in] localIds Element local IDs to construct vertices
     * \param[out] vertices Output array that will be sized (<code>localIds.size()</code>,#Vertices,#Dim)
     *
     * \note If not all elements have the same number of vertices an exception is thrown.
     *       If the size of <code>localIds</code> is 0, the function will silently return
     */
   template <typename ArrayT>
   void getElementVertices(std::vector<std::size_t> & localIds, ArrayT & vertices) const;

   // const stk::mesh::fem::FEMInterface & getFEMInterface() const 
   // { return *femPtr_; }

   const stk::mesh::EntityRank getElementRank() const { return metaData_->element_rank(); }
   const stk::mesh::EntityRank getSideRank() const { return metaData_->side_rank(); }
   const stk::mesh::EntityRank getFaceRank() const { return metaData_->face_rank(); }
   const stk::mesh::EntityRank getEdgeRank() const { return metaData_->edge_rank(); }
   const stk::mesh::EntityRank getNodeRank() const { return metaData_->node_rank(); }

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

   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > 
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
                             stk::mesh::EntityRank rank,bool setupIO);

   /** Build a safely handled Teuchos MPI communicator from a parallel machine.
     * This object asserts ownership of the communicator so that we can gurantee
     * existence so the outer user can do whatever they'd like with the original.
     */
   Teuchos::RCP<Teuchos::MpiComm<int> > getSafeCommunicator(stk::ParallelMachine parallelMach) const;

   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > periodicBCs_;

   Teuchos::RCP<stk::mesh::fem::FEMMetaData> metaData_;
   Teuchos::RCP<stk::mesh::BulkData> bulkData_;

   std::map<std::string, stk::mesh::Part*> elementBlocks_;   // Element blocks
   std::map<std::string, stk::mesh::Part*> sidesets_; // Side sets 
   std::map<std::string, stk::mesh::Part*> nodesets_; // Node sets 
   std::map<std::string, Teuchos::RCP<shards::CellTopology> > elementBlockCT_;

   // for storing/accessing nodes
   stk::mesh::Part * nodesPart_;
   std::vector<stk::mesh::Part*> nodesPartVec_;
   stk::mesh::Part * edgesPart_;
   std::vector<stk::mesh::Part*> edgesPartVec_;

   VectorFieldType * coordinatesField_;
   ProcIdFieldType * processorIdField_;
   LocalIdFieldType * localIdField_;
   
   // maps field names to solution field stk mesh handles
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToSolution_;
   std::map<std::pair<std::string,std::string>,SolutionFieldType*> fieldNameToCellField_;

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

#ifdef HAVE_IOSS
   // I/O support
   Teuchos::RCP<stk::io::MeshData> meshData_;
#endif

   // uses lazy evaluation
   mutable Teuchos::RCP<std::vector<stk::mesh::Entity*> > orderedElementVector_;

   // Object describing how to sort a vector of elements using
   // local ID as the key, very short lived object
   class LocalIdCompare {
   public:
     LocalIdCompare(const STK_Interface * mesh) : mesh_(mesh) {}
   
     // Compares two stk mesh entities based on local ID
     bool operator() (stk::mesh::Entity * a,stk::mesh::Entity * b) 
     { return mesh_->elementLocalId(a) < mesh_->elementLocalId(b);}
   
   private:
     const STK_Interface * mesh_;
   };
};

template <typename ArrayT>
void STK_Interface::setSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                                         const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues)
{
   const std::vector<stk::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity * element = elements[localId];

      // loop over nodes set solution values
      stk::mesh::PairIterRelation relations = element->relations(getNodeRank());
      for(std::size_t i=0;i<relations.size();++i) {
         stk::mesh::Entity * node = relations[i].entity();

         double * solnData = stk::mesh::field_data(*field,*node);
         // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
         solnData[0] = solutionValues(cell,i);
      }
   }
}

template <typename ArrayT>
void STK_Interface::getSolutionFieldData(const std::string & fieldName,const std::string & blockId,
                                         const std::vector<std::size_t> & localElementIds,ArrayT & solutionValues) const
{
   const std::vector<stk::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   solutionValues.resize(localElementIds.size(),elements[localElementIds[0]]->relations(getNodeRank()).size());

   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity * element = elements[localId];

      // loop over nodes set solution values
      stk::mesh::PairIterRelation relations = element->relations(getNodeRank());
      for(std::size_t i=0;i<relations.size();++i) {
         stk::mesh::Entity * node = relations[i].entity();

         double * solnData = stk::mesh::field_data(*field,*node);
         // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
         solutionValues(cell,i) = solnData[0]; 
      }
   }
}

template <typename ArrayT>
void STK_Interface::setCellFieldData(const std::string & fieldName,const std::string & blockId,
                                     const std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues)
{
   const std::vector<stk::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   SolutionFieldType * field = this->getCellField(fieldName,blockId);

   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity * element = elements[localId];

      double * solnData = stk::mesh::field_data(*field,*element);
      TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
      solnData[0] = solutionValues[cell];
   }
}

template <typename ArrayT>
void STK_Interface::getElementVertices(std::vector<std::size_t> & localElementIds, ArrayT & vertices) const
{

   // nothing to do! silently return
   if(localElementIds.size()==0) {
      vertices.resize(0,0,0);
      return;
   }

   const std::vector<stk::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   // get *master* cell toplogy...(belongs to first element)
   unsigned masterVertexCount 
      = stk::mesh::fem::get_cell_topology(elements[localElementIds[0]]->bucket()).getCellTopologyData()->vertex_count;

   // allocate space
   vertices.resize(localElementIds.size(),masterVertexCount,getDimension());

   // loop over each requested element
   unsigned dim = getDimension();
   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      stk::mesh::Entity * element = elements[localElementIds[cell]];
      TEUCHOS_ASSERT(element!=0);
 
      unsigned vertexCount 
         = stk::mesh::fem::get_cell_topology(element->bucket()).getCellTopologyData()->vertex_count;
      TEUCHOS_TEST_FOR_EXCEPTION(vertexCount!=masterVertexCount,std::runtime_error,
                         "In call to STK_Interface::getElementVertices all elements "
                         "must have the same vertex count!");

      // loop over all element nodes
      stk::mesh::PairIterRelation nodes = element->relations(getNodeRank());
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
 
}

#endif

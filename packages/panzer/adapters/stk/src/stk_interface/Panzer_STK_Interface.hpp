#ifndef __Panzer_STK_Interface_hpp__
#define __Panzer_STK_Interface_hpp__

#include <Teuchos_RCP.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <Shards_CellTopology.hpp>

#include <Panzer_STK_config.hpp>

namespace panzer_stk {

/** Pure virtial base class that builds a basic element. To be
  * overidden with several types of elements.
  */ 
class ElementDescriptor {
public:
   ElementDescriptor(stk::mesh::EntityId gid,const std::vector<stk::mesh::EntityId> & nodes);
   virtual ~ElementDescriptor();

   /** Function adds element and its relations to the bulk data. This function
     * will not be overriden by other classes. Its function is common to all elements.
     */ 
   void addOutlineToBulkData(stk::mesh::BulkData & bulkData,const std::vector<stk::mesh::Part*> & parts);

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
   /** Default constructor
     */
   STK_Interface();

   // functions called before initialize
   //////////////////////////////////////////

   /** Set the dimension for this mesh
     */
   void setDimension(unsigned dim)
   { dimension_ = dim; }

   /** Add an element block with a string name
     */
   template <typename TopologyTYpe>
   inline void addElementBlock(const std::string & name);

   /** Add a side set with a string name
     */
   void addSideset(const std::string & name);

   /** Add a solution field
     */ 
   void addSolutionField(const std::string & fieldName,const std::string & blockId);

   //////////////////////////////////////////

   /** Initialize the mesh with the current dimension This also calls
     * commit on the meta data causing it to be frozen. Information
     * about elements blocks has to be commited before this.
     */
   void initialize(stk::ParallelMachine parallelMach);

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

   void addElement(Teuchos::RCP<ElementDescriptor> & ed,stk::mesh::Part * block);

   /** Addes an entity to a specified side set.
     */
   void addEntityToSideset(stk::mesh::Entity & entity,stk::mesh::Part * sideset);

   // Methods to interrogate the mesh topology and structure
   //////////////////////////////////////////

   /** Look up a global node and get the coordinate.
     */
   const double * getNodeCoordinates(stk::mesh::EntityId nodeId) const;

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

   // Utility functions
   //////////////////////////////////////////

   /** Write this mesh to exodus
     */
   void writeToExodus(const std::string & filename);

   // Accessor functions
   //////////////////////////////////////////

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

   //! get the block count
   stk::mesh::Part * getElementBlockPart(const std::string & name) const
   { 
      std::map<std::string, stk::mesh::Part*>::const_iterator itr = elementBlocks_.find(name);   // Element blocks
      if(itr==elementBlocks_.end()) return 0;
      return elementBlocks_.find(name)->second; 
   }

   //! get the side set count
   std::size_t getNumSidesets() const
   { return sidesets_.size(); }

   stk::mesh::Part * getSideset(const std::string & name) const
   { return sidesets_.find(name)->second; }

   //! get the global counts for the entity of specified rank
   std::size_t getEntityCounts(unsigned entityRank) const;

   //! get max entity ID of type entityRank
   stk::mesh::EntityId getMaxEntityId(unsigned entityRank) const;

   // Utilities
   //////////////////////////////////////////

   //! get a set of elements sharing a single node
   void getElementsSharingNode(stk::mesh::EntityId nodeId,std::vector<stk::mesh::Entity *> & elements) const;

   //! get a set of elements sharing multiple nodes
   void getElementsSharingNodes(const std::vector<stk::mesh::EntityId> nodeId,std::vector<stk::mesh::Entity *> & elements) const;

   //! force the mesh to build the subcells of a particular rank
   void buildSubcells(unsigned subcellRank);

   /** Use a formula to determine the edge ID.  This formula is kind of bad
     * and should be viewed as a short term fix.
     */
   stk::mesh::EntityId getEdgeId(stk::mesh::EntityId n0,stk::mesh::EntityId n1) const;

   /** Get an elements local index
     */
   std::size_t elementLocalId(stk::mesh::Entity * elmt) const;

   /**  Get the containing block ID of this element.
     */ 
   std::string containingBlockId(stk::mesh::Entity * elmt);

   /** Get the stk mesh field pointer associated with a particular solution value
     * Assumes there is a field associated with "fieldName,blockId" pair. If none
     * is found an exception (std::runtime_error) is raised.
     */
   stk::mesh::Field<double> * getSolutionField(const std::string & fieldName,
                                               const std::string & blockId) const;

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
                             std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues);

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

   /** Add an edge and its relations to this cell
     */
   void addEdges_local(stk::mesh::Entity * cell);

   /** Add an edge with a given identifier to connect two nodes. If the edge
     * is already in the mesh then the entity is simply returned.
     * 
     * \param[in] n0 Node defining edge
     * \param[in] n1 Node defining edge
     * \param[in] edgeId Global identifier for the edge
     */
   stk::mesh::Entity * addEdge(stk::mesh::Entity * n0,stk::mesh::Entity * n1, stk::mesh::EntityId edgeId);

   typedef stk::mesh::Field<double> SolutionFieldType;
   typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;
   typedef stk::mesh::Field<int> ProcIdFieldType;
   typedef stk::mesh::Field<std::size_t> LocalIdFieldType;

   Teuchos::RCP<stk::mesh::MetaData> metaData_;
   Teuchos::RCP<stk::mesh::BulkData> bulkData_;

   std::map<std::string, stk::mesh::Part*> elementBlocks_;   // Element blocks
   std::map<std::string, stk::mesh::Part*> sidesets_; // Side sets 
   stk::mesh::Part * sidesetsPart_;
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

   unsigned dimension_;

   bool initialized_;

   // how many elements, faces, edges, and nodes are there globally
   std::vector<std::size_t> entityCounts_;

   // what is maximum entity ID
   std::vector<stk::mesh::EntityId> maxEntityId_;

   int procRank_;
   std::size_t currentLocalId_;

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

template <typename TopologyType> 
void STK_Interface::addElementBlock(const std::string & name)
{
   TEUCHOS_ASSERT(not initialized_);

   stk::mesh::Part * block = &metaData_->declare_part(name,stk::mesh::Element);
   stk::mesh::set_cell_topology(*block,shards::getCellTopologyData<TopologyType>());

   // construct cell topology object for this block
   Teuchos::RCP<shards::CellTopology> ct
         = Teuchos::rcp(new shards::CellTopology(stk::mesh::get_cell_topology(*block)));

   // add element block part and cell topology
   elementBlocks_.insert(std::make_pair(name,block));
   elementBlockCT_.insert(std::make_pair(name,ct));
}

template <typename ArrayT>
void STK_Interface::setSolutionFieldData(const std::string & fieldName,const std::string & blockId,std::vector<std::size_t> & localElementIds,const ArrayT & solutionValues)
{
   const std::vector<stk::mesh::Entity*> & elements = *(this->getElementsOrderedByLID());

   // SolutionFieldType * field = metaData_->get_field<SolutionFieldType>(fieldName); // if no blockId is specified you can get the field like this!
   SolutionFieldType * field = this->getSolutionField(fieldName,blockId);

   std::vector<std::size_t>::const_iterator itr;
   for(std::size_t cell=0;cell<localElementIds.size();cell++) {
      std::size_t localId = localElementIds[cell];
      stk::mesh::Entity * element = elements[localId];

      // loop over nodes set solution values
      stk::mesh::PairIterRelation relations = element->relations();
      for(std::size_t i=0;i<relations.size();++i) {
         stk::mesh::Entity * node = relations[i].entity();

         double * solnData = stk::mesh::field_data(*field,*node);
         // TEUCHOS_ASSERT(solnData!=0); // only needed if blockId is not specified
         solnData[0] = solutionValues(cell,i);
      }
   }
}
 
}

#endif

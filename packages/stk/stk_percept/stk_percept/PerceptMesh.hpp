#ifndef stk_percept_PerceptMesh_hpp
#define stk_percept_PerceptMesh_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <set>

#include <stk_percept/function/Function.hpp>
#include <stk_percept/Name.hpp>

#include "ShardsInterfaceTable.hpp"

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <stk_io/IossBridge.hpp>

#include <Intrepid_FieldContainer.hpp>


#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>

#include "Teuchos_RCP.hpp"

//#include "PerceptMeshReadWrite.hpp"
#include <stk_percept/function/ElementOp.hpp>
#include <stk_percept/function/BucketOp.hpp>

#include <stk_percept/SameRankRelation.hpp>

#include <stk_percept/function/internal/SimpleSearcher.hpp>

// if this is set, use stk_mesh relations to hold parent/child information, else use special data structures for this
#define PERCEPT_USE_FAMILY_TREE 1

// if this is set, extra high-rank entities are allowed to be created for use in generating shared nodes (see stk_adapt/NodeRegistry)
#define PERCEPT_USE_PSEUDO_ELEMENTS 0

//using namespace shards;

namespace Intrepid {
	template<class Scalar, class ArrayScalar>
	class Basis;
}

namespace stk {
  namespace percept {

    typedef mesh::Field<double>                          ScalarFieldType ;
    typedef mesh::Field<double, stk::mesh::Cartesian>    VectorFieldType ;

    static const unsigned EntityRankEnd = 6;

    enum FamiltyTreeLevel {
      FAMILY_TREE_LEVEL_0 = 0,
      FAMILY_TREE_LEVEL_1 = 1
    };

    enum FamiltyTreeParentIndex {
      FAMILY_TREE_PARENT = 0,
      FAMILY_TREE_CHILD_START_INDEX = 1
    };

    //using namespace interface_table;

    class GMeshSpec : public Name
    {
    public:
      explicit GMeshSpec(const std::string& name) : Name(name) {}
    };


    struct FieldCreateOrder
    {
      const std::string m_name;
      const unsigned m_entity_rank;
      const std::vector<int> m_dimensions;
      const mesh::Part* m_part;
      FieldCreateOrder();
      FieldCreateOrder(const std::string name, const unsigned entity_rank,
                      const std::vector<int> dimensions, const mesh::Part* part);
    };
    typedef std::vector<FieldCreateOrder> FieldCreateOrderVec;

    class PerceptMesh
    {
    public:
      typedef Intrepid::Basis<double, MDArray > BasisType;
      typedef Teuchos::RCP<BasisType>           BasisTypeRCP;
      typedef std::map<unsigned, BasisTypeRCP > BasisTableMap;

      static std::string s_omit_part;
      
    public:

      //========================================================================================================================
      /// high-level interface

      // ctor constructor
      /// Create a Mesh object that owns its constituent FEMMetaData and BulkData (which are created by this object)
      //PerceptMesh( stk::ParallelMachine comm =  MPI_COMM_WORLD );
      PerceptMesh(size_t spatialDimension = 3u, stk::ParallelMachine comm =  MPI_COMM_WORLD );

      /// Create a Mesh object that doesn't own its constituent FEMMetaData and BulkData, pointers to which are adopted
      /// by this constructor.
      PerceptMesh(const stk::mesh::fem::FEMMetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted=true);

      /// reads and commits mesh, editing disabled
      void
      openReadOnly(const std::string& in_filename);

      /// reads but doesn't commit mesh, enabling edit
      void
      open(const std::string& in_filename);

      /// opens an empty mesh, with a commit
      void
      openEmpty();

      /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec, Read Only mode, no edits allowed
      void
      newMeshReadOnly(const GMeshSpec gmesh_spec);

      /// creates a new mesh using the GeneratedMesh fixture with spec @param gmesh_spec
      void
      newMesh(const GMeshSpec gmesh_spec);

      /// add a field to the mesh
      stk::mesh::FieldBase *
      addField(const std::string& name, const unsigned entity_rank, int vectorDimension=0, const std::string part_name="universal_part");

      stk::mesh::FieldBase *
      getField(const std::string& name);

      /// commits mesh  - any operations done on a non-committed mesh, except to add fields will throw an exception
      void
      commit();

      /// reopens the mesh for editing - warning, this operation writes the mesh to a temp file then re-reads it and
      /// thus recreates the internal FEMMetaData and BulkData
      void
      reopen(const std::string temp_file_name="percept_tmp.e");

      /// commits mesh if not committed and saves it in new file
      void
      saveAs(const std::string& out_filename);

      /// closes this mesh to further changes
      void
      close();

      /// print number of parts and fields, and info on each
      void
      printInfo(std::ostream& stream, std::string header="", int print_level=0, bool do_endl=true);

      /// print number of parts and fields, and info on each
      void
      printInfo(std::string header="", int print_level = 0, bool do_endl=true);

      void
      printFields(std::string header="");

      int
      getSpatialDim();

      int
      getNumberElements();

      int
      getNumberNodes();

      int
      getNumberEdges();

      int
      getNumberElementsLocallyOwned();

      unsigned getRank() { return getBulkData()->parallel_rank(); }
      unsigned getParallelRank() { return getBulkData()->parallel_rank(); }
      unsigned getParallelSize() { return getBulkData()->parallel_size(); }

      void printEntity(const stk::mesh::Entity& entity, stk::mesh::FieldBase* field=0) { printEntity(std::cout, entity, field); };
      std::string printEntityCompact(const stk::mesh::Entity& entity, stk::mesh::FieldBase* field=0);

      stk::mesh::BulkData * getBulkData();
      stk::mesh::fem::FEMMetaData * getFEM_meta_data();

      mesh::Part* getNonConstPart(const std::string& part_name);

      // entity data setter/getters
      double get_field_data(const stk::mesh::FieldBase *field, const mesh::Entity* entity, unsigned ordinal=0)
      {
        double *val= field_data_entity(field,*entity);
        return val ? val[ordinal] : 0.0;
      }
      void set_field_data(double value, const stk::mesh::FieldBase *field, const mesh::Entity* entity, unsigned ordinal=0 )
      {
        double *val= field_data_entity(field,*entity);
        if( val) val[ordinal] = value;
      }

      double get_node_field_data(stk::mesh::FieldBase *field, const mesh::EntityId node_id, unsigned ordinal=0)
      {
        double *val = node_field_data(field, node_id);
        return val ? val[ordinal] : 0.0;
      }
      void set_node_field_data(double value, stk::mesh::FieldBase *field, const mesh::EntityId node_id, unsigned ordinal=0)
      {
        double *val = node_field_data(field, node_id);
        if( val) val[ordinal] = value;
      }

      stk::mesh::Entity *get_node(const mesh::EntityId node_id) 
      {
        return getBulkData()->get_entity(node_rank(), node_id);
      }

      stk::mesh::Entity *get_element(const mesh::EntityId element_id) 
      {
        return getBulkData()->get_entity(element_rank(), element_id);
      }

      stk::mesh::Entity *get_entity( mesh::EntityRank rank, const mesh::EntityId id) 
      {
        return getBulkData()->get_entity(rank, id);
      }

      /// find node closest to given point
      stk::mesh::Entity *get_node(double x, double y, double z=0, double t=0) ;

      /// find element that contains given point
      stk::mesh::Entity *get_element(double x, double y, double z=0, double t=0) ;

      //========================================================================================================================
      /// low-level interfaces

      //PerceptMesh(const stk::mesh::MetaData* metaData, stk::mesh::BulkData* bulkData, bool isCommitted=true);

      ~PerceptMesh() ;
      void init( stk::ParallelMachine comm  =  MPI_COMM_WORLD, bool no_alloc=false );      // FIXME - make private
      void destroy();       // FIXME - make private

      const mesh::Part*
      getPart(const std::string& part_name) ;

      void printEntity(std::ostream& out, const stk::mesh::Entity& entity, stk::mesh::FieldBase* field=0);

      // allow setting spatial dim after creation (for compatability with new FEMMetaData)
      void setSpatialDim(int sd);

      // streaming refine mesh
      void setStreamingSize(int streaming_size) { m_streaming_size= streaming_size; }
      int getStreamingSize() { return m_streaming_size; }

      /// reads the given file into a temporary model and prints info about it
      void dump(const std::string& file="");
      void dumpElements(const std::string& partName = "");
      void dumpElementsCompact();

      bool isGhostElement(const stk::mesh::Entity& element)
      {
        //throw std::runtime_error("not impl"); // FIXME
        bool isGhost = element.owner_rank() != getRank();
        return isGhost;
      }

      // checks if this entity has a duplicate (ie all nodes are the same)
      bool
      check_entity_duplicate(stk::mesh::Entity& entity);

      void delete_side_sets();

      /// add some fields that are useful for debugging or for exporting meshes to Mesquite - must be
      ///   done before commit()
      void addParallelInfoFields(bool elemental, bool nodal, 
                                 std::string elemental_proc_rank_name = "proc_rank",
                                 std::string nodal_fixed_flag="fixed", // boundary flag for telling Mesquite these nodes shouldn't be moved
                                 std::string nodal_global_id_name="GLOBAL_ID", 
                                 std::string nodal_proc_id_name="PROCESSOR_ID", 
                                 std::string nodal_local_id_name="LOCAL_ID");

      /// fill the fields from addParallelInfoFields with data from stk_mesh database
      void populateParallelInfoFields(bool elemental, bool nodal, 
                                      stk::mesh::Selector* fixed_node_selector=0,
                                      std::string elemental_proc_rank_name = "proc_rank",
                                      std::string nodal_fixed_flag="fixed", // boundary flag for telling Mesquite these nodes shouldn't be moved
                                      std::string nodal_global_id_name="GLOBAL_ID", 
                                      std::string nodal_proc_id_name="PROCESSOR_ID", 
                                      std::string nodal_local_id_name="LOCAL_ID");

      /**
       * A family tree relation holds the parent/child relations for a refined mesh.  
       * 
       * Case 0: a single refinement of a parent P_0 and its children C_0_0, C_0_1,...,C_0_N leads to a new
       *  family tree entity FT_0 that has down relations to {P_0, C_0_0, C_0_1,...,C_0_N}
       *  The back pointers from P_0, C_0_0, ... are initially stored as the 0'th index of their relations,
       *    i.e.: P_0.relations(FAMILY_TREE_RANK)[0] --> FT_0, 
       *          C_0_0.relations(FAMILY_TREE_RANK)[0] --> FT_0, etc.
       * Case 1: a previously refined child, say C_0_1, renamed to P_0_1, gets further refined leading to 
       *  a new family tree entity, FT_1 
       *  pointing to:  {P_0_1, C_0_1_0, C_0_1_1,... }
       *  but, now the relations indexing changes (actually, we can't predict what it will be, thus the
       *  need for this function getFamilyTreeRelationIndex):
       *     P_0_1.relations(FAMILY_TREE_RANK)[0] --> FT_1
       *     P_0_1.relations(FAMILY_TREE_RANK)[1] --> FT_0
       *     etc.
       *  So, we use this function to look for the family tree corresponding to if we are looking for the first
       *    level (if there's only one level, or we are looking for the family tree associated with the element
       *    when it was a child for the first time), orthe "level 1" family tree (corresponding to Case 1
       *    where we are looking for the family tree of the element associated with it being a parent).
       * 
       */
      unsigned getFamilyTreeRelationIndex(FamiltyTreeLevel level, const stk::mesh::Entity& element);

      /// the element is not a parent of the 0'th family_tree relation
      bool isChildElement( const stk::mesh::Entity& element, bool check_for_family_tree=true);

      // either has no family tree or is a child
      bool isLeafElement( const stk::mesh::Entity& element);

      /// the element is not a parent of any family tree relation
      bool isChildElementLeaf( const stk::mesh::Entity& element, bool check_for_family_tree=true);

      bool hasFamilyTree(const stk::mesh::Entity& element);

      /// if the element is a parent at any level, return true
      bool isParentElement( const stk::mesh::Entity& element, bool check_for_family_tree=true);

      /// is element a parent at the leaf level (either there is only one level, and it's a parent, or
      ///    if more than one, the element is a child and a parent and its children have no children)
      bool isParentElementLeaf( const stk::mesh::Entity& element, bool check_for_family_tree=true);

      /// is element a parent at level 2 (meaning that it is both a child and a parent)
      bool isParentElementLevel2( const stk::mesh::Entity& element, bool check_for_family_tree=true);

      /// is element a child with siblings with no nieces or nephews (siblings with children)
      ///  (alternative would be "is child and is parent not a grandparent")
      bool isChildWithoutNieces( const stk::mesh::Entity& element, bool check_for_family_tree=true);

      // return false if we couldn't get the children
      bool getChildren( const stk::mesh::Entity& element, std::vector<stk::mesh::Entity*>& children, bool check_for_family_tree=true, bool only_if_element_is_parent_leaf=false);

      void printParentChildInfo(const stk::mesh::Entity& element, bool check_for_family_tree=true);

      static inline
      stk::mesh::EntityRank fem_entity_rank( unsigned int t ) {
        return  t < EntityRankEnd ? stk::mesh::EntityRank(t) : stk::mesh::InvalidEntityRank ;
      }

      stk::mesh::EntityRank node_rank() const
      {
        return m_metaData->node_rank();
      }

      /** \brief Returns the edge rank which changes depending on spatial dimension
       */
      stk::mesh::EntityRank edge_rank() const
      {
        return m_metaData->edge_rank();
      }

      /** \brief Returns the face rank which changes depending on spatial dimension
       */
      stk::mesh::EntityRank face_rank() const
      {
        return m_metaData->face_rank();
      }

      /** \brief Returns the side rank which changes depending on spatial dimension
       */
      stk::mesh::EntityRank side_rank() const
      {
        return m_metaData->side_rank();
      }

      /** \brief Returns the element rank which is always equal to spatial dimension
       */
      stk::mesh::EntityRank element_rank() const
      {
        return m_metaData->element_rank();
      }

      stk::mesh::Entity & createOrGetNode(stk::mesh::EntityId nid, double* x=0);

      void createEntities(stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity *>& requested_entities);


      static double * field_data(const stk::mesh::FieldBase *field, const stk::mesh::Bucket & bucket, unsigned *stride=0);
      static double * field_data(const stk::mesh::FieldBase *field, const mesh::Entity& node, unsigned *stride=0);
      static double * field_data_entity(const stk::mesh::FieldBase *field, const mesh::Entity& entity, unsigned *stride=0);

      static inline double *
      field_data_inlined(const mesh::FieldBase *field, const mesh::Entity& node)
      {
        return field_data(field, node);
//         return
//           field->rank() == 0 ?
//           stk::mesh::field_data( *static_cast<const ScalarFieldType *>(field) , node )
//           :
//           stk::mesh::field_data( *static_cast<const VectorFieldType *>(field) , node );
      }


      double * node_field_data(stk::mesh::FieldBase *field, const mesh::EntityId node_id);

      static BasisTypeRCP getBasis(shards::CellTopology& topo);
      static void setupBasisTable();

      void nodalOpLoop(GenericFunction& nodalOp, stk::mesh::FieldBase *field=0);
      void elementOpLoop(ElementOp& elementOp, stk::mesh::FieldBase *field=0, stk::mesh::Part *part = 0);
      void bucketOpLoop(BucketOp& bucketOp, stk::mesh::FieldBase *field=0, stk::mesh::Part *part = 0);

      static unsigned size1(const stk::mesh::Bucket& bucket) { return bucket.size(); }
      static unsigned size1(const stk::mesh::Entity& element) { return 1; }

      /// \brief Fill the array cellNodes(numCells, numNodesPerCell, nDof) with DOF values from the given Field
      /// The stride of the data (last dimension in cellNodes) is taken to be that of the field's stride; however,
      /// in some cases you may want to pass in an array where nDof is less than the stride (e.g., pull out 2
      /// coordinates from a 3D coordinate field).  In that case, the dataStride argument can be set (to e.g. "2").
      template<class ArrayType>
      static void fillCellNodes( const stk::mesh::Bucket &bucket,
                                 //stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                 //VectorFieldType& coord_field,
                                 mesh::FieldBase* field,
                                 ArrayType& cellNodes, unsigned dataStride=0 );

      /// \brief see comment for fillCellNodes(Bucket& ...)
      template<class ArrayType>
      static void fillCellNodes( const stk::mesh::Entity &element,
                                 //stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                 //VectorFieldType& coord_field,
                                 mesh::FieldBase* field,
                                 ArrayType& cellNodes, unsigned dataStride=0 );

      double edge_length_ave(const stk::mesh::Entity &entity);

      static void findMinMaxEdgeLength(const mesh::Bucket &bucket,  stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                       Intrepid::FieldContainer<double>& elem_min_edge_length, Intrepid::FieldContainer<double>& elem_max_edge_length);

      VectorFieldType* getCoordinatesField() {
        // this should have been set by a previous internal call to setCoordinatesField
        return m_coordinatesField;
      }

      static void
      element_side_nodes( const mesh::Entity & elem , int local_side_id, stk::mesh::EntityRank side_entity_rank, std::vector<mesh::Entity *>& side_node_entities );

      static void
      element_side_permutation(const mesh::Entity& element, const mesh::Entity& side, unsigned iSubDimOrd, int& returnedIndex, int& returnedPolarity);

      // FIXME
      SameRankRelation& adapt_parent_to_child_relations() { return m_adapt_parent_to_child_relations; }

      bool
      isBoundarySurface(mesh::Part& block, mesh::Part& surface);

      /// here @param thing is a Part, Bucket, Entity, or Field or BulkData
      template<class T>
      static
      const stk::mesh::fem::FEMMetaData& get_fem_meta_data(const T& thing)
      {
        //const stk::mesh::fem::FEMMetaData& meta = stk::mesh::fem::FEMMetaData::get(thing);
        const stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get ( thing );
        return fem_meta;
      }

      /// here @param thing is a Part, Bucket, Entity
      template<class T>
      static
      const CellTopologyData * get_cell_topology(const T& thing)
      {
        const CellTopologyData * cell_topo_data = mesh::fem::get_cell_topology(thing).getCellTopologyData();
        return cell_topo_data;
      }

      const stk::mesh::PartVector& get_io_omitted_parts() { return m_io_omitted_parts; }
      void set_io_omitted_parts(stk::mesh::PartVector& io_omitted_parts) { m_io_omitted_parts = io_omitted_parts; }

      static bool mesh_difference(PerceptMesh& mesh1, PerceptMesh& mesh2,
                                  std::string msg,
                                  bool print=true, bool print_all_field_diffs=false);

      static bool mesh_difference(stk::mesh::fem::FEMMetaData& metaData_1,
                                  stk::mesh::fem::FEMMetaData& metaData_2,
                                  stk::mesh::BulkData& bulkData_1,
                                  stk::mesh::BulkData& bulkData_2,
                                  std::string msg,
                                  bool print=true, bool print_all_field_diffs=false);


    private:

      /// reads meta data, commits it, reads bulk data
      void readModel( const std::string& in_filename );

      /// read with no commit
      void read_metaDataNoCommit( const std::string& in_filename );

      /// create with no commit
      void create_metaDataNoCommit( const std::string& gmesh_spec);

      void commit_metaData();

      /// read the bulk data (no op in create mode)
      void readBulkData();

      /// Convenience method to read a model's meta data, create some new fields, commit meta data then read the bulk data
      /// deprecated
      void readModelAndCreateOptionalFields(const std::string file, bool print,  FieldCreateOrderVec create_field);

      /// after the meta data is read or created, create some fields using this method, then you can commit and read bulk data(if in
      /// reading mode, else no need to read bulk data in create mode)
      // deprecated
      void createFields(bool print, FieldCreateOrderVec create_field = FieldCreateOrderVec());

      /// Cache internal pointer to coordinate field
      void setCoordinatesField();

      // look for omitted parts
      void checkForPartsToAvoidWriting();

      // write in exodus format to given file
      void writeModel( const std::string& out_filename );

      stk::mesh::FieldBase * createField(const std::string& name, const unsigned entity_rank, const std::vector<int>& dimensions,
                                         const stk::mesh::Part* arg_part=0);

      //static void transformMesh(GenericFunction& coordinate_transform);


    private:
      //stk::mesh::fem::FEMMetaData *         m_fem_meta_data;
      stk::mesh::fem::FEMMetaData *                 m_metaData;
      stk::mesh::BulkData *                 m_bulkData;
      stk::io::util::Gmesh_STKmesh_Fixture* m_fixture;
      //Teuchos::RCP<Ioss::Region>            m_iossRegion;
      Ioss::Region*                         m_iossRegion;
      Teuchos::RCP<stk::io::MeshData>       m_iossMeshData;
      VectorFieldType*                      m_coordinatesField;
      int                                   m_spatialDim;
      bool                                  m_ownData;
      bool                                  m_isCommitted;
      bool                                  m_isOpen;
      bool                                  m_isInitialized;
      bool                                  m_isAdopted;
      bool                                  m_dontCheckState;
      std::string                           m_filename;
      stk::ParallelMachine                  m_comm;

      //static std::map<unsigned, BasisType *> m_basisTable;
      static BasisTableMap m_basisTable;

      SameRankRelation m_adapt_parent_to_child_relations;

      stk::mesh::PartVector                 m_io_omitted_parts;

      // normally 0, unless using streaming refine
      int                                   m_streaming_size;

      Searcher *                            m_searcher;

      void checkStateSpec(const std::string& function, bool cond1=true, bool cond2=true, bool cond3=true);

      void checkState(const std::string& function) {
        return checkStateSpec(function, m_isOpen, m_isInitialized, m_isCommitted);
      }

    }; // class PerceptMesh


    template<>
    const CellTopologyData *
    PerceptMesh::get_cell_topology(const mesh::Part& part) ;



#if 0
		inline
		std::string &operator<<(std::string& out, const char *str)
		{
			return out.append(str);
		}
#endif

    // static
    template<class ArrayType>
    void PerceptMesh::fillCellNodes( const mesh::Bucket &bucket,
                                  //stk::mesh::Cartesian>& coord_field,
                                  //VectorFieldType& coord_field,
                                  mesh::FieldBase* field,
                                  ArrayType& cellNodes, unsigned dataStrideArg)
    {
      unsigned number_elems = bucket.size();
      //const CellTopologyData * const bucket_cell_topo_data = stk::mesh::fem::get_cell_topology(bucket).getCellTopologyData();
      const CellTopologyData * const bucket_cell_topo_data = get_cell_topology(bucket);

      shards::CellTopology cell_topo(bucket_cell_topo_data);
      //unsigned numCells = number_elems;
      unsigned numNodes = cell_topo.getNodeCount();
      //unsigned spaceDim = cell_topo.getDimension();

      unsigned dataStride = dataStrideArg;
      if (!dataStrideArg)
        {
          const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, mesh::fem::FEMMetaData::get(*field).universal_part());
          dataStride = r.dimension() ;
        }
      //std::cout << "bucket dataStride= " << dataStride << std::endl;

      for ( unsigned iElemInBucketOrd = 0 ; iElemInBucketOrd < number_elems ; ++iElemInBucketOrd)
        {
          mesh::Entity & elem = bucket[iElemInBucketOrd] ;

          if (0) std::cout << "elemOfBucket= " << elem << std::endl;
          const mesh::PairIterRelation elem_nodes = elem.relations( mesh::fem::FEMMetaData::NODE_RANK );

          // FIXME: fill field data (node coordinates)
          for (unsigned iNodeOrd = 0; iNodeOrd < numNodes; iNodeOrd++)
            {
              mesh::Entity& node = *elem_nodes[iNodeOrd].entity();
              double * node_coord_data = PerceptMesh::field_data( field , node);

              for (unsigned iDOFOrd = 0; iDOFOrd < dataStride; iDOFOrd++)
                {
                  cellNodes(iElemInBucketOrd, iNodeOrd, iDOFOrd) = node_coord_data[iDOFOrd];
                }
            }


        }

    }

    // static
    template<class ArrayType>
    void PerceptMesh::fillCellNodes( const stk::mesh::Entity &element,
                                  //stk::mesh::Field<double, stk::mesh::Cartesian>& coord_field,
                                  //VectorFieldType& coord_field,
                                  mesh::FieldBase* field,
                                  ArrayType& cellNodes,
                                  unsigned dataStrideArg )
    {
      unsigned dataStride = dataStrideArg;
      if (!dataStrideArg)
        {
          const stk::mesh::FieldBase::Restriction & r = field->restriction(stk::mesh::fem::FEMMetaData::NODE_RANK, mesh::fem::FEMMetaData::get(*field).universal_part());
          dataStride = r.dimension() ;
        }
      //std::cout << "element dataStride= " << dataStride << std::endl;
      const mesh::PairIterRelation element_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
      unsigned numNodes = element_nodes.size();

      unsigned iCell = 0;
      for (unsigned iNodeOrd = 0; iNodeOrd < numNodes; iNodeOrd++)
        {
          mesh::Entity& node = *element_nodes[iNodeOrd].entity();
          double * node_coord_data = PerceptMesh::field_data( field , node);

          for (unsigned iDOFOrd = 0; iDOFOrd < dataStride; iDOFOrd++)
            {
              cellNodes(iCell, iNodeOrd, iDOFOrd) = node_coord_data[iDOFOrd];
            }
        }
    }

  }
}

#endif

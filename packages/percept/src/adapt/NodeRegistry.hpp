// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/**
 * NodeRegistry: class that defines a map of key->value where key is a SubDimCell
 *   (e.g. an edge or face bewtween two elements), and value contains various data
 *   about the SubDimCell, such as which element owns the SubDimCell, an array of
 *   new nodes needed on it for refinement, etc.
 */

#ifndef adapt_NodeRegistry_hpp
#define adapt_NodeRegistry_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_mesh/base/EntityKey.hpp>

#include <percept/stk_mesh.hpp>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelComm.hpp>

#include <percept/NoMallocArray.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>

#include <percept/PerceptBoostArray.hpp>

#define DEBUG_PRINT_11 0
#define NR_PRINT(a) do { if (DEBUG_PRINT_11) std::cout << #a << " = " << a ; } while(0)
#define NR_PRINT_OUT(a,out) do { if (DEBUG_PRINT_11) out << #a << " = " << a << std::endl; } while(0)

#define STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO 0
#define STK_ADAPT_NODEREGISTRY_DO_REHASH 1

#define DEBUG_NR_UNREF 0
#define DEBUG_NR_DEEP 0

#include <adapt/NodeRegistryType.hpp>

// use old PerceptMesh/BulkData create entities if set to 1 - if 0, use PerceptMesh ID server which is much faster (doesn't use DistributedIndex)
#define USE_CREATE_ENTITIES 0

  namespace percept {

    using std::vector;
    using std::map;
    using std::set;

    class Refiner;

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    class NodeRegistry
    {
    public:

      // this switches new method of determining ownership of subDimEntitys - to be removed
      static const bool s_use_new_ownership_check = true;

      static const unsigned NR_MARK_NONE             = 1u << 0;
      static const unsigned NR_MARK                  = 1u << 1;

      //========================================================================================================================
      // high-level interface

      NodeRegistry(percept::PerceptMesh& eMesh,
                   Refiner *refiner=0,
                   bool useCustomGhosting = true) : m_eMesh(eMesh),
                                                    m_refiner(refiner),
                                                    m_comm_all( new stk::CommSparse(eMesh.parallel()) ),
                                                    m_pseudo_entities(*eMesh.get_bulk_data()),
                                                    m_useCustomGhosting(useCustomGhosting),
                                                    m_useAddNodeSharing(false),
                                                    m_checkForGhostedNodes(false),
                                                    m_gee_cnt(0), m_gen_cnt(0),
                                                    m_debug(false),
                                                    m_state(NRS_NONE)
      {
        m_useCustomGhosting = false;
        m_useAddNodeSharing = !m_useCustomGhosting;
      }

      ~NodeRegistry() {
        if (m_comm_all)
          delete m_comm_all;
      }

      void init_comm_all();
      void clear_dangling_nodes(SetOfEntities* nodes_to_be_deleted);
      void initialize();

      // avoid adding subDimEntity if it has a ghosted node
      void setCheckForGhostedNodes(bool val) { m_checkForGhostedNodes = val; }
      bool getCheckForGhostedNodes() { return m_checkForGhostedNodes; }

#define CHECK_DB 0

      void beginRegistration(int ireg=0, int nreg=1);
      void endRegistration(int ireg=0, int nreg=1);

      void beginCheckForRemote(int ireg=0, int nreg=1);
      void endCheckForRemote(int ireg=0, int nreg=1);

      void beginGetFromRemote(int ireg=0, int nreg=1);
      void endGetFromRemote(int ireg=0, int nreg=1);


      /// Register the need for a new node on the sub-dimensional entity @param subDimEntity on element @param element.
      /// If the element is a ghost element, the entity is still registered: the locality/ownership of the new entity
      /// can be determined by the locality of the element (ghost or not).
      bool registerNeedNewNode(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes,const CellTopologyData * const bucket_topo_data);

      /// check the newly registered node from the registry, which does one of three things, depending on what mode we are in:
      ///   1. counts buffer in prep for sending (just does a pack)
      ///   2. packs the buffer (after buffers are alloc'd)
      ///   3. returns the new node after all communications are done
      bool checkForRemote(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes_notUsed,const CellTopologyData * const bucket_topo_data);
      bool getFromRemote(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes_notUsed,const CellTopologyData * const bucket_topo_data);


      /// util
      /// used to find which proc "owns" a subDimEntity - see more comments in implementation
      int proc_owning_subdim_entity(const SubDimCell_SDCEntityType& subDimEntity, std::vector<int>& other_procs, bool& all_shared);

      /// makes coordinates of this new node be the centroid of its sub entity
      void prolongateCoords(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd);
      void prolongateField(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd, stk::mesh::FieldBase *field);

      /// do interpolation for all fields, for the given nodes
      void prolongateFieldNodeVector(std::vector<stk::mesh::Entity> & nodes);

      double spacing_edge(std::vector<stk::mesh::Entity>& nodes,
                          unsigned iv0, unsigned iv1, unsigned nsz, unsigned nsp,  double lspc[8][3], double den_xyz[3], double *coord[8]);

      void normalize_spacing(stk::mesh::Entity element, std::vector<stk::mesh::Entity> &nodes,
                             unsigned nsz, unsigned nsp, double spc[8][3], double den_xyz[3], double *coord[8]);

      void prolongate(stk::mesh::FieldBase *field, unsigned *subDimSize_in=0, bool useFindValidCentroid = true);

      /// do interpolation for all fields
      void prolongateFields(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd);

      /// do interpolation for all fields
      void prolongateFields();

      /// check for adding new nodes to existing parts based on sub-entity part ownership
      void addToExistingParts(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd);

      void add_rbars(std::vector<std::vector<std::string> >& rbar_types );

      /// Check for adding new nodes to existing parts based on sub-entity part ownership.
      /// This version does it in bulk and thus avoids repeats on shared sub-dim entities.
      void addToExistingPartsNew();

      inline SubDimCellData& getNewNodeAndOwningElement(SubDimCell_SDCEntityType& subDimEntity)
      {
        return m_cell_2_data_map[subDimEntity];
      }

      inline SubDimCellData * getFromMapPtr(const SubDimCell_SDCEntityType& subDimEntity) const
      {
        const SubDimCellToDataMap::const_iterator i = m_cell_2_data_map.find( subDimEntity );
        return i != m_cell_2_data_map.end() ? (const_cast<SubDimCellData *>(&(i->second))) : 0 ;
      }

      inline SubDimCellData& getFromMap(const SubDimCell_SDCEntityType& subDimEntity) const
      {
        const SubDimCellToDataMap::const_iterator i = m_cell_2_data_map.find( subDimEntity );
        return * (const_cast<SubDimCellData *>(&(i->second)));
      }

      // allow client code to insert into map (e.g. for nodes created on-the-fly - see quad transition code)
      void forceInMap(std::vector<stk::mesh::Entity>& vecNodes, unsigned mark, stk::mesh::Entity element, stk::mesh::EntityRank needed_rank, unsigned iSubDimOrd );

      inline void putInMap(SubDimCell_SDCEntityType& subDimEntity, SubDimCellData& data)
      {
        bool debug = false;
        if (debug)
          {
            SubDimCellData& nodeId_elementOwnderId = data;

            NodeIdsOnSubDimEntityType& nodeIds_onSE = std::get<SDC_DATA_GLOBAL_NODE_IDS>(nodeId_elementOwnderId);
            stk::mesh::EntityId owning_elementId = std::get<SDC_DATA_OWNING_ELEMENT_KEY>(nodeId_elementOwnderId).id();

            if (1)
              std::cout << "put in map: nodeIds_onSE.size= " << (nodeIds_onSE.size())
                        << " owning_elementId= " << owning_elementId
                        << " subDimEntity.size= " << subDimEntity.size()
                        << "\n" << m_eMesh.demangled_stacktrace(20)
                        << std::endl;
          }
        m_cell_2_data_map[subDimEntity] = data;
      }

      /// @param needNodes should be true in general; it's used by registerNeedNewNode to generate actual data or not on the subDimEntity
      ///   For local refinement, subDimEntity's needs are not always known uniquely by the pair {elementId, iSubDimOrd}; for example, in
      ///   an element-based marking scheme, the shared face between two elements may be viewed differently.  So, we need the ability to
      ///   override the default behavior of always creating new nodes on the subDimEntity, but still allow the entity to be created in
      ///   the NodeRegistry databse.

      typedef bool (NodeRegistry::*ElementFunctionPrototype)( const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes,
                                                              const CellTopologyData * const bucket_topo_data
                                                              );

      /// this is a helper method that loops over all sub-dimensional entities whose rank matches on of those in @param needed_entity_ranks
      ///    and registers that sub-dimensional entity as needing a new node.
      /// @param isGhost should be true if this element is a ghost, in which case this will call the appropriate method to set up for
      //     communications

      void doForAllSubEntities(ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data);

      /// fill
      ///    @param subDimEntity with the stk::mesh::EntityId's of
      ///    the ordinal @param iSubDimOrd sub-dimensional entity of
      ///    @param element of rank
      ///    @param needed_entity_rank
      ///
      bool getSubDimEntity(SubDimCell_SDCEntityType& subDimEntity, const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd,
                           const CellTopologyData * const bucket_topo_data_element_0 =0
                           );

      //========================================================================================================================
      // low-level interface

      unsigned total_size();
      unsigned local_size();

      /// Replace element ownership
      /// When remeshing during unrefinement, replace ownership of sub-dim entities by non-deleted elements
      bool replaceElementOwnership(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes,
                                   const CellTopologyData * const bucket_topo_data
                                   );

      /// communicates pending deletes and removes subDimEntity's as requested
      bool replaceElementOwnership();

      bool initializeEmpty(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes_notUsed,
                           const CellTopologyData * const bucket_topo_data);

      void setAllReceivedNodeData();

      /// when a sub-dim entity is visited during node registration but is flagged as not being marked, and thus not requiring
      ///   any new nodes, we flag it with NR_MARK_NONE, then remove it here
      void removeUnmarkedSubDimEntities();

      void communicate_marks();

    private:

      void communicate_marks_pack(stk::CommSparse& commAll, int stage);
      void communicate_marks_unpack(stk::CommSparse& commAll);

      void do_add_node_sharing_comm();

      void prolongateUsingValidCentroid(std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity> >& findValidCentroidElementNodePairs,
              stk::mesh::Part* new_nodes_part, stk::mesh::FieldBase *field);
    public:

      bool is_empty( const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd);

      static std::string print_string(PerceptMesh& eMesh, const SubDimCell_SDCEntityType& subDimEntity)
      {
        std::ostringstream ost;
        ost << "SubDimCell:";
        for (unsigned ii=0; ii < subDimEntity.size(); ++ii)
          {
            ost << " " << eMesh.identifier(subDimEntity[ii]);
          }
        return ost.str();
      }

      static void print(PerceptMesh& eMesh, const SubDimCell_SDCEntityType& subDimEntity)
      {
        std::cout << print_string(eMesh, subDimEntity);
      }

      inline stk::mesh::Entity get_entity(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {
        return bulk.get_entity(rank, id);
      }

      inline stk::mesh::Entity get_entity_element(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id)
      {
        return get_entity(bulk, rank, id);
      }

      NodeIdsOnSubDimEntityType* getNewNodesOnSubDimEntity(const stk::mesh::Entity element,  stk::mesh::EntityRank& needed_entity_rank, unsigned iSubDimOrd);

      /// allocate the send/recv buffers for all-to-all communication
      bool allocateBuffers();
      void communicate();
      void unpack();

      /// after registering all needed nodes, this method is used to request new nodes on this processor
      void createNewNodesInParallel();

      /// unpacks the incoming information in @param buffer_entry and adds that information to my local node registry
      /// (i.e. the map of sub-dimensional entity to global node id is updated)
      void createNodeAndConnect(CommDataType& buffer_entry, NodeIdsOnSubDimEntityType& nodeIds_onSE, unsigned from_proc, vector<stk::mesh::EntityProc>& nodes_to_ghost);

    public:
      inline SubDimCellToDataMap& getMap() { return  m_cell_2_data_map; }
      inline PerceptMesh& getMesh() { return m_eMesh; }
      inline bool getUseCustomGhosting() { return m_useCustomGhosting; }

      /// remove any sub-dim entities from the map that have a node in deleted_nodes - it is assumed
      ///   that the list of deleted_nodes is parallel-consistent (a node being deleted on proc A, owned
      ///   by proc B, must also be in proc B's deleted_nodes list)
      void cleanDeletedNodes(SetOfEntities& deleted_nodes,
                             SetOfEntities& kept_nodes_orig_minus_kept_nodes,
                             SubDimCellToDataMap& to_save,
                             bool debug=false);

      void cleanInvalidNodes(bool debug=false);

      // further cleanup of the NodeRegistry db - some elements get deleted on some procs but the ghost elements
      //   are still in the db - the easiest way to detect this is as used here: during modification_begin(),
      //   cleanup of cached transactions are performed, then we find these by seeing if our element id's are
      //   no longer in the stk_mesh db.
      void clear_element_owner_data_phase_2(bool resize_nodeId_data=true, bool mod_begin_end=true, SetOfEntities * elemsToBeDeleted=0);
      void dumpDB(std::string msg="");

      // estimate of memory used by this object
      size_t get_memory_usage();

      void mod_begin();
      void mod_end(const std::string& msg="");

    private:
      percept::PerceptMesh& m_eMesh;
      Refiner *m_refiner;
      stk::CommSparse * m_comm_all;
      SubDimCellToDataMap m_cell_2_data_map;

      vector<stk::mesh::EntityProc> m_nodes_to_ghost;
      SetOfEntities m_pseudo_entities;
      bool m_useCustomGhosting;
      bool m_useAddNodeSharing;
      bool m_checkForGhostedNodes;

    public:
      int m_gee_cnt;
      int m_gen_cnt;

      bool m_debug;

      NodeRegistryState m_state;
    };
  } // namespace percept

#endif

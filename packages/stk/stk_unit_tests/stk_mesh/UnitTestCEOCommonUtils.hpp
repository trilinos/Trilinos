#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <gtest/gtest.h>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityRank, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_topology/topology.hpp"    // for topology, etc
#include "BulkDataTester.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/Types.hpp"

namespace CEOUtils
{

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::PartVector;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityRank;
using stk::mesh::PairIterEntityComm;

const EntityRank NODE_RANK = stk::topology::NODE_RANK;
const EntityRank EDGE_RANK = stk::topology::EDGE_RANK;
const EntityRank FACE_RANK = stk::topology::FACE_RANK;
const EntityRank ELEM_RANK = stk::topology::ELEM_RANK;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline bool isEntityValidOnCommList(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), entityKey);
    return iter != stkMeshBulkData.comm_list().end() && entityKey == iter->key;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline bool isEntityInGhostingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    bool is_entity_in_aura_comm_map = !stkMeshBulkData.entity_comm_map(entityKey, stkMeshBulkData.aura_ghosting()).empty();
    return is_entity_in_aura_comm_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline bool isEntityInSharingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)

{
    EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    bool is_entity_in_shared_comm_map = !stkMeshBulkData.entity_comm_map(entityKey, stkMeshBulkData.shared_ghosting()).empty();
    return is_entity_in_shared_comm_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void addSharingInfo(BulkDataTester& bulkData, stk::mesh::Entity entity, stk::mesh::BulkData::GHOSTING_ID ghostingId, int sharingProc)
{
    EXPECT_TRUE(bulkData.my_entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(ghostingId, sharingProc)));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void eraseSharingInfoUsingKey(BulkDataTester &bulkData, EntityKey key, stk::mesh::BulkData::GHOSTING_ID ghostingId)
{
    bulkData.my_entity_comm_map_erase(key, *bulkData.ghostings()[ghostingId]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void add_nodes_to_move(stk::mesh::BulkData& bulk,
                       stk::mesh::Entity elem,
                       int dest_proc,
                       std::vector<stk::mesh::EntityProc>& entities_to_move)
{
    const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
    for(unsigned i = 0; i < bulk.num_nodes(elem); ++i)
    {
        if(bulk.parallel_owner_rank(nodes[i]) == bulk.parallel_rank())
        {
            entities_to_move.push_back(stk::mesh::EntityProc(nodes[i], dest_proc));
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum EntityStates {
  STATE_VALID,
  STATE_NOT_VALID,
  STATE_OWNED,
  STATE_SHARED,
  STATE_NOT_SHARED,
  STATE_GHOSTED_FROM,
  STATE_NOT_GHOSTED_FROM,
  STATE_GHOSTED_TO,
  STATE_NOT_GHOSTED_TO,
  STATE_MESH_MODIFIED, // state mesh state checks
  STATE_MESH_UNCHANGED,
  STATE_MESH_CREATED,
  STATE_MESH_DELETED
};

inline bool check_state(const BulkDataTester & mesh, const EntityKey & entityKey, EntityStates state,
                 int p0 = -1, int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1, int p5 = -1)
{
  // Check to see if the state is as expected for the provided EntityKey.
  //
  // Meaning of the optional processor arguments for various states:
  //     STATE_OWNED: Processor that owns the Entity
  //    STATE_SHARED: List of Processors that we share this entity with
  //   STATE_GHOSTED_FROM: Processor that we ghost the Entity from
  //
  std::vector<int> expectedProcs;
  if (p0 >= 0) {
    expectedProcs.push_back(p0);
  }
  if (p1 >= 0) {
    expectedProcs.push_back(p1);
  }
  if (p2 >= 0) {
    expectedProcs.push_back(p2);
  }
  if (p3 >= 0) {
    expectedProcs.push_back(p3);
  }
  if (p4 >= 0) {
    expectedProcs.push_back(p4);
  }
  if (p5 >= 0) {
    expectedProcs.push_back(p5);
  }
  std::sort(expectedProcs.begin(), expectedProcs.end());

  const int p_rank = mesh.parallel_rank();
  stk::mesh::Entity entity = mesh.get_entity(entityKey);
  std::ostringstream oss;

  switch (state) {
    case STATE_VALID:
    {
      if (!expectedProcs.empty()) {
        oss << "check_state(): Cannot provide processors with validity check." << std::endl;
      }
      if (!mesh.is_valid(entity)) {
        oss << "check_state(): Entity " << entityKey << " is not valid when it should have been." << std::endl;
      }
      break;
    }
    case STATE_NOT_VALID:
    {
      if (!expectedProcs.empty()) {
        oss << "check_state(): Cannot provide processors with STATE_NOT_VALID check." << std::endl;
      }
      if (mesh.is_valid(entity)) {
        oss << "check_state(): Entity " << entityKey << " is valid when it shouldn't have been." << std::endl;
      }
      break;
    }
    case STATE_OWNED:
    {
      if (expectedProcs.size() != 1u) {
        oss << "check_state(): Entities can have only one owner." << std::endl;
      }
      if (mesh.is_valid(entity)) {
        if (expectedProcs[0] != mesh.parallel_owner_rank(entity) ) {
          oss << "check_state(): Owner of entity " << entityKey << " was proc " << mesh.parallel_owner_rank(entity)
              << " and not proc " << expectedProcs[0] << std::endl;
        }
      }
      else {
        oss << "check_state(): Can't check ownership of locally-invalid entity." << std::endl;
      }
      break;
    }
    case STATE_SHARED:
    {
      if (expectedProcs.empty()) {
        oss << "check_state(): Must provide processor(s) with STATE_SHARED check." << std::endl;
      }
      stk::mesh::PairIterEntityComm comm_it = mesh.entity_comm_map_shared(entityKey);
      std::vector<int>::const_iterator expected_procs_it = expectedProcs.begin();
      bool lists_match = true;

      if (comm_it.size() != expectedProcs.size()) {
        lists_match = false;
      }
      else {
        for ( ; expected_procs_it != expectedProcs.end(); ++comm_it, ++expected_procs_it) {
          int comm_proc = comm_it.first->proc;
          int user_proc = *expected_procs_it;
          if (comm_proc != user_proc) {
            lists_match = false;
            break;
          }
        }
      }

      if (!lists_match) {
        oss << "check_state(): Entity " << entityKey << " was shared with procs (";
        comm_it = mesh.entity_comm_map_shared(entityKey);
        for ( ; comm_it.first != comm_it.second; ++comm_it) {
          int proc = comm_it.first->proc;
          oss << proc << " ";
        }
        oss << ")" << std::endl
            << "               when it was expected to be shared with procs (";
        expected_procs_it = expectedProcs.begin();
        for ( ; expected_procs_it != expectedProcs.end(); ++expected_procs_it) {
          oss << *expected_procs_it << " ";
        }
        oss << ")" << std::endl;
      }

      break;
    }
    case STATE_NOT_SHARED:
    {
      if (!expectedProcs.empty()) {
        oss << "check_state(): Cannot provide processors with STATE_NOT_SHARED check." << std::endl;
      }
      if (!mesh.entity_comm_map_shared(entityKey).empty()) {
        oss << "check_state(): Entity " << entityKey << " was shared with procs (";
        stk::mesh::PairIterEntityComm comm_pit = mesh.entity_comm_map_shared(entityKey);
        for ( ; comm_pit.first != comm_pit.second; ++comm_pit) {
          int proc = comm_pit.first->proc;
          oss << proc << " ";
        }
        oss << ") when it shouldn't have been shared." << std::endl;
      }
      break;
    }
    case STATE_GHOSTED_FROM:
    {
      if (expectedProcs.size() != 1) {
        oss << "check_state(): Must provide one processor with STATE_GHOSTED_FROM_FROM check." << std::endl;
        break;   //need to break otherwise following call can segfault
      }
      bool entityIsInvalid = mesh.get_entity(entityKey) == Entity();
      bool inGhost = mesh.in_ghost(mesh.aura_ghosting() , entityKey , expectedProcs[0] );
      const int owner_rank_directly_from_comm_map = mesh.my_entity_comm_map().owner_rank(entityKey);

      if ( entityIsInvalid && inGhost )
      {
          break;
      }
      else if ( entityIsInvalid && !inGhost )
      {
        oss << "check_state(): Entity " << entityKey << " was not ghosted from proc " << owner_rank_directly_from_comm_map << "." << std::endl;
        break;
      }
      if (!mesh.in_receive_ghost( mesh.aura_ghosting() , entityKey )) {
        oss << "check_state(): Entity " << entityKey << " was not ghosted from any proc when it should have" << std::endl
            << "               been ghosted from proc " << expectedProcs[0] << "." << std::endl;
      }
      else {
        const int owner_rank = mesh.entity_comm_map_owner(entityKey);
        if (owner_rank != expectedProcs[0]) {
          oss << "check_state(): Entity " << entityKey << " was ghosted from proc " << owner_rank << std::endl
              << "               when it should have been ghosted from proc " << expectedProcs[0] << "." << std::endl;
        }
      }
      break;
    }
    case STATE_NOT_GHOSTED_FROM:
    {
      bool entityIsInvalid = mesh.get_entity(entityKey) == Entity();
      const int owner_rank_directly_from_comm_map = mesh.my_entity_comm_map().owner_rank(entityKey);

      if ( entityIsInvalid && owner_rank_directly_from_comm_map == stk::mesh::InvalidProcessRank)
      {
        break;
      }
      else if ( entityIsInvalid )
      {
          for (int i=0;i<mesh.parallel_size();i++)
          {
              bool inGhost = mesh.in_ghost(mesh.aura_ghosting() , entityKey , i );
              if ( inGhost )
              {
                  oss << "check_state(): Entity " << entityKey << " was ghosted from proc " <<  i << "." << std::endl;
              }
          }
        break;
      }
      if (!expectedProcs.empty()) {
        oss << "check_state(): Cannot provide processors with STATE_NOT_GHOSTED_FROM_FROM check." << std::endl;
      }
      if (mesh.in_receive_ghost( mesh.aura_ghosting() , entityKey )) {
        const int owner_rank = mesh.entity_comm_map_owner(entityKey);
        oss << "check_state(): Entity " << entityKey << " was ghosted from proc " << owner_rank << std::endl
            << "               when it shouldn't have been ghosted." << std::endl;
      }
      break;
    }
    case STATE_GHOSTED_TO:
    {
      if (expectedProcs.empty()) {
        oss << "check_state(): Must provide at least one processor with STATE_GHOSTED_FROM_TO check." << std::endl;
      }

      if (mesh.parallel_owner_rank(entity) != p_rank) {
        oss << "check_state(): Cannot check STATE_GHOSTED_FROM_TO for an Entity (" << entityKey << ")" << std::endl
            << "that we do not own" << std::endl;
      }
      std::vector<int> meshProcs;
      for ( PairIterEntityComm ec = mesh.entity_comm_map(entityKey); ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == 1 ) {
          meshProcs.push_back(ec->proc);
        }
      }
      std::sort(meshProcs.begin(), meshProcs.end());

      bool lists_match = true;
      if (meshProcs.size() != expectedProcs.size()) {
        lists_match = false;
      }
      else {
        for ( size_t i = 0; i < expectedProcs.size(); ++i ) {
          if (meshProcs[i] != expectedProcs[i]) {
            lists_match = false;
            break;
          }
        }
      }

      if (!lists_match) {
        oss << "check_state(): Entity " << entityKey << " was ghosted to procs (";
        for ( size_t i = 0; i < meshProcs.size(); ++i ) {
          oss << meshProcs[i] << " ";
        }
        oss << ")" << std::endl
            << "               when it was expected to be ghosted to procs (";
        for ( size_t i = 0; i < expectedProcs.size(); ++i ) {
          oss << expectedProcs[i] << " ";
        }
        oss << ")" << std::endl;
      }
      break;
    }
    case STATE_NOT_GHOSTED_TO:
    {
      if (!expectedProcs.empty()) {
        oss << "check_state(): Cannot provide processors with STATE_NOT_GHOSTED_TO check." << std::endl;
      }
      std::vector<int> meshProcs;
      for ( PairIterEntityComm ec = mesh.entity_comm_map(entityKey); ! ec.empty() ; ++ec ) {
        if ( (ec->ghost_id == 1) &&
             (mesh.parallel_owner_rank(entity) == p_rank) ) {
          meshProcs.push_back(ec->proc);
        }
      }

      if (!meshProcs.empty()) {
        std::sort(meshProcs.begin(), meshProcs.end());
        oss << "check_state(): Entity " << entityKey << " was ghosted to procs (";
        for ( size_t i = 0; i < meshProcs.size(); ++i ) {
          oss << meshProcs[i] << " ";
        }
        oss << ")" << std::endl << "               when it shouldn't have been ghosted to anyone." << std::endl;
      }
      break;
    }
    case STATE_MESH_MODIFIED:
    {
        if (mesh.state(entity) != stk::mesh::Modified) {
            oss << "check_state(): Entity" << entityKey << "has mesh state " <<
                    mesh.state(entity) << " when it should have been Modified." << std::endl;
        }
        break;
    }
    case STATE_MESH_CREATED:
    {
        if (mesh.state(entity) != stk::mesh::Created) {
            oss << "check_state(): Entity" << entityKey << "has mesh state " <<
                    mesh.state(entity) << " when it should have been Created." << std::endl;
        }
        break;
    }
    case STATE_MESH_UNCHANGED:
    {
        if (mesh.state(entity) != stk::mesh::Unchanged) {
            oss << "check_state(): Entity" << entityKey << "has mesh state " <<
                    mesh.state(entity) << " when it should have been Unchanged." << std::endl;
        }
        break;
    }
    case STATE_MESH_DELETED:
    {
        if (mesh.state(entity) != stk::mesh::Deleted) {
            oss << "check_state(): Entity" << entityKey << "has mesh state " <<
                    mesh.state(entity) << " when it should have been Deleted." << std::endl;
        }
        break;
    }
  }

  if (oss.str().size() > 0u) {
    std::cout << oss.str();
    return false;
  }
  else {
    return true;
  }

}

inline bool check_parts(const stk::mesh::BulkData & mesh, const EntityKey & entityKey,
                 stk::mesh::Part *p0 = NULL, stk::mesh::Part *p1 = NULL, stk::mesh::Part *p2 = NULL,
                 stk::mesh::Part *p3 = NULL, stk::mesh::Part *p4 = NULL, stk::mesh::Part *p5 = NULL,
                 stk::mesh::Part *p6 = NULL, stk::mesh::Part *p7 = NULL, stk::mesh::Part *p8 = NULL,
                 stk::mesh::Part *p9 = NULL, stk::mesh::Part *p10 = NULL, stk::mesh::Part *p11 = NULL)
{
  // Check to see if the state is as expected for the provided EntityKey.
  //
  // Meaning of the optional processor arguments for various states:
  //     STATE_OWNED: Processor that owns the Entity
  //    STATE_SHARED: List of Processors that we share this entity with
  //   STATE_GHOSTED_FROM: Processor that we ghost the Entity from
  //
  PartVector expectedParts;
  if ( p0 != NULL) { expectedParts.push_back(p0 ); }
  if ( p1 != NULL) { expectedParts.push_back(p1 ); }
  if ( p2 != NULL) { expectedParts.push_back(p2 ); }
  if ( p3 != NULL) { expectedParts.push_back(p3 ); }
  if ( p4 != NULL) { expectedParts.push_back(p4 ); }
  if ( p5 != NULL) { expectedParts.push_back(p5 ); }
  if ( p6 != NULL) { expectedParts.push_back(p6 ); }
  if ( p7 != NULL) { expectedParts.push_back(p7 ); }
  if ( p8 != NULL) { expectedParts.push_back(p8 ); }
  if ( p9 != NULL) { expectedParts.push_back(p9 ); }
  if (p10 != NULL) { expectedParts.push_back(p10); }
  if (p11 != NULL) { expectedParts.push_back(p11); }

  std::sort(expectedParts.begin(), expectedParts.end());

  Entity entity = mesh.get_entity(entityKey);
  std::ostringstream oss;

  if (expectedParts.empty()) {
    oss << "check_parts(): Must provide at least one Part." << std::endl;
  }
  if (!mesh.is_valid(entity)) {
    oss << "check_parts(): Entity " << entityKey << " is not valid, call not completed." << std::endl;
  }
  if (oss.str().size() > 0u) {
    std::cout << oss.str();
    return false; // get out early because doing the following line on an invalid entity crashes
  }
  const PartVector & unsortedMeshParts = mesh.bucket(entity).supersets();
  PartVector meshParts = unsortedMeshParts;
  std::sort(meshParts.begin(), meshParts.end());

  PartVector::const_iterator mesh_parts_it = meshParts.begin();
  PartVector::const_iterator expected_parts_it = expectedParts.begin();
  bool lists_match = true;

  if (expectedParts.size() != meshParts.size()) {
    lists_match = false;
  }
  else {
    for ( ; expected_parts_it != expectedParts.end(); ++expected_parts_it, ++mesh_parts_it) {
      if ((*expected_parts_it) != (*mesh_parts_it)) {
        lists_match = false;
        break;
      }
    }
  }

  if (!lists_match) {
    oss << "check_state(): Entity " << entityKey << " existed on Parts (";
    mesh_parts_it = meshParts.begin();
    for ( ; mesh_parts_it != meshParts.end(); ++mesh_parts_it) {
      oss << (*mesh_parts_it)->name() << " ";
    }
    oss << ")" << std::endl
        << "               when it was expected to exist on Parts (";
    expected_parts_it = expectedParts.begin();
    for ( ; expected_parts_it != expectedParts.end(); ++expected_parts_it) {
      oss << (*expected_parts_it)->name() << " ";
    }
    oss << ")" << std::endl;
  }


  if (oss.str().size() > 0u) {
    std::cout << oss.str();
    return false;
  }
  else {
    return true;
  }

}

inline bool check_relns(const stk::mesh::BulkData & mesh, const EntityKey & entityKey, stk::mesh::EntityRank to_rank,
                        EntityId id0 = EntityKey::MIN_ID, EntityId id1 = EntityKey::MIN_ID,
                        EntityId id2 = EntityKey::MIN_ID, EntityId id3 = EntityKey::MIN_ID,
                        EntityId id4 = EntityKey::MIN_ID, EntityId id5 = EntityKey::MIN_ID,
                        EntityId id6 = EntityKey::MIN_ID, EntityId id7 = EntityKey::MIN_ID,
                        EntityId id8 = EntityKey::MIN_ID, EntityId id9 = EntityKey::MIN_ID,
                        EntityId id10 = EntityKey::MIN_ID, EntityId id11 = EntityKey::MIN_ID)
{
  std::vector<EntityId> expectedIds;
  if (id0 != EntityKey::MIN_ID) {
    expectedIds.push_back(id0);
  }
  if (id1 != EntityKey::MIN_ID) {
    expectedIds.push_back(id1);
  }
  if (id2 != EntityKey::MIN_ID) {
    expectedIds.push_back(id2);
  }
  if (id3 != EntityKey::MIN_ID) {
    expectedIds.push_back(id3);
  }
  if (id4 != EntityKey::MIN_ID) {
    expectedIds.push_back(id4);
  }
  if (id5 != EntityKey::MIN_ID) {
    expectedIds.push_back(id5);
  }
  if (id6 != EntityKey::MIN_ID) {
    expectedIds.push_back(id6);
  }
  if (id7 != EntityKey::MIN_ID) {
    expectedIds.push_back(id7);
  }
  if (id8 != EntityKey::MIN_ID) {
    expectedIds.push_back(id8);
  }
  if (id9 != EntityKey::MIN_ID) {
    expectedIds.push_back(id9);
  }
  if (id10 != EntityKey::MIN_ID) {
    expectedIds.push_back(id10);
  }
  if (id11 != EntityKey::MIN_ID) {
    expectedIds.push_back(id11);
  }

  const int p_rank = mesh.parallel_rank();
  stk::mesh::Entity entity = mesh.get_entity(entityKey);
  std::ostringstream oss;

  if (!mesh.is_valid(entity)) {
    oss << "check_relns(): Entity " << entityKey << " is not valid, call not completed." << std::endl;
  }
  if (oss.str().size() > 0u) {
    std::cout << oss.str();
    return false; // get out early because doing the following line on an invalid entity crashes
  }
  Entity const * relations = mesh.begin(entity, to_rank);

  bool lists_match = true;

  unsigned numRels = mesh.num_connectivity(entity, to_rank);
  std::vector<EntityId> meshIds(numRels);
  for (unsigned i = 0; i < numRels; ++i) {
    meshIds[i] = mesh.identifier(relations[i]);
  }

  // If we're checking upward relations, then sort the lists of
  // ids since they are in an unpredictable order.
  if (to_rank > entityKey.rank()) {
    std::sort(meshIds.begin(), meshIds.end());
    std::sort(expectedIds.begin(), expectedIds.end());
  }


  if (meshIds.size() != expectedIds.size()) {
    lists_match = false;
  }
  else {
    for (unsigned i = 0; i < numRels; ++i) {
      if (expectedIds[i] != meshIds[i]) {
        lists_match = false;
        break;
      }
    }
  }

  if (!lists_match) {
    oss << "[p" << p_rank << "] check_state(): Entity " << entityKey << " relations to " << to_rank << " are (";
    for (unsigned i = 0; i < numRels; ++i) {
      oss << meshIds[i] << " ";
    }
    oss << ")" << std::endl
        << "                   but were expected to be (";
    for (unsigned i = 0; i < expectedIds.size(); ++i) {
      oss << expectedIds[i] << " ";
    }
    oss << ")" << std::endl;
  }

  if (oss.str().size() > 0u) {
    std::cout << oss.str();
    return false;
  }
  else {
    return true;
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Mesh Creation Functions With Tests              /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////// 2Elem2ProcMove //////////////////////////////////////
inline void fillMeshfor2Elem2ProcMoveAndTest(BulkDataTester& bulk, stk::mesh::MetaData &meta, std::vector<stk::mesh::Entity>& elems)
{
    //   id/owner_proc
    //
    //   1/0---4/0---5/0      1/0---4/1---5/1
    //    |     |     |        |     |     |
    //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
    //    |     |     |        |     |     |
    //   2/0---3/0---6/0      2/0---3/0---6/1

    stk::mesh::EntityId element_ids[2] = {1, 2};
    stk::mesh::EntityId elem_node_ids[][4] = { {1, 2, 3, 4}, {4, 3, 6, 5}};

    stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    meta.commit();

// Start with all entities on proc 0

    bulk.modification_begin();
    int p_rank = bulk.parallel_rank();
    if(p_rank == 0)
    {
        elems.push_back(stk::mesh::declare_element(bulk, *elem_part, element_ids[0], elem_node_ids[0]));
        elems.push_back(stk::mesh::declare_element(bulk, *elem_part, element_ids[1], elem_node_ids[1]));
    }
    bulk.modification_end();

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
//    Part * shared_part    = &meta.globally_shared_part();
//    Part * aura_part      = &meta.aura_part();

// Check initial state
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));
 
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));
 
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));
 
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));
 
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));
 
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    }
}

inline void checkStatesAfterCEO_2Elem2ProcMove(BulkDataTester &bulk)
{
    //   id/owner_proc
    //
    //   1/0---4/0---5/0      1/0---4/1---5/1
    //    |     |     |        |     |     |
    //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
    //    |     |     |        |     |     |
    //   2/0---3/0---6/0      2/0---3/0---6/1

    stk::mesh::MetaData & meta = bulk.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    // Part * shared_part    = &meta.globally_shared_part();
    // Part * aura_part      = &meta.aura_part();

    int p_rank = bulk.parallel_rank();
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));

  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));;
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 1)
    {
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
}

// valid vs not valid, owned.... shared vs not shared, ghosted_from vs not ghosted_from
// ghosted_to vs not_ghosted_to

inline void checkStatesAfterCEOME_2Elem2ProcMove(BulkDataTester &bulk)
{
    //   id/owner_proc
    //
    //   1/0---4/0---5/0      1/0---4/1---5/1
    //    |     |     |        |     |     |
    //    | 1/0 | 2/0 |   =>   | 1/0 | 2/1 |
    //    |     |     |        |     |     |
    //   2/0---3/0---6/0      2/0---3/0---6/1

    stk::mesh::MetaData & meta = bulk.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

    int p_rank = bulk.parallel_rank();
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(bulk, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0 ));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(bulk, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0 ));
        EXPECT_TRUE(check_parts(bulk, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(bulk, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
}


//////////////////////////////////// 2Elem2ProcFlip //////////////////////////////////////

inline void fillMeshfor2Elem2ProcFlipAndTest(BulkDataTester& mesh, stk::mesh::MetaData &meta)
{
    //   id/owner_proc
    //
    //   1/0---4/0---5/1        1/1---4/0---5/0
    //    |     |     |          |     |     |
    //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
    //    |     |     |          |     |     |
    //   2/0---3/0---6/1        2/1---3/0---6/0

    stk::mesh::EntityId element_ids[2] = {1, 2};
    stk::mesh::EntityId elem_node_ids[][4] = { {1, 2, 3, 4}, {4, 3, 6, 5}};

    stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    meta.commit();

// Start with all entities on proc 0
    std::vector<stk::mesh::Entity> elems;
    int p_rank = mesh.parallel_rank();

    mesh.modification_begin();
    if(p_rank == 0)
    {
        elems.push_back(stk::mesh::declare_element(mesh, *elem_part, element_ids[0], elem_node_ids[0]));
        mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3 )), 1);
        mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 1);
    }
    else if(p_rank == 1)
    {
        elems.push_back(stk::mesh::declare_element(mesh, *elem_part, element_ids[1], elem_node_ids[1]));
        mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 3 )), 0);
        mesh.add_node_sharing(mesh.get_entity(EntityKey(NODE_RANK, 4)), 0);
    }

    mesh.modification_end();

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, elem_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
}

inline void checkStatesAfterCEOME_2Elem2ProcFlip(BulkDataTester& mesh)
{
    //   id/owner_proc
    //
    //   1/0---4/0---5/1        1/1---4/0---5/0
    //    |     |     |          |     |     |
    //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
    //    |     |     |          |     |     |
    //   2/0---3/0---6/1        2/1---3/0---6/0

    stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

    int p_rank = mesh.parallel_rank();

    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
}

inline void checkStatesAfterCEO_2Elem2ProcFlip(BulkDataTester& mesh)
{
    //   id/owner_proc
    //
    //   1/0---4/0---5/1        1/1---4/0---5/0
    //    |     |     |          |     |     |
    //    | 1/0 | 2/1 |     =>   | 1/1 | 2/0 |
    //    |     |     |          |     |     |
    //   2/0---3/0---6/1        2/1---3/0---6/0

    stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
//    Part * aura_part      = &meta.aura_part();

    int p_rank = mesh.parallel_rank();

    if(p_rank == 0)
    {
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 1 ));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, topo_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 4, 3, 6, 5));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 1 ));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, topo_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 1 ));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, topo_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));  // Was (1, 2)

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));  // Was (1, 2)

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));  // Was (1, 2)

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));  // Was (1, 2)

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
    }
}

//////////////////////////////////// 3Elem2ProcMoveRight //////////////////////////////////////

inline void fillMeshfor3Elem2ProcMoveRightAndTest(BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector& elements)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

    stk::mesh::Part* elem_part = &meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
    stk::mesh::Part* node_part = &meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
    meta_data.commit();

    int p_rank = mesh.parallel_rank();
// Begin modification cycle so we can create the entities and relations
    mesh.modification_begin();

    if(p_rank == 0)
    {
        nodes.push_back(mesh.declare_entity(NODE_RANK, 1, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 2, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 3, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 4, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 5, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 6, *node_part));
        elements.push_back(mesh.declare_entity(ELEM_RANK, 1, *elem_part));
        elements.push_back(mesh.declare_entity(ELEM_RANK, 2, *elem_part));

        mesh.declare_relation(elements[0], nodes[0], 0);
        mesh.declare_relation(elements[0], nodes[1], 1);
        mesh.declare_relation(elements[0], nodes[2], 2);
        mesh.declare_relation(elements[0], nodes[3], 3);

        mesh.declare_relation(elements[1], nodes[2], 0);
        mesh.declare_relation(elements[1], nodes[3], 1);
        mesh.declare_relation(elements[1], nodes[4], 2);
        mesh.declare_relation(elements[1], nodes[5], 3);

        mesh.add_node_sharing(nodes[4], 1);
        mesh.add_node_sharing(nodes[5], 1);
    }
    else if(p_rank == 1)
    {
        nodes.push_back(mesh.declare_entity(NODE_RANK, 5, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 6, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 7, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 8, *node_part));
        elements.push_back(mesh.declare_entity(ELEM_RANK, 3, *elem_part));

        mesh.declare_relation(elements[0], nodes[0], 0);
        mesh.declare_relation(elements[0], nodes[1], 1);
        mesh.declare_relation(elements[0], nodes[2], 2);
        mesh.declare_relation(elements[0], nodes[3], 3);

        mesh.add_node_sharing(nodes[0], 0);
        mesh.add_node_sharing(nodes[1], 0);
    }

    mesh.modification_end();

    stk::mesh::Part * quad_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    stk::mesh::Part * node_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part * universal_part = &meta_data.universal_part();
    Part * owned_part     = &meta_data.locally_owned_part();
    Part * shared_part    = &meta_data.globally_shared_part();
    Part * aura_part      = &meta_data.aura_part();

    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
}

inline void checkStatesAfterCEO_3Elem2ProcMoveRight(BulkDataTester &mesh)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

    stk::mesh::MetaData & meta_data = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta_data.get_part("elem_part");
    stk::mesh::Part * node_part = meta_data.get_part("node_part");
    stk::mesh::Part * quad_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    stk::mesh::Part * node_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part * universal_part = &meta_data.universal_part();
    Part * owned_part     = &meta_data.locally_owned_part();
    Part * shared_part    = &meta_data.globally_shared_part();

    int p_rank = mesh.parallel_rank();
    if(p_rank == 0)
    {
        // /**/ is the delta from after modification_end to before modification_end

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 1)
    {

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
}

inline void checkStatesAfterCEOME_3Elem2ProcMoveRight(BulkDataTester &mesh)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/0---7/1         1/0---3/0---5/1---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/0 | 3/1 |     =>    | 1/0 | 2/1 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/0---8/1         2/0---4/0---6/1---8/1

    stk::mesh::MetaData & meta_data = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta_data.get_part("elem_part");
    stk::mesh::Part * node_part = meta_data.get_part("node_part");
    stk::mesh::Part * quad_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    stk::mesh::Part * node_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part * universal_part = &meta_data.universal_part();
    Part * owned_part     = &meta_data.locally_owned_part();
    Part * shared_part    = &meta_data.globally_shared_part();
    Part * aura_part      = &meta_data.aura_part();

    int p_rank = mesh.parallel_rank();
    if(p_rank == 0)
    {

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 1)
    {

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
}

//////////////////////////////////// 3Elem2ProcMoveLeft //////////////////////////////////////

inline void fillMeshfor3Elem2ProcMoveLeftAndTest(BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector &elements)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

    int p_rank = mesh.parallel_rank();

    stk::mesh::Part* elem_part = & meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
    stk::mesh::Part* node_part = & meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
    meta_data.commit();

    // Begin modification cycle so we can create the entities and relations
    mesh.modification_begin();

    if(p_rank == 0)
    {
        nodes.push_back(mesh.declare_entity(NODE_RANK, 1, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 2, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 3, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 4, *node_part));
        elements.push_back(mesh.declare_entity(ELEM_RANK, 1, *elem_part));

        mesh.declare_relation(elements[0], nodes[0], 0);
        mesh.declare_relation(elements[0], nodes[1], 1);
        mesh.declare_relation(elements[0], nodes[2], 2);
        mesh.declare_relation(elements[0], nodes[3], 3);

        mesh.add_node_sharing(nodes[2], 1);
        mesh.add_node_sharing(nodes[3], 1);
    }
    else if(p_rank == 1)
    {
        nodes.push_back(mesh.declare_entity(NODE_RANK, 3, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 4, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 5, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 6, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 7, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 8, *node_part));
        elements.push_back(mesh.declare_entity(ELEM_RANK, 2, *elem_part));
        elements.push_back(mesh.declare_entity(ELEM_RANK, 3, *elem_part));

        mesh.declare_relation(elements[0], nodes[0], 0);
        mesh.declare_relation(elements[0], nodes[1], 1);
        mesh.declare_relation(elements[0], nodes[2], 2);
        mesh.declare_relation(elements[0], nodes[3], 3);

        mesh.declare_relation(elements[1], nodes[2], 0);
        mesh.declare_relation(elements[1], nodes[3], 1);
        mesh.declare_relation(elements[1], nodes[4], 2);
        mesh.declare_relation(elements[1], nodes[5], 3);

        mesh.add_node_sharing(nodes[0], 0);
        mesh.add_node_sharing(nodes[1], 0);
    }

    mesh.modification_end();

    stk::mesh::Part * quad_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    stk::mesh::Part * node_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part * universal_part = &meta_data.universal_part();
    Part * owned_part     = &meta_data.locally_owned_part();
    Part * shared_part    = &meta_data.globally_shared_part();
    Part * aura_part      = &meta_data.aura_part();

    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
}

inline void checkStatesAfterCEO_3Elem2ProcMoveLeft(BulkDataTester &mesh)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

    stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * node_part = meta.get_part("node_part");
    stk::mesh::Part * quad_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    stk::mesh::Part * node_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    // Part * aura_part      = &meta.aura_part();

    int p_rank = mesh.parallel_rank();
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM)); ///////
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM)); ///////
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
}

inline void checkStatesAfterCEOME_3Elem2ProcMoveLeft(BulkDataTester &mesh)
{
    //   id/owner_proc
    //
    //   1/0---3/0---5/1---7/1         1/0---3/0---5/0---7/1
    //    |     |     |     |           |     |     |     |
    //    | 1/0 | 2/1 | 3/1 |     =>    | 1/0 | 2/0 | 3/1 |
    //    |     |     |     |           |     |     |     |
    //   2/0---4/0---6/1---8/1         2/0---4/0---6/0---8/1

    stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * node_part = meta.get_part("node_part");
    stk::mesh::Part * quad_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    stk::mesh::Part * node_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

    int p_rank = mesh.parallel_rank();
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, elem_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, aura_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, elem_part, node_part, node_topo_part, quad_topo_part, shared_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, elem_part, node_part, node_topo_part, quad_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));
    }
}

//////////////////////////////////// 4Elem4ProcEdge //////////////////////////////////////

inline void fillMeshfor4Elem4ProcEdgeAndTest(BulkDataTester &mesh, stk::mesh::MetaData &meta_data,
        EntityKey &elem_key_chg_own)
{
    // This unit-test is designed to test the conditions that results that
    // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
    // it will test the changing-of-ownership of a shared edge to a proc that
    // either ghosted it or did not know about it.
    //
    //         id/proc                             id/proc
    //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
    //        |      |     |    ||     |          |      |     |    ||     |
    //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
    //        |      |     |    ||     |          |      |     |    ||     |
    //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
    //  this edge moves to p0 --^
    //  element 3 moves to proc 0.
    //  nodes 7&8 move to proc 0.
    //  proc 2 forgets everything.
    //
    // To test this, we use the mesh above, with each elem going on a separate
    // proc, one elem per proc. We will take the edge shared by the last
    // two (rightmost) elements and change the ownership to proc 0.

    int p_rank = mesh.parallel_rank();

    Part* elem_part = &meta_data.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);
    Part* edge_part = &meta_data.declare_part_with_topology("edge_part", stk::topology::LINE_2);
    Part* node_part = &meta_data.declare_part_with_topology("node_part", stk::topology::NODE);
    Part* elem_block = &meta_data.declare_part("elem_block", stk::topology::ELEMENT_RANK);
    typedef stk::mesh::Field<double> rank_zero_field;
    rank_zero_field  & f0 = meta_data.declare_field< rank_zero_field >( NODE_RANK, "elem_field" );
    stk::mesh::put_field( f0 , *elem_block );
    meta_data.commit();

// Begin modification cycle so we can create the entities and relations
    mesh.modification_begin();

    elem_key_chg_own= EntityKey(ELEM_RANK, 3 /*id*/);
    EntityKey edge_key_chg_own(EDGE_RANK, 1 /*id*/);
    EntityKey node_A_key_chg_own(NODE_RANK, 7 /*id*/);
    EntityKey node_B_key_chg_own(NODE_RANK, 8 /*id*/);
    EntityKey node_C_key(NODE_RANK, 5 /*id*/);
    EntityKey node_D_key(NODE_RANK, 6 /*id*/);

// Create element
    stk::mesh::Entity elem = mesh.declare_entity(ELEM_RANK, p_rank + 1, //elem_id
                                                 *elem_part);

// If it is 2nd to last element, it is the one changing
    if(p_rank == 2)
    {
        EXPECT_TRUE(elem_key_chg_own == mesh.entity_key(elem));
    }

// Create nodes
    stk::mesh::EntityVector nodes;
    if(p_rank == 0)
    {
        nodes.push_back(mesh.declare_entity(NODE_RANK, 1, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 2, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 3, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 4, *node_part));
    }
    else if(p_rank == 1)
    {
        nodes.push_back(mesh.declare_entity(NODE_RANK, 3, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 4, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 5, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 6, *node_part));
    }
    else if(p_rank == 2)
    {
        PartVector add_parts;
        PartVector remove_parts;
        add_parts.push_back(elem_block);
        mesh.change_entity_parts(elem,add_parts,remove_parts);
        nodes.push_back(mesh.declare_entity(NODE_RANK, 5, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 6, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 7, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 8, *node_part));
    }
    else
    { // p_rank == 3
        nodes.push_back(mesh.declare_entity(NODE_RANK, 7, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 8, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 9, *node_part));
        nodes.push_back(mesh.declare_entity(NODE_RANK, 10,*node_part));
    }

// Add element relations to nodes
    unsigned rel_id = 0;
    for(stk::mesh::EntityVector::iterator itr = nodes.begin(); itr != nodes.end(); ++itr, ++rel_id)
    {
        mesh.declare_relation(elem, *itr, rel_id);
    }

// Create edge on last two procs

    if(p_rank >= 2)
    {
       stk::mesh::Entity edge = mesh.declare_entity(EDGE_RANK, 1, // id
                                                    *edge_part);
        EXPECT_TRUE(mesh.entity_key(edge) == edge_key_chg_own);

        // Add element relation to edge
        mesh.declare_relation(elem, edge, 1 /*rel-id*/);

        // Add edge relations to nodes
        if(p_rank == 2)
        {
            mesh.declare_relation(edge, nodes[2], 0);
            mesh.declare_relation(edge, nodes[3], 1);
        }
        else
        { // p_rank == 3
            mesh.declare_relation(edge, nodes[0], 0);
            mesh.declare_relation(edge, nodes[1], 1);
        }
    }
    if(p_rank == 0)
    {
        mesh.add_node_sharing(nodes[2], 1);
        mesh.add_node_sharing(nodes[3], 1);
    }
    else if((p_rank == 1) || (p_rank == 2))
    {
        mesh.add_node_sharing(nodes[0], p_rank - 1);
        mesh.add_node_sharing(nodes[1], p_rank - 1);
        mesh.add_node_sharing(nodes[2], p_rank + 1);
        mesh.add_node_sharing(nodes[3], p_rank + 1);
    }
    else
    { // p_rank ==3
        mesh.add_node_sharing(nodes[0], p_rank - 1);
        mesh.add_node_sharing(nodes[1], p_rank - 1);
    }

    mesh.modification_end();

    Part* elem_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    Part* edge_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::LINE_2));
    Part* node_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part* universal_part = &meta_data.universal_part();
    Part* owned_part     = &meta_data.locally_owned_part();
    Part* shared_part    = &meta_data.globally_shared_part();
    Part* aura_part      = &meta_data.aura_part();

//test pre-conditions
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED,  0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED,  1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 2 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 2 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 2  ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_GHOSTED_FROM, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, aura_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, elem_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 2   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 3  ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, owned_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 2   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 2  ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
    }
}


inline void checkStatesAfterCEO_4Elem4ProcEdge(BulkDataTester &mesh)
{
    // This unit-test is designed to test the conditions that results that
    // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
    // it will test the changing-of-ownership of a shared edge to a proc that
    // either ghosted it or did not know about it.
    //
    //         id/proc                             id/proc
    //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
    //        |      |     |    ||     |          |      |     |    ||     |
    //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
    //        |      |     |    ||     |          |      |     |    ||     |
    //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
    //  this edge moves to p0 --^
    //  element 3 moves to proc 0.
    //  nodes 7&8 move to proc 0.
    //  proc 2 forgets everything.
    //
    // To test this, we use the mesh above, with each elem going on a separate
    // proc, one elem per proc. We will take the edge shared by the last
    // two (rightmost) elements and change the ownership to proc 0.

    int p_rank = mesh.parallel_rank();

    const stk::mesh::MetaData &meta_data = mesh.mesh_meta_data();

    Part* elem_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    Part* edge_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::LINE_2));
    Part* node_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part* elem_part =      meta_data.get_part("elem_part");
    Part* elem_block =     meta_data.get_part("elem_block");
    Part* edge_part =      meta_data.get_part("edge_part");
    Part* node_part =      meta_data.get_part("node_part");
    Part* universal_part = &meta_data.universal_part();
    Part* owned_part     = &meta_data.locally_owned_part();
    Part* shared_part    = &meta_data.globally_shared_part();
    Part* aura_part      = &meta_data.aura_part();

//test post condition
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_block, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0   ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED           /* STATE_SHARED, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1 */));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, owned_part, /* shared_part, */ elem_block, elem_part, elem_topo_part, edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3 /*, 4 */ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED            /* STATE_SHARED, 1 */ ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1       /* STATE_NOT_GHOSTED_FROM */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, /* shared_part, */ aura_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED           /* STATE_SHARED, 1 */ ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1      /* STATE_NOT_GHOSTED_FROM */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, /* shared_part, */ aura_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED           /* STATE_SHARED, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1 */ ));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, /* shared_part,*/ elem_block, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3 /*, 4 */ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED           /* STATE_SHARED, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 1 */ ));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, /* shared_part, */ elem_block, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3 /*, 4 */));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM    /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, edge_part, edge_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 2 /* 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 0 /* 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2 /* , 3 */ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2 /* 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0/*, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, shared_part, elem_block, elem_part, elem_topo_part, node_part, node_topo_part));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2 /* , 3 */));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM /*STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, aura_part, elem_part, elem_topo_part,
        //                                                        edge_part, edge_topo_part, node_part, node_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, 1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, elem_part, elem_topo_part,
        //                                                        edge_part, edge_topo_part, node_part, node_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, 1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 2)
    { //amnesia
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 3 /*STATE_NOT_SHARED */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1 /* STATE_NOT_SHARED */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 /* STATE_NOT_SHARED */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 3 /* STATE_NOT_SHARED */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 /* STATE_NOT_SHARED */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO       /* STATE_GHOSTED_TO, 0 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0   ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 2 /*, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_block, elem_part, elem_topo_part, edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, /* 3, */ 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM     /* STATE_GHOSTED_FROM, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 2 /* 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, shared_part, elem_block, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, 1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, /* 3, */ 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2 /* 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, shared_part, elem_block, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, 1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, /* 3, */ 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0 */));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
    }
}


inline void checkStatesAfterCEOME_4Elem4ProcEdge(BulkDataTester &mesh)
{
    // This unit-test is designed to test the conditions that results that
    // resulted in the difficult-to-fix rebalance use-case bug. Specifically,
    // it will test the changing-of-ownership of a shared edge to a proc that
    // either ghosted it or did not know about it.
    //
    //         id/proc                             id/proc
    //        1/0---3/0---5/1---7/2---9/3         1/0---3/0---5/1---7/0---9/3
    //        |      |     |    ||     |          |      |     |    ||     |
    //        | 1/0  | 2/1 | 3/2|| 4/3 |          | 1/0  | 2/1 | 3/0|| 4/3 |
    //        |      |     |    ||     |          |      |     |    ||     |
    //        2/0---4/0---6/1---8/2---10/3        2/0---4/0---6/1---8/0---10/3
    //  this edge moves to p0 --^
    //  element 3 moves to proc 0.
    //  nodes 7&8 move to proc 0.
    //  proc 2 forgets everything.
    //
    // To test this, we use the mesh above, with each elem going on a separate
    // proc, one elem per proc. We will take the edge shared by the last
    // two (rightmost) elements and change the ownership to proc 0.

    int p_rank = mesh.parallel_rank();

    const stk::mesh::MetaData &meta_data = mesh.mesh_meta_data();

    Part* elem_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));
    Part* edge_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::LINE_2));
    Part* node_topo_part = &meta_data.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
    Part* elem_part =      meta_data.get_part("elem_part");
    Part* edge_part =      meta_data.get_part("edge_part");
    Part* node_part =      meta_data.get_part("node_part");
    Part* elem_block =     meta_data.get_part("elem_block");
    Part* universal_part = &meta_data.universal_part();
    Part* owned_part     = &meta_data.locally_owned_part();
    Part* shared_part    = &meta_data.globally_shared_part();
    Part* aura_part      = &meta_data.aura_part();

//test post condition
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, elem_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 3  ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, owned_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, shared_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 1 ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 3, 4, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_TO, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, shared_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, aura_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, aura_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 2)
    { //amnesia
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, elem_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, elem_part, elem_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 7, 8, 9, 10));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), EDGE_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 0  ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_part, elem_topo_part, edge_part, edge_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 7, 8));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, aura_part, elem_part, elem_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, shared_part, elem_part, elem_topo_part,
                                                               edge_part, edge_topo_part, node_part, node_topo_part, elem_block));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), EDGE_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, elem_part, elem_topo_part, node_part, node_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4));
    }
}

//////////////////////////////////// 8Elem4ProcMoveTop //////////////////////////////////////

inline void fillMeshfor8Elem4ProcMoveTopAndTest(BulkDataTester &mesh, stk::mesh::MetaData &meta)
{
    //
    //     id/proc                           id/proc
    //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
    //
    // This test moves ownership of elements 6 and 7 (as well as their locally-owned
    // nodes) to procs 3 and 0, respectively.
    //

    const int p_rank = mesh.parallel_rank();

    //we don't need this because the setup8Quad4ProcMesh2D declares this part already as block_1
    //stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::QUAD_4_2D);


    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    setup8Quad4ProcMesh2D(mesh);

    stk::mesh::Part * block_1 = meta.get_part("block_1");

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

// Check the initial state
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 2, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 3, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
}


inline void checkStatesAfterCEO_8Elem4ProcMoveTop(BulkDataTester &mesh)
{
    //
    //     id/proc                           id/proc
    //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
    //
    // This test moves ownership of elements 6 and 7 (as well as their locally-owned
    // nodes) to procs 3 and 0, respectively.
    //

    stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    stk::mesh::Part * block_1 = meta.get_part("block_1");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    // Part * aura_part      = &meta.aura_part();

    int p_rank = mesh.parallel_rank();

// Check the final state
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 2, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 2, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, owned_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1 /*, 2 */));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 1 /* , 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 2 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, owned_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 5 /* 1, 2, 5, 6*/ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_SHARED /* STATE_SHARED, 1, 2, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, /* shared_part, */  topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, /* 2, 3, 6, */  7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED /* STATE_SHARED, 2, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_TO, 1 */ ));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 7 /* 3, 4, 7, 8 */ ));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 1 /* 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 2 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, owned_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5 /* , 6 */ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED /* STATE_SHARED, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, /* shared_part, */ topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, /* 6, */  7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED /* STATE_SHARED, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_TO, 1, 2 */ ));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7 /* , 8 */ ));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
    else if(p_rank == 1)
    {
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
         // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, /* 1, */ 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 3 */));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, owned_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2 /* , 3 */ ));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 0 /* , 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 2 /* 1, 2, 5, 6 */ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 2 /* 0, 2, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2 /* , 3, 6, 7 */ ));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0 /* STATE_NOT_SHARED */ ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 2 /* STATE_NOT_SHARED */ ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 1, 3 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, /* 2, */  3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 1 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, owned_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3 /*, 4 */ ));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 2, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 1 /* 0, 1, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 3 /* 2, 3, 6, 7 */ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, /* 0, */ 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3 /* , 4, 7, 8 */ ));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 6));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 1 /* STATE_NOT_SHARED */ ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3 /* STATE_NOT_SHARED */ ));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 3 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
    else if(p_rank == 3)
    {
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM /** STATE_GHOSTED_FROM, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 2 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 1, 2 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM  /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM /** STATE_GHOSTED_FROM, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, /* 3, */ 4 ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED /* STATE_SHARED, 0, 1 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, /* shared_part, */ topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, /* 1, 2, 5, */ 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_SHARED /* STATE_SHARED, 0, 1, 2 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, /* 2, 3, */  6 /*, 7 */));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, /* 0, */  2   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, /*3,*/ 4, /*7,*/ 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        // EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM /* STATE_GHOSTED_FROM, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
        // EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
        // EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED /* STATE_SHARED, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, /* shared_part, */ topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, /* 5, */ 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED /* STATE_SHARED, 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 1, 2 */ ));
   /**/ EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, owned_part, /* shared_part, */ topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6 /*, 7*/ ));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 2 /* 0 */ ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, shared_part, topo_part, block_1));
   /**/ EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, /* 7, */ 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
   /**/ EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO /* STATE_GHOSTED_TO, 0, 2 */ ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
}


inline void checkStatesAfterCEOME_8Elem4ProcMoveTop(BulkDataTester &mesh)
{
    //
    //     id/proc                           id/proc
    //     11/0--12/0--13/1--14/2--15/3      11/0--12/0--13/3--14/0--15/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 5/0 | 6/1 | 7/2 | 8/3 |         | 5/0 | 6/3 | 7/0 | 8/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      6/0---7/0---8/1---9/2--10/3  -->  6/0---7/0---8/3---9/0--10/3
    //       |     |     |     |     |         |     |     |     |     |
    //       | 1/0 | 2/1 | 3/2 | 4/3 |         | 1/0 | 2/1 | 3/2 | 4/3 |
    //       |     |     |     |     |         |     |     |     |     |
    //      1/0---2/0---3/1---4/2---5/3       1/0---2/0---3/1---4/2---5/3
    //
    // This test moves ownership of elements 6 and 7 (as well as their locally-owned
    // nodes) to procs 3 and 0, respectively.
    //

    stk::mesh::MetaData & meta = mesh.mesh_meta_data();
    stk::mesh::Part * block_1 = meta.get_part("block_1");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

    int p_rank = mesh.parallel_rank();

// Check the final state
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_TO, 1, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_TO, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 1, 3    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_GHOSTED_TO, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 1, 2, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 2, 3    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_GHOSTED_TO, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_TO, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_TO, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_TO, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 0, 3   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 0, 2, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_GHOSTED_TO, 0, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 2, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 0, 1, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 0, 3   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 7, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 3, 4, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 10, 9));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 5), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 5), NODE_RANK, 6, 7, 12, 11));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 6), STATE_GHOSTED_TO, 0, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 6), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 6), NODE_RANK, 7, 8, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 7), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 7), NODE_RANK, 8, 9, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 8), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 8), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 8), NODE_RANK, 9, 10, 15, 14));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  1), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  2), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  3), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  3), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  4), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  4), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  5), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  5), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  6), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  6), ELEM_RANK, 1, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_SHARED, 0, 1   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  7), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  7), ELEM_RANK, 1, 2, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_SHARED, 0, 1, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  8), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  8), ELEM_RANK, 2, 3, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_SHARED, 0, 2   ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK,  9), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK,  9), ELEM_RANK, 3, 4, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 10), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 10), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 10), ELEM_RANK, 4, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 11), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 11), universal_part, aura_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 11), ELEM_RANK, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 12), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 12), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 12), ELEM_RANK, 5, 6));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_GHOSTED_TO, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 13), universal_part, owned_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 13), ELEM_RANK, 6, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 14), universal_part, shared_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 14), ELEM_RANK, 7, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 15), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 15), universal_part, owned_part, topo_part, block_1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 15), ELEM_RANK, 8));
    }
}

//////////////////////////////////// 4Elem4ProcRotate //////////////////////////////////////

inline void fillMeshfor4Elem4ProcRotateAndTest(BulkDataTester &mesh, stk::mesh::MetaData &meta)
{
    //
    //     id/proc                id/proc
    //      7/3---8/2---9/2        7/2---8/1---9/1
    //       |     |     |          |     |     |
    //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
    //       |     |     |          |     |     |
    //      4/0---5/0---6/1  -->   4/3---5/3---6/0
    //       |     |     |          |     |     |
    //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
    //       |     |     |          |     |     |
    //      1/0---2/0---3/1        1/3---2/3---3/0

    const int p_rank = mesh.parallel_rank();


    stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::QUAD_4_2D);
    stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::QUAD_4_2D);
    stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::QUAD_4_2D);
    stk::mesh::Part& block_4 = meta.declare_part_with_topology("block_4", stk::topology::QUAD_4_2D);
    meta.commit();

// 1 elem-id for each proc
    stk::mesh::EntityId proc_elemIDs[] = {1, 2, 3, 4};

// list of node-ids for each element
    const int nodesPerElem = 4;
    stk::mesh::EntityId elem_nodeIDs[][nodesPerElem] = {
            {1, 2, 5, 4},
            {2, 3, 6, 5},
            {5, 6, 9, 8},
            {4, 5, 8, 7}
    };

// list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    int shared_nodeIDs_and_procs[][3] = {
            {0, 2, 1}, {0, 5, 1}, {0, 5, 2}, {0, 5, 3}, {0, 4, 3}, // proc 0
            {1, 2, 0}, {1, 6, 2}, {1, 5, 0}, {1, 5, 2}, {1, 5, 3}, // proc 1
            {2, 5, 0}, {2, 5, 1}, {2, 5, 3}, {2, 6, 1}, {2, 8, 3}, // proc 2
            {3, 4, 0}, {3, 5, 0}, {3, 5, 1}, {3, 5, 2}, {3, 8, 2}  // proc 3
    };
    int numSharedNodeTriples = 20;

    mesh.modification_begin();

    stk::mesh::EntityId elemID = proc_elemIDs[p_rank];
    if (0 == p_rank) {
        stk::mesh::declare_element(mesh, block_1, elemID, elem_nodeIDs[p_rank]);
    }
    else if (1 == p_rank) {
        stk::mesh::declare_element(mesh, block_2, elemID, elem_nodeIDs[p_rank]);
    }
    else if (2 == p_rank) {
        stk::mesh::declare_element(mesh, block_3, elemID, elem_nodeIDs[p_rank]);
    }
    else if (3 == p_rank) {
        stk::mesh::declare_element(mesh, block_4, elemID, elem_nodeIDs[p_rank]);
    }

    for(int proc = 0; proc < numSharedNodeTriples; ++proc)
    {
        if(p_rank == shared_nodeIDs_and_procs[proc][0])
        {
            int nodeID = shared_nodeIDs_and_procs[proc][1];
            int sharingProc = shared_nodeIDs_and_procs[proc][2];
            stk::mesh::Entity node = mesh.get_entity(NODE_RANK, nodeID);
            mesh.add_node_sharing(node, sharingProc);
        }
    }

    mesh.modification_end();

    stk::mesh::Part * block1_part = meta.get_part("block_1");
    stk::mesh::Part * block2_part = meta.get_part("block_2");
    stk::mesh::Part * block3_part = meta.get_part("block_3");
    stk::mesh::Part * block4_part = meta.get_part("block_4");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

// Check the initial state

    //these are the same everywhere, everything is valid, always
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));

    if (p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 1, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, shared_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, owned_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, owned_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, shared_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 0, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, owned_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 2)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, owned_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 0, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, shared_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
}

inline void checkStatesAfterCEO_4Elem4ProcRotate(BulkDataTester &mesh, stk::mesh::MetaData &meta)
{
    // check intermediate state after change_entity_owner() but before internal_modification_end()
    //     id/proc                id/proc
    //      7/3---8/2---9/2        7/2---8/1---9/1
    //       |     |     |          |     |     |
    //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
    //       |     |     |          |     |     |
    //      4/0---5/0---6/1  -->   4/3---5/3---6/0
    //       |     |     |          |     |     |
    //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
    //       |     |     |          |     |     |
    //      1/0---2/0---3/1        1/3---2/3---3/0
    //
    //NOTES: sharing for nodes that have changed owner has not been deleted, still in previous state
    //       ghosting deleted entirely
    //       shared_part updated to match the commmap/list version of sharing
    //       aura_part never present (no ghosting)

    const int p_rank = mesh.parallel_rank();

    stk::mesh::Part * block1_part = meta.get_part("block_1");
    stk::mesh::Part * block2_part = meta.get_part("block_2");
    stk::mesh::Part * block3_part = meta.get_part("block_3");
    stk::mesh::Part * block4_part = meta.get_part("block_4");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
//    Part * aura_part      = &meta.aura_part(); //gone

    if (p_rank == 0)
    {
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 1));  //shared with 3 at end
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, shared_part, topo_part, block1_part, block2_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 2));  //connected to elem 1 & 2 at end

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));  //ghosted to 1,2,3 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 2));  //related to 1,2,3,4 at end

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED)); //shared with 1 at end
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO)); //ghosted to 2&3 at end
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), universal_part, owned_part, topo_part, block2_part, block3_part)); //shared_part, universal_part, owned_part, topo_part, block2_part, block3_part)); at end
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2));  //2 and 3 at end

  /**/    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  /**/    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  /**/    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 1)
    {
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO)); //0,2&3 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), shared_part, universal_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 3));  //1,2,3&4 at end

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 2)); // 0 at end, not 2
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, topo_part, block2_part, block3_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 3)); //2&3 at end

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2 ));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED)); //2 at end
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO)); //0&3 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), universal_part, owned_part, topo_part, block3_part, block4_part));//add shared_part at end
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3)); //3&4 at end

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO)); //0,2&3 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 2)
    {
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO)); //0,1&3 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));  //3 at end
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), universal_part, topo_part, block1_part, block4_part));  //shared_part as well at end
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 4)); //1 &4 at end

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 4)); //1,2,3&4 at end

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO)); //0,1&3 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 3));  //1 at end, not 3
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, topo_part, block3_part, block4_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 4)); //3&4 at end

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO)); //0,1&2 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO)); //0,1&2 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED)); //0 at end
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO)); //1&2 at end
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, owned_part, topo_part, block1_part, block2_part));  //add shares_part by end
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1));  //1 & 2 at end

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 0));  //shared to 2 at end
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));  //0&1 at end
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, owned_part, topo_part, block1_part, block4_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1)); // 1 & 4 at end

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), shared_part, universal_part, owned_part, topo_part, block1_part, block2_part, block3_part, block4_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1));  // 1,2,3 & 4 at end

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_VALID));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
//        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
//        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
//        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
}

inline void checkStatesAfterCEOME_4Elem4ProcRotate(BulkDataTester &mesh, stk::mesh::MetaData &meta)
{

    //     id/proc                id/proc
    //      7/3---8/2---9/2        7/2---8/1---9/1
    //       |     |     |          |     |     |
    //       | 4/3 | 3/2 |          | 4/2 | 3/1 |
    //       |     |     |          |     |     |
    //      4/0---5/0---6/1  -->   4/3---5/3---6/0
    //       |     |     |          |     |     |
    //       | 1/0 | 2/1 |          | 1/3 | 2/0 |
    //       |     |     |          |     |     |
    //      1/0---2/0---3/1        1/3---2/3---3/0

    const int p_rank = mesh.parallel_rank();

    stk::mesh::Part * block1_part = meta.get_part("block_1");
    stk::mesh::Part * block2_part = meta.get_part("block_2");
    stk::mesh::Part * block3_part = meta.get_part("block_3");
    stk::mesh::Part * block4_part = meta.get_part("block_4");
    stk::mesh::Part * topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::QUAD_4_2D));

    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

    //these are the same everywhere, everything is valid, always
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_VALID));
    EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_VALID));

    if (p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 1, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), universal_part, shared_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_TO, 1, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), universal_part, owned_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 1, 2, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_TO, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, owned_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), aura_part, universal_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 2, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_SHARED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), shared_part, universal_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_TO, 0, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, owned_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_TO, 0, 2, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), universal_part, owned_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 2)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_TO, 0, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), universal_part, owned_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), aura_part, universal_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), aura_part, universal_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), universal_part, shared_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_TO, 0, 1, 3));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), universal_part, owned_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), shared_part, universal_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
    else if (p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 0, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 9, 8));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_GHOSTED_FROM, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 4), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 4), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 4), NODE_RANK, 4, 5, 8, 7));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1), STATE_GHOSTED_TO, 0, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1), universal_part, owned_part, topo_part, block1_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1), ELEM_RANK, 1));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_SHARED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2), STATE_GHOSTED_TO, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2), shared_part, universal_part, owned_part, topo_part, block1_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2), ELEM_RANK, 1, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3), aura_part, universal_part, topo_part, block2_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3), ELEM_RANK, 2));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_SHARED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4), STATE_GHOSTED_TO, 0, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4), shared_part, universal_part, owned_part, topo_part, block1_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4), ELEM_RANK, 1, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_SHARED, 0, 1, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5), shared_part, universal_part, owned_part, topo_part, block1_part, block2_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5), ELEM_RANK, 1, 2, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6), aura_part, universal_part, topo_part, block2_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6), ELEM_RANK, 2, 3));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_GHOSTED_FROM, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7), aura_part, universal_part, topo_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7), ELEM_RANK, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8), aura_part, universal_part, topo_part, block3_part, block4_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8), ELEM_RANK, 3, 4));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 9), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 9), aura_part, universal_part, topo_part, block3_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 9), ELEM_RANK, 3));
    }
}

//////////////////////////////////// 3Elem4Proc1Edge3D //////////////////////////////////////

inline void fillMeshfor3Elem4Proc1Edge3DAndTest(BulkDataTester &mesh, stk::mesh::MetaData &meta)
{
    //  ID.proc
    //                    15.2--------16.2                      15.1--------16.1
    //                     /|          /|                        /|          /|
    //                    / |         / |                       / |         / |
    //                  7.2---------8.2 |                     7.1---------8.1 |
    //                   |  |  3.2   |  |                      |  |  3.1   |  |
    //                   |  |        |  |                      |  |        |  |
    //        12.0-------|13.0-------|14.1          12.3-------|13.3-------|14.0
    //         /|        | *|        | /|   -->      /|        | *|        | /|
    //        / |        |* |        |/ |           / |        |* |        |/ |
    //      4.0---------5.0---------6.1 |         4.3---------5.3---------6.0 |
    //       |  |  1.0   |  |  2.1   |  |          |  |  1.3   |  |  2.0   |  |
    //       |  |        |  |        |  |          |  |        |  |        |  |
    //       | 9.0-------|10.0-------|11.1         | 9.3-------|10.3-------|11.0
    //       | /         | /         | /           | /         | /         | /
    //       |/          |/          |/            |/          |/          |/
    //      1.0---------2.0---------3.1           1.3---------2.3---------3.0
    //
    //      (*)edge: 1.0                          (*)edge: 1.1

    const int p_rank = mesh.parallel_rank();

    stk::mesh::Part * elem_part = &meta.declare_part_with_topology("elem_part", stk::topology::HEX_8);
    stk::mesh::Part * edge_part = &meta.declare_part_with_topology("edge_part", stk::topology::LINE_2);
    meta.commit();

// 1 elem-id for each proc
    stk::mesh::EntityId proc_elemIDs[] =
    {   1, 2, 3};

// list of node-ids for each element
    const int nodesPerElem = 8;
    stk::mesh::EntityId elem_nodeIDs[][nodesPerElem] =
    {
        {   1, 2, 5, 4, 9, 10, 13, 12},
        {   2, 3, 6, 5, 10, 11, 14, 13},
        {   5, 6, 8, 7, 13, 14, 16, 15}
    };

// list of triplets: (owner-proc, shared-nodeID, sharing-proc)
    int shared_nodeIDs_and_procs[][3] =
    {
        {   0, 2, 1},
        {   0, 5, 1},
        {   0, 5, 2},
        {   0, 10, 1},
        {   0, 13, 1},
        {   0, 13, 2}, // proc 0
        {   1, 2, 0},
        {   1, 6, 2},
        {   1, 5, 0},
        {   1, 5, 2},
        {   1, 10, 0},
        {   1, 14, 2},
        {   1, 13, 0},
        {   1, 13, 2}, // proc 1
        {   2, 5, 0},
        {   2, 5, 1},
        {   2, 6, 1},
        {   2, 13, 0},
        {   2, 13, 1},
        {   2, 14, 1} // proc 2
    };
    int numSharedNodeTriples = 20;

    mesh.modification_begin();

    if(p_rank < 3)
    {
        stk::mesh::EntityId elemID = proc_elemIDs[p_rank];
        stk::mesh::Entity elem = stk::mesh::declare_element(mesh, *elem_part, elemID, elem_nodeIDs[p_rank]);
        stk::mesh::Entity edge = mesh.declare_entity(EDGE_RANK, 1, *edge_part);
        std::vector<stk::mesh::Entity> nodes;
        nodes.push_back(mesh.get_entity(NODE_RANK, 5));
        nodes.push_back(mesh.get_entity(NODE_RANK, 13 ));
        stk::mesh::impl::connectEntityToEdge(mesh, elem, edge, nodes);
        mesh.declare_relation(edge, nodes[0], 0);
        mesh.declare_relation(edge, nodes[1], 1);
    }

    for(int proc = 0; proc < numSharedNodeTriples; ++proc)
    {
        if(p_rank == shared_nodeIDs_and_procs[proc][0])
        {
            int nodeID = shared_nodeIDs_and_procs[proc][1];
            int sharingProc = shared_nodeIDs_and_procs[proc][2];
            stk::mesh::Entity node = mesh.get_entity(NODE_RANK, nodeID);
            mesh.add_node_sharing(node, sharingProc);
        }
    }

    mesh.modification_end();

    stk::mesh::Part * hex_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::HEXAHEDRON_8));
    stk::mesh::Part * edge_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::LINE_2));
    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

    // Check the initial state
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 1, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4, 9, 10, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5, 10, 11, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 8, 7, 13, 14, 16, 15));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 1, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, owned_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 1, 2, 3));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_GHOSTED_TO, 1, 2));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_SHARED, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_GHOSTED_TO, 2));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, owned_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 1, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_GHOSTED_FROM, 1 ));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3+i), ELEM_RANK, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_GHOSTED_TO, 1, 2));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 1, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, owned_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 1, 2, 3));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_GHOSTED_FROM, 1 ));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 2, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_OWNED, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_GHOSTED_FROM, 2 ));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_OWNED, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_GHOSTED_FROM, 2 ));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8+i), ELEM_RANK, 3));
        }
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4, 9, 10, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 0, 2));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5, 10, 11, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 2 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 8, 7, 13, 14, 16, 15));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 0, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 1, 2, 3));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_GHOSTED_FROM, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_SHARED, 0 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 1, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_GHOSTED_TO, 0, 2));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3+i), ELEM_RANK, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_GHOSTED_FROM, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 0, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 1, 2, 3));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_SHARED, 2 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_GHOSTED_TO, 0));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, owned_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 2, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_OWNED, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_GHOSTED_FROM, 2 ));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_OWNED, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_GHOSTED_FROM, 2 ));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8+i), ELEM_RANK, 3));
        }
    }
    else if(p_rank == 2)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4, 9, 10, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5, 10, 11, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 0, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 8, 7, 13, 14, 16, 15));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 0    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 0, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 1, 2, 3));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_GHOSTED_FROM, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_GHOSTED_FROM, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 1, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_GHOSTED_FROM, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3+i), ELEM_RANK, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_GHOSTED_FROM, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 0, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 1, 2, 3));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_SHARED, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 2, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_OWNED, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_GHOSTED_TO, 0, 1));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_OWNED, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_GHOSTED_TO, 0, 1));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8+i), ELEM_RANK, 3));
        }
    }
    else if(p_rank == 3)
    { //knows nothing
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));

        for (int i=1;i<=16;i++)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_TO));
        }
    }
}

inline void checkStatesAfterCEO_3Elem4Proc1Edge3D(BulkDataTester &mesh)
{
    //  ID.proc
    //                    15.2--------16.2                      15.1--------16.1
    //                     /|          /|                        /|          /|
    //                    / |         / |                       / |         / |
    //                  7.2---------8.2 |                     7.1---------8.1 |
    //                   |  |  3.2   |  |                      |  |  3.1   |  |
    //                   |  |        |  |                      |  |        |  |
    //        12.0-------|13.0-------|14.1          12.3-------|13.3-------|14.0
    //         /|        | *|        | /|   -->      /|        | *|        | /|
    //        / |        |* |        |/ |           / |        |* |        |/ |
    //      4.0---------5.0---------6.1 |         4.3---------5.3---------6.0 |
    //       |  |  1.0   |  |  2.1   |  |          |  |  1.3   |  |  2.0   |  |
    //       |  |        |  |        |  |          |  |        |  |        |  |
    //       | 9.0-------|10.0-------|11.1         | 9.3-------|10.3-------|11.0
    //       | /         | /         | /           | /         | /         | /
    //       |/          |/          |/            |/          |/          |/
    //      1.0---------2.0---------3.1           1.3---------2.3---------3.0
    //
    //      (*)edge: 1.0                          (*)edge: 1.1

    const int p_rank = mesh.parallel_rank();
    stk::mesh::MetaData &meta = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * edge_part = meta.get_part("edge_part");
    stk::mesh::Part * hex_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::HEXAHEDRON_8));
    stk::mesh::Part * edge_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::LINE_2));
    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();

// Check the final state
    if(p_rank == 0)
    {
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5, 10, 11, 14, 13));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 1    ));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 1, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 2));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 3));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_SHARED, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, shared_part, elem_part, hex_topo_part));
  /**/      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_FROM));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3+i), ELEM_RANK, 2));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 3));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 1, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
  /**/      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 2));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 0));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_FROM));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_TO));
  /**/      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, owned_part, elem_part, hex_topo_part));
  /**/      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 2));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_TO));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_TO));
        }
    }
    else if(p_rank == 1)
    {
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 8, 7, 13, 14, 16, 15));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 1    ));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 0, 2));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, owned_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 3));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_VALID));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_SHARED, 0));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_TO));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_TO));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 3));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 0, 2));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
  /**/      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 3));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 0));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_SHARED, 2 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, shared_part, elem_part, hex_topo_part));
  /**/      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_FROM));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_FROM));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8+i), ELEM_RANK, 3));
        }
    }
    else if(p_rank == 2)
    { //knows nothing
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 0, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));

        for (int i=1;i<=4;i++)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_TO));
        }
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_SHARED, 0, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  5), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  6), STATE_NOT_GHOSTED_TO));
        for (int i=7;i<=12;i++)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_TO));
        }
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_SHARED, 0, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 13), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_SHARED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 14), STATE_NOT_GHOSTED_TO));
        for (int i=15;i<=16;i++)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_TO));
        }
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4, 9, 10, 13, 12));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 1    ));
  /**/  EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
  /**/  EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
  /**/  EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 1));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_FROM));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 3));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_FROM));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_TO));
  /**/      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, owned_part, elem_part, hex_topo_part));
  /**/      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 1));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_FROM));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_GHOSTED_FROM));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 3));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
  /**/      EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, owned_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
  /**/      EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 1));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_TO));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_TO));

  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
  /**/      EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_TO));
        }
    }
}

inline void checkStatesAfterCEOME_3Elem4Proc1Edge3D(BulkDataTester &mesh)
{
    //  ID.proc
    //                    15.2--------16.2                      15.1--------16.1
    //                     /|          /|                        /|          /|
    //                    / |         / |                       / |         / |
    //                  7.2---------8.2 |                     7.1---------8.1 |
    //                   |  |  3.2   |  |                      |  |  3.1   |  |
    //                   |  |        |  |                      |  |        |  |
    //        12.0-------|13.0-------|14.1          12.3-------|13.3-------|14.0
    //         /|        | *|        | /|   -->      /|        | *|        | /|
    //        / |        |* |        |/ |           / |        |* |        |/ |
    //      4.0---------5.0---------6.1 |         4.3---------5.3---------6.0 |
    //       |  |  1.0   |  |  2.1   |  |          |  |  1.3   |  |  2.0   |  |
    //       |  |        |  |        |  |          |  |        |  |        |  |
    //       | 9.0-------|10.0-------|11.1         | 9.3-------|10.3-------|11.0
    //       | /         | /         | /           | /         | /         | /
    //       |/          |/          |/            |/          |/          |/
    //      1.0---------2.0---------3.1           1.3---------2.3---------3.0
    //
    //      (*)edge: 1.0                          (*)edge: 1.1

    const int p_rank = mesh.parallel_rank();
    stk::mesh::MetaData &meta = mesh.mesh_meta_data();
    stk::mesh::Part * elem_part = meta.get_part("elem_part");
    stk::mesh::Part * edge_part = meta.get_part("edge_part");
    stk::mesh::Part * hex_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::HEXAHEDRON_8));
    stk::mesh::Part * edge_topo_part = &meta.get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::LINE_2));
    Part * universal_part = &meta.universal_part();
    Part * owned_part     = &meta.locally_owned_part();
    Part * shared_part    = &meta.globally_shared_part();
    Part * aura_part      = &meta.aura_part();

// Check the final state
    if(p_rank == 0)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4, 9, 10, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_TO, 3, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5, 10, 11, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 8, 7, 13, 14, 16, 15));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 1    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 1, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 1, 2, 3));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_GHOSTED_FROM, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_SHARED, 3 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 1, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_GHOSTED_TO, 3, 1));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3+i), ELEM_RANK, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_GHOSTED_FROM, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 3, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 1, 2, 3));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_SHARED, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_GHOSTED_TO, 3));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, owned_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 2, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_GHOSTED_FROM, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_GHOSTED_FROM, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8+i), ELEM_RANK, 3));
        }
    }
    else if(p_rank == 1)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_FROM, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4, 9, 10, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5, 10, 11, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_TO, 3, 0));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 8, 7, 13, 14, 16, 15));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 1    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 0, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, owned_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 1, 2, 3));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_GHOSTED_FROM, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_GHOSTED_FROM, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 1, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_GHOSTED_FROM, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3+i), ELEM_RANK, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_GHOSTED_FROM, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 3, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 1, 2, 3));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_SHARED, 0 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 2, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_GHOSTED_TO, 3, 0));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_GHOSTED_TO, 3, 0));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8+i), ELEM_RANK, 3));
        }
    }
    else if(p_rank == 2)
    { //knows nothing
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_SHARED));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO));

        for (int i=1;i<=16;i++)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_SHARED));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK,  i), STATE_NOT_GHOSTED_TO));
        }
    }
    else if(p_rank == 3)
    {
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_OWNED, 3));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_NOT_GHOSTED_FROM));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 1), STATE_GHOSTED_TO, 0, 1));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 1), universal_part, owned_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 1), NODE_RANK, 1, 2, 5, 4, 9, 10, 13, 12));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_OWNED, 0));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_GHOSTED_FROM, 0 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 2), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 2), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 2), NODE_RANK, 2, 3, 6, 5, 10, 11, 14, 13));

        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_OWNED, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_GHOSTED_FROM, 1 ));
        EXPECT_TRUE(check_state(mesh, EntityKey(ELEM_RANK, 3), STATE_NOT_GHOSTED_TO));
        EXPECT_TRUE(check_parts(mesh, EntityKey(ELEM_RANK, 3), universal_part, aura_part, elem_part, hex_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(ELEM_RANK, 3), NODE_RANK, 5, 6, 8, 7, 13, 14, 16, 15));

        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_VALID));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_OWNED, 1    ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_SHARED, 0, 1));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_FROM ));
        EXPECT_TRUE(check_state(mesh, EntityKey(EDGE_RANK, 1), STATE_NOT_GHOSTED_TO ));
        EXPECT_TRUE(check_parts(mesh, EntityKey(EDGE_RANK, 1), universal_part, shared_part, elem_part, hex_topo_part,
                                                               edge_part, edge_topo_part));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), ELEM_RANK, 1, 2, 3));
        EXPECT_TRUE(check_relns(mesh, EntityKey(EDGE_RANK, 1), NODE_RANK, 5, 13));

        for (int i=0;i<=8;i+=8)
        {
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 1+i), STATE_GHOSTED_TO, 0, 1));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 1+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 1+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_SHARED, 0 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 2+i), STATE_GHOSTED_TO, 1));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 2+i), universal_part, owned_part, shared_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 2+i), ELEM_RANK, 1, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 3+i), STATE_GHOSTED_FROM, 0 ));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 3+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 3+i), ELEM_RANK, 2));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 4+i), STATE_GHOSTED_TO, 0, 1));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 4+i), universal_part, owned_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 4+i), ELEM_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_OWNED, 3));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_SHARED, 0, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_FROM));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 5+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 5+i), universal_part, owned_part, shared_part, elem_part, hex_topo_part,
                    edge_part, edge_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), ELEM_RANK, 1, 2, 3));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 5+i), EDGE_RANK, 1));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_OWNED, 0));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_GHOSTED_FROM, 0 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 6+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 6+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 6+i), ELEM_RANK, 2, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_GHOSTED_FROM, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 7+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 7+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 7+i), ELEM_RANK, 3));

            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_VALID));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_OWNED, 1));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_SHARED ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_GHOSTED_FROM, 1 ));
            EXPECT_TRUE(check_state(mesh, EntityKey(NODE_RANK, 8+i), STATE_NOT_GHOSTED_TO));
            EXPECT_TRUE(check_parts(mesh, EntityKey(NODE_RANK, 8+i), universal_part, aura_part, elem_part, hex_topo_part));
            EXPECT_TRUE(check_relns(mesh, EntityKey(NODE_RANK, 8+i), ELEM_RANK, 3));
        }
    }
}


}

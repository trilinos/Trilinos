#include "BulkDataTester.hpp"

namespace CEOUtils
{

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline bool isEntityValidOnCommList(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), entityKey);
    return iter != stkMeshBulkData.comm_list().end() && entityKey == iter->key;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline bool isEntityInGhostingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    bool is_entity_in_aura_comm_map = !stkMeshBulkData.entity_comm_map(entityKey, stkMeshBulkData.aura_ghosting()).empty();
    return is_entity_in_aura_comm_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline bool isEntityInSharingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    bool is_entity_in_shared_comm_map = !stkMeshBulkData.entity_comm_map(entityKey, stkMeshBulkData.shared_ghosting()).empty();
    return is_entity_in_shared_comm_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void addSharingInfo(BulkDataTester& bulkData, stk::mesh::Entity entity, stk::mesh::BulkData::GHOSTING_ID ghostingId, int sharingProc)
{
    EXPECT_TRUE(bulkData.my_entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(ghostingId, sharingProc)));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void eraseSharingInfoUsingKey(BulkDataTester &bulkData, stk::mesh::EntityKey key, stk::mesh::BulkData::GHOSTING_ID ghostingId)
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
  STATE_GHOSTED,
  STATE_NOT_GHOSTED
};

inline bool check_state(const stk::mesh::BulkData & mesh, const stk::mesh::EntityKey & entityKey, EntityStates state,
                 int p0 = -1, int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1, int p5 = -1)
{
  // Check to see if the state is as expected for the provided EntityKey.
  //
  // Meaning of the optional processor arguments for various states:
  //     STATE_OWNED: Processor that owns the Entity
  //    STATE_SHARED: List of Processors that we share this entity with
  //   STATE_GHOSTED: Processor that we ghost the Entity from
  //
  std::vector<int> procs;
  if (p0 >= 0) {
    procs.push_back(p0);
  }
  if (p1 >= 0) {
    procs.push_back(p1);
  }
  if (p2 >= 0) {
    procs.push_back(p2);
  }
  if (p3 >= 0) {
    procs.push_back(p3);
  }
  if (p4 >= 0) {
    procs.push_back(p4);
  }
  if (p5 >= 0) {
    procs.push_back(p5);
  }
  std::sort(procs.begin(), procs.end());

  stk::mesh::Entity entity = mesh.get_entity(entityKey);
  std::ostringstream oss;

  switch (state) {
    case STATE_VALID:
    {
      if (!procs.empty()) {
        oss << "check_state(): Cannot provide processors with validity check." << std::endl;
      }
      if (!mesh.is_valid(entity)) {
        oss << "check_state(): Entity " << entityKey << " is not valid when it should have been." << std::endl;
      }
      break;
    }
    case STATE_NOT_VALID:
    {
      if (!procs.empty()) {
        oss << "check_state(): Cannot provide processors with STATE_NOT_VALID check." << std::endl;
      }
      if (mesh.is_valid(entity)) {
        oss << "check_state(): Entity " << entityKey << " is valid when it shouldn't have been." << std::endl;
      }
      break;
    }
    case STATE_OWNED:
    {
      if (procs.size() != 1u) {
        oss << "check_state(): Entities can have only one owner." << std::endl;
      }
      if (mesh.is_valid(entity)) {
        if (procs[0] != mesh.parallel_owner_rank(entity) ) {
          oss << "check_state(): Owner of entity " << entityKey << " was proc " << mesh.parallel_owner_rank(entity)
              << " and not proc " << procs[0] << std::endl;
        }
      }
      else {
        oss << "check_state(): Can't check ownership of locally-invalid entity." << std::endl;
      }
      break;
    }
    case STATE_SHARED:
    {
      if (procs.empty()) {
        oss << "check_state(): Must provide processor(s) with STATE_SHARED check." << std::endl;
      }
      stk::mesh::PairIterEntityComm comm_it = mesh.entity_comm_map_shared(entityKey);
      std::vector<int>::const_iterator procs_it = procs.begin();
      bool lists_match = true;

      if (comm_it.size() != procs.size()) {
        lists_match = false;
      }
      else {
        for ( ; procs_it != procs.end(); ++comm_it, ++procs_it) {
          int comm_proc = comm_it.first->proc;
          int user_proc = *procs_it;
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
        procs_it = procs.begin();
        for ( ; procs_it != procs.end(); ++procs_it) {
          oss << *procs_it << " ";
        }
        oss << ")" << std::endl;
      }

      break;
    }
    case STATE_NOT_SHARED:
    {
      if (!procs.empty()) {
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
    case STATE_GHOSTED:
    {
      if (procs.size() != 1) {
        oss << "check_state(): Must provide one processor with STATE_GHOSTED check." << std::endl;
      }
      if (!mesh.in_receive_ghost( mesh.aura_ghosting() , entityKey )) {
        oss << "check_state(): Entity " << entityKey << " was not ghosted from any proc when it should have" << std::endl
            << "               been ghosted from proc " << procs[0] << "." << std::endl;
      }
      else {
        const int owner_rank = mesh.entity_comm_map_owner(entityKey);
        if (owner_rank != procs[0]) {
          oss << "check_state(): Entity " << entityKey << " was ghosted from proc " << owner_rank << std::endl
              << "               when it should have been ghosted from proc " << procs[0] << "." << std::endl;
        }
      }
      break;
    }
    case STATE_NOT_GHOSTED:
    {
      if (!procs.empty()) {
        oss << "check_state(): Cannot provide processors with STATE_NOT_GHOSTED check." << std::endl;
      }
      if (mesh.in_receive_ghost( mesh.aura_ghosting() , entityKey )) {
        const int owner_rank = mesh.entity_comm_map_owner(entityKey);
        oss << "check_state(): Entity " << entityKey << " was ghosted from proc " << owner_rank << std::endl
            << "               when it shouldn't have been ghosted." << std::endl;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


}

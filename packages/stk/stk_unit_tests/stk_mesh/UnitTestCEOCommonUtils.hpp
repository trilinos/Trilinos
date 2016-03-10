#include <stddef.h>                     // for NULL
#include <algorithm>                    // for lower_bound
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <vector>                       // for vector
#include <stk_unit_test_utils/BulkDataTester.hpp>           // for BulkDataTester
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityCommListInfo.hpp"
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, etc
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Part; } }

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

inline bool isEntityValidOnCommList(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::Entity entity)
{
    EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), entityKey);
    return iter != stkMeshBulkData.my_internal_comm_list().end() && entityKey == iter->key;
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

bool check_state(const stk::unit_test_util::BulkDataTester & mesh, const EntityKey & entityKey, EntityStates state,
                 int p0 = -1, int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1, int p5 = -1);

bool check_parts(const stk::mesh::BulkData & mesh, const EntityKey & entityKey,
                 stk::mesh::Part *p0 = NULL, stk::mesh::Part *p1 = NULL, stk::mesh::Part *p2 = NULL,
                 stk::mesh::Part *p3 = NULL, stk::mesh::Part *p4 = NULL, stk::mesh::Part *p5 = NULL,
                 stk::mesh::Part *p6 = NULL, stk::mesh::Part *p7 = NULL, stk::mesh::Part *p8 = NULL,
                 stk::mesh::Part *p9 = NULL, stk::mesh::Part *p10 = NULL, stk::mesh::Part *p11 = NULL);

bool check_relns(const stk::mesh::BulkData & mesh, const EntityKey & entityKey, stk::mesh::EntityRank to_rank,
                        EntityId id0 = EntityKey::MIN_ID, EntityId id1 = EntityKey::MIN_ID,
                        EntityId id2 = EntityKey::MIN_ID, EntityId id3 = EntityKey::MIN_ID,
                        EntityId id4 = EntityKey::MIN_ID, EntityId id5 = EntityKey::MIN_ID,
                        EntityId id6 = EntityKey::MIN_ID, EntityId id7 = EntityKey::MIN_ID,
                        EntityId id8 = EntityKey::MIN_ID, EntityId id9 = EntityKey::MIN_ID,
                        EntityId id10 = EntityKey::MIN_ID, EntityId id11 = EntityKey::MIN_ID);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Mesh Creation Functions With Tests              /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////// 2Elem2ProcMove //////////////////////////////////////
void fillMeshfor2Elem2ProcMoveAndTest(stk::unit_test_util::BulkDataTester& bulk, stk::mesh::MetaData &meta, std::vector<stk::mesh::Entity>& elems);

void checkStatesAfterCEO_2Elem2ProcMove(stk::unit_test_util::BulkDataTester &bulk);

// valid vs not valid, owned.... shared vs not shared, ghosted_from vs not ghosted_from
// ghosted_to vs not_ghosted_to

void checkStatesAfterCEOME_2Elem2ProcMove(stk::unit_test_util::BulkDataTester &bulk);

//////////////////////////////////// 2Elem2ProcFlip //////////////////////////////////////

void fillMeshfor2Elem2ProcFlipAndTest(stk::unit_test_util::BulkDataTester& mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEOME_2Elem2ProcFlip(stk::unit_test_util::BulkDataTester& mesh);

void checkStatesAfterCEO_2Elem2ProcFlip(stk::unit_test_util::BulkDataTester& mesh);

//////////////////////////////////// 3Elem2ProcMoveRight //////////////////////////////////////

void fillMeshfor3Elem2ProcMoveRightAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector& elements);

void checkStatesAfterCEO_3Elem2ProcMoveRight(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_3Elem2ProcMoveRight(stk::unit_test_util::BulkDataTester &mesh);

//////////////////////////////////// 3Elem2ProcMoveLeft //////////////////////////////////////

void fillMeshfor3Elem2ProcMoveLeftAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data, stk::mesh::EntityVector &nodes, stk::mesh::EntityVector &elements);

void checkStatesAfterCEO_3Elem2ProcMoveLeft(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_3Elem2ProcMoveLeft(stk::unit_test_util::BulkDataTester &mesh);

//////////////////////////////////// 4Elem4ProcEdge //////////////////////////////////////

void fillMeshfor4Elem4ProcEdgeAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta_data,
        EntityKey &elem_key_chg_own);

void checkStatesAfterCEO_4Elem4ProcEdge(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_4Elem4ProcEdge(stk::unit_test_util::BulkDataTester &mesh);

//////////////////////////////////// 8Elem4ProcMoveTop //////////////////////////////////////

void fillMeshfor8Elem4ProcMoveTopAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEO_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_8Elem4ProcMoveTop(stk::unit_test_util::BulkDataTester &mesh);

//////////////////////////////////// 4Elem4ProcRotate //////////////////////////////////////

void fillMeshfor4Elem4ProcRotateAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEO_4Elem4ProcRotate(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEOME_4Elem4ProcRotate(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

//////////////////////////////////// 3Elem4Proc1Edge3D //////////////////////////////////////

void fillMeshfor3Elem4Proc1Edge3DAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEO_3Elem4Proc1Edge3D(stk::unit_test_util::BulkDataTester &mesh);

void checkStatesAfterCEOME_3Elem4Proc1Edge3D(stk::unit_test_util::BulkDataTester &mesh);

//these tests are for turning regenerate_aura off in various places

void checkStatesAfterCEOME_2Elem2ProcMove_no_ghost(stk::unit_test_util::BulkDataTester &bulk);

void fillMeshfor2Elem2ProcFlipAndTest_no_ghost(stk::unit_test_util::BulkDataTester& mesh, stk::mesh::MetaData &meta);

void checkStatesAfterCEOME_2Elem2ProcFlip_no_ghost(stk::unit_test_util::BulkDataTester& mesh);

} //namespace CEOUtils

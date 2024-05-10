#include "UnitTestCEOCommonUtils.hpp"
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <sstream>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <string>                       // for char_traits, operator<<
#include <stk_unit_test_utils/BulkDataTester.hpp>           // for BulkDataTester
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityCommDatabase.hpp"  // for EntityCommDatabase
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, operator<<, etc
#include "stk_mesh/base/FEMHelpers.hpp"  // for declare_element
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, get_cell_topology, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for EntityVector, EntityId, etc
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc

namespace CEOUtils
{

bool check_state(const stk::unit_test_util::BulkDataTester & mesh, const EntityKey & entityKey, EntityStates state,
                 int p0, int p1, int p2, int p3, int p4, int p5)
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
    std::vector<int> shared_procs;
    mesh.comm_shared_procs(entityKey,shared_procs);
    std::vector<int>::const_iterator expected_procs_it = expectedProcs.begin();
    bool lists_match = true;

    if (shared_procs.size() != expectedProcs.size()) {
      lists_match = false;
    }
    else {
      size_t shared_procs_i = 0;
      for ( ; expected_procs_it != expectedProcs.end(); ++shared_procs_i, ++expected_procs_it) {
        int comm_proc = shared_procs[shared_procs_i];
        int user_proc = *expected_procs_it;
        if (comm_proc != user_proc) {
          lists_match = false;
          break;
        }
      }
    }

    if (!lists_match) {
      oss << "check_state(): Entity " << entityKey << " was shared with procs (";
      for (size_t i=0 ; i<shared_procs.size() ; ++i) {
        int proc = shared_procs[i];
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
    std::vector<int> shared_procs;
    mesh.comm_shared_procs(entityKey,shared_procs);
    if (!shared_procs.empty()) {
      oss << "check_state(): Entity " << entityKey << " was shared with procs (";
      for (size_t i=0 ; i<shared_procs.size() ; ++i) {
        int proc = shared_procs[i];
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
    stk::mesh::Entity thisEntity = mesh.get_entity(entityKey);
    bool entityIsInvalid = !mesh.is_valid(thisEntity);
    bool inGhost = mesh.in_ghost(mesh.aura_ghosting() , entityKey , expectedProcs[0] );
    const int owner_rank = mesh.parallel_owner_rank(thisEntity);

    if ( entityIsInvalid && inGhost )
    {
      break;
    }
    else if ( entityIsInvalid && !inGhost )
    {
      oss << "check_state(): Entity " << entityKey << " was not ghosted from proc " << owner_rank << "." << std::endl;
      break;
    }
    if (!mesh.in_receive_ghost( mesh.aura_ghosting() , entityKey )) {
      oss << "check_state(): Entity " << entityKey << " was not ghosted from any proc when it should have" << std::endl
          << "               been ghosted from proc " << expectedProcs[0] << "." << std::endl;
    }
    else {
      const int this_owner_rank = mesh.my_internal_entity_comm_map_owner(entityKey);
      if (this_owner_rank != expectedProcs[0]) {
        oss << "check_state(): Entity " << entityKey << " was ghosted from proc " << this_owner_rank << std::endl
            << "               when it should have been ghosted from proc " << expectedProcs[0] << "." << std::endl;
      }
    }
    break;
  }
  case STATE_NOT_GHOSTED_FROM:
  {
    stk::mesh::Entity thisEntity = mesh.get_entity(entityKey);
    bool entityIsInvalid = !mesh.is_valid(thisEntity);
    const int owner_rank = mesh.parallel_owner_rank(thisEntity);

    if ( entityIsInvalid && owner_rank == stk::mesh::InvalidProcessRank)
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
      const int this_owner_rank = mesh.my_internal_entity_comm_map_owner(entityKey);
      oss << "check_state(): Entity " << entityKey << " was ghosted from proc " << this_owner_rank << std::endl
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
    mesh.comm_procs(mesh.aura_ghosting(), entityKey, meshProcs);

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
    std::vector<int> auraProcs;
    mesh.comm_procs(mesh.aura_ghosting(), entityKey, auraProcs);
    std::vector<int> meshProcs;
    for ( size_t i=0; i<auraProcs.size(); i++ ) {
      if (mesh.parallel_owner_rank(entity) == p_rank) {
        meshProcs.push_back(auraProcs[i]);
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

bool check_parts(const stk::mesh::BulkData & mesh, const EntityKey & entityKey,
                 stk::mesh::Part *p0, stk::mesh::Part *p1, stk::mesh::Part *p2,
                 stk::mesh::Part *p3, stk::mesh::Part *p4, stk::mesh::Part *p5,
                 stk::mesh::Part *p6, stk::mesh::Part *p7, stk::mesh::Part *p8,
                 stk::mesh::Part *p9, stk::mesh::Part *p10, stk::mesh::Part *p11)
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

bool check_relns(const stk::mesh::BulkData & mesh, const EntityKey & entityKey, stk::mesh::EntityRank to_rank,
                 EntityId id0, EntityId id1,
                 EntityId id2, EntityId id3,
                 EntityId id4, EntityId id5,
                 EntityId id6, EntityId id7,
                 EntityId id8, EntityId id9,
                 EntityId id10, EntityId id11)
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


//////////////////////////////////// 3Elem4Proc1Edge3D //////////////////////////////////////

void fillMeshfor3Elem4Proc1Edge3DAndTest(stk::unit_test_util::BulkDataTester &mesh, stk::mesh::MetaData &meta)
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
  stk::mesh::EntityIdVector elem_nodeIDs[] {
    {1, 2, 5, 4, 9, 10, 13, 12},
    {2, 3, 6, 5, 10, 11, 14, 13},
    {5, 6, 8, 7, 13, 14, 16, 15}
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
    stk::mesh::Entity edge = mesh.declare_edge(1, stk::mesh::ConstPartVector{edge_part});
    std::vector<stk::mesh::Entity> nodes;
    nodes.push_back(mesh.get_entity(NODE_RANK, 5));
    nodes.push_back(mesh.get_entity(NODE_RANK, 13 ));
    mesh.declare_relation(edge, nodes[0], 0);
    mesh.declare_relation(edge, nodes[1], 1);
    stk::mesh::impl::connectUpwardEntityToEntity(mesh, elem, edge, nodes.data());
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

  stk::mesh::Part * hex_topo_part  = &meta.get_topology_root_part(stk::topology::HEXAHEDRON_8);
  stk::mesh::Part * edge_topo_part = &meta.get_topology_root_part(stk::topology::LINE_2);
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

void checkStatesAfterCEO_3Elem4Proc1Edge3D(stk::unit_test_util::BulkDataTester &mesh)
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
  stk::mesh::Part * hex_topo_part  = &meta.get_topology_root_part(stk::topology::HEXAHEDRON_8);
  stk::mesh::Part * edge_topo_part = &meta.get_topology_root_part(stk::topology::LINE_2);
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

void checkStatesAfterCEOME_3Elem4Proc1Edge3D(stk::unit_test_util::BulkDataTester &mesh)
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
  stk::mesh::Part * hex_topo_part  = &meta.get_topology_root_part(stk::topology::HEXAHEDRON_8);
  stk::mesh::Part * edge_topo_part = &meta.get_topology_root_part(stk::topology::LINE_2);
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

} //namespace CEOUtils

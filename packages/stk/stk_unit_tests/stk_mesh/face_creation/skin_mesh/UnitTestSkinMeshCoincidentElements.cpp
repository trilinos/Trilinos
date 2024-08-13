#include "mpi.h"                        // for MPI_Comm, etc
#include <stddef.h>                     // for size_t
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for TEST_F, ASSERT_EQ, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityIdVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

namespace stk { namespace mesh { class Part; } }

namespace
{

int get_other_proc(MPI_Comm comm)
{
  return stk::parallel_machine_size(comm) - stk::parallel_machine_rank(comm) - 1;
}

// void create_exposed_boundary_sides(BulkData &bulkData, const Selector& blocksToSkin, Part& partToPutSidesInto)

class CoincidentElements: public stk::unit_test_util::MeshTestFixture
{
protected:
  void make_coincident_element_mesh(unsigned numElemsToCreate, const stk::mesh::EntityIdVector &nodes, stk::mesh::Part &part)
  {
    get_bulk().modification_begin();
    declare_coincident_elements_round_robin(numElemsToCreate, nodes, part);
    make_nodes_shared(nodes);
    get_bulk().modification_end();
  }
  void run_skin_mesh()
  {
    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().universal_part(), stk::mesh::PartVector());
  }
  void expect_faces_connected_to_num_elements_locally(const std::vector<unsigned> &goldNumConnectedElems)
  {
    stk::mesh::EntityVector faces = get_all_faces();
    ASSERT_EQ(goldNumConnectedElems.size(), faces.size());
    for(size_t i=0; i<faces.size(); i++)
      EXPECT_EQ(goldNumConnectedElems[i], get_bulk().num_elements(faces[i])) << i;
  }
  stk::mesh::EntityVector get_all_faces()
  {
    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::FACE_RANK), faces);
    return faces;
  }
  void make_nodes_shared(const stk::mesh::EntityIdVector &nodes)
  {
    for(const stk::mesh::EntityId nodeId : nodes)
      get_bulk().add_node_sharing(get_bulk().get_entity(stk::topology::NODE_RANK, nodeId), get_other_proc(get_comm()));
  }
private:
  bool is_element_created_on_this_proc(int elementIndex)
  {
    return (elementIndex % get_bulk().parallel_size() == get_bulk().parallel_rank());
  }
  void declare_coincident_elements_round_robin(unsigned numElemsToCreate, const stk::mesh::EntityIdVector &nodes, stk::mesh::Part &part)
  {
    for(unsigned i = 0; i < numElemsToCreate; i++)
      if(is_element_created_on_this_proc(i))
        stk::mesh::declare_element(get_bulk(), part, i+1, nodes);
  }
};

class CoincidentQuad4Shells : public CoincidentElements
{
protected:
  void skin_num_coincident_shells(unsigned numElemsToCreate, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    create_stacked_shells(numElemsToCreate);
    run_skin_mesh();
  }
  void create_stacked_shells(unsigned numElemsToCreate)
  {
    stk::mesh::Part &shellPart = get_meta().get_topology_root_part(stk::topology::SHELL_QUADRILATERAL_4);
    stk::mesh::EntityIdVector nodes = {1, 2, 3, 4};
    make_coincident_element_mesh(numElemsToCreate, nodes, shellPart);
  }
  void expect_local_face_shared_to_other_proc(const std::vector<int> &goldSharedProcs, const stk::mesh::EntityIdVector &goldFaceIds)
  {
    stk::mesh::EntityVector faces = get_all_faces();
    ASSERT_EQ(goldSharedProcs.size(), faces.size());
    for(size_t i=0; i<faces.size(); i++)
      expect_face_identifier_and_sharing_proc(faces[i], goldFaceIds[i], goldSharedProcs[i]);
  }
  void expect_face_identifier_and_sharing_proc(stk::mesh::Entity face, stk::mesh::EntityId goldFaceId, int goldSharedProc)
  {
    EXPECT_EQ(goldFaceId, get_bulk().identifier(face));
    expect_sharing_proc(face, goldSharedProc);
  }
  void expect_sharing_proc(stk::mesh::Entity face, int goldSharedProc)
  {
    std::vector<int> procs;
    get_bulk().comm_shared_procs(get_bulk().entity_key(face), procs);
    ASSERT_EQ(1u, procs.size());
    EXPECT_EQ(goldSharedProc, procs[0]);
  }
};

class TwoCoincidentQuad4Shells : public CoincidentQuad4Shells
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    unsigned numElemsToCreate = 2;
    skin_num_coincident_shells(numElemsToCreate, auraOption);
    expect_faces_connected_to_num_elements_locally({numElemsToCreate, numElemsToCreate});
  }
};
TEST_F(TwoCoincidentQuad4Shells, Skin)
{
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}

class ThreeCoincidentQuad4Shells : public CoincidentQuad4Shells
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    unsigned numElemsToCreate = 3;
    skin_num_coincident_shells(numElemsToCreate, auraOption);
    expect_faces_connected_to_num_elements_locally({numElemsToCreate, numElemsToCreate});
  }
};
TEST_F(ThreeCoincidentQuad4Shells, Skin)
{
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}

class TwoCoincidentQuad4ShellsInParallel : public CoincidentQuad4Shells
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    skin_num_coincident_shells(2, auraOption);
    expect_faces_connected_to_num_elements_locally({1, 1});
    int otherProc = get_other_proc(get_comm());
    expect_local_face_shared_to_other_proc({otherProc, otherProc}, {1, 2});
  }
};
// disabled due to split coincident elements
TEST_F(TwoCoincidentQuad4ShellsInParallel, DISABLED_Skin)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}


class CoincidentHex8s : public CoincidentElements
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    create_two_coincident_hexes();
    run_skin_mesh();
    expect_faces_connected_to_num_elements_locally({2, 2, 2, 2, 2, 2});
  }
  void create_two_coincident_hexes()
  {
    block1 = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::mesh::EntityIdVector nodes = {1, 2, 3, 4, 5, 6, 7, 8};
    declare_num_coincident_elements(2, nodes, *block1);
  }
  void declare_num_coincident_elements(unsigned numElemsToCreate, const stk::mesh::EntityIdVector &nodes, stk::mesh::Part &part)
  {
    get_bulk().modification_begin();
    for(unsigned i = 0; i < numElemsToCreate; i++)
      stk::mesh::declare_element(get_bulk(), part, i+1, nodes);
    get_bulk().modification_end();
  }
  stk::mesh::Part *block1;
};
TEST_F(CoincidentHex8s, SkinMesh)
{
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}

class CoincidentHex8sWithAdjacentHex : public CoincidentHex8s
{
protected:
  void create_coincident_hex8s_with_adjacent_hex(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    block2 = &get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    create_two_coincident_hexes();
    create_adjacent_hex();
  }
private:
  void create_adjacent_hex()
  {
    stk::mesh::EntityIdVector nodes = {5, 6, 7, 8, 9, 10, 11, 12};
    get_bulk().modification_begin();
    stk::mesh::declare_element(get_bulk(), *block2, 3, nodes);
    get_bulk().modification_end();
  }
protected:
  stk::mesh::Part *block2;
};

class CoincidentHex8sWithAdjacentHexSerial : public CoincidentHex8sWithAdjacentHex
{
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_coincident_hex8s_with_adjacent_hex(auraOption);
    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_meta().universal_part(), {});
    expect_faces_connected_to_num_elements_locally({2, 2, 2, 2, 2, 1, 1, 1, 1, 1});
  }
};
TEST_F(CoincidentHex8sWithAdjacentHexSerial, Skin)
{
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}

class CoincidentHex8sWithAdjacentAirHex : public CoincidentHex8sWithAdjacentHex
{
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_coincident_hex8s_with_adjacent_hex(auraOption);
    stk::mesh::Selector air = *block2;
    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), *block1, {}, air);

    expect_faces_connected_to_num_elements_locally({2, 2, 2, 2, 2, 3});
  }
};
TEST_F(CoincidentHex8sWithAdjacentAirHex, Skin)
{
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}

class Hex8WithAdjacentCoincidentAirHex8s : public CoincidentHex8sWithAdjacentHex
{
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_coincident_hex8s_with_adjacent_hex(auraOption);

    stk::mesh::Selector air = *block1;
    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), *block2, {}, air);
    expect_faces_connected_to_num_elements_locally({3, 1, 1, 1, 1, 1});
  }
};
TEST_F(Hex8WithAdjacentCoincidentAirHex8s, Skin)
{
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}

class CoincidentHex8sWithAdjacentHexInParallel : public CoincidentHex8sWithAdjacentHex
{
protected:
  void create_coincident_hex8s_with_adjacent_hex_on_2_procs(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                                            const stk::mesh::EntityIdVector &ids)
  {
    setup_empty_mesh(auraOption);
    create_block_parts();
    fill_mesh_with_coincident_hex8s_and_adjacent_hex_in_parallel(ids);
  }
  void fill_mesh_with_coincident_hex8s_and_adjacent_hex_in_parallel(const stk::mesh::EntityIdVector &ids)
  {
    get_bulk().modification_begin();
    create_coincident_hexes_on_proc0_and_adjacent_hex_on_proc1(ids);
    make_nodes_shared({5, 6, 7, 8});
    get_bulk().modification_end();
  }
  void create_block_parts()
  {
    block1 = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    block2 = &get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
  }
  void create_coincident_hexes_on_proc0_and_adjacent_hex_on_proc1(const stk::mesh::EntityIdVector &ids)
  {
    if(get_bulk().parallel_rank() == 0)
      create_coincident_hexes(ids[0], ids[1]);
    else
      stk::mesh::declare_element(get_bulk(), *block2, ids[2], {5, 6, 7, 8, 9, 10, 11, 12});
  }
  void create_coincident_hexes(stk::mesh::EntityId id1, stk::mesh::EntityId id2)
  {
    stk::mesh::EntityIdVector nodes = {1, 2, 3, 4, 5, 6, 7, 8};
    stk::mesh::declare_element(get_bulk(), *block1, id1, nodes);
    stk::mesh::declare_element(get_bulk(), *block1, id2, nodes);
  }
  void expect_faces_connected_to_num_elements_locally_per_proc(const std::vector<unsigned> &goldNumConnectedElemsProc0,
                                                               const std::vector<unsigned> &goldNumConnectedElemsProc1)
  {
    if(get_bulk().parallel_rank() == 0)
      expect_faces_connected_to_num_elements_locally(goldNumConnectedElemsProc0);
    else
      expect_faces_connected_to_num_elements_locally(goldNumConnectedElemsProc1);
  }
  void skin_part_with_part2_as_air(stk::mesh::Selector partToSkin, stk::mesh::Selector partToConsiderAsAir)
  {
    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), partToSkin, {}, partToConsiderAsAir);
  }
  stk::mesh::Part *block2;
};

class CoincidentHex8sWithAdjacentAirHexInParallel : public CoincidentHex8sWithAdjacentHexInParallel
{
protected:
  void run_with_element_ids(stk::mesh::BulkData::AutomaticAuraOption auraOption, const stk::mesh::EntityIdVector &ids)
  {
    create_coincident_hex8s_with_adjacent_hex_on_2_procs(auraOption, ids);
    skin_part_with_part2_as_air(*block1, *block2);
    expect_faces_connected_to_num_elements_locally_per_proc( {2, 2, 2, 2, 2, 2}, {1});
  }
};
TEST_F(CoincidentHex8sWithAdjacentAirHexInParallel, SkinHex1Hex2)
{
  if(stk::parallel_machine_size(get_comm()) == 2)
    run_with_element_ids(stk::mesh::BulkData::NO_AUTO_AURA, {1, 2, 3});
}
TEST_F(CoincidentHex8sWithAdjacentAirHexInParallel, SkinHex2Hex1)
{
  if(stk::parallel_machine_size(get_comm()) == 2)
    run_with_element_ids(stk::mesh::BulkData::NO_AUTO_AURA, {2, 1, 3});
}

class Hex8WithAdjacentCoincidentAirHex8sInParallel : public CoincidentHex8sWithAdjacentHexInParallel
{
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    create_coincident_hex8s_with_adjacent_hex_on_2_procs(auraOption, {1, 2, 3});
    skin_part_with_part2_as_air(*block2, *block1);
    expect_faces_connected_to_num_elements_locally_per_proc({2}, {1, 1, 1, 1, 1, 1});
  }
};
TEST_F(Hex8WithAdjacentCoincidentAirHex8sInParallel, Skin)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

}

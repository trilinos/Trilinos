#include <gtest/gtest.h>                // for AssertHelper, EXPECT_TRUE, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <vector>                       // for vector, vector<>::reference
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_io/FillMesh.hpp"
#include "stk_mesh/baseImpl/MeshModLogObserver.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size
namespace stk { namespace mesh { class Ghosting; } }

namespace
{

using stk::unit_test_util::build_mesh;

class MeshModLogTest : public ::testing::Test
{
protected:
  MeshModLogTest()
    : bulkPtr(build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::AUTO_AURA)),
      meta(bulkPtr->mesh_meta_data()),
      bulk(*bulkPtr),
      ostrm()
  {
  }

  void setup(const stk::mesh::EntityKey& key,
             stk::mesh::impl::PatchType patchType=stk::mesh::impl::NO_PATCH,
             unsigned startAtModCycle = 1)
  {
    STK_ThrowRequireMsg(bulk.parallel_size() <= 4, "Must use no more than 4 procs for this test.");

    observer = std::make_shared<stk::mesh::impl::MeshModLogObserver>(bulk, key, ostrm, patchType, startAtModCycle);
    bulk.register_observer(observer);

    stk::io::fill_mesh("generated:1x1x4|sideset:xXyYzZ", bulk);
  }

  virtual ~MeshModLogTest()
  {
  }

  void delete_elem1()
  {
    bulk.modification_begin();

    stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
    if (bulk.is_valid(elem1) && bulk.bucket(elem1).owned()) {
      EXPECT_TRUE(bulk.destroy_entity(elem1));
    }

    bulk.modification_end();
  }

  void replace_node1()
  {
    bulk.modification_begin();

    if (bulk.parallel_rank() == 0) {
      stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
      stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
      stk::mesh::Entity face11 = bulk.get_entity(stk::topology::FACE_RANK, 11);
      ASSERT_TRUE(bulk.is_valid(node1));
      ASSERT_TRUE(bulk.is_valid(elem1));
      ASSERT_TRUE(bulk.is_valid(face11));

      stk::mesh::Entity node100 = bulk.declare_node(100);

      stk::mesh::ConnectivityOrdinal ordinal = 0;

      EXPECT_TRUE(bulk.destroy_relation(elem1, node1, ordinal));
      EXPECT_TRUE(bulk.destroy_relation(face11, node1, ordinal));

      bulk.declare_relation(elem1, node100, ordinal);
      bulk.declare_relation(face11, node100, ordinal);
    }

    bulk.modification_end();
  }

protected:
  std::shared_ptr<stk::mesh::BulkData> bulkPtr;
  stk::mesh::MetaData& meta;
  stk::mesh::BulkData& bulk;
  std::shared_ptr<stk::mesh::impl::MeshModLogObserver> observer;
  std::ostringstream ostrm;
};

TEST_F(MeshModLogTest, creation)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  setup(stk::mesh::EntityKey(stk::topology::FACE_RANK,11));

  std::string str = ostrm.str();

  if (bulk.parallel_rank() == 0) {
    const std::string expected(
      "modification_begin mod-cycle=0\n"
      "P0 mod-cycle=1, declare_entity (FACE_RANK,11)\n"
      "P0 mod-cycle=1, declare_relation (FACE_RANK,11) -> (NODE_RANK,1) ordinal=0\n"
      "P0 mod-cycle=1, declare_relation (NODE_RANK,1) -> (FACE_RANK,11) ordinal=0\n"
      "P0 mod-cycle=1, declare_relation (FACE_RANK,11) -> (NODE_RANK,2) ordinal=1\n"
      "P0 mod-cycle=1, declare_relation (NODE_RANK,2) -> (FACE_RANK,11) ordinal=1\n"
      "P0 mod-cycle=1, declare_relation (FACE_RANK,11) -> (NODE_RANK,6) ordinal=2\n"
      "P0 mod-cycle=1, declare_relation (NODE_RANK,6) -> (FACE_RANK,11) ordinal=2\n"
      "P0 mod-cycle=1, declare_relation (FACE_RANK,11) -> (NODE_RANK,5) ordinal=3\n"
      "P0 mod-cycle=1, declare_relation (NODE_RANK,5) -> (FACE_RANK,11) ordinal=3\n"
      "P0 mod-cycle=1, declare_relation (ELEMENT_RANK,1) -> (FACE_RANK,11) ordinal=0\n"
      "P0 mod-cycle=1, declare_relation (FACE_RANK,11) -> (ELEMENT_RANK,1) ordinal=0\n"
      "start modification_end mod-cycle=1\n"
      "P0 (FACE_RANK,11) {P0Created}, elems{1} edges{} nodes{1 2 6 5} \n"
      "finish modification_end mod-cycle=1\n"
      "P0 (FACE_RANK,11) {P0Created}, elems{1} edges{} nodes{1 2 6 5} \n"
    );
    EXPECT_EQ(expected, str);
  }
}

TEST_F(MeshModLogTest, deletion)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  setup(stk::mesh::EntityKey(stk::topology::ELEM_RANK,1));
  observer->add_key_to_watch(stk::mesh::EntityKey(stk::topology::ELEM_RANK,2));

  ostrm.str("");

  delete_elem1();

  std::string str = ostrm.str();

  if (bulk.parallel_rank() == 0) {
    const std::string expected(
      "modification_begin mod-cycle=1\n"
      "P0 (ELEMENT_RANK,2) {P0Created}, faces{21 22 23 24} edges{} nodes{5 6 8 7 9 10 12 11} \n"
      "P0 (ELEMENT_RANK,1) {P0Created}, faces{11 12 13 14 15} edges{} nodes{1 2 4 3 5 6 8 7} \n"
      "P0 mod-cycle=2, destroy_entity (ELEMENT_RANK,1)\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,15) ordinal=4\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,14) ordinal=3\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,13) ordinal=2\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,12) ordinal=1\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,11) ordinal=0\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,7) ordinal=7\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,8) ordinal=6\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,6) ordinal=5\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,5) ordinal=4\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,3) ordinal=3\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,4) ordinal=2\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,2) ordinal=1\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,1) ordinal=0\n"
      "start modification_end mod-cycle=2\n"
      "P0 (ELEMENT_RANK,2) {P0Modified}, faces{21 22 23 24} edges{} nodes{5 6 8 7 9 10 12 11} \n"
      "finish modification_end mod-cycle=2\n"
      "P0 (ELEMENT_RANK,2) {P0Modified}, faces{21 22 23 24} edges{} nodes{5 6 8 7 9 10 12 11} \n"
    );
    EXPECT_EQ(expected, str);
  }
}

TEST_F(MeshModLogTest, deletion_modcycle)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  const unsigned modCycle = 2;
  setup(stk::mesh::EntityKey(stk::topology::ELEM_RANK,1), stk::mesh::impl::NO_PATCH, modCycle);

  delete_elem1();

  std::string str = ostrm.str();

  if (bulk.parallel_rank() == 0) {
    const std::string expected(
      "modification_begin mod-cycle=1\n"
      "P0 (ELEMENT_RANK,1) {P0Created}, faces{11 12 13 14 15} edges{} nodes{1 2 4 3 5 6 8 7} \n"
      "P0 mod-cycle=2, destroy_entity (ELEMENT_RANK,1)\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,15) ordinal=4\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,14) ordinal=3\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,13) ordinal=2\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,12) ordinal=1\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (FACE_RANK,11) ordinal=0\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,7) ordinal=7\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,8) ordinal=6\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,6) ordinal=5\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,5) ordinal=4\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,3) ordinal=3\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,4) ordinal=2\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,2) ordinal=1\n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,1) ordinal=0\n"
      "start modification_end mod-cycle=2\n"
      "finish modification_end mod-cycle=2\n"
    );
    EXPECT_EQ(expected, str);
  }
}

TEST_F(MeshModLogTest, relations)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  setup(stk::mesh::EntityKey(stk::topology::ELEM_RANK,1));

  ostrm.str("");

  replace_node1();

  std::string str = ostrm.str();

  if (bulk.parallel_rank() == 0) {
    const std::string expected(
      "modification_begin mod-cycle=1\n"
      "P0 (ELEMENT_RANK,1) {P0Created}, faces{11 12 13 14 15} edges{} nodes{1 2 4 3 5 6 8 7} \n"
      "P0 mod-cycle=2, destroy_relation (ELEMENT_RANK,1) -> (NODE_RANK,1) ordinal=0\n"
      "P0 mod-cycle=2, declare_relation (ELEMENT_RANK,1) -> (NODE_RANK,100) ordinal=0\n"
      "P0 mod-cycle=2, declare_relation (NODE_RANK,100) -> (ELEMENT_RANK,1) ordinal=0\n"
      "start modification_end mod-cycle=2\n"
      "P0 (ELEMENT_RANK,1) {P0Modified}, faces{11 12 13 14 15} edges{} nodes{100 2 4 3 5 6 8 7} \n"
      "finish modification_end mod-cycle=2\n"
      "P0 (ELEMENT_RANK,1) {P0Modified}, faces{11 12 13 14 15} edges{} nodes{100 2 4 3 5 6 8 7} \n"
    );
    EXPECT_EQ(expected, str);
  }
}

TEST_F(MeshModLogTest, upDownPatch)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  setup(stk::mesh::EntityKey(stk::topology::FACE_RANK,11), stk::mesh::impl::UP_DOWN_PATCH);

  stk::mesh::EntityVector entityPatch = observer->get_entity_patch();

  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ(14u, entityPatch.size());
  }
  else {
    EXPECT_EQ(0u, entityPatch.size());
  }
}

TEST_F(MeshModLogTest, upwardPatch_face)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  setup(stk::mesh::EntityKey(stk::topology::FACE_RANK,11), stk::mesh::impl::UPWARD_PATCH);

  stk::mesh::EntityVector entityPatch = observer->get_entity_patch();

  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ(2u, entityPatch.size());
  }
  else {
    EXPECT_EQ(0u, entityPatch.size());
  }
}

TEST_F(MeshModLogTest, upwardPatch_node)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  setup(stk::mesh::EntityKey(stk::topology::NODE_RANK,1), stk::mesh::impl::UPWARD_PATCH);

  stk::mesh::EntityVector entityPatch = observer->get_entity_patch();

  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ(5u, entityPatch.size());
  }
  else {
    EXPECT_EQ(0u, entityPatch.size());
  }
}

TEST_F(MeshModLogTest, auraClosurePatch)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  setup(stk::mesh::EntityKey(stk::topology::FACE_RANK,11), stk::mesh::impl::AURA_CLOSURE_PATCH);

  stk::mesh::EntityVector entityPatch = observer->get_entity_patch();

  if (bulk.parallel_rank() == 0) {
    EXPECT_EQ(23u, entityPatch.size());
  }
  else {
    EXPECT_EQ(0u, entityPatch.size());
  }
}

}

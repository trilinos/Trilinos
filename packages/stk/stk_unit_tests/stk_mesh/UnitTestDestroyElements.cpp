#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/DestroyElements.hpp>

namespace
{

void expect_valid(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities, int line)
{
  for(stk::mesh::Entity entity : entities) {
    EXPECT_TRUE(bulk.is_valid(entity)) << "from line " << line;
  }
}

void expect_invalid(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities, int line)
{
  for(stk::mesh::Entity entity : entities) {
    EXPECT_FALSE(bulk.is_valid(entity)) << "from line " << line;
  }
}

void expect_not_shared(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities, int line)
{
  for(stk::mesh::Entity entity : entities) {
    EXPECT_TRUE(bulk.is_valid(entity)) << "from line " << line;
    EXPECT_FALSE(bulk.bucket(entity).shared()) << "from line " << line;
    std::vector<int> procs;
    bulk.comm_shared_procs(bulk.entity_key(entity), procs);
    EXPECT_TRUE(procs.empty()) << "from line " << line;
  }
}

void expect_shared(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities, int line)
{
  for(stk::mesh::Entity entity : entities) {
    EXPECT_TRUE(bulk.is_valid(entity)) << "from line " << line;
    EXPECT_TRUE(bulk.bucket(entity).shared()) << "from line " << line;
    std::vector<int> procs;
    bulk.comm_shared_procs(bulk.entity_key(entity), procs);
    EXPECT_TRUE(!procs.empty()) << "from line " << line;
  }
}

stk::mesh::EntityVector get_faces_for_entity(const stk::mesh::BulkData &bulk, const stk::mesh::Entity entity)
{
  unsigned numConnected = bulk.num_connectivity(entity, stk::topology::FACE_RANK);
  const stk::mesh::Entity* connectedEntities = bulk.begin(entity, stk::topology::FACE_RANK);
  stk::mesh::EntityVector entityFaces(connectedEntities, connectedEntities+numConnected);
  return entityFaces;
}

class HexMesh : public stk::unit_test_util::MeshTestFixture
{
protected:
  HexMesh()
  {
  }

  void initialize_my_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    set_spatial_dimension(3);
    setup_empty_mesh(auraOption);
    std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
        0,2,HEX_8,2,9,10,3,6,11,12,7";
        stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    initialize_my_mesh(auraOption);

    stk::mesh::Selector ownedNotShared = get_meta().locally_owned_part() & !get_meta().globally_shared_part();
    stk::mesh::EntityVector orphanedNodes;
    if (get_bulk().parallel_size() == 2 && get_bulk().parallel_rank() == procWithDelete)
      stk::mesh::get_selected_entities(ownedNotShared, get_bulk().buckets(stk::topology::NODE_RANK), orphanedNodes);

    stk::mesh::EntityVector sharedNodes;
    stk::mesh::get_selected_entities(get_meta().globally_shared_part(), get_bulk().buckets(stk::topology::NODE_RANK), sharedNodes);

    if (get_bulk().parallel_size() ==2)
      expect_shared(get_bulk(), sharedNodes, __LINE__);

    stk::mesh::EntityVector elementToDestroy;
    if (get_bulk().parallel_size() == 2 && get_bulk().parallel_rank() == procWithDelete)
      stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), elementToDestroy);

    if (get_bulk().parallel_size() == 2 && get_bulk().parallel_rank() == procWithDelete)
    {
      expect_valid(get_bulk(), orphanedNodes, __LINE__);
      expect_valid(get_bulk(), elementToDestroy, __LINE__);
    }

    stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

    if (get_bulk().parallel_size() == 2 && get_bulk().parallel_rank() == procWithDelete)
    {
      expect_invalid(get_bulk(), orphanedNodes, __LINE__);
      expect_invalid(get_bulk(), elementToDestroy, __LINE__);
      expect_invalid(get_bulk(), sharedNodes, __LINE__);
    }
    if (get_bulk().parallel_size() == 2 && get_bulk().parallel_rank() != procWithDelete)
      expect_not_shared(get_bulk(), sharedNodes, __LINE__);
  }
  int procWithDelete = -1;
};

TEST_F(HexMesh, DeleteOneElement)
{
  procWithDelete = 0;
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(HexMesh, DeleteOnProcZeroWithSharedNodes_withAura)
{
  procWithDelete = 0;
  run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(HexMesh, DeleteOnProcZeroWithSharedNodes_NoAura)
{
  procWithDelete = 0;
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(HexMesh, DeleteOnProcOneWithSharedNodes_withAura)
{
  procWithDelete = 1;
  run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(HexMesh, DeleteOnProcOneWithSharedNodes_NoAura)
{
  procWithDelete = 1;
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

class TetMesh : public stk::unit_test_util::MeshFixture
{
protected:
  TetMesh()
  {
  }

  void initialize_my_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    set_spatial_dimension(3);
    setup_empty_mesh(auraOption);
    std::string meshDesc =
        "0,1,TET_4,1,2,3,4\n\
        0,2,TET_4,2,5,3,4";

        if(get_bulk().parallel_size() == 2)
    {
      meshDesc =  "0,1,TET_4,1,2,3,4\n\
          1,2,TET_4,2,5,3,4";
    }
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }
};

TEST_F(TetMesh, DeleteOneElement)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    initialize_my_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);

    stk::mesh::EntityVector orphanedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 5)
    };

    stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

    stk::mesh::EntityVector facesOfDestroyedElement = get_faces_for_entity(get_bulk(), elementToDestroy[0]);
    EXPECT_EQ(facesOfDestroyedElement.size(), 4u);
    for(stk::mesh::Entity face : facesOfDestroyedElement)
      EXPECT_TRUE(get_bulk().is_valid(face));

    expect_valid(get_bulk(), orphanedNodes, __LINE__);
    expect_valid(get_bulk(), elementToDestroy, __LINE__);

    stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

    expect_invalid(get_bulk(), orphanedNodes, __LINE__);
    expect_invalid(get_bulk(), elementToDestroy, __LINE__);

    unsigned numValid = 0;
    for(stk::mesh::Entity face : facesOfDestroyedElement)
      if(get_bulk().is_valid(face))
        numValid++;

    EXPECT_EQ(1u, numValid);
  }
}

TEST_F(TetMesh, DeleteOneElement_DontDeleteAnyOrphans)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    initialize_my_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);

    stk::mesh::EntityVector orphanedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 5)
    };

    stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

    stk::mesh::EntityVector facesOfDestroyedElement = get_faces_for_entity(get_bulk(), elementToDestroy[0]);
    EXPECT_EQ(facesOfDestroyedElement.size(), 4u);
    for(stk::mesh::Entity face : facesOfDestroyedElement)
      EXPECT_TRUE(get_bulk().is_valid(face));

    expect_valid(get_bulk(), orphanedNodes, __LINE__);
    expect_valid(get_bulk(), elementToDestroy, __LINE__);

    stk::mesh::Selector selectNoOrphans;
    stk::mesh::destroy_elements(get_bulk(), elementToDestroy, selectNoOrphans);

    expect_valid(get_bulk(), orphanedNodes, __LINE__);
    expect_invalid(get_bulk(), elementToDestroy, __LINE__);

    unsigned numValid = 0;
    for(stk::mesh::Entity face : facesOfDestroyedElement)
      if(get_bulk().is_valid(face))
        numValid++;

    EXPECT_EQ(4u, numValid);
  }
}

TEST_F(TetMesh, DeleteElement_afterChangeParts)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
  {
    initialize_my_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);

    stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

    expect_valid(get_bulk(), elementToDestroy, __LINE__);

    stk::mesh::Part* block1 = get_meta().get_part("block_TETRAHEDRON_4"); //huh?? why not 'block_1'?
    ASSERT_TRUE(block1 != nullptr);
    stk::mesh::PartVector empty;
    stk::mesh::PartVector rmParts = {block1};

    get_bulk().modification_begin();
    get_bulk().change_entity_parts(elementToDestroy, empty, rmParts);
    get_bulk().modification_end(stk::mesh::ModEndOptimizationFlag::MOD_END_NO_SORT);

    stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

    expect_invalid(get_bulk(), elementToDestroy, __LINE__);
  }
}

TEST_F(TetMesh, DeleteElementOnProcBoundaryWithOwnedFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    initialize_my_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);

    stk::mesh::EntityVector orphanedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 1)
    };

    stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1)};
    stk::mesh::EntityVector facesOfDestroyedElement;

    if(get_parallel_rank() == 0)
    {
      facesOfDestroyedElement = get_faces_for_entity(get_bulk(), elementToDestroy[0]);
      EXPECT_EQ(facesOfDestroyedElement.size(), 4u);

      for(stk::mesh::Entity face : facesOfDestroyedElement)
        EXPECT_TRUE(get_bulk().is_valid(face));

      expect_valid(get_bulk(), orphanedNodes, __LINE__);
      expect_valid(get_bulk(), elementToDestroy, __LINE__);
    }
    else if(get_parallel_rank() == 1)
    {
      stk::mesh::Entity element2 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2);
      EXPECT_TRUE(get_bulk().is_valid(element2));
    }
    stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

    expect_invalid(get_bulk(), orphanedNodes, __LINE__);
    expect_invalid(get_bulk(), elementToDestroy, __LINE__);

    if(get_parallel_rank() == 0)
    {
      for(stk::mesh::Entity face : facesOfDestroyedElement)
        EXPECT_FALSE(get_bulk().is_valid(face));
    }
    else if(get_parallel_rank() == 1)
    {
      stk::mesh::EntityVector faces = get_faces_for_entity(get_bulk(), get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2));
      expect_not_shared(get_bulk(), faces, __LINE__);
    }
  }
}

TEST_F(TetMesh, DeleteGhostedElement)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    initialize_my_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::EntityVector orphanedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 1)
    };

    stk::mesh::EntityVector sharedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 2),
          get_bulk().get_entity(stk::topology::NODE_RANK, 3),
          get_bulk().get_entity(stk::topology::NODE_RANK, 4)
    };

    get_bulk().modification_begin();
    stk::mesh::Ghosting &ghosting = get_bulk().create_ghosting("ghost");
    stk::mesh::EntityProcVec ghostedElement;
    if(get_parallel_rank() == 0)
    {
      stk::mesh::Entity element1 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
      ghostedElement = {stk::mesh::EntityProc(element1, 1)};
    }
    get_bulk().change_ghosting(ghosting, ghostedElement);
    get_bulk().modification_end();

    stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1)};
    expect_valid(get_bulk(), elementToDestroy, __LINE__);

    if(get_parallel_rank() == 0)
    {
      expect_valid(get_bulk(), orphanedNodes, __LINE__);
    }
    else if(get_parallel_rank() == 1)
    {
      stk::mesh::Entity element2 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2);
      EXPECT_TRUE(get_bulk().is_valid(element2));
    }
    stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

    expect_invalid(get_bulk(), orphanedNodes, __LINE__);
    expect_invalid(get_bulk(), elementToDestroy, __LINE__);

    if(get_parallel_rank() == 0)
    {
      expect_invalid(get_bulk(), sharedNodes, __LINE__);
    }
    else if(get_parallel_rank() == 1)
    {
      expect_not_shared(get_bulk(), sharedNodes, __LINE__);
    }
  }
}

class BeamMesh : public stk::unit_test_util::MeshTestFixture
{
protected:
  BeamMesh()
  {
  }

  void initialize_my_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    set_spatial_dimension(2);
    setup_empty_mesh(auraOption);
    std::string meshDesc =
        "0,1,BEAM_2,1,2\n\
        0,2,BEAM_2,2,3";
        stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    initialize_my_mesh(auraOption);

    stk::mesh::EntityVector orphanedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 3)
    };

    stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

    expect_valid(get_bulk(), orphanedNodes, __LINE__);
    expect_valid(get_bulk(), elementToDestroy, __LINE__);

    stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

    expect_invalid(get_bulk(), orphanedNodes, __LINE__);
    expect_invalid(get_bulk(), elementToDestroy, __LINE__);
  }
};

TEST_F(BeamMesh, DeleteOneElement)
{
  run_test_on_num_procs(1, stk::mesh::BulkData::NO_AUTO_AURA);
}

class QuadMesh : public stk::unit_test_util::MeshTestFixture
{
protected:
  QuadMesh()
  {
  }

  void initialize_my_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    set_spatial_dimension(2);
    setup_empty_mesh(auraOption);
    std::string meshDesc =
        "0,1,QUAD_4_2D,1,2,3,4\n\
        1,2,QUAD_4_2D,2,5,6,3";
        stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    initialize_my_mesh(auraOption);

    stk::mesh::EntityVector orphanedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 5),
          get_bulk().get_entity(stk::topology::NODE_RANK, 6)
    };

    stk::mesh::EntityVector sharedNodes{
      get_bulk().get_entity(stk::topology::NODE_RANK, 2),
          get_bulk().get_entity(stk::topology::NODE_RANK, 3)
    };

    stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

    expect_shared(get_bulk(), sharedNodes, __LINE__);

    if(get_parallel_rank() == 1)
    {
      expect_valid(get_bulk(), orphanedNodes, __LINE__);
      expect_valid(get_bulk(), elementToDestroy, __LINE__);
    }

    stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

    expect_invalid(get_bulk(), orphanedNodes, __LINE__);
    expect_invalid(get_bulk(), elementToDestroy, __LINE__);

    if(get_parallel_rank() == 1)
    {
      expect_invalid(get_bulk(), sharedNodes, __LINE__);
    }
    else if(get_parallel_rank() == 0)
    {
      expect_not_shared(get_bulk(), sharedNodes, __LINE__);
    }
  }
};

TEST_F(QuadMesh, DeleteProcBoundaryElementWithoutAura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(QuadMesh, DeleteProcBoundaryElementWithAura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}

class HexShellHex : public stk::unit_test_util::MeshFixture
{
public:
  void make_hex_shell_hex_mesh()
  {
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

    std::vector<stk::mesh::EntityId> nodeIds = {5, 6, 7, 8};
    stk::mesh::EntityId elementId = 8;
    get_bulk().modification_begin();
    stk::mesh::Part& partWithTopology = get_meta().get_topology_root_part(stk::topology::SHELL_QUADRILATERAL_4);
    stk::mesh::declare_element(get_bulk(), partWithTopology, elementId, nodeIds);
    get_bulk().modification_end();
  }

  void delete_elements_in_ascending_order(unsigned num_elements, const stk::mesh::Entity* elements)
  {
    for(unsigned i=0;i<num_elements;++i)
    {
      if(i==0 || i == 1)
        EXPECT_TRUE(get_bulk().is_valid(elements[i]));
      else
        EXPECT_FALSE(get_bulk().is_valid(elements[i]));

      get_bulk().destroy_entity(elements[i]);
    }
  }

  void delete_elements_in_descending_order(unsigned num_elements, const stk::mesh::Entity* elements)
  {
    ASSERT_TRUE(num_elements==3u);

    for(unsigned i=0;i<num_elements;++i)
    {
      unsigned reverse_index = num_elements - i - 1;
      EXPECT_TRUE(get_bulk().is_valid(elements[reverse_index]));
      get_bulk().destroy_entity(elements[reverse_index]);
    }
  }

  void try_to_delete_incorrectly_all_elements_on_node_5()
  {
    stk::mesh::Entity node5 = get_bulk().get_entity(stk::topology::NODE_RANK, 5);
    get_bulk().modification_begin();
    delete_elements_in_ascending_order(get_bulk().num_elements(node5),  get_bulk().begin_elements(node5));
    get_bulk().modification_end();
  }

  void delete_correctly_all_elements_on_node_5()
  {
    stk::mesh::Entity node5 = get_bulk().get_entity(stk::topology::NODE_RANK, 5);
    get_bulk().modification_begin();
    delete_elements_in_descending_order(get_bulk().num_elements(node5),  get_bulk().begin_elements(node5));
    get_bulk().modification_end();
  }
};

TEST_F(HexShellHex, testDeletionInAscendingOrder_Wrong)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    make_hex_shell_hex_mesh();
    try_to_delete_incorrectly_all_elements_on_node_5();
    unsigned numElements = stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(1u, numElements); // WRONG ANSWER
  }
}

TEST_F(HexShellHex, testDeletionInDescendingOrder_Correct)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    make_hex_shell_hex_mesh();
    delete_correctly_all_elements_on_node_5();
    unsigned numElements = stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(0u, numElements); // Right Answer
  }
}

TEST(DestroyElements, emptyAndDeleteBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_initial_bucket_capacity(2);
  builder.set_maximum_bucket_capacity(2);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

  stk::io::fill_mesh("generated:1x1x4", *bulkPtr);

  stk::mesh::Selector all = bulkPtr->mesh_meta_data().universal_part();
  {
    stk::mesh::BucketVector elemBuckets = bulkPtr->get_buckets(stk::topology::ELEM_RANK, all);
    EXPECT_EQ(2u, elemBuckets.size());
    stk::mesh::BucketVector nodeBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, all);
    EXPECT_EQ(10u, nodeBuckets.size());
  }

  stk::mesh::EntityVector elemsToDestroy = {
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 1),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 3)
  };

  stk::mesh::destroy_elements(*bulkPtr, elemsToDestroy);

  {
    stk::mesh::BucketVector elemBuckets = bulkPtr->get_buckets(stk::topology::ELEM_RANK, all);
    EXPECT_EQ(1u, elemBuckets.size());
    for(const stk::mesh::Bucket* bptr : elemBuckets) {
      for(unsigned i=0; i<bptr->size(); ++i) {
        stk::mesh::Entity ent = (*bptr)[i];
        EXPECT_TRUE(bulkPtr->is_valid(ent))<<"elem bkt-id "<<bptr->bucket_id()<<", ord "<<i;
      }
    }
    stk::mesh::BucketVector nodeBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, all);
    EXPECT_EQ(8u, nodeBuckets.size());
    for(const stk::mesh::Bucket* bptr : nodeBuckets) {
      for(unsigned i=0; i<bptr->size(); ++i) {
        stk::mesh::Entity ent = (*bptr)[i];
        EXPECT_TRUE(bulkPtr->is_valid(ent))<<"node bkt-id "<<bptr->bucket_id()<<", ord "<<i;
      }
    }
  }
}

TEST(DestroyElements, emptyAndDeleteBucketsOnBlockBoundary)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_initial_bucket_capacity(2);
  builder.set_maximum_bucket_capacity(2);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

  stk::io::fill_mesh("generated:1x1x4", *bulkPtr);
  stk::mesh::Part* block1 = bulkPtr->mesh_meta_data().get_part("block_1");
  stk::mesh::Part& block2 = bulkPtr->mesh_meta_data().declare_part("block_2", stk::topology::ELEM_RANK);
  stk::mesh::EntityVector elemsToMove = {
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 1),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 2)
  };

  bulkPtr->modification_begin();
  bulkPtr->change_entity_parts(elemsToMove, stk::mesh::PartVector{&block2}, stk::mesh::PartVector{block1});
  bulkPtr->modification_end();

  stk::mesh::Selector all = bulkPtr->mesh_meta_data().universal_part();
  {
    stk::mesh::BucketVector elemBuckets = bulkPtr->get_buckets(stk::topology::ELEM_RANK, all);
    EXPECT_EQ(2u, elemBuckets.size());
    stk::mesh::BucketVector nodeBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, all);
    EXPECT_EQ(10u, nodeBuckets.size());
  }

  stk::mesh::EntityVector elemsToDestroy = {
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 1),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 2),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 3)
  };

  stk::mesh::destroy_elements(*bulkPtr, elemsToDestroy);

  {
    stk::mesh::BucketVector elemBuckets = bulkPtr->get_buckets(stk::topology::ELEM_RANK, all);
    EXPECT_EQ(1u, elemBuckets.size());
    for(const stk::mesh::Bucket* bptr : elemBuckets) {
      for(unsigned i=0; i<bptr->size(); ++i) {
        stk::mesh::Entity ent = (*bptr)[i];
        EXPECT_TRUE(bulkPtr->is_valid(ent))<<"elem bkt-id "<<bptr->bucket_id()<<", ord "<<i;
      }
    }
    stk::mesh::BucketVector nodeBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, all);
    EXPECT_EQ(4u, nodeBuckets.size());
    for(const stk::mesh::Bucket* bptr : nodeBuckets) {
      for(unsigned i=0; i<bptr->size(); ++i) {
        stk::mesh::Entity ent = (*bptr)[i];
        EXPECT_TRUE(bulkPtr->is_valid(ent))<<"node bkt-id "<<bptr->bucket_id()<<", ord "<<i;
      }
    }
  }
}

TEST(DestroyElements, destroyAll)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_initial_bucket_capacity(2);
  builder.set_maximum_bucket_capacity(2);
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();

  stk::io::fill_mesh("generated:1x1x4", *bulkPtr);
  stk::mesh::Part* block1 = bulkPtr->mesh_meta_data().get_part("block_1");
  stk::mesh::Part& block2 = bulkPtr->mesh_meta_data().declare_part("block_2", stk::topology::ELEM_RANK);
  stk::mesh::EntityVector elemsToMove = {
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 1),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 2)
  };

  bulkPtr->modification_begin();
  bulkPtr->change_entity_parts(elemsToMove, stk::mesh::PartVector{&block2}, stk::mesh::PartVector{block1});
  bulkPtr->modification_end();

  stk::mesh::Selector all = bulkPtr->mesh_meta_data().universal_part();
  {
    stk::mesh::BucketVector elemBuckets = bulkPtr->get_buckets(stk::topology::ELEM_RANK, all);
    EXPECT_EQ(2u, elemBuckets.size());
    stk::mesh::BucketVector nodeBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, all);
    EXPECT_EQ(10u, nodeBuckets.size());
  }

  stk::mesh::EntityVector elemsToDestroy = {
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 1),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 2),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 3),
    bulkPtr->get_entity(stk::topology::ELEM_RANK, 4)
  };

  stk::mesh::destroy_elements(*bulkPtr, elemsToDestroy);

  {
    stk::mesh::BucketVector elemBuckets = bulkPtr->get_buckets(stk::topology::ELEM_RANK, all);
    EXPECT_EQ(0u, elemBuckets.size());
    stk::mesh::BucketVector nodeBuckets = bulkPtr->get_buckets(stk::topology::NODE_RANK, all);
    EXPECT_EQ(0u, nodeBuckets.size());
  }
}

}

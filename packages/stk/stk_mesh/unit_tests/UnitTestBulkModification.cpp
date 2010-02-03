
#include <unit_tests/stk_utest_macros.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <unit_tests/UnitTestRingMeshFixture.hpp>

#include <algorithm>

class UnitTestStkMeshBulkModification {
 public:
  UnitTestStkMeshBulkModification() : m_num_procs(0), m_rank(0)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &m_num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  }

  void test_bulkdata_not_syncronized();
  void test_closure_of_non_locally_used_entities();
  void test_all_local_nodes();
  void test_all_local_edges();
  void test_parallel_consistency();

  int m_num_procs;
  int m_rank;
};

namespace {

STKUNIT_UNIT_TEST( UnitTestBulkDataNotSyrncronized , testUnit )
{
  UnitTestStkMeshBulkModification unit;
  unit.test_bulkdata_not_syncronized(); 
}

STKUNIT_UNIT_TEST( UnitTestClosureOfNonLocallyUsedEntities , testUnit )
{
  UnitTestStkMeshBulkModification unit;
  unit.test_closure_of_non_locally_used_entities();
}

STKUNIT_UNIT_TEST( UnitTestAllLocalNodes , testUnit )
{
  UnitTestStkMeshBulkModification unit;
  unit.test_all_local_nodes();
}

STKUNIT_UNIT_TEST( UnitTestAllLocalEdges , testUnit )
{
  UnitTestStkMeshBulkModification unit;
  unit.test_all_local_edges();
}

STKUNIT_UNIT_TEST( UnitTestParallelConsistency , testUnit )
{
  UnitTestStkMeshBulkModification unit;
  unit.test_parallel_consistency();
}

} //end namespace

void UnitTestStkMeshBulkModification::test_bulkdata_not_syncronized()
{
  RingMeshFixture ring_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = ring_mesh.bulk_data();

  bulk_data.modification_begin(); // Intentially make things unsynced

  std::vector< stk::mesh::Entity *> entities;
  std::vector< stk::mesh::Entity *> entities_closure;
  STKUNIT_ASSERT_THROW(stk::mesh::find_closure(bulk_data, entities, entities_closure), std::runtime_error);
}

void UnitTestStkMeshBulkModification::test_closure_of_non_locally_used_entities()
{
  RingMeshFixture ring_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = ring_mesh.bulk_data();

  const std::vector<stk::mesh::Ghosting *> & ghostings = bulk_data.ghostings();

  if (!ghostings.empty()) {
    const stk::mesh::Ghosting& ghost = *ghostings.front();
    const std::vector< stk::mesh::Entity* > & ghost_receive = ghost.receive();

    if (!ghost_receive.empty()) {
      std::vector< stk::mesh::Entity *> entities;
      std::vector< stk::mesh::Entity *> entities_closure;

      entities.push_back(ghost_receive.front());

      STKUNIT_ASSERT_THROW(stk::mesh::find_closure(bulk_data, entities, entities_closure), std::runtime_error);
    }
  }
}

void UnitTestStkMeshBulkModification::test_all_local_nodes()
{
  RingMeshFixture ring_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = ring_mesh.bulk_data();

  {
    std::vector< stk::mesh::Entity *> entities;
    std::vector< stk::mesh::Entity *> entities_closure;
    find_closure(bulk_data, entities, entities_closure);

    // the closure of the an empty set of entities on all procs should be empty
    STKUNIT_EXPECT_TRUE(entities_closure.empty());
  }

  {
    // Get a selector for the univeral part (contains local, shared, and ghosted)
    const stk::mesh::Part& universal = ring_mesh.meta_data().universal_part();
    stk::mesh::Selector universal_selector(universal);

    // Get the buckets that will give us the universal nodes
    const std::vector<stk::mesh::Bucket*>& node_buckets = bulk_data.buckets(stk::mesh::Node);
    std::vector<stk::mesh::Bucket*> buckets;
    stk::mesh::get_buckets(universal_selector, node_buckets, buckets);

    // Get the universal nodes
    std::vector< stk::mesh::Entity *> universal_entities;
    for (std::vector<stk::mesh::Bucket*>::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      stk::mesh::Bucket& b = **itr;
      for (stk::mesh::BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(&(*bitr));
      }
    }
    buckets.clear();

    // sort and unique the universal nodes
    std::sort(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityLess());
    std::vector<stk::mesh::Entity*>::iterator new_end = std::unique(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityEqual());
    universal_entities.erase(new_end, universal_entities.end());

    // Get the buckets that will give us the locally used nodes
    const stk::mesh::Part& local = ring_mesh.meta_data().locally_used_part();
    stk::mesh::Selector locally_used_selector(local);
    stk::mesh::get_buckets(locally_used_selector, node_buckets, buckets);

    // Get the locally used nodes
    std::vector< stk::mesh::Entity *> entities;
    for (std::vector<stk::mesh::Bucket*>::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      stk::mesh::Bucket& b = **itr;
      for (stk::mesh::BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(&(*bitr));
      }
    }

    // Get the closure, passing in the locally used nodes on each proc
    std::vector< stk::mesh::Entity *> entities_closure;
    stk::mesh::find_closure(bulk_data, entities, entities_closure);

    // The ghosted nodes on this part will be locally used on one of the other
    // procs, so we expect that they will be part of the closure. In other
    // words, the set of nodes returned by find_closure should exactly match
    // the set of universal nodes.
    STKUNIT_ASSERT_TRUE(universal_entities.size() == entities_closure.size());
    stk::mesh::EntityEqual ee;
    for (size_t i = 0; i < entities_closure.size(); ++i) {
      STKUNIT_EXPECT_TRUE(ee(universal_entities[i], entities_closure[i]));
    }
  }
}

void UnitTestStkMeshBulkModification::test_all_local_edges()
{
  RingMeshFixture ring_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = ring_mesh.bulk_data();

  {
    const stk::mesh::Part& universal = ring_mesh.meta_data().universal_part();
    stk::mesh::Selector universal_selector(universal);

    const std::vector<stk::mesh::Bucket*>& node_buckets = bulk_data.buckets(stk::mesh::Node);
    const std::vector<stk::mesh::Bucket*>& edge_buckets = bulk_data.buckets(stk::mesh::Edge);
    std::vector<stk::mesh::Bucket*> buckets;

    stk::mesh::get_buckets(universal_selector, node_buckets, buckets);

    // get all the nodes that this process knows about
    std::vector< stk::mesh::Entity *> universal_entities;
    for (std::vector<stk::mesh::Bucket*>::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      stk::mesh::Bucket& b = **itr;
      for (stk::mesh::BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(&(*bitr));
      }
    }
    buckets.clear();

    stk::mesh::get_buckets(universal_selector, edge_buckets, buckets);

    // get all the edges that this process knows about
    for (std::vector<stk::mesh::Bucket*>::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      stk::mesh::Bucket& b = **itr;
      for (stk::mesh::BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(&(*bitr));
      }
    }
    buckets.clear();

    // universal entities should now have all the universal nodes and edges
    // sort and uniq the universal nodes/edges
    std::sort(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityLess());
    std::vector<stk::mesh::Entity*>::iterator new_end = std::unique(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityEqual());
    universal_entities.erase(new_end, universal_entities.end());

    // get the buckets that we need to traverse to get the locally used edges
    const stk::mesh::Part& local = ring_mesh.meta_data().locally_used_part();
    stk::mesh::Selector locally_used_selector(local);
    stk::mesh::get_buckets(locally_used_selector, edge_buckets, buckets);

    // get the locally used edges and store them in entities
    std::vector< stk::mesh::Entity *> entities;
    for (std::vector<stk::mesh::Bucket*>::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      stk::mesh::Bucket& b = **itr;
      for (stk::mesh::BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(&(*bitr));
      }
    }

    // call find_closure, passing in the locally used edges
    std::vector< stk::mesh::Entity *> entities_closure;
    stk::mesh::find_closure(bulk_data, entities, entities_closure);

    // The ghosted entities on this proc (edge or node) should be contained
    // in the closure of the locally-used edge on some other proc, so we
    // expect that they will be part of the closure. In other
    // words, the set of entities returned by find_closure should exactly match
    // the set of universal entities (nodes and edges).
    STKUNIT_ASSERT_TRUE(universal_entities.size() == entities_closure.size());
    stk::mesh::EntityEqual ee;
    for (size_t i = 0; i < entities_closure.size(); ++i) {
      STKUNIT_EXPECT_TRUE(ee(universal_entities[i], entities_closure[i]));
    }
  }
}

void UnitTestStkMeshBulkModification::test_parallel_consistency()
{
  RingMeshFixture ring_mesh(MPI_COMM_WORLD);

  stk::mesh::BulkData& bulk_data = ring_mesh.bulk_data();

  stk::CommBroadcast all(bulk_data.parallel(), 0);

  std::vector< stk::mesh::Entity *> entities;
  std::vector< stk::mesh::Entity *> entities_closure;

  // For proc 0 only, add locally used nodes to entities, for all other
  // procs, leave entities empty.
  if (m_rank == 0) {
    const std::vector<stk::mesh::Bucket*>& node_buckets = bulk_data.buckets(stk::mesh::Node);

    const stk::mesh::Part& local = ring_mesh.meta_data().locally_used_part();
    stk::mesh::Selector locally_used_selector(local);

    std::vector<stk::mesh::Bucket*> buckets;
    stk::mesh::get_buckets(locally_used_selector, node_buckets, buckets);

    for (std::vector<stk::mesh::Bucket*>::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      stk::mesh::Bucket& b = **itr;
      for (stk::mesh::BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(&(*bitr));
      }
    }
  }

  // Call find_closure with proc 0 passing in locally-used nodes
  stk::mesh::find_closure(bulk_data, entities, entities_closure);

  // Proc 0 will broadcast the global ids of the nodes it passed to
  // find_closure

  // pack entities for sizing
  for (std::vector<stk::mesh::Entity*>::const_iterator
          ep = entities.begin() ; ep != entities.end() ; ++ep ) {
    all.send_buffer().pack<stk::mesh::EntityKey>((*ep)->key());
  }

  all.allocate_buffer();

  // pack for real
  for (std::vector<stk::mesh::Entity*>::const_iterator
         ep = entities.begin() ; ep != entities.end() ; ++ep ) {
    all.send_buffer().pack<stk::mesh::EntityKey>((*ep)->key());
  }

  all.communicate();

  // clear-out entities and put the nodes that correspond to the keys
  // broadcast by proc 0 into entities.
  entities.clear();
  stk::CommBuffer& buf = all.recv_buffer();
  stk::mesh::EntityKey k ;
  while ( buf.remaining() ) {
    buf.unpack<stk::mesh::EntityKey>(k);
    stk::mesh::Entity * e = bulk_data.get_entity(k);
    // If a proc is not aware of a key, that means it has no relationship 
    // with that entity, so it can ignore it.
    if (e != NULL) {
      entities.push_back(e);
    }
  }

  // sort and unique entities
  std::sort(entities.begin(), entities.end(), stk::mesh::EntityLess());
  std::vector<stk::mesh::Entity*>::iterator new_end = std::unique(entities.begin(), entities.end(), stk::mesh::EntityEqual());
  entities.erase(new_end, entities.end());

  // If any processor had ghosted nodes that were local to proc 0, those
  // nodes should be in the closure because proc 0 passed them in to 
  // find_closure.
  STKUNIT_ASSERT_TRUE(entities.size() == entities_closure.size());
  stk::mesh::EntityEqual ee;
  for (size_t i = 0; i < entities_closure.size(); ++i) {
    STKUNIT_EXPECT_TRUE(ee(entities[i], entities_closure[i]));
  }
}

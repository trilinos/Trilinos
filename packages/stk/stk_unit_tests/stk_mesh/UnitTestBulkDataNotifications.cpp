#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/ModificationObserver.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_max
#include <string>                       // for string
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <vector>                       // for vector
#include "mpi.h"                        // for ompi_communicator_t, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Types.hpp"      // for EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }

namespace {

void reduce_max(std::vector<size_t> &localVector, stk::ParallelMachine comm)
{
  std::vector<size_t> globalVector(localVector.size());
  stk::all_reduce_max(comm, localVector.data(), globalVector.data(), globalVector.size());
  localVector.swap(globalVector);
}

class TestListener : public stk::mesh::ModificationObserver {
public:
  TestListener(stk::ParallelMachine communicator)
    : stk::mesh::ModificationObserver(stk::mesh::ModificationObserverPriority::APPLICATION),
      comm(communicator),
      buckets_changed(5, 0)
  {}

  virtual void local_buckets_changed_notification(stk::mesh::EntityRank rank)
  {
    buckets_changed[rank]++;
  }

  virtual void finished_modification_end_notification()
  {
    reduce_max(buckets_changed, comm);
  }

  size_t get_buckets_changed(stk::mesh::EntityRank rank) const
  { return buckets_changed[rank]; }

private:
  stk::ParallelMachine comm;
  std::vector<size_t> buckets_changed;
};

TEST(BulkDataMod, synchronized_count_basic_fill_mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) > 2) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();

  const std::string generatedMeshSpec = "generated:1x1x2";
  stk::io::fill_mesh(generatedMeshSpec, *bulk);

  EXPECT_EQ(1u, bulk->synchronized_count());
}

TEST(BulkDataMod, synchronized_count_empty_mod_cycle)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(comm) > 2) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();

  const std::string generatedMeshSpec = "generated:1x1x2";
  stk::io::fill_mesh(generatedMeshSpec, *bulk);

  unsigned modCount = bulk->synchronized_count();

  bulk->modification_begin();
  bulk->modification_end();

  EXPECT_EQ(modCount, bulk->synchronized_count());
}

TEST(BulkDataNotifications, test_listener_buckets_changed)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(comm);
  if (numProcs==2)
  {
    stk::mesh::MetaData meta(3);
    stk::unit_test_util::BulkDataTester mesh(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);

    const std::string generatedMeshSpec = "generated:1x1x2";
    stk::io::fill_mesh(generatedMeshSpec, mesh);

    std::shared_ptr<TestListener> listener = std::make_shared<TestListener>(comm);
    mesh.register_observer(listener);

    int procId = stk::parallel_machine_rank(comm);

    stk::mesh::Part& newPart = meta.declare_part("new part", stk::topology::NODE_RANK);
    stk::mesh::Entity node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);

    if (procId == 0) {
      EXPECT_TRUE(mesh.is_valid(node1));
    }
    else {
      EXPECT_FALSE(mesh.is_valid(node1));
    }

    mesh.modification_begin();
    if (procId == 0) {
      mesh.change_entity_parts(node1, stk::mesh::ConstPartVector{&newPart});
    }
    mesh.modification_end();

    EXPECT_EQ(1u, listener->get_buckets_changed(stk::topology::NODE_RANK));
  }
}

class TestOrderingObserver : public stk::mesh::ModificationObserver {
public:
  TestOrderingObserver(int id, std::vector<int> & evaluationOrder)
    : stk::mesh::ModificationObserver(stk::mesh::ModificationObserverPriority::APPLICATION),
      m_id(id),
      m_evaluationOrder(evaluationOrder)
  {}

  ~TestOrderingObserver() = default;

  virtual void finished_modification_end_notification()
  {
    m_evaluationOrder.push_back(m_id);
  }

private:
  int m_id;
  std::vector<int> & m_evaluationOrder;
};

TEST(BulkDataNotifications, test_ordering_different_priorities)
{
  stk::mesh::ModificationNotifier notifier;

  std::vector<int> evaluationOrder;
  std::shared_ptr<TestOrderingObserver> observer1 = std::make_shared<TestOrderingObserver>(1, evaluationOrder);
  std::shared_ptr<TestOrderingObserver> observer2 = std::make_shared<TestOrderingObserver>(2, evaluationOrder);
  std::shared_ptr<TestOrderingObserver> observer3 = std::make_shared<TestOrderingObserver>(3, evaluationOrder);

  observer3->set_priority(stk::mesh::ModificationObserverPriority::APPLICATION);
  observer1->set_priority(stk::mesh::ModificationObserverPriority::STK_INTERNAL_LOW_PRIORITY);
  observer2->set_priority(stk::mesh::ModificationObserverPriority::STK_TRANSITION);

  notifier.register_observer(observer3);
  notifier.register_observer(observer1);
  notifier.register_observer(observer2);

  notifier.notify_finished_modification_end(MPI_COMM_WORLD);

  EXPECT_EQ(evaluationOrder, (std::vector<int>{1, 2, 3}));
}

TEST(BulkDataNotifications, test_ordering_duplicate_priorities)
{
  stk::mesh::ModificationNotifier notifier;

  std::vector<int> evaluationOrder;
  std::shared_ptr<TestOrderingObserver> observer1 = std::make_shared<TestOrderingObserver>(1, evaluationOrder);
  std::shared_ptr<TestOrderingObserver> observer2 = std::make_shared<TestOrderingObserver>(2, evaluationOrder);
  std::shared_ptr<TestOrderingObserver> observer3 = std::make_shared<TestOrderingObserver>(3, evaluationOrder);

  observer3->set_priority(stk::mesh::ModificationObserverPriority::STK_INTERNAL_LOW_PRIORITY);
  observer1->set_priority(stk::mesh::ModificationObserverPriority::STK_INTERNAL_LOW_PRIORITY);
  observer2->set_priority(stk::mesh::ModificationObserverPriority::STK_INTERNAL_LOW_PRIORITY);

  notifier.register_observer(observer3);
  notifier.register_observer(observer1);
  notifier.register_observer(observer2);

  notifier.notify_finished_modification_end(MPI_COMM_WORLD);

  EXPECT_EQ(evaluationOrder, (std::vector<int>{3, 1, 2}));
}

}

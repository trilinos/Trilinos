
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/baseImpl/elementGraph/ProcessKilledElements.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <vector>                       // for vector
#include "BulkDataElementGraphTester.hpp"
#include "ElementGraphTester.hpp"       // for ElemElemGraphTester
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for EntityVector, PartVector
#include "stk_unit_test_utils/ioUtils.hpp"  // for fill_mesh_using_stk_io
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_mesh/base/ModificationObserver.hpp>

namespace stk { namespace mesh { class Part; } }

namespace
{

unsigned count_elements(const stk::mesh::BulkData &bulk, stk::mesh::Part &activePart)
{
  std::vector<size_t> countsPerEntityType;
  stk::mesh::count_entities(activePart, bulk, countsPerEntityType);
  return countsPerEntityType[stk::topology::ELEMENT_RANK];
}

unsigned get_gold_num_elements_by_proc(int procRank, int procRankToKillElementsOn)
{
  unsigned goldNumElements = 2u;
  if (procRank == procRankToKillElementsOn)
  {
    goldNumElements = 1u;
  }
  return goldNumElements;
}

// Test cases: (all cases must involve a "element death" starter)
//  * create neighbor then deactivate => "refinement"
//  * move processor boundary (change_entity_owner) => rebalance
//  ? framework change global id ?

class ElementDeathMock
{
public:
  ElementDeathMock(BulkDataElementGraphTester &bulk_, stk::mesh::Part &activePart_, stk::mesh::Part &boundaryPart_) :
    bulk(bulk_),
    activePart(activePart_),
    boundaryPart(boundaryPart_)
  { }

  void kill_elements(stk::mesh::EntityVector &elementsToKill)
  {
    deactivate_elements(elementsToKill);

    stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(bulk, bulk.get_face_adjacent_element_graph(), activePart, remoteActiveSelector);

    stk::mesh::PartVector newlyCreatedBoundaryFacesGetAddedToTheseParts = {&boundaryPart};
    stk::mesh::PartVector exposedButExistingBoundaryFacesGetAddedToTheseParts = {&boundaryPart};
    process_killed_elements(bulk,
                            elementsToKill,
                            activePart,
                            remoteActiveSelector,
                            newlyCreatedBoundaryFacesGetAddedToTheseParts,
                            &exposedButExistingBoundaryFacesGetAddedToTheseParts);
  }
private:
  void deactivate_elements(stk::mesh::EntityVector elementsToDeactivate)
  {
    bulk.modification_begin();
    for(stk::mesh::Entity element : elementsToDeactivate)
    {
      deactivate_element(element);
    }
    bulk.modification_end();
  }

  void deactivate_element(stk::mesh::Entity elementToDeactivate)
  {
    bulk.change_entity_parts(elementToDeactivate, stk::mesh::ConstPartVector{}, stk::mesh::ConstPartVector{&activePart});
  }
private:
  BulkDataElementGraphTester &bulk;
  stk::mesh::Part &activePart;
  stk::mesh::Part &boundaryPart;
};

class MeshRefinementMock
{
public:
  MeshRefinementMock(stk::mesh::BulkData &bulk_, stk::mesh::Part &activePart_)
    : bulk(bulk_), activePart(activePart_)
  {
  }

  void create_element()
  {
    stk::mesh::EntityId newElementId = 5;
    stk::mesh::EntityIdVector nodesForNewElement = {17, 18, 20, 19, 21, 22, 24, 23};
    stk::mesh::Part *block1Part = bulk.mesh_meta_data().get_part("block_1");
    stk::mesh::PartVector elementParts = {block1Part, &activePart};
    stk::mesh::declare_element(bulk, elementParts, newElementId, nodesForNewElement);
  }

  void create_element_on_proc1()
  {
    bulk.modification_begin();
    if(bulk.parallel_rank() == 1)
    {
      create_element();
    }
    bulk.modification_end();
  }

private:
  stk::mesh::BulkData &bulk;
  stk::mesh::Part &activePart;
};

class UpdateElemElemGraphTest : public ::testing::Test
{
public:
  UpdateElemElemGraphTest() :
    meta(3),
    bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA),
    activePart(meta.declare_part("active")),
    boundaryPart(meta.declare_part("boundary"))
  {
    stk::io::fill_mesh("generated:1x1x4", bulk);
    stk::unit_test_util::put_mesh_into_part(bulk, activePart);
  }

  void add_element_to_deactivate_on_proc(int elementId, int procRankToKillElementsOn, stk::mesh::EntityVector &elementsDeactivated)
  {
    int myProcRank = stk::parallel_machine_rank(bulk.parallel());
    if(myProcRank == procRankToKillElementsOn)
    {
      stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elementId);
      elementsDeactivated.push_back(element);
    }
  }

  void kill_element_1_on_proc_0(ElementDeathMock &elementKiller)
  {
    stk::mesh::EntityVector elementsDeactivated;
    int element1Id = 1;
    int procRankToKillElementsOn = 0;
    add_element_to_deactivate_on_proc(element1Id, procRankToKillElementsOn, elementsDeactivated);

    elementKiller.kill_elements(elementsDeactivated);

    unsigned goldNumElements = get_gold_num_elements_by_proc(bulk.parallel_rank(), procRankToKillElementsOn);
    EXPECT_EQ(goldNumElements, count_elements(bulk, activePart));
  }

  void kill_element_4_on_proc_1(ElementDeathMock &elementKiller)
  {
    stk::mesh::EntityVector elementsDeactivated;
    int element4Id = 4;
    int procRankToKillElementsOn = 1;
    add_element_to_deactivate_on_proc(element4Id, procRankToKillElementsOn, elementsDeactivated);

    elementKiller.kill_elements(elementsDeactivated);
  }

  void move_element_3_to_proc_0()
  {
    stk::mesh::EntityProcVec entityProcs;
    if(bulk.parallel_rank() == 1)
    {
      stk::mesh::Entity element3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3);
      entityProcs.push_back(stk::mesh::EntityProc(element3, 0));
    }
    bulk.change_entity_owner(entityProcs);
  }

  void delete_element_3_on_proc_1()
  {
    bulk.modification_begin();
    if(bulk.parallel_rank() == 1)
    {
      stk::mesh::Entity element3 = bulk.get_entity(stk::topology::ELEMENT_RANK, 3);
      bulk.destroy_entity(element3);
    }
    bulk.modification_end();
  }

protected:
  stk::mesh::MetaData meta;
  BulkDataElementGraphTester bulk;
  stk::mesh::Part &activePart;
  stk::mesh::Part &boundaryPart;
};

TEST_F(UpdateElemElemGraphTest, killRefineUpdateGraphKill)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) == 2)
  {
    ElementDeathMock elementKiller(bulk, activePart, boundaryPart);
    kill_element_1_on_proc_0(elementKiller);

    MeshRefinementMock meshRefinement(bulk, activePart);
    meshRefinement.create_element_on_proc1();

    kill_element_4_on_proc_1(elementKiller);

    unsigned goldNumFaces = 1;
    unsigned goldNumElems = 1;
    if(bulk.parallel_rank() == 1)
    {
      goldNumFaces = 2;
      goldNumElems = 2;
    }
    EXPECT_EQ(goldNumFaces, stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(goldNumElems, stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::ELEM_RANK)));
  }
}

TEST_F(UpdateElemElemGraphTest, killChangeOwnerUpdateGraphKill)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) == 2)
  {
    ElementDeathMock elementKiller(bulk, activePart, boundaryPart);
    kill_element_1_on_proc_0(elementKiller);

    move_element_3_to_proc_0();

    kill_element_4_on_proc_1(elementKiller);

    unsigned goldNumFaces = 2;
    unsigned goldNumElems = 2;
    if(bulk.parallel_rank() == 1)
    {
      goldNumFaces = 1;
      goldNumElems = 0;
    }
    EXPECT_EQ(goldNumFaces, stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(goldNumElems, stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::ELEM_RANK)));
  }
}

TEST_F(UpdateElemElemGraphTest, killDeleteUpdateGraphKill)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) == 2)
  {
    ElementDeathMock elementKiller(bulk, activePart, boundaryPart);
    kill_element_1_on_proc_0(elementKiller);

    delete_element_3_on_proc_1();

    kill_element_4_on_proc_1(elementKiller);

    unsigned goldNumFaces = 1;
    unsigned goldNumElems = 1;
    if(bulk.parallel_rank() == 1)
    {
      goldNumFaces = 0;
      goldNumElems = 0;
    }
    EXPECT_EQ(goldNumFaces, stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::FACE_RANK)));
    EXPECT_EQ(goldNumElems, stk::mesh::count_selected_entities(activePart, bulk.buckets(stk::topology::ELEM_RANK)));
  }
}


class ElemElemGraphUpdaterMock : public stk::mesh::ModificationObserver
{
public:
  ElemElemGraphUpdaterMock()
  : stk::mesh::ModificationObserver(stk::mesh::ModificationObserverPriority::APPLICATION),
    numAdded(0)
  {

  }

  virtual void entity_added(stk::mesh::Entity entity)
  {
    numAdded++;
  }

  int get_num_entities_added()
  {
    return numAdded;
  }
private:
  int numAdded;
};

TEST_F(UpdateElemElemGraphTest, NewEntityNotification)
{
  if(bulk.parallel_size() == 2)
  {
    std::shared_ptr<ElemElemGraphUpdaterMock> observer = std::make_shared<ElemElemGraphUpdaterMock>();
    bulk.register_observer(observer);

    MeshRefinementMock meshRefinement(bulk, activePart);
    meshRefinement.create_element_on_proc1();

    if(bulk.parallel_rank() == 1)
    {
      EXPECT_EQ(5, observer->get_num_entities_added());
    }
  }
}

TEST_F(UpdateElemElemGraphTest, MultipleObservers)
{
  if(bulk.parallel_size() == 2)
  {
    std::shared_ptr<ElemElemGraphUpdaterMock> observer1 = std::make_shared<ElemElemGraphUpdaterMock>();
    bulk.register_observer(observer1);
    std::shared_ptr<ElemElemGraphUpdaterMock> observer2 = std::make_shared<ElemElemGraphUpdaterMock>();
    bulk.register_observer(observer2);

    MeshRefinementMock meshRefinement(bulk, activePart);
    meshRefinement.create_element_on_proc1();

    if(bulk.parallel_rank() == 1)
    {
      EXPECT_EQ(5, observer1->get_num_entities_added());
      EXPECT_EQ(5, observer2->get_num_entities_added());
    }
  }
}

TEST(Observer, NewEntityNotification)
{
  std::shared_ptr<ElemElemGraphUpdaterMock> observer = std::make_shared<ElemElemGraphUpdaterMock>();
  stk::mesh::Entity goldNewEntity(2);
  observer->entity_added(goldNewEntity);

  EXPECT_EQ(1, observer->get_num_entities_added());
}

} // namespace

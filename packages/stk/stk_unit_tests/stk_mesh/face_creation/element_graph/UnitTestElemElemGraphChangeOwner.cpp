
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <ostream>                      // for basic_ostream::operator<<, etc
#include <stdexcept>                    // for logic_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include <stk_util/environment/WallTime.hpp>  // for wall_time
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <string>                       // for string
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include "ElementGraphTester.hpp"       // for ElemElemGraphTester
#include "mpi.h"                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>  // for change_entity_owner, etc
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>  // for parallel_info
#include <stk_mesh/baseImpl/elementGraph/ProcessKilledElements.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>

namespace stk { namespace mesh { class Part; } }

namespace {

void expect_elem_connected_to_local_elem_id_via_side(const stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph& elemGraph,
                                                     stk::mesh::Entity elem, stk::mesh::EntityId connectedId, int viaSide)
{
  stk::mesh::EntityId elemId = bulkData.identifier(elem);
  size_t numConnected = elemGraph.get_num_connected_elems(elem);
  bool foundLocallyConnectedId = false;
  for(size_t i=0; i<numConnected; ++i) {
    if (elemGraph.is_connected_elem_locally_owned(elem, i)) {
      stk::mesh::impl::ElementViaSidePair elemViaSidePair = elemGraph.get_connected_element_and_via_side(elem, i);
      if (bulkData.identifier(elemViaSidePair.element) == connectedId && elemViaSidePair.side == viaSide) {
        foundLocallyConnectedId = true;
      }
    }
  }
  ASSERT_TRUE(foundLocallyConnectedId) << "elem " << elemId << " expected local elem " << connectedId;
}

void expect_elem_connected_to_remote_elem_id_via_side(const stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph& elemGraph,
                                                      stk::mesh::Entity elem, stk::mesh::EntityId connectedId, int viaSide)
{
  stk::mesh::EntityId elemId = bulkData.identifier(elem);
  size_t numConnected = elemGraph.get_num_connected_elems(elem);
  bool foundRemotelyConnectedId = false;
  for(size_t i=0; i<numConnected; ++i) {
    if (!elemGraph.is_connected_elem_locally_owned(elem, i)) {
      stk::mesh::impl::IdViaSidePair idViaSidePair = elemGraph.get_connected_remote_id_and_via_side(elem, i);
      if (idViaSidePair.id == connectedId && idViaSidePair.side == viaSide) {
        foundRemotelyConnectedId = true;
      }
    }
  }
  ASSERT_TRUE(foundRemotelyConnectedId) << "elem " << elemId << " expected remote elem " << connectedId;
}

class ElemGraphChangeOwner : public stk::unit_test_util::MeshTestFixture
{
protected:
  void expect_initial_graph_correct()
  {
    if(get_bulk().parallel_rank() == 0)
      check_element2_connected_to_element1_locally_and_element3_remotely();
    else if(get_bulk().parallel_rank() == 1)
      check_element3_conencted_to_element4_locally_and_element2_remotely();
  }

  void check_element2_connected_to_element1_locally_and_element3_remotely()
  {
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem2));
    expect_connected_to_local_elem_id_via_side(elem2, 1, 4);
    expect_element2_connected_to_3_remotely_via_side_5();
  }

  void check_element3_conencted_to_element4_locally_and_element2_remotely()
  {
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
    ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem3));
    expect_connected_to_local_elem_id_via_side(elem3, 4, 5);
    expect_element3_connected_to_2_remotely_via_side_4();
  }

  void expect_element2_connected_to_3_remotely_via_side_5()
  {
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    expect_connected_to_remote_elem_id_via_side(elem2, 3, 5);
    expect_parallel_info_from_elem2_to_3();
  }

  void expect_element3_connected_to_2_remotely_via_side_4()
  {
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
    expect_connected_to_remote_elem_id_via_side(elem3, 2, 4);
    expect_parallel_info_from_elem3_to_2();
  }

  void expect_parallel_info_from_elem2_to_3()
  {
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    const stk::mesh::impl::ParallelInfo &parInfo = get_elem_graph().get_parallel_edge_info(elem2, 5, 3, 4);
    expect_otherProc_permutation_chosenId(parInfo, 1, 4);
  }

  void expect_parallel_info_from_elem3_to_2()
  {
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::impl::ParallelInfo &parInfo = get_elem_graph().get_parallel_edge_info(elem3, 4, 2, 5);
    expect_otherProc_permutation_chosenId(parInfo, 0, 4);
  }

  void expect_connected_to_local_elem_id_via_side(stk::mesh::Entity elem, stk::mesh::EntityId connectedId, int viaSide)
  {
    expect_elem_connected_to_local_elem_id_via_side(get_bulk(), get_elem_graph(), elem, connectedId, viaSide);
  }

  void expect_connected_to_remote_elem_id_via_side(stk::mesh::Entity elem, stk::mesh::EntityId connectedId, int viaSide)
  {
    expect_elem_connected_to_remote_elem_id_via_side(get_bulk(), get_elem_graph(), elem, connectedId, viaSide);
  }

  void expect_connected_to_remote_elem_id(stk::mesh::Entity elem,
                                          size_t connectedIndex,
                                          stk::mesh::EntityId connectedId)
  {
    stk::mesh::EntityId elemId = get_bulk().identifier(elem);
    ASSERT_TRUE(!get_elem_graph().is_connected_elem_locally_owned(elem, connectedIndex))
        << "elem " << elemId << " expected remote elem " << connectedId;
    EXPECT_EQ(connectedId, get_elem_graph().get_connected_remote_id_and_via_side(elem, connectedIndex).id) << "elem " << elemId;
  }

  void expect_otherProc_permutation_chosenId(const stk::mesh::impl::ParallelInfo &parInfo,
                                             int otherProc,
                                             int perm)
  {
    EXPECT_EQ(otherProc, parInfo.get_proc_rank_of_neighbor());
    EXPECT_EQ(perm, parInfo.m_permutation);
  }

  void move_elements(const stk::mesh::EntityIdProcVec &elementIdProcsToMove)
  {
    stk::mesh::EntityProcVec elemProcPairsToMove;
    for(const stk::mesh::EntityIdProc &entityIdProc : elementIdProcsToMove)
      append_element_if_owned(entityIdProc.first, entityIdProc.second, elemProcPairsToMove);
    get_bulk().change_entity_owner(elemProcPairsToMove);
  }

  void append_element_if_owned(stk::mesh::EntityId elementId, int destProc, stk::mesh::EntityProcVec &elemProcPairsToMove)
  {
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elementId);
    if(is_owned_on_this_proc(element))
      elemProcPairsToMove.push_back(stk::mesh::EntityProc(element, destProc));
  }

  bool is_owned_on_this_proc(stk::mesh::Entity element)
  {
    return (get_bulk().is_valid(element) && get_bulk().bucket(element).owned());
  }

  void create_elem_graph()
  {
    elemElemGraph = new ElemElemGraphTester(get_bulk());
    elemElemGraphUpdater = std::make_shared<stk::mesh::ElemElemGraphUpdater>(get_bulk(), *elemElemGraph);
    get_bulk().register_observer(elemElemGraphUpdater);
  }

  ElemElemGraphTester &get_elem_graph()
  {
    return *elemElemGraph;
  }

  ElemGraphChangeOwner()
    : elemElemGraph(nullptr), elemElemGraphUpdater()
  {
  }

  ~ElemGraphChangeOwner()
  {
    delete elemElemGraph;
  }

protected:
  ElemElemGraphTester *elemElemGraph;
  std::shared_ptr<stk::mesh::ElemElemGraphUpdater> elemElemGraphUpdater;
};


class ElemGraphChangeOwnerMoveFrom1To0 : public ElemGraphChangeOwner
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x1x4", auraOption);
    test_graph_updated_with_elem_3_moving_to_proc_0();
  }

  void test_graph_updated_with_elem_3_moving_to_proc_0()
  {
    create_elem_graph();
    expect_initial_graph_correct();
    move_elements({stk::mesh::EntityIdProc(3, 0)});
    expect_graph_updated_after_elem_3_moved_to_0();
  }

  void expect_graph_updated_after_elem_3_moved_to_0()
  {
    if(get_bulk().parallel_rank() == 0)
    {
      check_element2_connected_to_element1_and_element3_locally();
      check_element3_connected_to_element2_locally_and_element4_remotely();
    }
    else
    {
      check_element4_connected_to_element3_remotely();
    }
  }

  void check_element2_connected_to_element1_and_element3_locally()
  {
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem2));
    expect_connected_to_local_elem_id_via_side(elem2, 1, 4);
    expect_connected_to_local_elem_id_via_side(elem2, 3, 5);
  }

  void check_element4_connected_to_element3_remotely()
  {
    expect_element4_connected_to_3_remotely_via_side_3();
    expect_parallel_info_from_elem4_to_3();
  }

  void check_element3_connected_to_element2_locally_and_element4_remotely()
  {
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
    ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem3));
    expect_element3_connected_to_4_remotely_via_side_5();
    expect_connected_to_local_elem_id_via_side(elem3, 2, 4);
  }

  void expect_element3_connected_to_4_remotely_via_side_5()
  {
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
    expect_connected_to_remote_elem_id_via_side(elem3, 4, 5);
    expect_parallel_info_from_elem3_to_4();
  }

  void expect_parallel_info_from_elem3_to_4()
  {
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::impl::ParallelInfo &parInfo = get_elem_graph().get_parallel_edge_info(elem3, 5, 4, 4);
    expect_otherProc_permutation_chosenId(parInfo, 1, 4);
  }

  void expect_element4_connected_to_3_remotely_via_side_3()
  {
    stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4);
    ASSERT_EQ(1u, get_elem_graph().get_num_connected_elems(elem4));
    expect_connected_to_remote_elem_id_via_side(elem4, 3, 4);
  }

  void expect_parallel_info_from_elem4_to_3()
  {
    stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4);
    const stk::mesh::impl::ParallelInfo &parInfo = get_elem_graph().get_parallel_edge_info(elem4, 4, 3, 5);
    expect_otherProc_permutation_chosenId(parInfo, 0, 4);
  }
};
TEST_F(ElemGraphChangeOwnerMoveFrom1To0, withAura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveFrom1To0, withoutAura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerMoveEverythingFromProc1 : public ElemGraphChangeOwner
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x1x4", auraOption);
    expect_graph_correct_after_moving_everything_to_proc0();
  }
  void expect_graph_correct_after_moving_everything_to_proc0()
  {
    create_elem_graph();
    expect_initial_graph_correct();
    move_elements({stk::mesh::EntityIdProc(3, 0), stk::mesh::EntityIdProc(4, 0)});

    if (get_bulk().parallel_rank() == 0) {
      stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
      ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(elem3));
      expect_connected_to_local_elem_id_via_side(elem3, 2, 4);
      expect_connected_to_local_elem_id_via_side(elem3, 4, 5);
    }
    else {
      EXPECT_EQ(0u, get_elem_graph().num_edges());
      EXPECT_EQ(0u, get_elem_graph().num_parallel_edges());
    }
  }
};
TEST_F(ElemGraphChangeOwnerMoveEverythingFromProc1, withAura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveEverythingFromProc1, withoutAura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerLeapFrog : public ElemGraphChangeOwner
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x1x6", auraOption);
    expect_graph_correct_after_leaps();
  }
  void expect_graph_correct_after_leaps()
  {
    create_elem_graph();
    expect_initial_graph_correct();
    move_elements({stk::mesh::EntityIdProc(2, 1), stk::mesh::EntityIdProc(3, 2)});

    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      expect_connected_to_remote_elem_id_via_side(elem2, 3, 5);
      const stk::mesh::impl::ParallelInfo& parInfo = get_elem_graph().get_const_parallel_edge_info(elem2, 5, 3, 4);
      expect_otherProc_permutation_chosenId(parInfo, 2, 4);
    }
    else if (get_bulk().parallel_rank() == 2) {
      stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
      expect_connected_to_remote_elem_id_via_side(elem3, 2, 4);
      const stk::mesh::impl::ParallelInfo& parInfo = get_elem_graph().get_const_parallel_edge_info(elem3, 4, 2, 5);
      expect_otherProc_permutation_chosenId(parInfo, 1, 4);
    }
  }
};
TEST_F(ElemGraphChangeOwnerLeapFrog, withAura)
{
  run_test_on_num_procs(3, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerLeapFrog, withoutAura)
{
  run_test_on_num_procs(3, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerMoveNeighborsToEnd : public ElemGraphChangeOwner
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x1x6", auraOption);
    expect_graph_correct_after_moving_neighbors_to_last_proc();
  }
  void expect_graph_correct_after_moving_neighbors_to_last_proc()
  {
    create_elem_graph();
    expect_initial_graph_correct();
    move_elements({stk::mesh::EntityIdProc(2, 2), stk::mesh::EntityIdProc(3, 2)});

    if (get_bulk().parallel_rank() == 2) {
      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
      stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
      expect_connected_to_local_elem_id_via_side(elem2, 3, 5);
      expect_connected_to_local_elem_id_via_side(elem3, 2, 4);
    }
  }
};
TEST_F(ElemGraphChangeOwnerMoveNeighborsToEnd, withAura)
{
  run_test_on_num_procs(3, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveNeighborsToEnd, withoutAura)
{
  run_test_on_num_procs(3, stk::mesh::BulkData::NO_AUTO_AURA);
}


class ElemGraphChangeOwnerMoveNeighborsToDifferentProcs : public ElemGraphChangeOwner
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh_with_cyclic_decomp("1x2x2", auraOption);
    expect_mesh_created_correctly();
    expect_graph_correct_after_moving_neighbors_to_different_procs();
  }
  void expect_mesh_created_correctly()
  {
    expect_element_on_proc(1, 0);
    expect_element_on_proc(2, 1);
    expect_element_on_proc(3, 2);
    expect_element_on_proc(4, 3);
  }
  void expect_element_on_proc(stk::mesh::EntityId id, int proc)
  {
    if(get_bulk().parallel_rank() == proc)
    {
      EXPECT_TRUE(is_owned_on_this_proc(get_bulk().get_entity(stk::topology::ELEMENT_RANK, id)));
    }
  }
  void expect_graph_correct_after_moving_neighbors_to_different_procs()
  {
    create_elem_graph();
    check_initial_graph();
    move_elements({stk::mesh::EntityIdProc(1, 1), stk::mesh::EntityIdProc(3, 3)});

    if (get_bulk().parallel_rank() == 1) {
      stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
      expect_connected_to_remote_elem_id_via_side(elem1, 3, 5);
      const stk::mesh::impl::ParallelInfo& parInfo = get_elem_graph().get_const_parallel_edge_info(elem1, 5, 3, 4);
      expect_otherProc_permutation_chosenId(parInfo, 3, 4);
    }
    else if (get_bulk().parallel_rank() == 3) {
      stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
      expect_connected_to_remote_elem_id_via_side(elem3, 1, 4);
      const stk::mesh::impl::ParallelInfo& parInfo = get_elem_graph().get_const_parallel_edge_info(elem3, 4, 1, 5);
      expect_otherProc_permutation_chosenId(parInfo, 1, 4);
    }
  }
  void check_initial_graph()
  {
    expect_element1_connected_correctly();
    expect_element2_connected_correctly();
  }

  void expect_element1_connected_correctly()
  {
    stk::mesh::Entity element1 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
    if(is_owned_on_this_proc(element1))
      expect_element1_connected_to_2_and_3(element1);
  }

  void expect_element1_connected_to_2_and_3(stk::mesh::Entity element1)
  {
    ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(element1));
    expect_connected_to_remote_elem_id_on_proc(element1, 0, 2, 0, 1);
    expect_connected_to_remote_elem_id_on_proc(element1, 1, 3, 4, 2);
  }

  void expect_element2_connected_correctly()
  {
    stk::mesh::Entity element2 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2);
    if(is_owned_on_this_proc(element2))
      expect_element2_connected_to_1_and_4(element2);
  }

  void expect_element2_connected_to_1_and_4(stk::mesh::Entity element2)
  {
    ASSERT_EQ(2u, get_elem_graph().get_num_connected_elems(element2));
    expect_connected_to_remote_elem_id_on_proc(element2, 0, 1, 2, 0);
    expect_connected_to_remote_elem_id_on_proc(element2, 1, 4, 4, 3);
  }

  void expect_connected_to_remote_elem_id_on_proc(stk::mesh::Entity elem,
                                                  size_t connectedIndex,
                                                  stk::mesh::EntityId connectedId,
                                                  int side2,
                                                  int expectedProc)
  {
    expect_connected_to_remote_elem_id(elem, connectedIndex, connectedId);
    int side1 = get_elem_graph().get_connected_remote_id_and_via_side(elem, connectedIndex).side;
    const stk::mesh::impl::ParallelInfo &parInfo = get_elem_graph().get_parallel_edge_info(elem, side1, connectedId, side2);
    EXPECT_EQ(expectedProc, parInfo.get_proc_rank_of_neighbor());
  }
};
TEST_F(ElemGraphChangeOwnerMoveNeighborsToDifferentProcs, withAura)
{
  run_test_on_num_procs(4, stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(ElemGraphChangeOwnerMoveNeighborsToDifferentProcs, withoutAura)
{
  run_test_on_num_procs(4, stk::mesh::BulkData::NO_AUTO_AURA);
}

std::shared_ptr<stk::mesh::BulkData> build_mesh(stk::ParallelMachine comm,
                                                stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                                unsigned spatialDim = 0,
                                                const std::vector<std::string>& entityRankNames = std::vector<std::string>())
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
  builder.set_entity_rank_names(entityRankNames);
  builder.set_aura_option(auraOption);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  return bulk;
}

void change_entity_owner_hex_test_2_procs(bool aura_on)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int proc = stk::parallel_machine_rank(comm);

  if(stk::parallel_machine_size(comm) == 2)
  {
    stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
    if (!aura_on)
    {
      aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
    }
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, aura_option);
    stk::mesh::BulkData& bulkData = *bulkPtr;

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
    EXPECT_EQ(2, numLocallyOwnedElems);

    ElemElemGraphTester elem_graph(bulkData);
    bulkData.register_observer(std::make_shared<stk::mesh::ElemElemGraphUpdater>(bulkData, elem_graph));

    // Create a vector of the elements to be moved
    std::vector <stk::mesh::Entity> elems_to_move;

    stk::mesh::EntityId elem_2_id = 2;
    std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
    stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

    if (proc == 0)
    {
      elems_to_move.push_back(elem_2);

      int other_proc = 1;
      for (unsigned i=0; i<elems_to_move.size(); i++)
      {
        EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
        EXPECT_EQ(0, bulkData.parallel_owner_rank(elems_to_move[i]));
        elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
      }
    }

    bulkData.change_entity_owner(elem_proc_pairs_to_move);

    elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

    if (proc == 0)
    {
      ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, 5, stk::mesh::EntityId(3), 4), std::logic_error);
    }

    if (proc == 1)
    {
      EXPECT_TRUE(bulkData.is_valid(elem_2));
      EXPECT_EQ(1, bulkData.parallel_owner_rank(elem_2));

      EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_2));

      expect_elem_connected_to_remote_elem_id_via_side(bulkData, elem_graph, elem_2, 1, 4);
      expect_elem_connected_to_local_elem_id_via_side(bulkData, elem_graph, elem_2, 3, 5);

      stk::mesh::Entity elem_3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
      ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_3, 4, stk::mesh::EntityId(2), 5), std::logic_error);
    }

    EXPECT_EQ(1u, elem_graph.num_parallel_edges());
  }
}

TEST(ElementGraph, test_change_entity_owner_2_procs_hex_mesh_with_aura)
{
  bool aura_on = true;
  change_entity_owner_hex_test_2_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_2_procs_hex_mesh_without_aura)
{
  bool aura_on = false;
  change_entity_owner_hex_test_2_procs(aura_on);
}

void change_entity_owner_then_death_hex_test_2_procs(bool aura_on)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int proc = stk::parallel_machine_rank(comm);

  if(stk::parallel_machine_size(comm) == 2)
  {
    unsigned spatial_dim = 3;
    stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
    if (!aura_on)
    {
      aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
    }
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, aura_option, spatial_dim, stk::mesh::entity_rank_names());
    stk::mesh::BulkData& bulkData = *bulkPtr;
    stk::mesh::MetaData& meta = bulkData.mesh_meta_data();

    stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
    stk::mesh::Part& active = meta.declare_part("active", stk::topology::ELEMENT_RANK);
    stk::mesh::PartVector boundary_mesh_parts { &faces_part };

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    stk::unit_test_util::put_mesh_into_part(bulkData, active);

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
    EXPECT_EQ(2, numLocallyOwnedElems);

    stk::mesh::ElemElemGraph &elem_graph = bulkData.get_face_adjacent_element_graph();

    // Create a vector of the elements to be moved
    std::vector <stk::mesh::Entity> elems_to_move;

    stk::mesh::EntityId elem_2_id = 2;
    std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
    stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

    if (proc == 0)
    {
      elems_to_move.push_back(elem_2);

      int other_proc = 1;
      for (unsigned i=0; i<elems_to_move.size(); i++)
      {
        EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
        EXPECT_EQ(0, bulkData.parallel_owner_rank(elems_to_move[i]));
        elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
      }
    }

    bulkData.change_entity_owner(elem_proc_pairs_to_move);

    elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem_2_id);

    stk::mesh::EntityVector killedElements;
    std::vector<stk::mesh::PartVector> add_parts, remove_parts;
    if (proc == 1)
    {
      killedElements.push_back(elem_2);
      add_parts.push_back(stk::mesh::PartVector());
      remove_parts.push_back(stk::mesh::PartVector{&active});
    }

    bulkData.batch_change_entity_parts(killedElements, add_parts, remove_parts);
    boundary_mesh_parts.push_back(&active);

    stk::mesh::impl::ParallelSelectedInfo remoteActiveSelector;
    stk::mesh::impl::populate_selected_value_for_remote_elements(bulkData, elem_graph, bulkData.mesh_meta_data().universal_part(), remoteActiveSelector);

    process_killed_elements(bulkData, killedElements, active, remoteActiveSelector, boundary_mesh_parts, &boundary_mesh_parts);

    if (proc == 1)
    {
      EXPECT_TRUE(bulkData.is_valid(elem_2));
      EXPECT_EQ(1, bulkData.parallel_owner_rank(elem_2));
      EXPECT_FALSE(bulkData.bucket(elem_2).member(active));

      EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_2));

      expect_elem_connected_to_remote_elem_id_via_side(bulkData, elem_graph, elem_2, 1, 4);
      expect_elem_connected_to_local_elem_id_via_side(bulkData, elem_graph, elem_2, 3, 5);
    }
    if (proc == 0)
    {
      EXPECT_FALSE(remoteActiveSelector[-2]);

      ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_2, 5, stk::mesh::EntityId(3), 4), std::logic_error);
    }

    EXPECT_EQ(1u, elem_graph.num_parallel_edges());
  }
}

TEST(ElementGraph, test_change_entity_owner_and_death_hex_mesh_2_procs_with_aura)
{
  bool aura_on = true;
  change_entity_owner_then_death_hex_test_2_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_and_death_hex_mesh_2_procs_without_aura)
{
  bool aura_on = false;
  change_entity_owner_then_death_hex_test_2_procs(aura_on);
}


void setup_hex_shell_hex_mesh(stk::mesh::BulkData& bulkData)
{
  //
  //                proc 0               proc 1           proc 2
  //
  //               block_1          |   block_2  |      block_3
  //
  //          3---------------7        7            7-------------11
  //          /|             /|       /|           /|             /|
  //         / |            / |      / |          / |            / |
  //        /  |           /  |     /  |         /  |           /  |
  //       4--------------8   |    8   |        8--------------12  |
  //       |   |          |   |    |   |        |   |          |   |
  //       |   |   1      |   |    | 2 |        |   |   3      |   |
  //       |   |          |   |    |   |        |   |          |   |
  //       |   2----------|---6    |   6        |   6----------|---10
  //       |  /           |  /     |  /         |  /           |  /
  //       | /            | /      | /          | /            | /
  //       |/             |/       |/           |/             |/
  //       1--------------5        5            5--------------9

  stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
  unsigned spatial_dimension = 3;
  meta.initialize(spatial_dimension, stk::mesh::entity_rank_names());

  stk::mesh::Field<double>& field = meta.declare_field<double>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), nullptr);

  stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  stk::mesh::Part& block_2 = meta.declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);
  stk::mesh::Part& block_3 = meta.declare_part_with_topology("block_3", stk::topology::HEX_8);
  meta.commit();

  bulkData.modification_begin();

  stk::mesh::EntityIdVector elem1_nodes {1, 2, 3, 4, 5, 6, 7, 8};
  stk::mesh::EntityIdVector elem2_nodes {5, 6, 7, 8};
  stk::mesh::EntityIdVector elem3_nodes {5, 6, 7, 8, 9, 10, 11, 12};

  stk::mesh::EntityId elemId = 1;
  if (bulkData.parallel_rank() == 0) {
    stk::mesh::declare_element(bulkData, block_1, elemId, elem1_nodes);
    stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
    stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
    stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
    stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
    bulkData.add_node_sharing(node5, 1);
    bulkData.add_node_sharing(node6, 1);
    bulkData.add_node_sharing(node7, 1);
    bulkData.add_node_sharing(node8, 1);
    bulkData.add_node_sharing(node5, 2);
    bulkData.add_node_sharing(node6, 2);
    bulkData.add_node_sharing(node7, 2);
    bulkData.add_node_sharing(node8, 2);
  }
  else if (bulkData.parallel_rank() == 1) {
    elemId = 2;
    stk::mesh::declare_element(bulkData, block_2, elemId, elem2_nodes);
    stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
    stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
    stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
    stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
    bulkData.add_node_sharing(node5, 0);
    bulkData.add_node_sharing(node6, 0);
    bulkData.add_node_sharing(node7, 0);
    bulkData.add_node_sharing(node8, 0);
    bulkData.add_node_sharing(node5, 2);
    bulkData.add_node_sharing(node6, 2);
    bulkData.add_node_sharing(node7, 2);
    bulkData.add_node_sharing(node8, 2);
  }
  else if (bulkData.parallel_rank() == 2) {
    elemId = 3;
    stk::mesh::declare_element(bulkData, block_3, elemId, elem3_nodes);
    stk::mesh::Entity node5 = bulkData.get_entity(stk::topology::NODE_RANK, 5);
    stk::mesh::Entity node6 = bulkData.get_entity(stk::topology::NODE_RANK, 6);
    stk::mesh::Entity node7 = bulkData.get_entity(stk::topology::NODE_RANK, 7);
    stk::mesh::Entity node8 = bulkData.get_entity(stk::topology::NODE_RANK, 8);
    bulkData.add_node_sharing(node5, 0);
    bulkData.add_node_sharing(node6, 0);
    bulkData.add_node_sharing(node7, 0);
    bulkData.add_node_sharing(node8, 0);
    bulkData.add_node_sharing(node5, 1);
    bulkData.add_node_sharing(node6, 1);
    bulkData.add_node_sharing(node7, 1);
    bulkData.add_node_sharing(node8, 1);
  }

  bulkData.modification_end();
}

void change_entity_owner_hex_shell_hex_test_3_procs(bool aura_on)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int proc = stk::parallel_machine_rank(comm);
  if(stk::parallel_machine_size(comm) == 3)
  {
    stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
    if (!aura_on)
    {
      aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
    }
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, aura_option);
    stk::mesh::BulkData& bulkData = *bulkPtr;

    setup_hex_shell_hex_mesh(bulkData);

    bulkData.initialize_face_adjacent_element_graph();
    stk::mesh::ElemElemGraph &elementGraph = bulkData.get_face_adjacent_element_graph();

    const stk::mesh::Entity hex1   = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
    const stk::mesh::Entity hex3   = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
    const stk::mesh::Entity shell2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);

    if (proc == 0) {
      // Connectivity for Hex Element 1
      EXPECT_EQ(1u, elementGraph.get_num_connected_elems(hex1));
      EXPECT_EQ(5,  elementGraph.get_connected_remote_id_and_via_side(hex1, 0).side);
      EXPECT_EQ(2u, elementGraph.get_connected_remote_id_and_via_side(hex1, 0).id);
      EXPECT_FALSE(elementGraph.is_connected_elem_locally_owned(hex1, 0));
      EXPECT_EQ(1u, elementGraph.num_edges());
      EXPECT_EQ(1u, elementGraph.num_parallel_edges());
    }
    else if (proc == 1) {
      // Connectivity for Shell Element 2
      unsigned hex1Index = 0;
      unsigned hex3Index = 1;
      EXPECT_EQ(2u, elementGraph.get_num_connected_elems(shell2));
      EXPECT_EQ(0,  elementGraph.get_connected_remote_id_and_via_side(shell2, hex3Index).side);
      EXPECT_EQ(1,  elementGraph.get_connected_remote_id_and_via_side(shell2, hex1Index).side);
      EXPECT_EQ(3u, elementGraph.get_connected_remote_id_and_via_side(shell2, hex3Index).id);
      EXPECT_EQ(1u, elementGraph.get_connected_remote_id_and_via_side(shell2, hex1Index).id);
      EXPECT_FALSE(elementGraph.is_connected_elem_locally_owned(shell2, hex3Index));
      EXPECT_FALSE(elementGraph.is_connected_elem_locally_owned(shell2, hex1Index));
      EXPECT_EQ(2u, elementGraph.num_edges());
      EXPECT_EQ(2u, elementGraph.num_parallel_edges());
    }
    else if (proc == 2) {
      // Connectivity for Hex Element 3
      EXPECT_EQ(1u, elementGraph.get_num_connected_elems(hex3));
      EXPECT_EQ(4,  elementGraph.get_connected_remote_id_and_via_side(hex3, 0).side);
      EXPECT_EQ(2u, elementGraph.get_connected_remote_id_and_via_side(hex3, 0).id);
      EXPECT_FALSE(elementGraph.is_connected_elem_locally_owned(hex3, 0));
      EXPECT_EQ(1u, elementGraph.num_edges());
      EXPECT_EQ(1u, elementGraph.num_parallel_edges());
    }

    stk::mesh::EntityId elem_to_move_global_id = 2;
    std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
    stk::mesh::Entity elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_to_move_global_id);

    if (proc == 1)
    {
      int destination_proc = 2;
      elem_proc_pairs_to_move.push_back(std::make_pair(elem_to_move, destination_proc));
    }

    bulkData.change_entity_owner(elem_proc_pairs_to_move);

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numLocallyOwnedElemsInMesh = counts[stk::topology::ELEM_RANK];

    size_t size_of_elem_graph = elementGraph.size();

    if (proc == 0)
    {
      EXPECT_EQ(1, numLocallyOwnedElemsInMesh);
      EXPECT_EQ(1u, size_of_elem_graph);

      stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
      stk::mesh::impl::ParallelInfo& elem1_to_elem2_info = elementGraph.get_parallel_edge_info(elem_1, 5, stk::mesh::EntityId(2), 1);
      EXPECT_EQ(2, elem1_to_elem2_info.get_proc_rank_of_neighbor());
      EXPECT_EQ(1u, elementGraph.num_edges());
      EXPECT_EQ(1u, elementGraph.num_parallel_edges());
    }
    if (proc == 1)
    {
      EXPECT_EQ(0, numLocallyOwnedElemsInMesh);
      EXPECT_EQ(0u, size_of_elem_graph);

      stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
      ASSERT_THROW(elementGraph.get_parallel_edge_info(elem_2, 1, stk::mesh::EntityId(1), 5), std::logic_error);

      ASSERT_THROW(elementGraph.get_parallel_edge_info(elem_2, 0, stk::mesh::EntityId(3), 4), std::logic_error);
      EXPECT_EQ(0u, elementGraph.num_edges());
      EXPECT_EQ(0u, elementGraph.num_parallel_edges());
    }
    if (proc == 2)
    {
      EXPECT_EQ(2, numLocallyOwnedElemsInMesh);
      EXPECT_EQ(2u, size_of_elem_graph);

      stk::mesh::Entity elem_2 = bulkData.get_entity(stk::topology::ELEM_RANK, 2);
      stk::mesh::impl::ParallelInfo& elem2_to_elem1_info = elementGraph.get_parallel_edge_info(elem_2, 1, stk::mesh::EntityId(1), 5);
      EXPECT_EQ(0, elem2_to_elem1_info.get_proc_rank_of_neighbor());
      EXPECT_EQ(3u, elementGraph.num_edges());
      EXPECT_EQ(1u, elementGraph.num_parallel_edges());
    }
  }
}

TEST(ElementGraph, test_change_entity_owner_3_procs_hex_shell_hex_mesh_with_aura)
{
  bool aura_on = true;
  change_entity_owner_hex_shell_hex_test_3_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_3_procs_hex_shell_hex_mesh_without_aura)
{
  bool aura_on = false;
  change_entity_owner_hex_shell_hex_test_3_procs(aura_on);
}

void change_entity_owner_hex_test_4_procs(bool aura_on)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int proc = stk::parallel_machine_rank(comm);
  std::vector<double> wall_times;
  wall_times.reserve(10);
  std::vector<std::string> msgs;
  msgs.reserve(10);

  std::vector<size_t> mem_usage;

  wall_times.push_back(stk::wall_time());
  msgs.push_back("program-start");
  mem_usage.push_back(stk::get_memory_usage_now());

  if(stk::parallel_machine_size(comm) == 4)
  {
    stk::mesh::BulkData::AutomaticAuraOption aura_option = stk::mesh::BulkData::AUTO_AURA;
    if (!aura_on)
    {
      aura_option = stk::mesh::BulkData::NO_AUTO_AURA;
    }
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, aura_option);
    stk::mesh::BulkData& bulkData = *bulkPtr;

    stk::io::fill_mesh("generated:1x1x4", bulkData);

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after mesh-read");
    mem_usage.push_back(stk::get_memory_usage_now());

    std::vector<size_t> counts;
    stk::mesh::count_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData, counts);
    int numLocallyOwnedElems = counts[stk::topology::ELEM_RANK];
    EXPECT_EQ(1, numLocallyOwnedElems);

    stk::mesh::ElemElemGraph &elem_graph = bulkData.get_face_adjacent_element_graph();

    wall_times.push_back(stk::wall_time());
    msgs.push_back("after fill-graph");
    mem_usage.push_back(stk::get_memory_usage_now());

    // Create a vector of the elements to be moved
    std::vector <stk::mesh::Entity> elems_to_move;

    stk::mesh::EntityId elem_global_id = 2;
    std::vector< std::pair< stk::mesh::Entity, int > > elem_proc_pairs_to_move;
    stk::mesh::Entity elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_global_id);
    if (proc == 1)
    {
      elems_to_move.push_back(elem_to_move);

      int other_proc = 2;
      for (unsigned i=0; i<elems_to_move.size(); i++)
      {
        EXPECT_TRUE(bulkData.is_valid(elems_to_move[i]));
        EXPECT_EQ(1, bulkData.parallel_owner_rank(elems_to_move[i]));
        elem_proc_pairs_to_move.push_back(std::make_pair(elems_to_move[i], other_proc));
      }
    }

    bulkData.change_entity_owner(elem_proc_pairs_to_move);

    elem_to_move = bulkData.get_entity(stk::topology::ELEM_RANK, elem_global_id);

    if (proc == 2)
    {
      EXPECT_TRUE(bulkData.is_valid(elem_to_move));
      EXPECT_EQ(2, bulkData.parallel_owner_rank(elem_to_move));

      EXPECT_EQ(2u, elem_graph.get_num_connected_elems(elem_to_move));

      expect_elem_connected_to_remote_elem_id_via_side(bulkData, elem_graph, elem_to_move, 1, 4);
      expect_elem_connected_to_local_elem_id_via_side(bulkData, elem_graph, elem_to_move, 3, 5);

      stk::mesh::Entity elem_3 = bulkData.get_entity(stk::topology::ELEM_RANK, 3);
      ASSERT_THROW(elem_graph.get_parallel_edge_info(elem_3, 4, stk::mesh::EntityId(2), 5), std::logic_error);

      EXPECT_EQ(2u, elem_graph.num_parallel_edges());
    }
    else if (proc == 0)
    {
      stk::mesh::Entity elem_1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1);
      stk::mesh::impl::ParallelInfo &elem_1_to_2_p_info = elem_graph.get_parallel_edge_info(elem_1, 5, stk::mesh::EntityId(2), 4);
      EXPECT_EQ(2, elem_1_to_2_p_info.get_proc_rank_of_neighbor());
      EXPECT_EQ(1u, elem_graph.num_parallel_edges());
    }
    else if (proc == 1)
    {
      EXPECT_EQ(0u, elem_graph.size());
      EXPECT_EQ(0u, elem_graph.num_edges());
      EXPECT_EQ(0u, elem_graph.num_parallel_edges());
    }
    else if (proc == 3)
    {
      EXPECT_EQ(1u, elem_graph.num_parallel_edges());
    }
  }
}

TEST(ElementGraph, test_change_entity_owner_4_procs_hex_mesh_with_aura)
{
  bool aura_on = true;
  change_entity_owner_hex_test_4_procs(aura_on);
}

TEST(ElementGraph, test_change_entity_owner_4_procs_hex_mesh_without_aura)
{
  bool aura_on = false;
  change_entity_owner_hex_test_4_procs(aura_on);
}

} //namespace

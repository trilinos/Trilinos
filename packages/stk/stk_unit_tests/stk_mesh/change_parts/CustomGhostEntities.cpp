#include "gtest/gtest.h"
#include "mpi.h"
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace
{

class CustomGhostEntities: public stk::unit_test_util::MeshFixture
{
protected:
  CustomGhostEntities()
    : stk::unit_test_util::MeshFixture(),
      myPart(nullptr)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    myPart = &(get_meta().declare_part("myPart",stk::topology::FACE_RANK));
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
      stk::io::fill_mesh("generated:2x1x2|sideset:z", get_bulk());
      get_face_entities();
      add_myPart_to_both_faces();
    }
  }

  void get_face_entities()
  {
    face15 = get_bulk().get_entity(stk::topology::FACE_RANK, 15);
    face25 = get_bulk().get_entity(stk::topology::FACE_RANK, 25);
  }

  void add_myPart_to_face15_and_face25()
  {
    if (get_bulk().is_valid(face15)) {
      get_bulk().change_entity_parts(face15, stk::mesh::ConstPartVector{myPart});
    }
    if (get_bulk().is_valid(face25)) {
      get_bulk().change_entity_parts(face25, stk::mesh::ConstPartVector{myPart});
    }
  }

  void add_myPart_to_both_faces()
  {
    get_bulk().modification_begin();
    add_myPart_to_face15_and_face25();
    get_bulk().modification_end();
  }

  void send_face_15_from_proc_0_to_proc_1(stk::mesh::Ghosting& ghosting)
  {
    stk::mesh::EntityProcVec entitiesToGhost;
    if (get_bulk().parallel_rank() == 0) {
      entitiesToGhost.push_back(stk::mesh::EntityProc(face15, proc1));
    }
    get_bulk().change_ghosting(ghosting, entitiesToGhost);
  }

  void custom_ghost_face_15_from_proc_0_to_proc_1()
  {
    get_bulk().modification_begin();
    stk::mesh::Ghosting& ghosting = get_bulk().create_ghosting("myCustomGhosting");
    send_face_15_from_proc_0_to_proc_1(ghosting);
    get_bulk().modification_end();
  }

  void remove_face_15_from_myPart_on_proc_0()
  {
    get_bulk().modification_begin();
    if (get_bulk().parallel_rank() == 0) {
      get_bulk().change_entity_parts(face15, stk::mesh::ConstPartVector{}, stk::mesh::ConstPartVector{myPart});
    }
    get_bulk().modification_end();
  }

  void expect_face15_not_in_myPart_on_any_proc()
  {
    face15 = get_bulk().get_entity(stk::topology::FACE_RANK, 15);
    EXPECT_FALSE(get_bulk().bucket(face15).member(*myPart));
  }

  void expect_proc_0_has_4_nodes_in_myPart()
  {
    if (get_bulk().parallel_rank() == 0) {
      EXPECT_EQ(4u, stk::mesh::count_selected_entities(*myPart, get_bulk().buckets(stk::topology::NODE_RANK)));
    }
  }

  void expect_proc_1_has_2_nodes_in_myPart()
  {
    if (get_bulk().parallel_rank() == 1) {
      EXPECT_EQ(2u, stk::mesh::count_selected_entities(*myPart, get_bulk().buckets(stk::topology::NODE_RANK)));
    }
  }

  void check_results()
  {
    expect_face15_not_in_myPart_on_any_proc();
    expect_proc_0_has_4_nodes_in_myPart();
    expect_proc_1_has_2_nodes_in_myPart();
  }

private:
  stk::mesh::Part* myPart;
  stk::mesh::Entity face15, face25;
  int proc1 = 1;
};

TEST_F(CustomGhostEntities, removePart)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
    custom_ghost_face_15_from_proc_0_to_proc_1();
    remove_face_15_from_myPart_on_proc_0();
    check_results();
  }
}

}

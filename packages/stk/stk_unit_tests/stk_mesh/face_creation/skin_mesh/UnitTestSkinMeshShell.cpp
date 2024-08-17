#include <gtest/gtest.h>

#include <stk_unit_test_utils/BuildMesh.hpp>
#include "UnitTestSkinMeshUseCaseUtils.hpp"

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;
using stk::unit_test_util::build_mesh;

TEST(ElementGraph, two_wedge_sandwich_with_quad_shell)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) <= 3)
  {
    std::string fileName("made_up.g");
    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);
      stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
      stk::mesh::BulkData& bulk_data = *bulkPtr;
      stk::mesh::fixtures::VectorFieldType & node_coord =
          meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates");
      stk::mesh::put_field_on_mesh(node_coord, meta_data.universal_part(), 3, nullptr);

      stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("wedges", stk::topology::WEDGE_6));
      stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("quad_shells", stk::topology::SHELL_QUAD_4));

      meta_data.commit();

      static const char method[] =
          "stk_mesh::fixtures::heterogenous_mesh_bulk_data";
      ////////////////////////////////

      bulk_data.modification_begin();

      stk::mesh::Part & wedge_block = *meta_data.get_part("wedges", method);
      stk::mesh::Part & quad_shell_block = *meta_data.get_part("quad_shells", method);

      unsigned elem_id = 1;

      unsigned number_wedge = 2;

      stk::mesh::EntityIdVector wedge_node_ids[2] = {
        {1, 7, 3, 2, 8, 4},
        {7, 5, 3, 8, 6, 4}
      };

      for(unsigned i = 0; i < number_wedge; ++i, ++elem_id)
      {
        stk::mesh::declare_element(bulk_data, wedge_block, elem_id, wedge_node_ids[i]);
      }

      stk::mesh::EntityIdVector shell_quad_node_ids = {3, 4, 8, 7};

      stk::mesh::declare_element(bulk_data, quad_shell_block, elem_id, shell_quad_node_ids);

      const unsigned node_count = 8;
      static const double node_coord_data[node_count][3] = {
        {0, 0, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 1},
        {1, 1, 0}, {1, 1, 1}, {0, 1, 0}, {0, 1, 1}};

      for(unsigned i = 0; i < node_count; ++i)
      {
        stk::mesh::Entity const node = bulk_data.get_entity(stk::topology::NODE_RANK, i + 1);

        double * const coord = stk::mesh::field_data(node_coord, node);

        coord[0] = node_coord_data[i][0];
        coord[1] = node_coord_data[i][1];
        coord[2] = node_coord_data[i][2];
      }

      bulk_data.modification_end();
      ////////////////////////////////

      if(stk::parallel_machine_rank(comm) == 0)
      {
        stk::io::write_mesh(fileName, bulk_data);
      }
    }
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, comm, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& bulk_data = *bulkPtr;
    stk::mesh::Part &skin = meta_data.declare_part("skin", meta_data.side_rank());
    stk::io::put_io_part_attribute(skin);
    stk::unit_test_util::read_from_serial_file_and_decompose(fileName, bulk_data, "RIB");
    unlink(fileName.c_str());
    EXPECT_NO_FATAL_FAILURE(ElemGraphTestUtils::skin_boundary(bulk_data, meta_data.locally_owned_part(), {&skin}));
    std::vector<size_t> mesh_counts;
    stk::mesh::comm_mesh_counts(bulk_data, mesh_counts);
    EXPECT_EQ(8u, mesh_counts[meta_data.side_rank()]);

    std::vector<std::pair<stk::mesh::EntityId, int>> id_and_num_faces = {
      {1, 4},
      {2, 4},
      {3, 0}
    };

    for(size_t i = 0; i < id_and_num_faces.size(); ++i)
    {
      stk::mesh::EntityId id = id_and_num_faces[i].first;
      int gold_num_faces = id_and_num_faces[i].second;
      stk::mesh::Entity elem = bulk_data.get_entity(stk::topology::ELEM_RANK, id);
      if(bulk_data.is_valid(elem))
      {
        int num_faces = bulk_data.num_sides(elem);
        EXPECT_EQ(gold_num_faces, num_faces)<< "element " << id << " has topology " << bulk_data.bucket(elem).topology() << " with num faces " << num_faces << " not same as gold value " << gold_num_faces << std::endl;
      }
    }
  }
}

}

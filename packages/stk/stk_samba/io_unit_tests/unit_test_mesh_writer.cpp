#include <gtest/gtest.h>

#include <Ioss_SubSystem.h>

#include <samba_io/mesh_reader.hpp>
#include <samba_io/mesh_reader_fixture.hpp>

#include <samba_io/mesh_writer.hpp>

#include <init/Ionit_Initializer.h>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

#include "unit_test_utils.hpp"

void check_same_nodes(samba::mesh mesh1, samba::io::ioss_mapper::Ptr mapper1,
                      samba::mesh mesh2, samba::io::ioss_mapper::Ptr mapper2)
{
  samba::io::entity_key_to_ioss_id_field to_global1 = mapper1->m_to_ioss_global_id;
  samba::io::entity_key_to_ioss_id_field to_global2 = mapper2->m_to_ioss_global_id;

  samba::io::entity_key_vector nodes1, nodes2;
  mesh1.get_entities(samba::entity_rank::node(), nodes1);
  mesh2.get_entities(samba::entity_rank::node(), nodes2);

  EXPECT_EQ(nodes1.size(), nodes2.size());
  size_t num_nodes = nodes1.size();

  std::map<int, samba::entity_key> global_to_nodes2;
  for (size_t i = 0; i < num_nodes; ++i)
  {
    int global_id = to_global2[nodes2[i]];
    global_to_nodes2[global_id] = nodes2[i];
  }

  samba::io::ioss_mapper &mapper_ref1 = *mapper1;
  samba::io::ioss_mapper &mapper_ref2 = *mapper2;

  bool nontrivial = false;
  for (size_t i = 0; i < num_nodes; ++i)
  {
    nontrivial = true;

    samba::entity_key node1 = nodes1[i];
    int global_id = to_global1[node1];
    samba::entity_key node2 = global_to_nodes2[global_id];

    // As long as the orginating mesh in the test has not had nodes deleted,
    // this can be expected.  Should the idea of trying to preserve local_ids
    // be ditched entirely, though?
    EXPECT_EQ(node1.local_id(), node2.local_id());

    bool nodes_equal = is_nodally_equal(node1, mesh1, mapper_ref1, node2, mesh2, mapper_ref2);
    EXPECT_TRUE(nodes_equal);
  }
  EXPECT_TRUE(nontrivial);
}


void check_same_elements(samba::mesh mesh1, samba::io::ioss_mapper::Ptr mapper1,
                         samba::mesh mesh2, samba::io::ioss_mapper::Ptr mapper2)
{
  samba::io::entity_key_to_ioss_id_field to_global1 = mapper1->m_to_ioss_global_id;
  samba::io::entity_key_to_ioss_id_field to_global2 = mapper2->m_to_ioss_global_id;

  samba::io::entity_key_vector elements1, elements2;
  mesh1.get_entities(samba::entity_rank::element(), elements1);
  mesh2.get_entities(samba::entity_rank::element(), elements2);

  EXPECT_EQ(elements1.size(), elements2.size());
  size_t num_elements = elements1.size();

  std::map<int, samba::entity_key> global_to_elements2;
  for (size_t i = 0; i < num_elements; ++i)
  {
    int global_id = to_global2[elements2[i]];
    global_to_elements2[global_id] = elements2[i];
  }

  samba::io::ioss_mapper &mapper_ref1 = *mapper1;
  samba::io::ioss_mapper &mapper_ref2 = *mapper2;

  bool nontrivial = false;
  for (size_t i = 0; i < num_elements; ++i)
  {
    nontrivial = true;

    samba::entity_key element1 = elements1[i];
    int global_id = to_global1[element1];
    samba::entity_key element2 = global_to_elements2[global_id];
    // EXPECT_EQ(element1.local_id(), element2.local_id());

    bool elements_equal = is_nodally_equal(element1, mesh1, mapper_ref1, element2, mesh2, mapper_ref2);
    EXPECT_TRUE(elements_equal);
  }
  EXPECT_TRUE(nontrivial);
}


void check_same_faces_on_elements(samba::mesh mesh1, samba::io::ioss_mapper::Ptr mapper1,
                                  samba::mesh mesh2, samba::io::ioss_mapper::Ptr mapper2)
{
  samba::io::entity_key_to_ioss_id_field to_global1 = mapper1->m_to_ioss_global_id;
  samba::io::entity_key_to_ioss_id_field to_global2 = mapper2->m_to_ioss_global_id;

  samba::io::entity_key_vector elements1, elements2;
  mesh1.get_entities(samba::entity_rank::element(), elements1);
  mesh2.get_entities(samba::entity_rank::element(), elements2);

  EXPECT_EQ(elements1.size(), elements2.size());
  size_t num_elements = elements1.size();

  std::map<int, samba::entity_key> global_to_elements2;
  for (size_t i = 0; i < num_elements; ++i)
  {
    int global_id = to_global2[elements2[i]];
    global_to_elements2[global_id] = elements2[i];
  }

  samba::io::ioss_mapper &mapper_ref1 = *mapper1;
  samba::io::ioss_mapper &mapper_ref2 = *mapper2;

  bool nontrivial = false;
  for (size_t i = 0; i < num_elements; ++i)
  {
    nontrivial = true;

    samba::entity_key element1 = elements1[i];
    int global_id = to_global1[element1];
    samba::entity_key element2 = global_to_elements2[global_id];

    samba::entity_proxy proxy1 = mesh1[element1];
    samba::entity_proxy proxy2 = mesh2[element2];

    samba::entity_key_iterator faces1_iter = proxy1.begin_faces<samba::entity_key>();
    samba::entity_key_iterator faces2_iter = proxy2.begin_faces<samba::entity_key>();
    for (; faces1_iter != proxy1.end_faces<samba::entity_key>(); ++faces1_iter, ++faces2_iter)
    {
      samba::entity_key face1 = *faces1_iter;
      samba::entity_key face2 = *faces2_iter;
      // EXPECT_EQ(face1, face2);

      if (face1 != samba::entity_key::invalid())
      {
        bool faces_equal = is_nodally_equal(face1, mesh1, mapper_ref1, face2, mesh2, mapper_ref2);
        EXPECT_TRUE(faces_equal);
      }
    }
  }
  EXPECT_TRUE(nontrivial);
}


TEST(samba_io, mesh_writer_generated_basic)
{
  samba::io::mesh_reader_fixture model1(0, "generated", "3x4x5|sideset:xXyYzZ");

  samba::mesh mesh1 = model1.m_mesh;
  samba::io::ioss_mapper::Ptr mapper1 = model1.m_ioss_mapper_ptr;

  samba::io::mesh_writer write1(0, "test_writer_1.exo", mesh1, mapper1);

  samba::io::mesh_reader_fixture model2(0, "exodus", "test_writer_1.exo");

  samba::mesh mesh2 = model2.m_mesh;

  size_t mesh2_num_nodes = mesh2.num_entities(samba::entity_rank::node());
  EXPECT_EQ(mesh2_num_nodes, 120u);

  size_t mesh2_num_elements = mesh2.num_entities(samba::entity_rank::element());
  EXPECT_EQ(mesh2_num_elements, 60u);

  size_t mesh2_num_faces = mesh2.num_entities(samba::entity_rank::face());
  EXPECT_EQ(mesh2_num_faces, 94u);
}


TEST(samba_io, mesh_writer_generated_detailed)
{
  samba::io::mesh_reader_fixture model1(0, "generated", "3x4x5|sideset:xXyYzZ");

  samba::mesh mesh1 = model1.m_mesh;
  samba::io::ioss_mapper::Ptr mapper1 = model1.m_ioss_mapper_ptr;

  samba::io::mesh_writer write1(0, "test_writer_1.exo", mesh1, mapper1);

  samba::io::mesh_reader_fixture model2(0, "exodus", "test_writer_1.exo");

  samba::mesh mesh2 = model2.m_mesh;
  samba::io::ioss_mapper::Ptr mapper2 = model2.m_ioss_mapper_ptr;

  //
  // Check correspondence of global exodus ids, samba entity keys, and coordinate values.
  //

  check_same_nodes(mesh1, mapper1, mesh2, mapper2);

  check_same_elements(mesh1, mapper1, mesh2, mapper2);

  check_same_faces_on_elements(mesh1, mapper1, mesh2, mapper2);
}



TEST(samba_io, mesh_writer_generated_nodal_field)
{
  typedef samba::field<double, samba::scalar_functor> real_field_type;

  samba::io::mesh_reader_fixture model1(MPI_COMM_NULL, "generated", "3x4x5|sideset:xXyYzZ");
  samba::mesh mesh1 = model1.m_mesh;
  samba::io::ioss_mapper::Ptr mapper1 = model1.m_ioss_mapper_ptr;
  real_field_type ids_as_reals1(mesh1, samba::entity_topology::node());

  samba::io::entity_key_vector nodes1;
  mesh1.get_entities(samba::entity_rank::node(), nodes1);
  for (size_t i = 0; i < nodes1.size(); ++i)
  {
    samba::entity_key curr_node = nodes1[i];
    samba::partition_index pi   = mesh1.convert(curr_node);
    ids_as_reals1[pi] = mapper1->m_to_ioss_global_id[pi];
  }

  samba::io::mesh_writer writer1(MPI_COMM_NULL, "test_writer_1.exo", mesh1, mapper1);

  const double time = 0.0;
  writer1.write_nodal_field(ids_as_reals1, "samba_saved_global_ids", time);
  // writer1.write_nodal_field(mapper1->m_node_coordinates, "samba_saved_node_coordinates", time);

  samba::io::mesh_reader_fixture model2(MPI_COMM_NULL, "exodus", "test_writer_1.exo");

  samba::mesh mesh2 = model2.m_mesh;
  samba::io::ioss_mapper::Ptr mapper2 = model2.m_ioss_mapper_ptr;
  real_field_type ids_as_reals2(mesh2, samba::entity_topology::node());

  model2.m_reader->read_nodal_field(ids_as_reals2, "samba_saved_global_ids", time);

  //
  // Check whether the global ids of the fields match.
  //
  bool entered = false;
  samba::io::entity_key_vector nodes2;
  mesh2.get_entities(samba::entity_rank::node(), nodes2);
  for (size_t i = 0; i < nodes2.size(); ++i)
  {
    samba::entity_key curr_node = nodes2[i];
    double gid_from_mapper = mapper2->m_to_ioss_global_id[curr_node];
    double gid_as_real = ids_as_reals2[curr_node];
    EXPECT_EQ(gid_from_mapper, gid_as_real);

    entered = true;
  }

  EXPECT_TRUE(entered);
}

#if 0
#include <gtest/gtest.h>

#include <samba/distributed_mesh.hpp>

template<typename SambaMeshT>
void connect_hex8_nodes_only(SambaMeshT mesh_arg, samba::entity_key elt,
                             const std::vector<samba::entity_key> &nodes)
{
  BOOST_ASSERT_MSG(nodes.size() == 8, "connect_hex8_nodes_only expected entity_keys for 8 nodes");

  for (size_t i = 0; i < 8; ++i)
  {
    mesh_arg.add_connectivity(elt, nodes[i], samba::connectivity_ordinal::create(i));
  }
}

TEST(samba_distributed, share_simplest)
{
  using namespace samba;

  process_id p0 = {0};
  process_id p1 = {1};

  // 4 hex_8s in a row (technically, a chain, since the test does not assign
  // nodal coordinates).  First two owned by p0, second two owned by p1.
  // p_0 owns the nodes where the middle two elements abut.
  // The middle two elements are shared

  //
  // Pair of adjacent hex_8s, all nodes owned by p0.
  //
  distributed_mesh mesh_p0(p0, 2, connectivity_map::default_map() );
  mesh_p0.begin_modification();
  entity_key_interval nodes_p0 = mesh_p0.add_entities(entity_topology::node(), 12);
  EXPECT_EQ(nodes_p0.size(), 12u);
  entity_key_interval hex_elts_p0 = mesh_p0.add_entities(entity_topology::hex_8(), 2);
  EXPECT_EQ(hex_elts_p0.size(), 2u);

  std::vector<entity_key> hex_nodes(8);

  hex_nodes[0] = nodes_p0[0];
  hex_nodes[1] = nodes_p0[1];
  hex_nodes[2] = nodes_p0[4];
  hex_nodes[3] = nodes_p0[3];
  hex_nodes[4] = nodes_p0[6];
  hex_nodes[5] = nodes_p0[7];
  hex_nodes[6] = nodes_p0[10];
  hex_nodes[7] = nodes_p0[9];
  connect_hex8_nodes_only(mesh_p0, hex_elts_p0[0], hex_nodes);

  hex_nodes[0] = nodes_p0[1];
  hex_nodes[1] = nodes_p0[2];
  hex_nodes[2] = nodes_p0[5];
  hex_nodes[3] = nodes_p0[4];
  hex_nodes[4] = nodes_p0[7];
  hex_nodes[5] = nodes_p0[8];
  hex_nodes[6] = nodes_p0[11];
  hex_nodes[7] = nodes_p0[10];
  connect_hex8_nodes_only(mesh_p0, hex_elts_p0[1], hex_nodes);

  std::vector<entity_key> shared_nodes;
  shared_nodes.push_back(nodes_p0[2]);
  shared_nodes.push_back(nodes_p0[5]);
  shared_nodes.push_back(nodes_p0[8]);
  shared_nodes.push_back(nodes_p0[11]);

  bool res = false;

  for (size_t i = 0, num_sns = shared_nodes.size(); i < num_sns; ++i)
  {
    res = mesh_p0.add_sharer(p1, shared_nodes[i]);
    EXPECT_TRUE(res);
  }
  res = mesh_p0.add_sharer(p1, hex_elts_p0[1]);
  EXPECT_TRUE(res);

  // Makes no sense for a mesh to share with itself.
  for (size_t i = 0, num_sns = shared_nodes.size(); i < num_sns; ++i)
  {
    res = mesh_p0.add_sharer(p0, shared_nodes[i]);
    EXPECT_FALSE(res);
  }
  res = mesh_p0.add_sharer(p0, hex_elts_p0[1]);
  EXPECT_FALSE(res);

  mesh_p0.end_modification();


  //
  // Pair of adjacent hex_8s. One shares 4 nodes that p0 owns.  p1 owns all nodes of the other.
  //
  distributed_mesh mesh_p1(p1, 2, connectivity_map::default_map() );
  mesh_p1.begin_modification();
  entity_key_interval nodes_p1 = mesh_p1.add_entities(entity_topology::node(), 8);
  EXPECT_EQ(nodes_p1.size(), 8u);
  entity_block_key none = entity_block_key::invalid();
  res =  mesh_p1.add_unowned_entities(shared_nodes.begin(), shared_nodes.end(),
                                      entity_state::shared(), &none, &none);
  EXPECT_TRUE(res);
  entity_key_interval hex_elts_p1 = mesh_p1.add_entities(entity_topology::hex_8(), 2);
  EXPECT_EQ(hex_elts_p1.size(), 2u);

  hex_nodes[0] = shared_nodes[0];
  hex_nodes[1] = nodes_p1[0];
  hex_nodes[2] = nodes_p1[2];
  hex_nodes[3] = shared_nodes[1];
  hex_nodes[4] = shared_nodes[2];
  hex_nodes[5] = nodes_p1[4];
  hex_nodes[6] = nodes_p1[6];
  hex_nodes[7] = shared_nodes[3];
  connect_hex8_nodes_only(mesh_p1, hex_elts_p1[0], hex_nodes);

  hex_nodes[0] = nodes_p1[0];
  hex_nodes[1] = nodes_p1[1];
  hex_nodes[2] = nodes_p1[3];
  hex_nodes[3] = nodes_p1[2];
  hex_nodes[4] = nodes_p1[4];
  hex_nodes[5] = nodes_p1[5];
  hex_nodes[6] = nodes_p1[7];
  hex_nodes[7] = nodes_p1[6];
  connect_hex8_nodes_only(mesh_p1, hex_elts_p1[1], hex_nodes);

  res = mesh_p1.add_sharer(p0, hex_elts_p1[0]);
  EXPECT_TRUE(res);

  mesh_p1.end_modification();

  //
  // Now tell meshes which entities they are ghosting.
  //
  mesh_p0.begin_modification();
  res =  mesh_p0.add_unowned_entities(hex_elts_p1.begin(), hex_elts_p1.begin() + 1,
                                      entity_state::ghosted(), &none, &none);
  EXPECT_TRUE(res);
  mesh_p0.end_modification();

  mesh_p1.begin_modification();
  res =  mesh_p1.add_unowned_entities(hex_elts_p0.begin() + 1, hex_elts_p0.begin() + 2,
                                      entity_state::ghosted(), &none, &none);
  EXPECT_TRUE(res);
  mesh_p1.end_modification();


  ////
  //// WHOSE REPSONSIBILITY IS IT TO UPDATE THE CONNECTIVITY INFORMATION FOR ENTITIES
  //// THAT ARE OWNED BY ANOTHER PROC AND SHARED/GHOSTED?  THE ALGORITHM?  THE MESH?
  ////

  ////
  //// NOTE NEED TO CHANGE THE FOLLOWING TWO SECTIONS IF GHOST PARTS AUTOMATICALLY
  //// DO INDUCED MEMBERSHIP.
  ////

  //
  // Check p0's comm list
  //
  distributed_mesh::comm_list comm_list_p0;
  mesh_p0.get_comm_list(comm_list_p0);
  EXPECT_EQ(comm_list_p0.size(), 1u);
  EXPECT_EQ((comm_list_p0[0].first)(), 1);
  std::vector<entity_key> &p0_to_p1 = *(comm_list_p0[0].second);
  EXPECT_EQ(p0_to_p1.size(), 5u);
  for (size_t i = 0, end_i = shared_nodes.size(); i < end_i; ++i)
  {
    EXPECT_TRUE(std::find(p0_to_p1.begin(), p0_to_p1.end(), shared_nodes[i]) != p0_to_p1.end());
  }
  EXPECT_TRUE(std::find(p0_to_p1.begin(), p0_to_p1.end(), hex_elts_p0[1]) != p0_to_p1.end());


  //
  // Check p1's comm list
  //
  distributed_mesh::comm_list comm_list_p1;
  mesh_p1.get_comm_list(comm_list_p1);
  EXPECT_EQ(comm_list_p1.size(), 1u);
  EXPECT_EQ((comm_list_p1[0].first)(), 0);
  std::vector<entity_key> &p1_to_p0 = *(comm_list_p1[0].second);
  EXPECT_EQ(p1_to_p0.size(), 1u);
  EXPECT_TRUE(std::find(p1_to_p0.begin(), p1_to_p0.end(), hex_elts_p1[0]) != p1_to_p0.end());

}

//// Need a test that mimics meshes being constructed separately and then merged.
//// Point is to exercise the mesh features that would be used by distributed algorithm
//// that discovers overlap, snips, and reattaches with sharing and ghosting.
#endif

#include <sierra/mesh/array_mesh/array_mesh.hpp>

#include <gtest/gtest.h>

#include <iostream>

//-------------------------------------------------------------------
TEST( array_mesh, topologies)
{
  using namespace sierra;
  using namespace sierra::mesh;

  int tet4_nodes = num_nodes<Tet4>::value;
  EXPECT_EQ(tet4_nodes, 4);

  int hex8_nodes = num_nodes<Hex8>::value;
  EXPECT_EQ(hex8_nodes, 8);
}

//-------------------------------------------------------------------
TEST( array_mesh, basic)
{
  using namespace sierra;
  using namespace sierra::mesh;

  array_mesh amesh;

  //This test will create a mesh with 1 hex8 element:

  const size_t num_elems = 1;
  const size_t num_nodes = 8;

  //create 1-based elem-id list and node-id list.

  std::vector<int> elem_ids(num_elems);
  for(size_t i=0; i<elem_ids.size(); ++i) elem_ids[i] = i+1;

  std::vector<int> node_ids(num_nodes);
  for(size_t i=0; i<node_ids.size(); ++i) node_ids[i] = i+1;

  //create 0-based node-index list
  std::vector<int> node_idx(num_nodes);
  for(size_t i=0; i<node_idx.size(); ++i) node_idx[i] = i;

  //give the node-ids and elem-ids to the mesh:
  //in exodus lingo, these are the node-number map and element-number map:

  amesh.add_node_ids(node_ids.begin(), node_ids.end());
  amesh.add_element_ids(elem_ids.begin(), elem_ids.end());

  EXPECT_EQ(num_elems, amesh.get_num_elements());
  EXPECT_EQ(num_nodes, amesh.get_num_nodes());

  int node_id = amesh.get_node_id(5);
  EXPECT_EQ(node_id, node_ids[5]);

  int block_id = 1;
  array_mesh::BlockIndex blk = amesh.add_block<Hex8>(array_mesh::Element, block_id, num_elems);

  const int nodes_per_elem = sierra::mesh::num_nodes<sierra::mesh::Hex8>::value;

  //add node connectivity for the first (and only) element,
  //using the first 'nodes_per_elem' node indices:
  int elem_num = 0;
  amesh.add_connectivity(blk, elem_num, node_idx.begin(), node_idx.begin()+nodes_per_elem);

  const std::vector<int>& blk_conn = amesh.get_block_connectivity(blk);

  //Now insist that the 5th node in the element's connectivity is the same
  //as the 5th node in our node_ids list.
  EXPECT_EQ(node_id, amesh.get_node_id(blk_conn[5]));
}


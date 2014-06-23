#include <sierra/mesh/array_mesh/array_mesh.hpp>

#include <gtest/gtest.h>

#include <iostream>

//-------------------------------------------------------------------
TEST( array_mesh, connections)
{
  using namespace sierra;
  using namespace sierra::mesh;

  bool create_upward_connectivity = true;
  array_mesh amesh(create_upward_connectivity);

  //This test will create a mesh with 2 hex8 elements:

  stk::topology elem_topo = stk::topology::HEX_8;
  const size_t num_elems = 2;
  const size_t num_nodes = 12;
  const int nodes_per_elem = elem_topo.num_nodes();

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

  //create an element-block. this tells the mesh how many elements will be in the
  //block, and what topology they will be, but doesn't create connectivity yet.
  int block_id = 1;
  array_mesh::BlockIndex blk = amesh.add_block(array_mesh::Element, block_id, num_elems,elem_topo);

  //Now specify the nodal connectivity for the elements in the block:

  //add node connectivity for the first element,
  //using the first 'nodes_per_elem' node indices:
  int elem_num = 0;
  amesh.add_connectivity(blk, elem_num, node_idx.begin(), node_idx.begin()+nodes_per_elem);

  //add node connectivity for the second element,
  //using the last 'nodes_per_elem' node indices:
  int elem_node_offset = num_nodes - nodes_per_elem;
  elem_num = 1;
  amesh.add_connectivity(blk, elem_num, node_idx.begin()+elem_node_offset, node_idx.end());

  array_mesh::ConstIntRange node_range = amesh.get_connected_nodes(elem_num);
  int num_connected_nodes = std::distance(node_range.first, node_range.second);
  EXPECT_EQ(num_connected_nodes, nodes_per_elem);

  int node_id = amesh.get_node_id(*node_range.first);
  EXPECT_EQ(node_id, node_ids[elem_node_offset]);

  //the 3rd node should have 1 connected element:
  int node_num1 = 3;
  array_mesh::ConstIntRange elem_range1 = amesh.get_connected_elements(node_num1);
  int num_elems1 = std::distance(elem_range1.first, elem_range1.second);
  EXPECT_EQ(num_elems1, 1);

  //the 4th node should have 2 connected element:
  int node_num2 = 4;
  array_mesh::ConstIntRange elem_range2 = amesh.get_connected_elements(node_num2);
  int num_elems2 = std::distance(elem_range2.first, elem_range2.second);
  EXPECT_EQ(num_elems2, 2);
}


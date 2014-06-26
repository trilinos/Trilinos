#include <gtest/gtest.h>

#include <vector>
#include <fstream>
#include <iostream>

#include <sierra/mesh/array_mesh/array_mesh.hpp>
#include <sierra/mesh/fixture/array_mesh_hex_fixture.hpp>

//-------------------------------------------------------------------
TEST( array_mesh, get_side_nodes )
{

  using namespace sierra;
  using namespace sierra::mesh;

  int Hex8_key[6][4]= { {0, 1, 5, 4},
                        {1, 2, 6, 5},
                        {2, 3, 7, 6},
                        {0, 4, 7, 3},
                        {0, 3, 2, 1},
                        {4, 5, 6, 7}};

  int Tet4_key[4][3] = {{8, 9, 11},
                        {9, 10, 11},
                        {8, 11, 10},
                        {8, 10, 9}};

  std::vector < std::vector < std::vector <int> > > gold_key(2);

  gold_key[0] = std::vector < std::vector < int > >(6,
                              std::vector < int >  (4));

  for(int i = 0; i < 6; ++i){
    for(int j = 0; j < 4; ++j){
      gold_key[0][i][j] = Hex8_key[i][j];
    }
  }

  gold_key[1] = std::vector < std::vector < int > >(4,
                              std::vector < int >  (3));
  for(int i = 0; i < 4; ++i){
    for(int j = 0; j < 3; ++j){
      gold_key[1][i][j] = Tet4_key[i][j];
    }
  }

  const bool create_upward_connectivity = false;
  array_mesh amesh(create_upward_connectivity);

  //This test will create a mesh with a HEX8 and TET4 elements:

  stk::topology Hex8_topo = stk::topology::HEX_8;
  stk::topology Tet4_topo = stk::topology::TET_4;
  const size_t num_elems_per_blk = 1;
  const size_t num_nodes = 12;
  const int Hex8_nodes_per_elem = Hex8_topo.num_nodes();
  //create 1-based elem-id list and node-id list.

  std::vector<int> elem_ids(2 * num_elems_per_blk);
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
  const int Hex8_block_id = 1;
  const int Tet4_block_id = 2;
  amesh.add_block(array_mesh::Element, Hex8_block_id, num_elems_per_blk, Hex8_topo);
  amesh.add_block(array_mesh::Element, Tet4_block_id, num_elems_per_blk, Tet4_topo);

  //add node connectivity for the Hex8 element,
  //using the first 'nodes_per_elem' node indices:
  int elem_num = 0;
  amesh.add_connectivity(amesh.get_block(Hex8_block_id), elem_num, node_idx.begin(), node_idx.begin()+Hex8_nodes_per_elem);

  //add node connectivity for the Tet4 element,
  //using the last 'nodes_per_elem' node indices:
  int elem_node_offset = Hex8_nodes_per_elem;

  amesh.add_connectivity(amesh.get_block(Tet4_block_id), elem_num, node_idx.begin()+elem_node_offset, node_idx.end());

  for(size_t ielem = 0; ielem < amesh.get_num_elements(); ++ielem) {
	  stk::topology topology = amesh.get_element_topology(ielem);
	  std::cout << topology << " " << __FILE__ << " " <<__LINE__ << std::endl;
	  size_t num_sides = topology.num_faces();

	  for(size_t iside = 0; iside < num_sides; ++iside) {
		  std::vector<int> side_node_vect;
		  amesh.get_side_nodes(ielem, iside, side_node_vect);

		  for(size_t inode = 0; inode < side_node_vect.size(); ++inode){
		    EXPECT_EQ ( gold_key[ielem][iside][inode] , side_node_vect[inode] );
		  }
	  }
  }
}

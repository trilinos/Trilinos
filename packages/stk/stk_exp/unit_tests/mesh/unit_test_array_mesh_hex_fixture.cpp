#include <sierra/mesh/fixture/array_mesh_hex_fixture.hpp>

#include <gtest/gtest.h>

#include <iostream>

//-------------------------------------------------------------------
TEST( array_mesh_hex_fixture, basic)
{
  using namespace sierra;
  using namespace sierra::mesh;

  const unsigned nx = 2;
  const unsigned ny = 3;
  const unsigned nz = 4;

  fixture::array_mesh_hex_fixture array_mesh_hex(nx, ny, nz);

  array_mesh_hex.generate_mesh();

  array_mesh& cmesh = array_mesh_hex.m_mesh;
  std::cout<<"num-nodes: "<<cmesh.get_num_nodes()<<std::endl;
  EXPECT_EQ(60u, cmesh.get_num_nodes());
}


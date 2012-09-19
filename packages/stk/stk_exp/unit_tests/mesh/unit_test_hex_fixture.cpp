#ifndef __IBMCPP__
#include <gtest/gtest.h>
#include <sierra/mesh/fixture/hex_fixture.hpp>
#include <sierra/mesh/details/selected_buckets.hpp>

#include <iostream>

TEST( hex_fixture, basic)
{
  using namespace sierra;
  using namespace sierra::mesh;
  using namespace sierra::mesh::fixture;
  using namespace sierra::mesh::details;

  unsigned nx=2, ny=3, nz=4;
  hex_fixture fixture(nx,ny,nz);

  fixture.generate_mesh();

  hex_fixture::CoordinateField& coord_field = fixture.m_coordinates;

  for(size_t i=0; i<fixture.m_num_nodes; ++i) {
    entity_key node = static_cast<entity_key>(i);
    double* coords = coord_field[fixture.m_mesh.get_bucket_location(node)];
    EXPECT_EQ(coords[0], (double)(i%(nx+1)));
    size_t tmpy = i/(nx+1);
    EXPECT_EQ(coords[1], (double)(tmpy%(ny+1)));
    size_t tmpz = tmpy/(ny+1);
    EXPECT_EQ(coords[2], (double)(tmpz%(nz+1)));
//    std::cout <<"node " << i << " coords: " << coords[0]<<", "<<coords[1]<<", "<<coords[2]<<std::endl;
  }
}
#endif

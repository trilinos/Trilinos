#include <gtest/gtest.h>

#include <samba/mesh.hpp>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <sstream>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/equal.hpp>

#include <samba_fixtures/tri_fixture.hpp>

TEST(samba, mesh_tri_fixture)
{
  samba::mesh mesh(samba::connectivity_map::default_map_2d());

  samba::fixtures::tri_fixture tf;
  tf.samba_mesh_create(mesh);

  EXPECT_TRUE(mesh.num_nodes() == 121);
  EXPECT_TRUE(mesh.num_elements() == 200);
}

#ifndef __IBMCPP__

#include <sierra/mesh/details/selected_buckets.hpp>
#include <sierra/mesh/modifiable_mesh_api.hpp>

#include <gtest/gtest.h>

using namespace sierra::mesh;
using namespace sierra::mesh::details;

TEST( selected_buckets , basic )
{
  boost::shared_ptr<modifiable_mesh> mesh(new modifiable_mesh);

  part_key part1 = mesh->declare_part("part1");
  part_key part2 = mesh->declare_part("part2");

  std::vector<entity_key> keys;
  for ( size_t i = 0; i<10; ++i) {
    keys.push_back( mesh->add_entity(entity_property()));
    mesh->change_entity_parts(keys.back(), &part1, &part1+1);
  }

  for ( size_t i = 0; i<10; ++i) {
    keys.push_back( mesh->add_entity(entity_property()));
    mesh->change_entity_parts(keys.back(), &part2, &part2+1);
  }

  //select with 1 part (should result in 1 bucket):
  selector select_part1 = part1;

  BOOST_FOREACH( bucket_key bucket, mesh->get_buckets(part1)) {
    EXPECT_TRUE( is_selected( bucket, select_part1, *mesh));
  }




  modifiable_mesh::selected_bucket_range selected_buckets = get_selected_buckets(select_part1, *mesh);

  EXPECT_FALSE(selected_buckets.first == selected_buckets.second);

  //incrementing the first once, should make it equal to the end.
  ++selected_buckets.first;
  EXPECT_TRUE(selected_buckets.first == selected_buckets.second);

  //select with two parts (should result in 2 buckets):
  selector select = part1 | part2;

  selected_buckets = get_selected_buckets(select, *mesh);

  size_t num_selected_buckets = std::distance(selected_buckets.first, selected_buckets.second);
  EXPECT_EQ(size_t(2), num_selected_buckets);
};
#endif

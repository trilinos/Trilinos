#ifndef __IBMCPP__

#include <gtest/gtest.h>

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>

#include <sierra/mesh/details/constant_size_field.hpp>

using namespace sierra::mesh;
using namespace sierra::mesh::details;

TEST( constant_size_field , basic )
{
  modifiable_mesh mesh;

  part_key part1 = mesh.declare_part("part 1");
  part_key part2 = mesh.declare_part("part 2");

  std::vector<entity_key> keys;
  for ( size_t i = 0; i<10; ++i) {
    keys.push_back( mesh.add_entity(entity_property(mesh.element_rank())));
    mesh.change_entity_parts(keys.back(), &part1, &part1+1);
  }

  for ( size_t i = 0; i<10; ++i) {
    keys.push_back( mesh.add_entity(entity_property(mesh.node_rank())));
    mesh.change_entity_parts(keys.back(), &part2, &part2+1);
  }

  selector select = part1 | part2;

  constant_size_field<double,3> coord_field;

  coord_field.update_from_mesh(select, mesh);

  size_t total_num_field_scalars = coord_field.field_data_flat_array().size();
  EXPECT_EQ(total_num_field_scalars, keys.size()*coord_field.field_length_per_entity);
};

#endif

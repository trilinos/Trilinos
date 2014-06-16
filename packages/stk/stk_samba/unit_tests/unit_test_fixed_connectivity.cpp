#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include <samba/mesh.hpp>



TEST(samba, fixed_node_connectivity)
{
  using namespace samba;

  typedef detail::partition_connectivity<entity_rank::node_type, connectivity_kind::fixed_type> fixed_node_type;

  fixed_node_type tet_2_nodes(entity_topology::tet_4(),spatial_dimension::create(3));

  uint32_t how_many = 16;
  const uint32_t num_node_connectivity = num_nodes(entity_topology::tet_4());

  tet_2_nodes.add_entities(how_many);

  //make sure space has been allocated
  ASSERT_EQ(num_node_connectivity*how_many, tet_2_nodes.end<entity_key>(partition_offset::create(how_many-1)) - tet_2_nodes.begin<entity_key>(partition_offset::create(0)));


  //assign initial connectivities
  {
    entity_local_id id = entity_local_id::create(0);

    for (partition_offset i = partition_offset::create(0); i < how_many; ++i) {
      for (connectivity_ordinal ord = connectivity_ordinal::create(0); ord < num_node_connectivity; ++ord) {
        entity_key to = entity_key::create(entity_topology::node(), process_id::invalid(), id++);
        tet_2_nodes.add_connectivity(i,to,ord,connectivity_orientation::invalid());
      }
    }
  }

  {
    entity_local_id id = entity_local_id::create(0);

    for (partition_offset i = partition_offset::create(0); i < how_many; ++i) {
      for (const entity_key * b = tet_2_nodes.begin<entity_key>(i), * e = tet_2_nodes.end<entity_key>(i); b<e; ++b) {
        entity_key to = entity_key::create(entity_topology::node(), process_id::invalid(), id++);
        EXPECT_EQ(*b,to);
      }
    }
  }

  //reverse connectivities by swapping first with last

  {
    for (partition_offset i = partition_offset::create(0), e = partition_offset::create(how_many-1); i < e; ++i, --e) {
      tet_2_nodes.swap(i,e);
    }
  }

  {
    for (partition_offset i = partition_offset::create(0); i < how_many; ++i) {

      entity_local_id id = entity_local_id::create((how_many-i()-1)*num_node_connectivity);

      for (const entity_key * b = tet_2_nodes.begin<entity_key>(i), * e = tet_2_nodes.end<entity_key>(i); b<e; ++b) {
        entity_key to = entity_key::create(entity_topology::node(), process_id::invalid(), id++);
        EXPECT_EQ(*b,to);
      }
    }
  }

  //move entities to new partition

  fixed_node_type tet_2_nodes_2(entity_topology::tet_4(),spatial_dimension::create(3));

  tet_2_nodes.move_entities(tet_2_nodes_2, how_many);

  {
    for (partition_offset i = partition_offset::create(0); i < how_many; ++i) {

      entity_local_id id = entity_local_id::create((how_many-i()-1)*num_node_connectivity);

      for (const entity_key * b = tet_2_nodes_2.begin<entity_key>(i), * e = tet_2_nodes_2.end<entity_key>(i); b<e; ++b) {
        entity_key to = entity_key::create(entity_topology::node(), process_id::invalid(), id++);
        EXPECT_EQ(*b,to);
      }
    }
  }

  //move entities back to original partition

  tet_2_nodes_2.move_entities(tet_2_nodes, how_many);

  {
    for (partition_offset i = partition_offset::create(0); i < how_many; ++i) {

      entity_local_id id = entity_local_id::create((how_many-i()-1)*num_node_connectivity);

      for (const entity_key * b = tet_2_nodes.begin<entity_key>(i), * e = tet_2_nodes.end<entity_key>(i); b<e; ++b) {
        entity_key to = entity_key::create(entity_topology::node(), process_id::invalid(), id++);
        EXPECT_EQ(*b,to);
      }
    }
  }

}

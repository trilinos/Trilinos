#include <gtest/gtest.h>

#include <samba/mesh.hpp>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <sstream>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/equal.hpp>

TEST(samba, mesh_impl_basic)
{
  using namespace samba;

  samba::mesh mesh;

  entity_block_key block_0 = mesh.add_entity_block("Block_0");
  entity_block_key block_1 = mesh.add_entity_block("Block_1");

  {
    entity_key_interval keys = mesh.add_entities( entity_topology::node()
                                                  ,2
                                                  ,&block_0
                                                  ,&block_0+1
                                                  );

    partition_proxy partition = mesh[keys[0]].partition();

    entity_part_vector parts;
    parts.push_back(entity_rank::node());
    parts.push_back(entity_topology::node());
    parts.push_back(entity_state::universe());
    parts.push_back(entity_state::owned());
    parts.push_back(block_0);

    EXPECT_TRUE(boost::equal(parts,partition.parts()));

    EXPECT_EQ(partition.size(),2U);

    for (partition_offset i={0},e=partition_offset::create(partition.size()); i<e; ++i) {
      EXPECT_EQ(keys[i()],partition[i].key());
    }

    //remove the first 1 keys
    mesh.remove_entities(keys.begin(),keys.begin()+1);

    EXPECT_EQ(partition.size(),1U);

    //move the remaining entities to block_1
    mesh.move_entities( keys.begin()+1,keys.end()
                       ,&block_1, &block_1+1  //add to block_1
                       ,&block_0, &block_0+1  //remove from block_0
                      );

    partition = mesh[partition_id::create(0)];

    EXPECT_EQ(partition.size(),0U);

    // get the new partition for the entities
    partition = mesh[keys[1]].partition();

    EXPECT_EQ(partition.size(),1U);

    EXPECT_EQ(entity_rank::invalid(), mesh.get_rank(block_0));
    EXPECT_EQ(entity_rank::invalid(), mesh.get_rank(block_1));
  }

  mesh.end_modification();
}

namespace {

const uint32_t NX = 10;
const uint32_t NY = 10;
const uint32_t NZ = 10;

inline size_t node_offset(uint32_t x, uint32_t y, uint32_t z)
{ return (x + ( NX + 1 ) * ( y + ( NY + 1 ) * z ));  }

inline size_t hex_offset(uint32_t x, uint32_t y, uint32_t z)
{ return (x + NX * ( y +  NY * z ));  }

} //unnamed namespace

TEST(samba, mesh_hex_fixture)
{

  const uint32_t num_nodes = (NX+1)*(NY+1)*(NZ+1);
  const uint32_t num_hexes = NX*NY*NZ;

  samba::mesh mesh;

  samba::entity_key_interval nodes = mesh.add_entities( samba::entity_topology::node(),  num_nodes);
  samba::entity_key_interval hexes = mesh.add_entities( samba::entity_topology::hex_8(), num_hexes);

  for (uint32_t x=0; x<NX; ++x) {
  for (uint32_t y=0; y<NY; ++y) {
  for (uint32_t z=0; z<NZ; ++z) {

    samba::connectivity_ordinal ordinal = {0};

    samba::entity_key hex = hexes[hex_offset(x,y,z)];
    samba::entity_key node;

    node = nodes[node_offset( x   , y   , z   )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

    node = nodes[node_offset( x+1 , y   , z   )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

    node = nodes[node_offset( x+1 , y   , z+1 )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

    node = nodes[node_offset( x   , y   , z+1 )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

    node = nodes[node_offset( x   , y+1 , z   )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

    node = nodes[node_offset( x+1 , y+1 , z   )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

    node = nodes[node_offset( x+1 , y+1 , z+1 )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

    node = nodes[node_offset( x   , y+1 , z+1 )];
    mesh.add_connectivity(hex,node,ordinal);
    mesh.add_connectivity(node,hex,ordinal);
    ++ordinal;

  }
  }
  }

  for (uint32_t x=0; x<NX; ++x) {
  for (uint32_t y=0; y<NY; ++y) {
  for (uint32_t z=0; z<NZ; ++z) {

    samba::connectivity_ordinal ordinal = {0};

    samba::entity_key hex = hexes[hex_offset(x,y,z)];
    samba::entity_key node;

    samba::entity_key_iterator  node_key_connectivity = mesh[hex].begin_nodes<samba::entity_key>();
    samba::ordinal_iterator node_ordinal_connectivity = mesh[hex].begin_nodes<samba::connectivity_ordinal>();

    node = nodes[node_offset( x   , y   , z   )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    node = nodes[node_offset( x+1 , y   , z   )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    node = nodes[node_offset( x+1 , y   , z+1 )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    node = nodes[node_offset( x   , y   , z+1 )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    node = nodes[node_offset( x   , y+1 , z   )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    node = nodes[node_offset( x+1 , y+1 , z   )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    node = nodes[node_offset( x+1 , y+1 , z+1 )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    node = nodes[node_offset( x   , y+1 , z+1 )];
    EXPECT_EQ( *node_key_connectivity, node);
    EXPECT_EQ( *node_ordinal_connectivity, ordinal);
    ++ordinal;
    ++node_key_connectivity;
    ++node_ordinal_connectivity;

    EXPECT_TRUE(node_key_connectivity == mesh[hex].end_nodes<samba::entity_key>());
    EXPECT_TRUE(node_ordinal_connectivity == mesh[hex].end_nodes<samba::connectivity_ordinal>());
  }
  }
  }

  // Check back-relations
  for (uint32_t x=0; x<=NX; ++x) {
    for (uint32_t y=0; y<=NY; ++y) {
      for (uint32_t z=0; z<=NZ; ++z) {
         const size_t num_dim_not_on_bounary = (x != 0 && x != NX) + (y != 0 && y != NY) + (z != 0 && z != NZ);
         const int expected_num_back_rels = 1 << num_dim_not_on_bounary; // 2^num
         samba::entity_key node = nodes[node_offset( x, y, z)];

         EXPECT_EQ( expected_num_back_rels, boost::size(mesh[node].elements<samba::entity_key>()) );
      }
    }
  }

  mesh.end_modification();
}

TEST(samba, mesh_impl_modification)
{
  //
  // We set up the following scenario:
  //        chunk1     chunk2     chunk3     chunk4     chunk5     chunk6     chunk7
  //      <--------- ---------- -block 0 - --------->
  //                 <--------- ---------- -block 1 - --------->
  //                            <--------- ---------- -block 2 - --------->
  //                                        <--------- ---------- -block 3 - --------->
  //
  // Each "chunk" of entities should get it's own partition.
  //
  // We then remove entities out of their blocks, 1 block per modification cycle.
  //

  using namespace samba;

  mesh fixture;

  // New mesh should start as modifiable
  EXPECT_TRUE(fixture.modifiable());

  // Try to begin modifcation while still in a modification cycle. This
  // is legal, but should return false.
  EXPECT_FALSE(fixture.begin_modification());

  const size_t entity_chunk_size = 2;
  const size_t num_blocks        = 4;
  const size_t entity_block_size = num_blocks * entity_chunk_size;
  const size_t num_entities      = entity_block_size + entity_chunk_size * (num_blocks - 1);

  // Add blocks
  std::vector<entity_block_key> entity_block_keys;
  for (size_t i = 0; i < num_blocks; ++i) {
    std::ostringstream oss;
    oss << "Block_" << i;
    entity_block_keys.push_back(fixture.add_entity_block(oss.str()));
  }

  // Add all entities
  entity_key_interval all_keys = fixture.add_entities( entity_topology::node(), num_entities );

  // Add entities to appropriate entity_blocks
  for (size_t block_id = 0; block_id < num_blocks; ++block_id) {
    fixture.move_entities(all_keys.begin() + block_id * entity_chunk_size,
                          all_keys.begin() + block_id * entity_chunk_size + entity_block_size,
                          entity_block_keys.begin() + block_id,
                          entity_block_keys.begin() + block_id + 1,
                          entity_block_keys.end(),
                          entity_block_keys.end());
  }

  EXPECT_TRUE(fixture.end_modification());

  size_t expected_num_partitions = num_entities / entity_chunk_size;

  ASSERT_EQ(expected_num_partitions, fixture.num_partitions());

  // Sanity check on partition size and back-and-forth conversions between
  // partition_indexs and entity_keys is working.
  for (size_t partition = 0; partition < expected_num_partitions; ++partition) {
    partition_id pd = partition_id::create(partition);
    partition_proxy      pp = fixture[pd];
    EXPECT_EQ(entity_chunk_size, pp.size());
    for (detail::partition_iterator pi = pp.begin(), end = pp.end(); pi != end; ++pi) {
      partition_index d = pi->descriptor();
      entity_key        k = fixture.convert(d);
      EXPECT_EQ(d, fixture.convert(k));
    }
  }

  // We remove an entity_block, one at a time, and then do a similar sanity check
  // as was done above.
  for (size_t to_remove = 0; to_remove < num_blocks; ++to_remove) {
    EXPECT_TRUE(fixture.begin_modification());

    fixture.move_entities(all_keys.begin() + to_remove * entity_chunk_size,
                          all_keys.begin() + to_remove * entity_chunk_size + entity_block_size,
                          entity_block_keys.end(),
                          entity_block_keys.end(),
                          entity_block_keys.begin() + to_remove,
                          entity_block_keys.begin() + to_remove + 1);

    EXPECT_TRUE(fixture.end_modification());

    // In general, we lose 2 partitions per entity_block removal, except for the first and last
    // removals, which only incur the loss of one partition.
    expected_num_partitions = expected_num_partitions - ( (to_remove == 0 || to_remove == num_blocks-1) ? 1 : 2);

    ASSERT_EQ(expected_num_partitions, fixture.num_partitions());

    for (size_t partition = 0; partition < expected_num_partitions; ++partition) {
      partition_id pd = partition_id::create(partition);
      partition_proxy      pp = fixture[pd];
      if (expected_num_partitions == 1) {
        EXPECT_EQ(num_entities, pp.size());
      }
      else {
        if (partition != 0) {
          if (partition == (num_blocks - (to_remove+1))) {
            EXPECT_EQ(entity_chunk_size * (to_remove+2), pp.size());
          }
          else {
            EXPECT_EQ(entity_chunk_size, pp.size());
          }
        }
        else {
          EXPECT_EQ( (to_remove+1) * entity_chunk_size, pp.size() );
        }
      }
      for (detail::partition_iterator pi = pp.begin(), end = pp.end(); pi != end; ++pi) {
        partition_index d = pi->descriptor();
        entity_key        k = fixture.convert(d);
        EXPECT_EQ(d, fixture.convert(k));
      }
    }
  }

  // Try to end modification when there is no modification cycle to end.
  // This is legal, but should return false.
  EXPECT_FALSE(fixture.end_modification());
}

TEST(samba, induced_parts)
{
  //
  // Set up a very simple mesh with one element, two sides, one node:
  // rels: E->(S1,S2)->N
  // initial parts: E in unranked_not_inducable, unranked_inducable, element_inducable, element_non_inducable
  //                S1 in side_inducable, side_non_inducable
  //                S2 in unranked_inducable, side_non_inducable
  //                N is not initially in any parts
  //

  using namespace samba;

  mesh fixture;

  entity_block_key unranked_not_inducable = fixture.add_entity_block("unranked-not-inducable");
  entity_block_key unranked_inducable     = fixture.add_entity_block("unranked-inducable", entity_rank::invalid(), true /*inducable*/);
  entity_block_key element_inducable      = fixture.add_entity_block("element-inducable", entity_rank::element(), true /*inducable*/);
  entity_block_key element_non_inducable  = fixture.add_entity_block("element-not-inducable", entity_rank::element());
  entity_block_key side_inducable         = fixture.add_entity_block("side-inducable", entity_rank::face(), true /*inducable*/);
  entity_block_key side_non_inducable     = fixture.add_entity_block("side-non-inducable", entity_rank::face());

  std::vector<entity_block_key> elem_parts, side1_parts, side2_parts;

  elem_parts.push_back(unranked_not_inducable);
  elem_parts.push_back(unranked_inducable);
  elem_parts.push_back(element_inducable);
  elem_parts.push_back(element_non_inducable);

  side1_parts.push_back(side_inducable);
  side1_parts.push_back(side_non_inducable);

  side2_parts.push_back(unranked_inducable);
  side2_parts.push_back(side_non_inducable);

  entity_key_interval elems = fixture.add_entities(entity_topology::hex_8(),
                                                   1 /*num to add*/,
                                                   elem_parts.begin(),
                                                   elem_parts.end());
  EXPECT_EQ(1u, elems.size());
  entity_key elem = elems[0];

  entity_key_interval sides1 = fixture.add_entities(entity_topology::quad_4(),
                                                   1 /*num to add*/,
                                                   side1_parts.begin(),
                                                   side1_parts.end());
  EXPECT_EQ(1u, sides1.size());
  entity_key side1 = sides1[0];

  entity_key_interval sides2 = fixture.add_entities(entity_topology::quad_4(),
                                                    1 /*num to add*/,
                                                    side2_parts.begin(),
                                                    side2_parts.end());
  EXPECT_EQ(1u, sides2.size());
  entity_key side2 = sides2[0];

  entity_key_interval nodes = fixture.add_entities(entity_topology::node(),
                                                   1 /*num to add*/);
  EXPECT_EQ(1u, nodes.size());
  entity_key node = nodes[0];

  EXPECT_TRUE(fixture.end_modification());

  // We have not specified any connectivity yet, so the entities' parts should
  // match the initial parts

  {
    entity_part_range parts = fixture[elem].parts();
    entity_part expected_parts[] = {entity_rank::element(),
                                    entity_topology::hex_8(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    unranked_not_inducable,
                                    unranked_inducable,
                                    element_inducable,
                                    element_non_inducable};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[side1].parts();
    entity_part expected_parts[] = {entity_rank::face(),
                                    entity_topology::quad_4(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    side_inducable,
                                    side_non_inducable};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[side2].parts();
    entity_part expected_parts[] = {entity_rank::face(),
                                    entity_topology::quad_4(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    unranked_inducable,
                                    side_non_inducable};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[node].parts();
    entity_part expected_parts[] = {entity_rank::node(),
                                    entity_topology::node(),
                                    entity_state::universe(),
                                    entity_state::owned()};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  // Add connectivity

  EXPECT_TRUE(fixture.begin_modification());

  fixture.add_connectivity(elem, side1, connectivity_ordinal::create(0));
  fixture.add_connectivity(elem, side2, connectivity_ordinal::create(1));
  fixture.add_connectivity(side1, node, connectivity_ordinal::create(0));
  fixture.add_connectivity(side2, node, connectivity_ordinal::create(0));

  EXPECT_TRUE(fixture.end_modification());

  // Check that correct parts were induced

  {
    entity_part_range parts = fixture[elem].parts();
    entity_part expected_parts[] = {entity_rank::element(),
                                    entity_topology::hex_8(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    unranked_not_inducable,
                                    unranked_inducable,
                                    element_inducable,
                                    element_non_inducable};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[side1].parts();
    entity_part expected_parts[] = {entity_rank::face(),
                                    entity_topology::quad_4(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    side_inducable,
                                    side_non_inducable,
                                    entity_part(entity_topology::hex_8(), true),
                                    entity_part(unranked_inducable, true),
                                    entity_part(element_inducable, true)};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[side2].parts();
    entity_part expected_parts[] = {entity_rank::face(),
                                    entity_topology::quad_4(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    unranked_inducable,
                                    side_non_inducable,
                                    entity_part(entity_topology::hex_8(), true),
                                    entity_part(unranked_inducable, true),
                                    entity_part(element_inducable, true)};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[node].parts();
    entity_part expected_parts[] = {entity_rank::node(),
                                    entity_topology::node(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    entity_part(entity_topology::quad_4(), true),
                                    entity_part(unranked_inducable, true),
                                    entity_part(side_inducable, true)};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  // Remove some connectivity

  EXPECT_TRUE(fixture.begin_modification());

  fixture.remove_connectivity(elem, side1, connectivity_ordinal::create(0));
  fixture.remove_connectivity(side1, node, connectivity_ordinal::create(0));
  fixture.remove_connectivity(side2, node, connectivity_ordinal::create(0));

  EXPECT_TRUE(fixture.end_modification());

  // Check that correct parts were removed

  {
    entity_part_range parts = fixture[elem].parts();
    entity_part expected_parts[] = {entity_rank::element(),
                                    entity_topology::hex_8(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    unranked_not_inducable,
                                    unranked_inducable,
                                    element_inducable,
                                    element_non_inducable};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[side1].parts();
    entity_part expected_parts[] = {entity_rank::face(),
                                    entity_topology::quad_4(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    side_inducable,
                                    side_non_inducable};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[side2].parts();
    entity_part expected_parts[] = {entity_rank::face(),
                                    entity_topology::quad_4(),
                                    entity_state::universe(),
                                    entity_state::owned(),
                                    unranked_inducable,
                                    side_non_inducable,
                                    entity_part(entity_topology::hex_8(), true),
                                    entity_part(unranked_inducable, true),
                                    entity_part(element_inducable, true)};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }

  {
    entity_part_range parts = fixture[node].parts();
    entity_part expected_parts[] = {entity_rank::node(),
                                    entity_topology::node(),
                                    entity_state::universe(),
                                    entity_state::owned()};
    EXPECT_EQ(static_cast<size_t>(parts.size()), sizeof(expected_parts) / sizeof(entity_part));
    for (int i = 0; i < parts.size(); ++i) {
      EXPECT_EQ(expected_parts[i], parts[i]);
    }
  }
}

TEST(samba, connectivity_to_deleted_entity)
{
  //
  // Set up a very simple mesh with one element, one node. Establish and relation/back-relation
  // between the two entities. Then remove the element; after end-modification, the node should
  // have no back relation.
  //

  using namespace samba;

  mesh fixture;

  entity_key_interval elems = fixture.add_entities(entity_topology::hex_8(), 1 /*just add 1*/);
  entity_key_interval nodes = fixture.add_entities(entity_topology::node() , 1 /*just add 1*/);

  EXPECT_EQ(1u, elems.size());
  EXPECT_EQ(1u, nodes.size());

  entity_key elem = elems[0];
  entity_key node = nodes[0];

  connectivity_ordinal ord = {0};
  fixture.add_connectivity(elem, node, ord);
  fixture.add_connectivity(node, elem, ord); //back rel

  // Check that node now had back relation
  EXPECT_EQ(1u, fixture[node].num_elements());

  fixture.remove_entities(&elem, &elem + 1);

  // Node should still have back relation since end_modification hasn't cleaned it up yet
  EXPECT_EQ(1u, fixture[node].num_elements());

  EXPECT_TRUE(fixture.end_modification());

  // Connectivity to removed element should be gone

  EXPECT_EQ(0u, fixture[node].num_elements());
}

TEST(samba, connectivity_to_deleted_entity_key_recycling)
{
  //
  // Set up a very simple mesh with one element, one node. Establish and relation/back-relation
  // between the two entities. Then remove the element and add an element with the same key.
  // Then add an element. The node should no longer have any relations.
  //

  using namespace samba;

  mesh fixture;

  entity_key_interval elems = fixture.add_entities(entity_topology::hex_8(), 1 /*just add 1*/);
  entity_key_interval nodes = fixture.add_entities(entity_topology::node() , 1 /*just add 1*/);

  EXPECT_EQ(1u, elems.size());
  EXPECT_EQ(1u, nodes.size());

  entity_key elem = elems[0];
  entity_key node = nodes[0];

  connectivity_ordinal ord = {0};
  fixture.add_connectivity(elem, node, ord);
  fixture.add_connectivity(node, elem, ord); //back rel

  // Check that node now had back relation
  EXPECT_EQ(1u, fixture[node].num_elements());

  fixture.remove_entities(&elem, &elem + 1);

  elems = fixture.add_entities(entity_topology::hex_8(), 1 /*just add 1*/);

  EXPECT_TRUE(fixture.end_modification());

  // Connectivity to removed element should be gone

  EXPECT_EQ(0u, fixture[node].num_elements());

  EXPECT_EQ(1u, fixture.num_elements());
}

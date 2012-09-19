#ifndef __IBMCPP__
#include <sierra/mesh/modifiable/modifiable_mesh.hpp>
#include <sierra/mesh/details/cell_topology.hpp>

#include <boost/range.hpp>

#include <gtest/gtest.h>

using namespace sierra::mesh;
using namespace sierra::mesh::details;

TEST( modifiable_mesh , basic )
{
  modifiable_mesh mesh;

  for ( size_t i = 0u; i<20u; ++i) {
    mesh.add_entity();
  }

  EXPECT_EQ( 20u, mesh.num_entities() );
  EXPECT_EQ(std::distance(mesh.get_entities().first, mesh.get_entities().second), static_cast<int>(mesh.num_entities()) );

  part_key odd = mesh.declare_part("odd");
  part_key even = mesh.declare_part("even");

  BOOST_FOREACH( entity_key key, mesh.get_entities() ) {
    if (static_cast<size_t>(key)&1) { //odd
      mesh.change_entity_parts(key, &odd, &odd+1);
    }
    else {
      mesh.change_entity_parts(key, &even, &even+1);
    }
  }

  BOOST_FOREACH( bucket_key b, mesh.get_buckets(odd)) {
    BOOST_FOREACH( entity_key key, mesh.get_entities(b)) {
      EXPECT_TRUE( 1 == key%2);
    }
  }
  BOOST_FOREACH( bucket_key b, mesh.get_buckets(even)) {
    BOOST_FOREACH( entity_key key, mesh.get_entities(b)) {
      EXPECT_TRUE( 0 == key%2);
    }
  }

  //  We expect all entities to be in the universal part:
  size_t num_entities_in_universal_part = 0;
  BOOST_FOREACH(bucket_key b, mesh.get_buckets(mesh.universal_part())) {
    BOOST_FOREACH(entity_key key, mesh.get_entities(b)) {
      (void)key; // avoid an unused but set variable compiler warning.
      ++num_entities_in_universal_part;
    }
  }

  EXPECT_EQ(num_entities_in_universal_part, mesh.num_entities());

  //For this serial test, we expect all entities to be in the locally-owned part:
  size_t num_entities_in_locally_owned_part = 0;
  BOOST_FOREACH(bucket_key b, mesh.get_buckets(mesh.locally_owned_part())) {
    BOOST_FOREACH(entity_key key, mesh.get_entities(b)) {
      (void)key; // avoid an unused but set variable compiler warning.
      ++num_entities_in_locally_owned_part;
    }
  }

  EXPECT_EQ(num_entities_in_locally_owned_part, mesh.num_entities());
}

TEST( modifiable_mesh , relation_order )
{
  modifiable_mesh mesh;

  for ( size_t i = 0u; i<9u; ++i) {
    mesh.add_entity();
  }

  EXPECT_EQ( 9u, mesh.num_entities() );
  EXPECT_EQ(std::distance(mesh.get_entities().first, mesh.get_entities().second), static_cast<int>(mesh.num_entities()) );

  part_key elem = mesh.declare_part("elem");
  part_key node = mesh.declare_part("node");

  entity_key elem_key;
  std::vector<entity_key> node_keys;
  BOOST_FOREACH( entity_key key, mesh.get_entities() ) {
    if (static_cast<size_t>(key) == 0) { //elem
      elem_key = key;
      mesh.change_entity_parts(key, &elem, &elem+1);
    }
    else {
      mesh.change_entity_parts(key, &node, &node+1);
      //store the node-keys in reverse order (highest key at the front of the vector)
      node_keys.insert(node_keys.begin(), key);
    }
  }

  //add element-to-node relations in reverse order, highest key first
  for(size_t i=0; i<node_keys.size(); ++i) {
    mesh.add_relation(elem_key, node_keys[i], relation_position(0, i));
  }

  out_relation_range node_range = mesh.get_out_relations(elem_key);

  //now verify that the node relations for the element are in the same order
  //that we added them:
  size_t i=0;
  for(out_relation_iterator iter = boost::begin(node_range), end = boost::end(node_range); iter!=end; ++iter, ++i) {
    entity_key node_key = mesh.target_entity(*iter);
    EXPECT_EQ(node_key, node_keys[i]);
  }
}

TEST ( modifiable_mesh, parts )
{
  modifiable_mesh mesh;

  const entity_rank element_rank(3);
  const entity_rank node_rank(0);

  part_key element_rank_part = mesh.declare_part("element_rank_part",element_rank);
  part_key node_rank_part = mesh.declare_part("node_rank_part",node_rank);

  const cell_topology hex_top = shards::getCellTopologyData<shards::Hexahedron<8> >();
  const cell_topology wedge_top = shards::getCellTopologyData<shards::Wedge<6> >();
  const cell_topology node_top = shards::getCellTopologyData<shards::Node >();

  part_key hex_top_part = mesh.declare_part("hex_top_part",hex_top);
  part_key node_top_part = mesh.declare_part("node_top_part",node_top);

  { // test get_rank
    EXPECT_TRUE( mesh[element_rank_part].has_property<entity_rank>());
    EXPECT_EQ( mesh[element_rank_part].get_property<entity_rank>(), element_rank);
  }

  //add element-rank to topology part:
  part_property& part_prop = mesh.get(hex_top_part);
  part_prop.add_property<entity_rank>(element_rank);

  { // test get topology part property
    EXPECT_TRUE( mesh[hex_top_part].has_property<cell_topology>());
    EXPECT_EQ( mesh[hex_top_part].get_property<cell_topology>(), hex_top);
  }

  {//test getting an entity-rank from the topology part:
    EXPECT_TRUE( mesh[hex_top_part].has_property<entity_rank>());
    EXPECT_EQ( mesh[hex_top_part].get_property<entity_rank>(), element_rank);
  }

  { // add an element to the mesh
    std::vector<part_key> parts;
    parts.push_back(element_rank_part);
    parts.push_back(hex_top_part);

    entity_key elem = mesh.add_entity(entity_property(element_rank));
    mesh.change_entity_parts(elem, parts.begin(), parts.end());

    std::vector<part_key> node_parts;
    node_parts.push_back(node_rank_part);
    node_parts.push_back(node_top_part);
    node_parts.push_back(hex_top_part);
    for (size_t i=0; i<8u; ++i) {
      entity_key node = mesh.add_entity(entity_property(node_rank));
      mesh.add_relation(elem,node, relation_position(node_rank,i));
      mesh.change_entity_parts(node,node_parts.begin(),node_parts.end());
    }
  }

  std::cout << "Element Rank Part: " << element_rank_part.base() << std::endl;
  std::cout << "Node Rank Part:    " << node_rank_part.base()  << std::endl;
  std::cout << "Hex Top Part:      " << hex_top_part.base() << std::endl;
  std::cout << "Node Top Part:     " << node_top_part.base() << std::endl;
  std::cout << std::endl;

  std::cout << "Num buckets: " << mesh.num_buckets() << std::endl;

//  BOOST_FOREACH( bucket_key bucket, mesh.get_buckets() ) {
//    std::cout << bucket << std::endl;
//    std::cout << "\tParts: ";
//    BOOST_FOREACH( part_key part, mesh.get_parts(bucket)) {
//      std::cout << part.base() << " ";
//    }
//    std::cout << std::endl;
//    std::cout << "\tNum Entities: " << mesh.num_entities(bucket) << std::endl;
//    BOOST_FOREACH( entity_key entity, mesh.get_entities(bucket)) {
//      std::cout << "\t" << mesh[entity].m_rank << std::endl;
//      std::cout << "\t   Num Out relations: " << mesh.num_out_relations(entity) << std::endl;
//      std::cout << "\t   Num In relations: " << mesh.num_in_relations(entity) << std::endl;
//      std::cout << std::endl;
//    }
//  }
}
#endif

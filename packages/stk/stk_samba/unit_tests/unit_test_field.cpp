#include <gtest/gtest.h>

#include <iostream>

#include <samba/mesh.hpp>
#include <samba/field.hpp>
#include <samba/rank_field.hpp>

TEST(samba, scalar_field_basic)
{
  typedef samba::field<double> field_type;
  typedef samba::field<const double> const_field_type;

  samba::mesh mesh;

  samba::entity_block_key block_1 = mesh.add_entity_block("block_1");

  //set the const field equal to the field
  const double default_value = 3;
  field_type field(mesh, block_1,default_value);
  const_field_type const_field = field;

  const size_t num_add = 100;

  samba::entity_key_interval keys = mesh.add_entities(samba::entity_topology::hex_8(),num_add,&block_1,&block_1+1);

  EXPECT_EQ(1u, mesh.num_partitions());

  //get the partition for the added keys
  samba::partition_index entity_desc = mesh.convert(keys[0]);
  samba::partition_id partition = entity_desc.partition();

  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    EXPECT_EQ(field[partition][offset],default_value);
  }

  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    field[partition][offset] = 0.0;
  }

  //remove the entities...this will NOT reset the data
  mesh.remove_entities(keys.begin(),keys.end());

  keys = mesh.add_entities(samba::entity_topology::hex_8(),num_add,&block_1,&block_1+1);

  //data should be reset to the default value
  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    EXPECT_EQ(field[partition][offset],default_value);
  }

  mesh.end_modification();
}

TEST(samba, scalar_rank_field_basic)
{
  typedef samba::rank_field<samba::entity_rank::element_type, double> field_type;
  typedef samba::rank_field<samba::entity_rank::element_type, const double> const_field_type;

  samba::mesh mesh;

  samba::entity_block_key block_1 = mesh.add_entity_block("block_1");

  //set the const field equal to the field
  const double default_value = 3.0;
  field_type field(mesh, default_value);
  const_field_type const_field = field;

  const size_t num_add = 100;

  samba::entity_key_interval keys = mesh.add_entities(samba::entity_topology::hex_8(),num_add,&block_1,&block_1+1);

  EXPECT_EQ(1u, mesh.num_partitions());

  //get the partition for the added keys
  samba::partition_index entity_desc = mesh.convert(keys[0]);
  samba::partition_id partition = entity_desc.partition();

  // access by partition_index
  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    samba::partition_index pi = samba::partition_index::create(samba::entity_rank::element(), partition, offset);
    field[pi] = static_cast<double>(offset());
  }

  mesh.end_modification();

  for (size_t i = 0; i < num_add; ++i) {
    EXPECT_EQ(static_cast<double>(i), field[i]);
  }

  mesh.begin_modification();

  //remove the entities...this will NOT reset the data
  mesh.remove_entities(keys.begin(),keys.end());

  keys = mesh.add_entities(samba::entity_topology::hex_8(),num_add,&block_1,&block_1+1);

  // access by partition_index
  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    samba::partition_index pi = samba::partition_index::create(samba::entity_rank::element(), partition, offset);
    EXPECT_EQ(default_value, field[pi]);
  }

  mesh.end_modification();

  //data should be reset to the default value; should be indexable by raw int
  for (size_t i = 0; i < num_add; ++i) {
    EXPECT_EQ(default_value, field[i]);
  }
}


TEST(samba, field_move_default)
{
  typedef samba::field<double, samba::num_nodes_functor> field_type;
  typedef samba::field<const double, samba::num_nodes_functor> const_field_type;

  samba::mesh mesh;

  samba::entity_block_key block_1 = mesh.add_entity_block("block_1");

  //set the const field equal to the field
  const double default_value = 3.0;
  field_type field(mesh, block_1,default_value);
  const_field_type const_field = field;

  const size_t num_add = 100;

  samba::entity_key_interval keys = mesh.add_entities(samba::entity_topology::hex_8(),num_add,&block_1,&block_1+1);

  EXPECT_EQ(1u, mesh.num_partitions());

  //get the partition for the added keys
  samba::partition_index entity_desc = mesh.convert(keys[0]);
  samba::partition_id partition = entity_desc.partition();

  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    for (size_t i = 0, ee = samba::num_nodes(mesh[partition].topology()); i<ee; ++i) {
      EXPECT_EQ(field[partition][offset][i],default_value);
      field[partition][offset][i] = 42.0; // any value other than default
    }
  }

  samba::entity_block_key none = samba::entity_block_key::invalid();

  //remove the entities from block_1
  mesh.move_entities(keys.begin(),keys.end(), &none, &none, &block_1, &block_1+1);

  EXPECT_EQ(2u, mesh.num_partitions());


  //move entities back to block_1
  mesh.move_entities(keys.begin(),keys.end(), &block_1, &block_1+1, &none, &none);

  EXPECT_EQ(2u, mesh.num_partitions());

  //should be back in original partition
  //NOTE: this is not a safe thing to do in general,
  //but partition_ids can only change when
  //end_modification is signaled
  EXPECT_EQ(num_add, mesh[partition].size());

  //check that the field has reset the values to the default
  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    for (size_t i = 0, ee = samba::num_nodes(mesh[partition].topology()); i<ee; ++i) {
      EXPECT_EQ(field[partition][offset][i],default_value);
    }
  }

  mesh.end_modification();
}

TEST(samba, init_mesh_before_field_during_mod)
{
  typedef samba::field<double> field_type;

  samba::mesh mesh;

  samba::entity_block_key block_1 = mesh.add_entity_block("block_1");

  const size_t num_add = 100;

  samba::entity_key_interval keys = mesh.add_entities(samba::entity_topology::hex_8(),num_add,&block_1,&block_1+1);

  EXPECT_EQ(1u, mesh.num_partitions());

  //set the const field equal to the field
  const double default_value = 3;
  field_type field(mesh, block_1, default_value);

  //get the partition for the added keys
  samba::partition_index entity_desc = mesh.convert(keys[0]);
  samba::partition_id partition = entity_desc.partition();

  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    EXPECT_EQ(field[partition][offset],default_value);
  }

  mesh.end_modification();
}

TEST(samba, init_mesh_before_field_after_mod)
{
  typedef samba::field<double> field_type;

  samba::mesh mesh;

  samba::entity_block_key block_1 = mesh.add_entity_block("block_1");

  const size_t num_add = 100;

  samba::entity_key_interval keys = mesh.add_entities(samba::entity_topology::hex_8(),num_add,&block_1,&block_1+1);

  mesh.end_modification();

  EXPECT_EQ(1u, mesh.num_partitions());

  //set the const field equal to the field
  const double default_value = 3;
  field_type field(mesh, block_1, default_value);

  //get the partition for the added keys
  samba::partition_index entity_desc = mesh.convert(keys[0]);
  samba::partition_id partition = entity_desc.partition();

  for (samba::partition_offset offset ={0}, e={mesh[partition].size()}; offset < e; ++offset) {
    EXPECT_EQ(field[partition][offset],default_value);
  }
}

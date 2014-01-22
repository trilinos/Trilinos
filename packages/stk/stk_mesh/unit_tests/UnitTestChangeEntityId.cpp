/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>

#include <boost/foreach.hpp>

STKUNIT_UNIT_TEST( UnitTestChangeEntityId, change_id_small )
{
  // change_entity_id is a broken concept in parallel:
  //   Only entity owners should have to call change-id. As it stands now,
  //   every process that might know about the entity, even as a ghost,
  //   must make the call. If any sharer/ghoster of the changing entity
  //   forgets to make the call, or changes to wrong id, awful errors ensue
  //   in the next modification_end. The only way to make this safe is to have
  //   change_entity_id be a parallel call similar to change_entity_owner; this
  //   would be very expensive. We need to implement the tracking of exodus-ids
  //   in a different fashion in the long run.

  // Demonstrate how change_entity_id should work

  using namespace stk::mesh;

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned spatial_dim = 2;
  MetaData meta_data(spatial_dim);
  meta_data.commit();

  BulkData mesh(meta_data, pm);
  int p_rank = mesh.parallel_rank();
  int p_numProcs = mesh.parallel_size();
  if (p_numProcs > 1) {
    //change-entity-id test only supported in serial
    return;
  }

  mesh.modification_begin();

  Entity elem                 = mesh.declare_entity(MetaData::ELEMENT_RANK, p_rank + 1 /*id*/);

  Entity node1_local_chg_id   = mesh.declare_entity(MetaData::NODE_RANK, p_rank*3 + 1 /*id*/);
  Entity node2_shared_chg_id  = mesh.declare_entity(MetaData::NODE_RANK,            2 /*id*/);
  Entity node3                = mesh.declare_entity(MetaData::NODE_RANK, p_rank*3 + 3 /*id*/);

  mesh.declare_relation(elem, node1_local_chg_id , 1 /*relation ordinal*/);
  mesh.declare_relation(elem, node2_shared_chg_id, 2 /*relation ordinal*/);
  mesh.declare_relation(elem, node3              , 3 /*relation ordinal*/);

  mesh.modification_end();

  mesh.modification_begin();

  const EntityId new_id_for_local_node = 42 + p_rank + 1;
  const EntityId new_id_for_shared_node = 666;

  mesh.change_entity_id(new_id_for_local_node, node1_local_chg_id);
  if (mesh.parallel_owner_rank(node2_shared_chg_id) == p_rank) {
    mesh.change_entity_id(new_id_for_shared_node, node2_shared_chg_id);
  }

  STKUNIT_EXPECT_EQ(new_id_for_local_node, mesh.identifier(node1_local_chg_id));

  mesh.modification_end();

  STKUNIT_EXPECT_EQ(new_id_for_local_node, mesh.identifier(node1_local_chg_id));
  STKUNIT_EXPECT_EQ(new_id_for_shared_node, mesh.identifier(node2_shared_chg_id));
}

STKUNIT_UNIT_TEST( UnitTestChangeEntityId, change_id_large )
{
  using namespace stk::mesh;

  const unsigned NX = 20;
  const unsigned NY = 20;
  const unsigned NZ = 20;
  const unsigned num_elems = NX * NY * NZ;

  fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);

  Field<int> & simple_nodal_field = hf.m_meta.declare_field<Field<int> >(stk::topology::NODE_RANK, "simple_nodal_field");

  put_field( simple_nodal_field, *hf.m_elem_parts[0]);

  //create nodal field on hex topo

  hf.m_meta.commit();

  hf.generate_mesh();

  stk::mesh::BulkData & mesh = hf.m_bulk_data;
  if (mesh.parallel_size() > 1) {
    return;
  }

  mesh.modification_begin();

  const BucketVector & node_buckets = mesh.buckets(MetaData::NODE_RANK);

  BOOST_FOREACH(Bucket * b, node_buckets) {
    int* nodal_field = mesh.field_data( simple_nodal_field, *b );
    for (size_t i =0; i<b->size(); ++i) {
      nodal_field[i] = 1;
    }
  }

  const BucketVector & elem_buckets = mesh.buckets(MetaData::ELEMENT_RANK);

  std::vector<EntityId> old_ids;
  old_ids.reserve(num_elems);
  BOOST_FOREACH(Bucket * b, elem_buckets) {
    for (size_t i =0; i<b->size(); ++i) {
      Entity e = (*b)[i];
      old_ids.push_back(mesh.identifier(e));
      mesh.change_entity_id( mesh.identifier(e)+num_elems, e);
    }
  }

  mesh.modification_end();

  mesh.modification_begin();
  mesh.modification_end();

  std::vector<EntityId> new_ids_minus_num_elems;
  new_ids_minus_num_elems.reserve(num_elems);
  BOOST_FOREACH(Bucket * b, elem_buckets) {
    for (size_t i =0; i<b->size(); ++i) {
      Entity e = (*b)[i];
      new_ids_minus_num_elems.push_back(mesh.identifier(e)-num_elems);
    }
  }

  STKUNIT_EXPECT_TRUE(old_ids == new_ids_minus_num_elems);

  BOOST_FOREACH(Bucket * b, node_buckets) {
    int* nodal_field = mesh.field_data( simple_nodal_field, *b );
    for (size_t i =0; i<b->size(); ++i) {
      STKUNIT_EXPECT_TRUE( nodal_field[i] == 1);
    }
  }
}

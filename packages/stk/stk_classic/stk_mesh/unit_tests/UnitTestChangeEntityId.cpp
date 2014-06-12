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
#include <stk_mesh/base/FieldData.hpp>

#include <boost/foreach.hpp>


STKUNIT_UNIT_TEST( UnitTestChangeEntityId, change_id )
{
  using namespace stk_classic::mesh;

  const unsigned NX = 50;
  const unsigned NY = 50;
  const unsigned NZ = 50;
  const unsigned num_elems = NX * NY * NZ;

  fixtures::HexFixture hf(MPI_COMM_WORLD,NX,NY,NZ);

  Field<int> & simple_nodal_field = hf.m_fem_meta.declare_field<Field<int> >("simple_nodal_field");

  put_field( simple_nodal_field,
             fem::FEMMetaData::NODE_RANK,
             hf.m_hex_part);


  //create nodal field on hex topo

  hf.m_fem_meta.commit();

  hf.generate_mesh();

  stk_classic::mesh::BulkData & mesh = hf.m_bulk_data;

  mesh.modification_begin();

  const BucketVector & nodes = mesh.buckets(fem::FEMMetaData::NODE_RANK);

  BOOST_FOREACH(Bucket * b, nodes) {
    BucketArray< Field<int> > nodal_field(simple_nodal_field,*b);
    for (int i =0; i<nodal_field.size(); ++i) {
      nodal_field[i] = 1;
    }
  }


  const BucketVector & elems = mesh.buckets(hf.m_fem_meta.element_rank());

  std::vector<EntityId> old_ids;
  old_ids.reserve(num_elems);
  BOOST_FOREACH(Bucket * b, elems) {
    for (size_t i =0; i<b->size(); ++i) {
      Entity & e = (*b)[i];
      old_ids.push_back(e.identifier());
      mesh.change_entity_id( e.identifier()+num_elems, e);
    }
  }

  mesh.modification_end();

  mesh.modification_begin();
  mesh.modification_end();

  std::vector<EntityId> new_ids_minus_num_elems;
  new_ids_minus_num_elems.reserve(num_elems);
  BOOST_FOREACH(Bucket * b, elems) {
    for (size_t i =0; i<b->size(); ++i) {
      Entity & e = (*b)[i];
      new_ids_minus_num_elems.push_back(e.identifier()-num_elems);
    }
  }

  STKUNIT_EXPECT_TRUE(old_ids == new_ids_minus_num_elems);

  BOOST_FOREACH(Bucket * b, nodes) {
    BucketArray< Field<int> > nodal_field(simple_nodal_field,*b);
    for (int i =0; i<nodal_field.size(); ++i) {
      STKUNIT_EXPECT_TRUE( nodal_field[i] == 1);
    }
  }

}


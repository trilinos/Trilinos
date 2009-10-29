
#include <unit_tests/stk_utest_macros.hpp>

#include <mpi.h>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_mesh/base/GetEntities.hpp>

#include <stk_linsys/DofMapper.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include <unit_tests/UnitTest_helpers.hpp>

namespace stk_linsys_unit_tests {


//------------- here is the DofMapper unit-test... -----------------------

void testDofMapper( MPI_Comm comm )
{
  //First create and fill MetaData and BulkData objects:

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_type_names() );

  fill_utest_mesh_meta_data( meta_data );
  meta_data.commit();

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk::mesh::BulkData bulk_data( meta_data, comm, bucket_size );
  fill_utest_mesh_bulk_data( bulk_data );

  bulk_data.modification_end();

  stk::mesh::Selector selector(meta_data.locally_used_part());
  std::vector<unsigned> count;
  stk::mesh::count_entities(selector, bulk_data, count);

  STKUNIT_ASSERT_EQUAL( count[stk::mesh::Element], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[stk::mesh::Node],    (unsigned)20 );

  std::vector<stk::mesh::Entity*> nodes;
  stk::mesh::get_entities(bulk_data, stk::mesh::Node, nodes);

  stk::mesh::ScalarField* temperature_field =
      meta_data.get_field<stk::mesh::ScalarField>("temperature");

  //Now we're ready to test the DofMapper:

  stk::linsys::DofMapper dof_mapper(comm);
 
  stk::mesh::PartVector part_intersection;
  part_intersection.push_back(&(meta_data.locally_used_part()));

  dof_mapper.add_dof_mappings(bulk_data, part_intersection,
                              stk::mesh::Node, *temperature_field);

  stk::mesh::EntityType ent_type;
  stk::mesh::EntityId ent_id;
  const stk::mesh::FieldBase* field = NULL;
  int offset_into_field;
  int index = 0;
  //DofMapper::get_dof can't be called until after DofMapper::finalize() has
  //been called.
  //We'll call it now to verify that an exception is thrown:
  std::cout << "Testing error condition: " << std::endl;
  STKUNIT_ASSERT_THROW(dof_mapper.get_dof(index, ent_type, ent_id, field, offset_into_field), std::runtime_error );
  std::cout << "...Completed testing error condition." << std::endl;

  dof_mapper.finalize();

  //find a node that is in the locally-used part:
  size_t i_node = 0;
  while(!nodes[i_node]->bucket().member(meta_data.locally_used_part()) && i_node<nodes.size()) {
    ++i_node;
  }

  //test the get_global_index function:
  stk::mesh::EntityId node_id = nodes[i_node]->identifier();
  index = dof_mapper.get_global_index(stk::mesh::Node, node_id, *temperature_field);
  STKUNIT_ASSERT_EQUAL( index, (int)(node_id-1) );

  std::cout << "Testing error condition: " << std::endl;
  //call DofMapper::get_global_index with a non-existent ID and verify that an
  //exception is thrown:
  STKUNIT_ASSERT_THROW(dof_mapper.get_global_index(stk::mesh::Node, (stk::mesh::EntityId)999999, *temperature_field), std::runtime_error);
  std::cout << "...Completed testing error condition." << std::endl;

  int numProcs = 1;
  MPI_Comm_size(comm, &numProcs);

  fei::SharedPtr<fei::VectorSpace> fei_vspace = dof_mapper.get_fei_VectorSpace();
  int numIndices = fei_vspace->getGlobalNumIndices();
  STKUNIT_ASSERT_EQUAL( numIndices, (int)(numProcs*20 - (numProcs-1)*4) );

  dof_mapper.get_dof(index, ent_type, ent_id, field, offset_into_field);

  STKUNIT_ASSERT_EQUAL( ent_type, nodes[i_node]->entity_type() );
  STKUNIT_ASSERT_EQUAL( ent_id,   nodes[i_node]->identifier() );
  STKUNIT_ASSERT_EQUAL( field->name(), temperature_field->name() );
  STKUNIT_ASSERT_EQUAL( offset_into_field, (int)0 );
}

} // namespace stk_linsys_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfDofMapper, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_linsys_unit_tests::testDofMapper ( MPI_COMM_WORLD );
}


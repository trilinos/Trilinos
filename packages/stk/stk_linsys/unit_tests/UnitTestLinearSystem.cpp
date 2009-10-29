
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
#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/LinsysFunctions.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include <unit_tests/UnitTest_helpers.hpp>

#include <fei_Factory_Trilinos.hpp>


namespace stk_linsys_unit_tests {

//------------- here is the LinearSystem unit-test... -----------------------

void testLinearSystem( MPI_Comm comm )
{
  //First create and fill MetaData and BulkData objects:

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_type_names() );

  fill_utest_mesh_meta_data( meta_data );

  //create a boundary-condition part for testing later:
  stk::mesh::Part& bcpart = meta_data.declare_part("bcpart");

  meta_data.commit();

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk::mesh::BulkData bulk_data( meta_data, comm, bucket_size );
  fill_utest_mesh_bulk_data( bulk_data );

  bulk_data.modification_end();

  //put a node in our boundary-condition part. arbitrarily choose the
  //first locally-owned node:

  bulk_data.modification_begin();

  std::vector<stk::mesh::Entity*> local_nodes;
  stk::mesh::Selector select_owned(meta_data.locally_owned_part());
  stk::mesh::get_selected_entities(select_owned,
                                   bulk_data.buckets(stk::mesh::Node),
                                   local_nodes);

  stk::mesh::EntityId bc_node_id = 0;

  if (local_nodes.size() > 0) {
    stk::mesh::PartVector partvector;
    partvector.push_back(&bcpart);
    bulk_data.change_entity_parts(*local_nodes[0], partvector);
    bc_node_id = stk::linsys::impl::entityid_to_int(local_nodes[0]->identifier());
  }

  bulk_data.modification_end();

  stk::mesh::Selector select_used(meta_data.locally_used_part());
  std::vector<unsigned> count;
  stk::mesh::count_entities(select_used, bulk_data, count);

  STKUNIT_ASSERT_EQUAL( count[stk::mesh::Element], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[stk::mesh::Node],    (unsigned)20 );

  std::vector<stk::mesh::Entity*> nodes;
  stk::mesh::get_entities(bulk_data, stk::mesh::Node, nodes);

  stk::mesh::ScalarField* temperature_field =
      meta_data.get_field<stk::mesh::ScalarField>("temperature");

  //Now we're ready to test the LinearSystem object:

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));
  stk::linsys::LinearSystem ls(comm,factory);
 
  stk::mesh::PartVector owned_parts;
  owned_parts.push_back(&(meta_data.locally_owned_part()));

  stk::linsys::add_connectivities(ls, stk::mesh::Element, stk::mesh::Node,
                              *temperature_field, owned_parts, bulk_data);

  ls.synchronize_mappings_and_structure();

  const stk::linsys::LinearSystem& const_ls = ls;
  const fei::SharedPtr<fei::MatrixGraph> const_mgraph = const_ls.get_fei_MatrixGraph();
  fei::SharedPtr<fei::MatrixGraph> mgraph = ls.get_fei_MatrixGraph();

  const fei::MatrixGraph* const_mgraph_ptr = const_mgraph.get();
  fei::MatrixGraph* mgraph_ptr = mgraph.get();

  STKUNIT_ASSERT_EQUAL( mgraph_ptr, (fei::MatrixGraph*)const_mgraph_ptr );
  
  ls.create_fei_LinearSystem();

  const fei::LinearSystem* const_linsys_ptr = const_ls.get_fei_LinearSystem().get();
  fei::LinearSystem* linsys_ptr = ls.get_fei_LinearSystem().get();

  STKUNIT_ASSERT_EQUAL( (fei::LinearSystem*)const_linsys_ptr, linsys_ptr);

  int numProcs = 1, myProc = 0;
  MPI_Comm_rank(comm, &myProc);
  MPI_Comm_size(comm, &numProcs);

  fei::SharedPtr<fei::Matrix> matrix = ls.get_fei_LinearSystem()->getMatrix();
  int num_local_rows = matrix->getLocalNumRows();
  int expected_num_local_rows = myProc==0 ? 20 : 16;

  STKUNIT_ASSERT_EQUAL( num_local_rows, expected_num_local_rows );

  stk::linsys::DofMapper& dof_mapper = ls.get_DofMapper();

  const stk::linsys::DofMapper& const_dof_mapper = const_ls.get_DofMapper();

  STKUNIT_ASSERT_EQUAL( (stk::linsys::DofMapper*)&dof_mapper , (stk::linsys::DofMapper*)&const_dof_mapper );
}

} // namespace stk_linsys_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfLinearSystem, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_linsys_unit_tests::testLinearSystem ( MPI_COMM_WORLD );
}


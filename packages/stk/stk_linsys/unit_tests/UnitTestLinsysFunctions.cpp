#include <unit_tests/stk_utest_macros.hpp>

#include <cmath>

#include <mpi.h>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SetOwners.hpp>

#include <stk_mesh/base/GetEntities.hpp>

#include <stk_linsys/DofMapper.hpp>
#include <stk_linsys/LinsysFunctions.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include <fei_Factory_Trilinos.hpp>
#include <Teuchos_ParameterList.hpp>


#include <unit_tests/UnitTest_helpers.hpp>

namespace stk_linsys_unit_tests {

//------------- here is the LinsysFunctions unit-test... -----------------------

void testLinsysFunctions( MPI_Comm comm )
{
  //First create and fill MetaData and BulkData objects:

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_type_names() );
  stk::mesh::MetaData meta_data2( stk::mesh::fem_entity_type_names() );

  fill_utest_mesh_meta_data( meta_data );

  bool use_temperature=false;
  fill_utest_mesh_meta_data( meta_data2, use_temperature );

  //create a boundary-condition part for testing later:
  stk::mesh::Part& bcpart = meta_data.declare_part("bcpart");

  meta_data.commit();
  meta_data2.commit();

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk::mesh::BulkData bulk_data( meta_data, comm, bucket_size );
  stk::mesh::BulkData bulk_data2( meta_data2, comm, bucket_size );
  fill_utest_mesh_bulk_data( bulk_data );
  fill_utest_mesh_bulk_data( bulk_data2 );

  bulk_data.modification_end();
  bulk_data2.modification_end();

  //set owner-processors to lowest-sharing (stk::mesh defaults to
  //highest-sharing) If highest-sharing owns, then it isn't correct for the
  //way the fei library sets ownership of shared nodes for vectors etc.
  stk::mesh::set_owners<stk::mesh::LowestRankSharingProcOwns>( bulk_data );

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

  stk::mesh::PartVector parts;
  parts.push_back(&(meta_data.locally_used_part()));
  parts.push_back(meta_data.get_part("block_1"));

  stk::mesh::Selector selector(parts);
  std::vector<unsigned> count;
  stk::mesh::count_entities(selector, bulk_data, count);

  STKUNIT_ASSERT_EQUAL( count[stk::mesh::Element], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[stk::mesh::Node],    (unsigned)20 );

  stk::mesh::ScalarField* temperature_field =
      meta_data.get_field<stk::mesh::ScalarField>("temperature");

  //Create a DofMapper:

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));

  stk::linsys::LinearSystem ls(comm, factory);

  stk::linsys::add_connectivities(ls, stk::mesh::Element, stk::mesh::Node,
                                  *temperature_field, parts, bulk_data);

  fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();
  int num_blocks = matgraph->getNumConnectivityBlocks();

  STKUNIT_ASSERT_EQUAL( num_blocks, (int)1 );

  ls.synchronize_mappings_and_structure();
  ls.create_fei_LinearSystem();

  //put 0 throughout the matrix and 3 throughout the rhs:
  fei::SharedPtr<fei::Matrix> mat = ls.get_fei_LinearSystem()->getMatrix();
  ls.get_fei_LinearSystem()->getMatrix()->putScalar(0);
  ls.get_fei_LinearSystem()->getRHS()->putScalar(3.0);

  //put 10 on the matrix diagonal to ensure it will be easy to solve later.
  fei::SharedPtr<fei::VectorSpace> vspace = ls.get_fei_LinearSystem()->getRHS()->getVectorSpace();
  int numLocalRows = vspace->getNumIndices_Owned();
  std::vector<int> local_rows(numLocalRows);
  vspace->getIndices_Owned(numLocalRows, &local_rows[0], numLocalRows);

  for(size_t i=0; i<local_rows.size(); ++i) {
    int col = local_rows[i];
    double coef = 10;
    double* coefPtr = &coef;
    mat->sumIn(1, &local_rows[i], 1, &col, &coefPtr);
  }

  //now we'll impose a dirichlet bc on our one-node bcpart:
  stk::linsys::dirichlet_bc(ls, bulk_data, bcpart, stk::mesh::Node,
                  *temperature_field, 0, 9.0);

  ls.finalize_assembly();

  //now confirm that the rhs value for the equation corresponding to our
  //bc node is 9.0:

  fei::SharedPtr<fei::Vector> rhsvec = ls.get_fei_LinearSystem()->getRHS();
  double rhs_bc_val = 0;
  int bc_eqn_index = ls.get_DofMapper().get_global_index(stk::mesh::Node,
           bc_node_id, *temperature_field);
  rhsvec->copyOut(1, &bc_eqn_index, &rhs_bc_val);

  bool bc_val_is_correct = std::abs(rhs_bc_val - 9.0) < 1.e-13;
  STKUNIT_ASSERT( bc_val_is_correct );

  stk::linsys::copy_vector_to_mesh( *rhsvec, ls.get_DofMapper(), bulk_data);

  stk::mesh::Entity* bc_node = bulk_data.get_entity(stk::mesh::Node, local_nodes[0]->identifier());

  stk::mesh::FieldTraits<stk::mesh::ScalarField>::data_type* bc_node_data =
    stk::mesh::field_data(*temperature_field, *bc_node);

  bool bc_node_data_is_correct = std::abs(bc_node_data[0] - 9.0) < 1.e-13;
  STKUNIT_ASSERT( bc_node_data_is_correct );

  //now make sure we get a throw if we use the wrong bulk-data (that doesn't have the
  //temperature field defined)
  STKUNIT_ASSERT_THROW(stk::linsys::copy_vector_to_mesh( *rhsvec, ls.get_DofMapper(), bulk_data2), std::runtime_error);

  //obtain and zero the solution vector
  fei::SharedPtr<fei::Vector> solnvec = ls.get_fei_LinearSystem()->getSolutionVector();
  solnvec->putScalar(0);

  //copy the vector of zeros into the mesh:
  stk::linsys::copy_vector_to_mesh( *solnvec, ls.get_DofMapper(), bulk_data);

  //assert that our bc node's data is now zero.
  bc_node_data_is_correct = std::abs(bc_node_data[0] - 0) < 1.e-13;
  STKUNIT_ASSERT( bc_node_data_is_correct );

  //call the linear-system solve function.
  //(note that when we add options to the solve method, we'll need to enhance this
  //testing to exercise various specific solves.)
  Teuchos::ParameterList params;

  int status = 0;
  ls.solve(status, params);

  //copy the solution-vector into the mesh:
  stk::linsys::copy_vector_to_mesh( *solnvec, ls.get_DofMapper(), bulk_data);

  //now assert that the value 9 (bc value) produced by the solve is in this
  //node's data.
  //note that we use a looser tolerance now, because the default solver tolerance
  //is only (I think) 1.e-6.
  bc_node_data_is_correct = std::abs(bc_node_data[0] - 9.0) < 1.e-6;
  STKUNIT_ASSERT( bc_node_data_is_correct );
}

} // namespace stk_linsys_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfLinsysFunctions, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_linsys_unit_tests::testLinsysFunctions ( MPI_COMM_WORLD );
}


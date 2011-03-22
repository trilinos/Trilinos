/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <cmath>

#include <stk_util/parallel/Parallel.hpp>
#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>
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
#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include <fei_Factory_Trilinos.hpp>
#include <Teuchos_ParameterList.hpp>

#include <unit_tests/UnitTest_helpers.hpp>

using stk::mesh::fem::NODE_RANK;

namespace stk_linsys_unit_tests {

typedef stk::mesh::Field<double>                          ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>    VectorField ;

bool confirm_vector_values(const fei::Vector& vec, double expected_value)
{
  std::vector<int> indices;
  fei::SharedPtr<fei::VectorSpace> vspace = vec.getVectorSpace();
  vspace->getIndices_Owned(indices);
  bool result = true;
  if (indices.size() > 0) {
    std::vector<double> coefs(indices.size());
    vec.copyOut(indices.size(), &indices[0], &coefs[0]);

    for(size_t i=0; i<indices.size(); ++i) {
      if (std::abs(coefs[i] - expected_value) > 1.e-13) {
        result = false;
        break;
      }
    }
  }
  return result;
}

bool confirm_matrix_values(const fei::Matrix& mat, double expected_value)
{
  std::vector<int> rows;
  fei::SharedPtr<fei::VectorSpace> rspace = mat.getMatrixGraph()->getRowSpace();
  rspace->getIndices_Owned(rows);
  bool result = true;
  std::vector<int> indices;
  std::vector<double> coefs;
  for(size_t i=0; i<rows.size(); ++i) {
    int rowlength = 0;
    mat.getRowLength(rows[i], rowlength);
    indices.resize(rowlength);
    coefs.resize(rowlength);
    mat.copyOutRow(rows[i], rowlength, &coefs[0], &indices[0]);

    for(size_t j=0; j<indices.size(); ++j) {
      if (std::abs(coefs[j] - expected_value) > 1.e-13) {
        result = false;
        break;
      }
    }
  }
  return result;
}

//------------- LinsysFunctions unit-tests... -----------------------

STKUNIT_UNIT_TEST(UnitTestLinsysFunctions, test1)
{
  static const size_t spatial_dimension = 3;

  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Comm comm = MPI_COMM_WORLD;
  //First create and fill MetaData and BulkData objects:

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk::mesh::fem::FEMMetaData fem_meta;
  stk::mesh::fem::FEMMetaData fem_meta2;
  fem_meta.FEM_initialize(spatial_dimension);
  fem_meta2.FEM_initialize(spatial_dimension);

  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::MetaData & meta_data2 = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta2);

  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  stk::mesh::BulkData bulk_data( meta_data, comm, bucket_size );
  stk::mesh::BulkData bulk_data2( meta_data2, comm, bucket_size );

  //create a boundary-condition part for testing later:
  stk::mesh::Part& bcpart = fem_meta.declare_part("bcpart");

  fill_utest_mesh_meta_data( fem_meta );

  bool use_temperature=false;
  fill_utest_mesh_meta_data( fem_meta2, use_temperature );

  fill_utest_mesh_bulk_data( bulk_data );
  fill_utest_mesh_bulk_data( bulk_data2 );

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
                                   bulk_data.buckets(NODE_RANK),
                                   local_nodes);

  stk::mesh::EntityId bc_node_id = 0;

  if (local_nodes.size() > 0) {
    stk::mesh::PartVector partvector;
    partvector.push_back(&bcpart);
    bulk_data.change_entity_parts(*local_nodes[0], partvector);
    bc_node_id = stk::linsys::impl::entityid_to_int(local_nodes[0]->identifier());
  }

  bulk_data.modification_end();

  stk::mesh::Selector selector = ( meta_data.locally_owned_part() | meta_data.globally_shared_part() ) & *meta_data.get_part("block_1");
  std::vector<unsigned> count;
  stk::mesh::count_entities(selector, bulk_data, count);

  STKUNIT_ASSERT_EQUAL( count[element_rank], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[NODE_RANK],     (unsigned)20 );

  ScalarField* temperature_field = meta_data.get_field<ScalarField>("temperature");

  //Create a fei Factory and stk::linsys::LinearSystem object:

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));

  stk::linsys::LinearSystem ls(comm, factory);

  stk::linsys::add_connectivities(ls, element_rank, NODE_RANK,
                                  *temperature_field, selector, bulk_data);

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
  stk::linsys::dirichlet_bc(ls, bulk_data, bcpart, NODE_RANK,
                            *temperature_field, 0, 9.0);

  ls.finalize_assembly();

  //now confirm that the rhs value for the equation corresponding to our
  //bc node is 9.0:

  fei::SharedPtr<fei::Vector> rhsvec = ls.get_fei_LinearSystem()->getRHS();
  double rhs_bc_val = 0;
  int bc_eqn_index = ls.get_DofMapper().get_global_index(NODE_RANK,
           bc_node_id, *temperature_field);
  rhsvec->copyOut(1, &bc_eqn_index, &rhs_bc_val);

  bool bc_val_is_correct = std::abs(rhs_bc_val - 9.0) < 1.e-13;
  STKUNIT_ASSERT( bc_val_is_correct );

  stk::linsys::copy_vector_to_mesh( *rhsvec, ls.get_DofMapper(), bulk_data);

  stk::mesh::Entity* bc_node = bulk_data.get_entity(NODE_RANK, local_nodes[0]->identifier());

  stk::mesh::FieldTraits<ScalarField>::data_type* bc_node_data = stk::mesh::field_data(*temperature_field, *bc_node);

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
  //note that we use a loose tolerance, because the default solver tolerance
  //is (I think) only 1.e-6.
  bc_node_data_is_correct = std::abs(bc_node_data[0] - 9.0) < 1.e-6;
  STKUNIT_ASSERT( bc_node_data_is_correct );
  STKUNIT_ASSERT(bc_node_data_is_correct);
}

STKUNIT_UNIT_TEST(UnitTestLinsysFunctions, test2)
{
  static const size_t spatial_dimension = 3;

  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Comm comm = MPI_COMM_WORLD;
  //First create and fill MetaData and BulkData objects:

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dimension);

  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);

  stk::mesh::BulkData bulk_data( meta_data, comm, bucket_size );

  fill_utest_mesh_meta_data( fem_meta );

  fill_utest_mesh_bulk_data( bulk_data );

  //set owner-processors to lowest-sharing (stk::mesh defaults to
  //highest-sharing) If highest-sharing owns, then it isn't correct for the
  //way the fei library sets ownership of shared nodes for vectors etc.
  stk::mesh::set_owners<stk::mesh::LowestRankSharingProcOwns>( bulk_data );

  stk::mesh::Selector selector = ( meta_data.locally_owned_part() | meta_data.globally_shared_part() ) & *meta_data.get_part("block_1");
  std::vector<unsigned> count;
  stk::mesh::count_entities(selector, bulk_data, count);
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  STKUNIT_ASSERT_EQUAL( count[element_rank], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[NODE_RANK],     (unsigned)20 );

  ScalarField* temperature_field = meta_data.get_field<ScalarField>("temperature");

  //Create a fei Factory and stk::linsys::LinearSystem object:

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));

  stk::linsys::LinearSystem ls(comm, factory);

  stk::linsys::add_connectivities(ls, element_rank, NODE_RANK,
                                  *temperature_field, selector, bulk_data);

  fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();
  int num_blocks = matgraph->getNumConnectivityBlocks();

  STKUNIT_ASSERT_EQUAL( num_blocks, (int)1 );

  ls.synchronize_mappings_and_structure();
  ls.create_fei_LinearSystem();

  //put 3 throughout the matrix and 3 throughout the rhs:
  fei::SharedPtr<fei::Matrix> mat = ls.get_fei_LinearSystem()->getMatrix();
  mat->putScalar(3.0);
  ls.get_fei_LinearSystem()->getRHS()->putScalar(3.0);

  fei::SharedPtr<fei::Vector> rhsvec = ls.get_fei_LinearSystem()->getRHS();

  stk::linsys::scale_vector(2, *rhsvec);
  stk::linsys::scale_matrix(2, *mat);

  //now confirm that the rhs and matrix contain 6:

  bool result = confirm_vector_values(*rhsvec, 6);
  STKUNIT_ASSERT( result );
  result = confirm_matrix_values(*mat, 6);
  STKUNIT_ASSERT( result );
}

STKUNIT_UNIT_TEST(UnitTestLinsysFunctions, test3)
{
  static const size_t spatial_dimension = 3;

  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Comm comm = MPI_COMM_WORLD;
  //First create and fill MetaData and BulkData objects:

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dimension);
  stk::mesh::MetaData & meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk::mesh::BulkData bulk_data( meta_data, comm, bucket_size );

  fill_utest_mesh_meta_data( fem_meta );

  fill_utest_mesh_bulk_data( bulk_data );

  //set owner-processors to lowest-sharing (stk::mesh defaults to
  //highest-sharing) If highest-sharing owns, then it isn't correct for the
  //way the fei library sets ownership of shared nodes for vectors etc.
  stk::mesh::set_owners<stk::mesh::LowestRankSharingProcOwns>( bulk_data );

  stk::mesh::Selector selector = ( meta_data.locally_owned_part() | meta_data.globally_shared_part() ) & *meta_data.get_part("block_1");
  std::vector<unsigned> count;
  stk::mesh::count_entities(selector, bulk_data, count);
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  STKUNIT_ASSERT_EQUAL( count[element_rank], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[NODE_RANK],     (unsigned)20 );

  ScalarField* temperature_field = meta_data.get_field<ScalarField>("temperature");

  //Create a fei Factory and stk::linsys::LinearSystem object:

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));

  stk::linsys::LinearSystem ls(comm, factory);

  stk::linsys::add_connectivities(ls, element_rank, NODE_RANK,
                                  *temperature_field, selector, bulk_data);

  fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();
  int num_blocks = matgraph->getNumConnectivityBlocks();

  STKUNIT_ASSERT_EQUAL( num_blocks, (int)1 );

  ls.synchronize_mappings_and_structure();
  ls.create_fei_LinearSystem();

  //put 3 throughout the matrix and 3 throughout the rhs:
  fei::SharedPtr<fei::Matrix> mat = ls.get_fei_LinearSystem()->getMatrix();
  mat->putScalar(3.0);
  ls.get_fei_LinearSystem()->getRHS()->putScalar(3.0);

  fei::SharedPtr<fei::Vector> rhsvec = ls.get_fei_LinearSystem()->getRHS();

  stk::linsys::scale_vector(2, *rhsvec);
  stk::linsys::scale_matrix(2, *mat);

  //now the rhs and matrix contain 6.

  //create another matrix and vector:
  fei::SharedPtr<fei::Matrix> mat2 = factory->createMatrix(matgraph);
  fei::SharedPtr<fei::Vector> vec2 = factory->createVector(matgraph);
  mat2->putScalar(3.0);
  vec2->putScalar(3.0);

  //add 3*mat to mat2
  stk::linsys::add_matrix_to_matrix(3.0, *mat, *mat2);

  //confirm that mat2 contains 21:
  bool result = confirm_matrix_values(*mat2, 21);
  STKUNIT_ASSERT(result);

  //add 3*rhsvec to vec2:
  stk::linsys::add_vector_to_vector(3.0, *rhsvec, *vec2);

  //confirm that vec2 contains 21:
  result = confirm_vector_values(*vec2, 21);
  STKUNIT_ASSERT(result);
}

} // namespace stk_linsys_unit_tests

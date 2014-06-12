/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_mesh/base/GetEntities.hpp>

#include <stk_linsys/DofMapper.hpp>
#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/AggregateLinearSystem.hpp>
#include <stk_linsys/LinsysFunctions.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include <unit_tests/UnitTest_helpers.hpp>

#include <fei_Factory_Trilinos.hpp>

static const stk_classic::mesh::EntityRank NODE_RANK = stk_classic::mesh::fem::FEMMetaData::NODE_RANK;

namespace stk_linsys_unit_tests {

typedef stk_classic::mesh::Field<double>                          ScalarField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian>    VectorField ;

void get_local_coefs(const fei::Vector& vec, std::vector<double>& coefs)
{
  fei::SharedPtr<fei::VectorSpace> vspace = vec.getVectorSpace();
  std::vector<int> indices;
  vspace->getIndices_Owned(indices);
  coefs.resize(indices.size());
  vec.copyOut(indices.size(), &indices[0], &coefs[0]);
}

void get_first_local_row_coefs(const fei::Matrix& mat, std::vector<double>& row_coefs)
{
  row_coefs.clear();
  fei::SharedPtr<fei::VectorSpace> vspace = mat.getMatrixGraph()->getRowSpace();
  std::vector<int> row_numbers;
  vspace->getIndices_Owned(row_numbers);
  std::vector<int> col_indices;
  int rowlength = 0;
  if (row_numbers.size() > 0) {
    mat.getRowLength(row_numbers[0], rowlength);
    if (rowlength > 0) {
      col_indices.resize(rowlength);
      row_coefs.resize(rowlength);
      mat.copyOutRow(row_numbers[0], rowlength, &row_coefs[0], &col_indices[0]);
    }
  }
}

STKUNIT_UNIT_TEST(UnitTestLinearSystem, test1)
{
  static const size_t spatial_dimension = 3;
  
  MPI_Comm comm = MPI_COMM_WORLD;

  //First create and fill MetaData and BulkData objects:

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk_classic::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dimension);
  stk_classic::mesh::MetaData & meta_data = stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk_classic::mesh::BulkData bulk_data( meta_data, comm, bucket_size );

  //create a boundary-condition part for testing later:
  stk_classic::mesh::Part& bcpart = fem_meta.declare_part("bcpart");

  fill_utest_mesh_meta_data( fem_meta );
  fill_utest_mesh_bulk_data( bulk_data );

  //------------------------------
  //put a node in our boundary-condition part. arbitrarily choose the
  //first locally-owned node:

  bulk_data.modification_begin();

  std::vector<stk_classic::mesh::Entity*> local_nodes;
  stk_classic::mesh::Selector select_owned(meta_data.locally_owned_part());
  stk_classic::mesh::get_selected_entities(select_owned,
                                   bulk_data.buckets(NODE_RANK),
                                   local_nodes);

  stk_classic::mesh::EntityId bc_node_id = 0;

  if (local_nodes.size() > 0) {
    stk_classic::mesh::PartVector partvector;
    partvector.push_back(&bcpart);
    bulk_data.change_entity_parts(*local_nodes[0], partvector);
    bc_node_id = stk_classic::linsys::impl::entityid_to_int(local_nodes[0]->identifier());
  }
  (void)bc_node_id; // to avoid a set but not used warning from gcc 4.6.3

  bulk_data.modification_end();
  //------------------------------

  stk_classic::mesh::Selector select_used = meta_data.locally_owned_part() | meta_data.globally_shared_part() ;
  std::vector<unsigned> count;
  stk_classic::mesh::count_entities(select_used, bulk_data, count);
  const stk_classic::mesh::EntityRank element_rank = fem_meta.element_rank();

  STKUNIT_ASSERT_EQUAL( count[element_rank], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[NODE_RANK],    (unsigned)20 );

  std::vector<stk_classic::mesh::Entity*> nodes;
  stk_classic::mesh::get_entities(bulk_data, NODE_RANK, nodes);

  ScalarField* temperature_field = meta_data.get_field<ScalarField>("temperature");

  //Now we're ready to test the LinearSystem object:

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));
  stk_classic::linsys::LinearSystem ls(comm,factory);

  stk_classic::linsys::add_connectivities(ls, element_rank, NODE_RANK,
                              *temperature_field, select_owned, bulk_data);

  ls.synchronize_mappings_and_structure();

  const stk_classic::linsys::LinearSystem& const_ls = ls;
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
  myProc = stk_classic::parallel_machine_rank( MPI_COMM_WORLD );
  numProcs = stk_classic::parallel_machine_size( MPI_COMM_WORLD );
  STKUNIT_ASSERT(numProcs>=1);

  fei::SharedPtr<fei::Matrix> matrix = ls.get_fei_LinearSystem()->getMatrix();
  int num_local_rows = matrix->getLocalNumRows();
  int expected_num_local_rows = myProc==0 ? 20 : 16;

  STKUNIT_ASSERT_EQUAL( num_local_rows, expected_num_local_rows );

  stk_classic::linsys::DofMapper& dof_mapper = ls.get_DofMapper();

  const stk_classic::linsys::DofMapper& const_dof_mapper = const_ls.get_DofMapper();

  STKUNIT_ASSERT_EQUAL( (stk_classic::linsys::DofMapper*)&dof_mapper , (stk_classic::linsys::DofMapper*)&const_dof_mapper );

  assemble_elem_matrices_and_vectors(bulk_data, *temperature_field, ls);

  Teuchos::ParameterList params;
  int status = 0, return_code = 0;
  return_code = ls.solve(status, params);
  STKUNIT_ASSERT_EQUAL(return_code, (int)0);
}

STKUNIT_UNIT_TEST(UnitTestAggregateLinearSystem, test1)
{
  static const size_t spatial_dimension = 3;

  MPI_Comm comm = MPI_COMM_WORLD;

  //First create and fill MetaData and BulkData objects:

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk_classic::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dimension);
  stk_classic::mesh::MetaData & meta_data = stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk_classic::mesh::BulkData bulk_data( meta_data, comm, bucket_size );

  //create a boundary-condition part for testing later:
  stk_classic::mesh::Part& bcpart = fem_meta.declare_part("bcpart");

  fill_utest_mesh_meta_data( fem_meta );
  fill_utest_mesh_bulk_data( bulk_data );

  //------------------------------
  //put a node in our boundary-condition part. arbitrarily choose the
  //first locally-owned node:

  bulk_data.modification_begin();

  std::vector<stk_classic::mesh::Entity*> local_nodes;
  stk_classic::mesh::Selector select_owned(meta_data.locally_owned_part());
  stk_classic::mesh::get_selected_entities(select_owned,
                                   bulk_data.buckets(NODE_RANK),
                                   local_nodes);

  stk_classic::mesh::EntityId bc_node_id = 0;

  if (local_nodes.size() > 0) {
    stk_classic::mesh::PartVector partvector;
    partvector.push_back(&bcpart);
    bulk_data.change_entity_parts(*local_nodes[0], partvector);
    bc_node_id = stk_classic::linsys::impl::entityid_to_int(local_nodes[0]->identifier());
  }
  (void)bc_node_id; // to avoid a set but not used warning from gcc 4.6.3

  bulk_data.modification_end();
  //------------------------------

  stk_classic::mesh::Selector select_used = meta_data.locally_owned_part() | meta_data.globally_shared_part() ;
  std::vector<unsigned> count;
  stk_classic::mesh::count_entities(select_used, bulk_data, count);
  const stk_classic::mesh::EntityRank element_rank = fem_meta.element_rank();

  STKUNIT_ASSERT_EQUAL( count[element_rank], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[NODE_RANK],    (unsigned)20 );

  std::vector<stk_classic::mesh::Entity*> nodes;
  stk_classic::mesh::get_entities(bulk_data, NODE_RANK, nodes);

  ScalarField* temperature_field = meta_data.get_field<ScalarField>("temperature");

  //Now we're ready to test the AggregateLinearSystem object:

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));
  size_t num_matrices = 2;
  size_t num_rhsvecs = 2;
  stk_classic::linsys::AggregateLinearSystem ls(comm,factory, num_matrices, num_rhsvecs);

  stk_classic::linsys::add_connectivities(ls, element_rank, NODE_RANK,
                              *temperature_field, select_owned, bulk_data);

  ls.synchronize_mappings_and_structure();

  const stk_classic::linsys::LinearSystemInterface& const_ls = ls;
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
  myProc = stk_classic::parallel_machine_rank( MPI_COMM_WORLD );
  numProcs = stk_classic::parallel_machine_size( MPI_COMM_WORLD );
  STKUNIT_ASSERT(numProcs >= 1);

  fei::SharedPtr<fei::Matrix> matrix = ls.get_fei_LinearSystem()->getMatrix();
  int num_local_rows = matrix->getLocalNumRows();
  int expected_num_local_rows = myProc==0 ? 20 : 16;

  STKUNIT_ASSERT_EQUAL( num_local_rows, expected_num_local_rows );

  stk_classic::linsys::DofMapper& dof_mapper = ls.get_DofMapper();

  const stk_classic::linsys::DofMapper& const_dof_mapper = const_ls.get_DofMapper();

  STKUNIT_ASSERT_EQUAL( (stk_classic::linsys::DofMapper*)&dof_mapper , (stk_classic::linsys::DofMapper*)&const_dof_mapper );

  fei::SharedPtr<fei::Matrix> mat0 = ls.get_matrix(0);
  fei::SharedPtr<fei::Vector> rhs0 = ls.get_rhsvec(0);
  assemble_elem_matrices_and_vectors(bulk_data, *temperature_field, dof_mapper, *mat0, *rhs0);
  fei::SharedPtr<fei::Matrix> mat1 = ls.get_matrix(1);
  fei::SharedPtr<fei::Vector> rhs1 = ls.get_rhsvec(1);
  assemble_elem_matrices_and_vectors(bulk_data, *temperature_field, dof_mapper, *mat1, *rhs1);

  //grab a copy of rhs1's coefficients for later comparison:
  std::vector<double> rhs1_coefs;
  get_local_coefs(*rhs1, rhs1_coefs);

  //grab a copy of mat1's first-row coefs for later comparison:
  std::vector<double> mat1_row_coefs;
  get_first_local_row_coefs(*mat1, mat1_row_coefs);

  std::vector<double> mat_scalars(num_matrices, 1.0);
  std::vector<double> rhs_scalars(num_rhsvecs, 1.0);
  ls.aggregate_system(mat_scalars, rhs_scalars);

  fei::SharedPtr<fei::Vector> rhs = ls.get_fei_LinearSystem()->getRHS();
  fei::SharedPtr<fei::Matrix> mat = ls.get_fei_LinearSystem()->getMatrix();

  //grab a copy of rhs's coefficients:
  std::vector<double> rhs_coefs;
  get_local_coefs(*rhs, rhs_coefs);

  //grab a copy of mat's first-row coefs:
  std::vector<double> mat_row_coefs;
  get_first_local_row_coefs(*mat, mat_row_coefs);

  //rhs and rhs1 should have the same length:
  bool rhs_test_passed = rhs_coefs.size() == rhs1_coefs.size();
  STKUNIT_ASSERT(rhs_test_passed);

  //mat and mat1 should have the same length for the first local row:
  bool mat_test_passed = mat_row_coefs.size() == mat1_row_coefs.size();
  STKUNIT_ASSERT(mat_test_passed);

  for(size_t i=0; i<rhs_coefs.size(); ++i) {
    //rhs coefs should be twice rhs1 coefs:
    if (std::abs(rhs1_coefs[i]) > 0.1) {
      rhs_test_passed = std::abs(rhs_coefs[0]/rhs1_coefs[0] - 2.0) < 1.e-8;
      STKUNIT_ASSERT(rhs_test_passed);
    }
  }

  for(size_t i=0; i<mat_row_coefs.size(); ++i) {
    //mat's row coefs should be twice mat1's row coefs:
    if (std::abs(mat1_row_coefs[i]) > 0.1) {
      mat_test_passed = std::abs(mat_row_coefs[0]/mat1_row_coefs[0] - 2.0) < 1.e-8;
      STKUNIT_ASSERT(mat_test_passed);
    }
  }

  Teuchos::ParameterList params;
  int status = 0, return_code = 0;
  return_code = ls.solve(status, params);
  STKUNIT_ASSERT_EQUAL(return_code, (int)0);
}

} // namespace stk_linsys_unit_tests


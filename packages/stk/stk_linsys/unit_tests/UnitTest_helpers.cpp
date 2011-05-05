/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_linsys/LinearSystemInterface.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

static const size_t NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

typedef stk::mesh::Field<double>                          ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>    VectorField ;

void fill_utest_mesh_meta_data(stk::mesh::fem::FEMMetaData& fem_meta, bool use_temperature)
{
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  stk::mesh::Part& elem_block = fem_meta.declare_part( "block_1" , element_rank );

  stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<8> >());
  stk::mesh::fem::set_cell_topology( elem_block, hex_top );

  const unsigned number_of_states = 1;

  if (use_temperature) {
    ScalarField& temperature_field = fem_meta.declare_field<ScalarField>( "temperature", number_of_states );
    stk::mesh::put_field( temperature_field, NODE_RANK,    elem_block );
  }

  ScalarField& pressure_field = fem_meta.declare_field<ScalarField>( "pressure", number_of_states );
  VectorField& velocity_field = fem_meta.declare_field<VectorField>( "velocity",    number_of_states );

  stk::mesh::put_field( pressure_field, element_rank, elem_block );
  stk::mesh::put_field( velocity_field, NODE_RANK,     elem_block );

  fem_meta.commit();
}

// Assume the following mesh of 4 hex8 elements on the first processor (proc 0).
//
// Global node and element numbering
//     3       7      11      15      19
//     +-------+-------+-------+-------+
//    /       /       /       /       /|
//  4/      8/     12/     16/     20/ |
//  +-------+-------+-------+-------+  |
//  |       |       |       |       |  +18
//  |  e1   |  e2   |  e3   |  e4   | /
//  |       |       |       |       |/
//  +-------+-------+-------+-------+
//  1       5      9       13      17
//
// Local node numbering
//     8       7
//     +-------+
//    /       /|
//  5/      6/ |
//  +-------+  |
//  |       |  +3
//  |  e1   | /
//  |       |/
//  +-------+
//  1       2
//
// A similar mesh will be created on the other processors, but using different elem-ids.
// i.e., proc 0 has elements 1 - 4, proc 1 has elements 5 - 8, and so on. 4 nodes will be
// shared between each neighboring pair of processors.
//----------------------------------------------------------------------

void elem_node_ids( stk::mesh::EntityId elem_id , stk::mesh::EntityId node_ids[] )
{
  if ( elem_id == 0 ) {
    std::cout << "use_case_1, elem_node_ids: ERROR, elem_id ("
        << elem_id << ") must be greater than 0." << std::endl;
    return;
  }

  const unsigned base = ( elem_id - 1 ) * 4 ;
  node_ids[0] = base + 1 ;
  node_ids[1] = base + 5 ;
  node_ids[2] = base + 6 ;
  node_ids[3] = base + 2 ;
  node_ids[4] = base + 4 ;
  node_ids[5] = base + 8 ;
  node_ids[6] = base + 7 ;
  node_ids[7] = base + 3 ;
}

void fill_utest_mesh_bulk_data(stk::mesh::BulkData& bulk_data)
{
  bulk_data.modification_begin();

  const unsigned num_elems = 4;

  int numProcs = 1, myProc = 0;
  numProcs = stk::parallel_machine_size( MPI_COMM_WORLD );
  myProc = stk::parallel_machine_rank( MPI_COMM_WORLD );

  stk::mesh::EntityId elem_id = 1 + myProc*num_elems;
  stk::mesh::EntityId node_ids[ shards::Hexahedron<8>::node_count ];

  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(bulk_data);
  stk::mesh::Part& elem_block = *(fem_meta.get_part( "block_1" ));

  for(unsigned i = 0; i<num_elems; ++i) {
    elem_node_ids( elem_id+i, node_ids );
    stk::mesh::fem::declare_element( bulk_data, elem_block, elem_id+i, node_ids );
  }

  bulk_data.modification_end();
}

void assemble_elem_matrices_and_vectors(stk::mesh::BulkData& mesh, ScalarField& field, stk::linsys::LinearSystemInterface& ls)
{
  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(mesh);
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  const std::vector<stk::mesh::Bucket*>& mesh_buckets = mesh.buckets(element_rank);

  std::vector<stk::mesh::Bucket*> part_buckets;
  stk::mesh::Selector select_owned(stk::mesh::MetaData::get(mesh).locally_owned_part());
  stk::mesh::get_buckets(select_owned, mesh_buckets, part_buckets);

  stk::linsys::DofMapper& dof_mapper = ls.get_DofMapper();

  int field_id = dof_mapper.get_field_id(field);

  stk::mesh::Entity& first_entity = *(part_buckets[0]->begin());
  stk::mesh::PairIterRelation rel = first_entity.relations(NODE_RANK);
  int num_nodes_per_elem = rel.second - rel.first;

  fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();
  int pattern_id = matgraph->definePattern(num_nodes_per_elem, NODE_RANK, field_id);

  std::vector<int> node_ids(num_nodes_per_elem);

  const int field_size = dof_mapper.get_fei_VectorSpace()->getFieldSize(field_id);
  const int matsize = num_nodes_per_elem*field_size*num_nodes_per_elem*field_size;
  const int vecsize = num_nodes_per_elem*field_size;

  std::vector<double> elem_matrix_1d(matsize, 0);
  std::vector<double*> elem_matrix_2d(vecsize);

  std::vector<double> elem_vector(vecsize, 0);

  for(size_t i=0; i<elem_matrix_2d.size(); ++i) {
    elem_matrix_2d[i] = &elem_matrix_1d[i*vecsize];
  }

  //fill our dummy elem-matrix:
  //This dummy matrix will be the same for every element. A real application
  //would form a different elem-matrix for each element.
  for(size_t i=0; i<elem_matrix_2d.size(); ++i) {
    double* row = elem_matrix_2d[i];
    if (i>=1) row[i-1] = -1;
    row[i] = 2;
    if (i<elem_matrix_2d.size()-1) row[i+1] = -1;

    elem_vector[i] = 1;
  }

  std::vector<int> eqn_indices(vecsize);
  fei::SharedPtr<fei::Matrix> matrix = ls.get_fei_LinearSystem()->getMatrix();
  fei::SharedPtr<fei::Vector> rhs = ls.get_fei_LinearSystem()->getRHS();

  for(size_t i=0; i<part_buckets.size(); ++i) {
    stk::mesh::Bucket::iterator
      b_iter = part_buckets[i]->begin(),
             b_end  = part_buckets[i]->end();
    for(; b_iter != b_end; ++b_iter) {
      stk::mesh::Entity& elem = *b_iter;
      rel = elem.relations(NODE_RANK);
      for(int j=0; rel.first != rel.second; ++rel.first, ++j) {
        node_ids[j] = rel.first->entity()->identifier();
      }

      matgraph->getPatternIndices(pattern_id, &node_ids[0], eqn_indices);

      matrix->sumIn(vecsize, &eqn_indices[0], vecsize, &eqn_indices[0],
                    &elem_matrix_2d[0]);
      rhs->sumIn(vecsize, &eqn_indices[0], &elem_vector[0]);
    }
  }
}

void assemble_elem_matrices_and_vectors(stk::mesh::BulkData& mesh, ScalarField& field, stk::linsys::DofMapper& dof_mapper, fei::Matrix& matrix, fei::Vector& rhs)
{
  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get(mesh);
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  const std::vector<stk::mesh::Bucket*>& mesh_buckets = mesh.buckets(element_rank);

  std::vector<stk::mesh::Bucket*> part_buckets;
  stk::mesh::Selector select_owned(stk::mesh::MetaData::get(mesh).locally_owned_part());
  stk::mesh::get_buckets(select_owned, mesh_buckets, part_buckets);

  int field_id = dof_mapper.get_field_id(field);

  stk::mesh::Entity& first_entity = *(part_buckets[0]->begin());
  stk::mesh::PairIterRelation rel = first_entity.relations(NODE_RANK);
  int num_nodes_per_elem = rel.second - rel.first;

  fei::SharedPtr<fei::MatrixGraph> matgraph = matrix.getMatrixGraph();
  int pattern_id = matgraph->definePattern(num_nodes_per_elem, NODE_RANK, field_id);

  std::vector<int> node_ids(num_nodes_per_elem);

  const int field_size = dof_mapper.get_fei_VectorSpace()->getFieldSize(field_id);
  const int matsize = num_nodes_per_elem*field_size*num_nodes_per_elem*field_size;
  const int vecsize = num_nodes_per_elem*field_size;

  std::vector<double> elem_matrix_1d(matsize, 0);
  std::vector<double*> elem_matrix_2d(vecsize);

  std::vector<double> elem_vector(vecsize, 0);

  for(size_t i=0; i<elem_matrix_2d.size(); ++i) {
    elem_matrix_2d[i] = &elem_matrix_1d[i*vecsize];
  }

  //fill our dummy elem-matrix:
  //This dummy matrix will be the same for every element. A real application
  //would form a different elem-matrix for each element.
  for(size_t i=0; i<elem_matrix_2d.size(); ++i) {
    double* row = elem_matrix_2d[i];
    if (i>=1) row[i-1] = -1;
    row[i] = 2;
    if (i<elem_matrix_2d.size()-1) row[i+1] = -1;

    elem_vector[i] = 1;
  }

  std::vector<int> eqn_indices(vecsize);

  for(size_t i=0; i<part_buckets.size(); ++i) {
    stk::mesh::Bucket::iterator
      b_iter = part_buckets[i]->begin(),
             b_end  = part_buckets[i]->end();
    for(; b_iter != b_end; ++b_iter) {
      stk::mesh::Entity& elem = *b_iter;
      rel = elem.relations(NODE_RANK);
      for(int j=0; rel.first != rel.second; ++rel.first, ++j) {
        node_ids[j] = rel.first->entity()->identifier();
      }

      matgraph->getPatternIndices(pattern_id, &node_ids[0], eqn_indices);

      matrix.sumIn(vecsize, &eqn_indices[0], vecsize, &eqn_indices[0],
                    &elem_matrix_2d[0]);
      rhs.sumIn(vecsize, &eqn_indices[0], &elem_vector[0]);
    }
  }
}


/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_linsys/FeiBaseIncludes.hpp>
#include <stk_linsys/LinsysFunctions.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <fei_trilinos_macros.hpp>
#include <fei_Solver_AztecOO.hpp>
#include <fei_Trilinos_Helpers.hpp>


namespace stk {
namespace linsys {

void add_connectivities(stk::linsys::LinearSystemInterface& ls,
                        stk::mesh::EntityRank entity_rank,
                        stk::mesh::EntityRank connected_entity_rank,
                        const stk::mesh::FieldBase& field,
                        const stk::mesh::Selector& selector,
                        const stk::mesh::BulkData& mesh_bulk)
{
  const std::vector<mesh::Bucket*>& mesh_buckets = mesh_bulk.buckets(entity_rank);
  std::vector<mesh::Bucket*> part_buckets;
  stk::mesh::get_buckets(selector, mesh_buckets, part_buckets);

  if (part_buckets.empty()) return;

  DofMapper& dof_mapper = ls.get_DofMapper();

  dof_mapper.add_dof_mappings(mesh_bulk, selector, connected_entity_rank, field);

  int field_id = dof_mapper.get_field_id(field);

  stk::mesh::Entity& first_entity = *(part_buckets[0]->begin());
  stk::mesh::PairIterRelation rel = first_entity.relations(connected_entity_rank);
  int num_connected = rel.second - rel.first;

  fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();

  int pattern_id = matgraph->definePattern(num_connected, connected_entity_rank, field_id);

  int num_entities = 0;
  for(size_t i=0; i<part_buckets.size(); ++i) {
    num_entities += part_buckets[i]->size();
  }

  int block_id = matgraph->initConnectivityBlock(num_entities, pattern_id);

  std::vector<int> connected_ids(num_connected);

  for(size_t i=0; i<part_buckets.size(); ++i) {
    stk::mesh::Bucket::iterator
      b_iter = part_buckets[i]->begin(),
      b_end  = part_buckets[i]->end();
    for(; b_iter != b_end; ++b_iter) {
      stk::mesh::Entity& entity = *b_iter;
      rel = entity.relations(connected_entity_rank);
      for(int j=0; rel.first != rel.second; ++rel.first, ++j) {
        connected_ids[j] = rel.first->entity()->identifier();
      }
      int conn_id = entity.identifier();
      matgraph->initConnectivity(block_id, conn_id, &connected_ids[0]);
    }
  }
}

void dirichlet_bc(stk::linsys::LinearSystemInterface& ls,
                  const stk::mesh::BulkData& mesh,
                  const stk::mesh::Part& bcpart,
                  stk::mesh::EntityRank entity_rank,
                  const stk::mesh::FieldBase& field,
                  unsigned field_component,
                  double prescribed_value)
{
  std::vector<stk::mesh::Bucket*> buckets;
  stk::mesh::Selector selector(bcpart);
  stk::mesh::get_buckets(selector, mesh.buckets(entity_rank), buckets);

  int field_id = ls.get_DofMapper().get_field_id(field);
  std::vector<int> entity_ids;

  for(size_t i=0; i<buckets.size(); ++i) {
    stk::mesh::Bucket::iterator
      iter = buckets[i]->begin(), iend = buckets[i]->end();
    for(; iter!=iend; ++iter) {
      const stk::mesh::Entity& entity = *iter;
      entity_ids.push_back(stk::linsys::impl::entityid_to_int(entity.identifier()));
    }
  }

  std::vector<double> prescribed_values(entity_ids.size(), prescribed_value);

  ls.get_fei_LinearSystem()->loadEssentialBCs(entity_ids.size(),
                                       &entity_ids[0],
                                       stk::linsys::impl::entitytype_to_int(entity_rank),
                                       field_id, field_component,
                                       &prescribed_values[0]);
}

int fei_solve(int & status, fei::LinearSystem &fei_ls, const Teuchos::ParameterList & params )
{
//Note: hard-coded to Trilinos/Aztec!!!
  Solver_AztecOO solver_az;

  fei::ParameterSet param_set;

  Trilinos_Helpers::copy_parameterlist(params, param_set);

  int iterationsTaken = 0;

  return solver_az.solve( & fei_ls,
                     NULL,
                     param_set,
                     iterationsTaken,
                     status
                  );
}

double compute_residual_norm2(fei::LinearSystem& fei_ls, fei::Vector& r)
{
  fei::SharedPtr<fei::Matrix> A = fei_ls.getMatrix();
  fei::SharedPtr<fei::Vector> x = fei_ls.getSolutionVector();
  fei::SharedPtr<fei::Vector> b = fei_ls.getRHS();

  //form r = A*x
  A->multiply(x.get(), &r);
  //form r = b - r   (i.e., r = b - A*x)
  r.update(1, b.get(), -1);

//!!!!!! fix this !!!!!!!!!
//terrible data copy. fei::Vector should provide a norm operation instead
//of making me roll my own here...

  std::vector<int> indices;
  r.getVectorSpace()->getIndices_Owned(indices);
  std::vector<double> coefs(indices.size());
  r.copyOut(indices.size(), &indices[0], &coefs[0]);
  double local_sum = 0;
  for(size_t i=0; i<indices.size(); ++i) {
    local_sum += coefs[i]*coefs[i];
  }
  MPI_Comm comm = r.getVectorSpace()->getCommunicator();
  double global_sum = 0;
  int num_doubles = 1;
  MPI_Allreduce(&local_sum, &global_sum, num_doubles, MPI_DOUBLE, MPI_SUM, comm);
  return std::sqrt(global_sum);
}

void copy_vector_to_mesh( fei::Vector & vec,
                          const DofMapper & dof,
                          stk::mesh::BulkData & mesh_bulk_data
                        )
{
  vec.scatterToOverlap();

  std::vector<int> shared_and_owned_indices;

  vec.getVectorSpace()->getIndices_SharedAndOwned(shared_and_owned_indices);

  size_t num_values = shared_and_owned_indices.size();

  if(num_values == 0) {
    return;
  }

  std::vector<double> values(num_values);
  vec.copyOut(num_values,&shared_and_owned_indices[0],&values[0]);

  stk::mesh::EntityRank ent_type;
  stk::mesh::EntityId ent_id;
  const stk::mesh::FieldBase * field;
  int offset_into_field;

  for(size_t i = 0; i < num_values; ++i)
  {

    dof.get_dof( shared_and_owned_indices[i],
                 ent_type,
                 ent_id,
                 field,
                 offset_into_field
                );

    stk::mesh::Entity & entity = *mesh_bulk_data.get_entity(ent_type, ent_id);

    void * data = stk::mesh::field_data(*field,entity);

    if(!(field->type_is<double>()) || data == NULL) {
      std::ostringstream oss;
      oss << "stk::linsys::copy_vector_to_mesh ERROR, bad data type, or ";
      oss << " field (" << field->name() << ") not found on entity with type " << entity.entity_rank();
      oss << " and ID " << entity.identifier();
      std::string str = oss.str();
      throw std::runtime_error(str.c_str());
    }

    double * double_data = reinterpret_cast<double *>(data);

    double_data[offset_into_field] = values[i];

  }
}

void scale_matrix(double scalar, fei::Matrix& matrix)
{
  fei::SharedPtr<fei::VectorSpace> vspace = matrix.getMatrixGraph()->getRowSpace();

  int numRows = vspace->getNumIndices_Owned();
  std::vector<int> rows(numRows);
  vspace->getIndices_Owned(numRows, &rows[0], numRows);

  std::vector<int> indices;
  std::vector<double> coefs;

  for(size_t i=0; i<rows.size(); ++i) {
    int rowlen = 0;
    matrix.getRowLength(rows[i], rowlen);

    if ((int)indices.size() < rowlen) {
      indices.resize(rowlen);
      coefs.resize(rowlen);
    }

    matrix.copyOutRow(rows[i], rowlen, &coefs[0], &indices[0]);

    for(int j=0; j<rowlen; ++j) {
      coefs[j] *= scalar;
    }

    double* coefPtr = &coefs[0];
    matrix.copyIn(1, &rows[i], rowlen, &indices[0], &coefPtr);
  }
}

void add_matrix_to_matrix(double scalar,
                          const fei::Matrix& src_matrix,
                          fei::Matrix& dest_matrix)
{
  fei::SharedPtr<fei::VectorSpace> vspace = src_matrix.getMatrixGraph()->getRowSpace();

  int numRows = vspace->getNumIndices_Owned();
  std::vector<int> rows(numRows);
  vspace->getIndices_Owned(numRows, &rows[0], numRows);

  std::vector<int> indices;
  std::vector<double> coefs;

  for(size_t i=0; i<rows.size(); ++i) {
    int rowlen = 0;
    src_matrix.getRowLength(rows[i], rowlen);

    if ((int)indices.size() < rowlen) {
      indices.resize(rowlen);
      coefs.resize(rowlen);
    }

    src_matrix.copyOutRow(rows[i], rowlen, &coefs[0], &indices[0]);

    for(int j=0; j<rowlen; ++j) {
      coefs[j] *= scalar;
    }

    double* coefPtr = &coefs[0];
    dest_matrix.sumIn(1, &rows[i], rowlen, &indices[0], &coefPtr);
  }
}

void scale_vector(double scalar,
                  fei::Vector& vec)
{
  fei::SharedPtr<fei::VectorSpace> vspace = vec.getVectorSpace();

  int numIndices = vspace->getNumIndices_Owned();
  std::vector<int> indices(numIndices);
  vspace->getIndices_Owned(numIndices, &indices[0], numIndices);

  std::vector<double> coefs(numIndices);

  vec.copyOut(numIndices, &indices[0], &coefs[0]);

  for(size_t j=0; j<coefs.size(); ++j) {
    coefs[j] *= scalar;
  }

  vec.copyIn(numIndices, &indices[0], &coefs[0]);
}

void add_vector_to_vector(double scalar,
                          const fei::Vector& src_vector,
                          fei::Vector& dest_vector)
{
  fei::SharedPtr<fei::VectorSpace> vspace = src_vector.getVectorSpace();

  int numIndices = vspace->getNumIndices_Owned();
  std::vector<int> indices(numIndices);
  vspace->getIndices_Owned(numIndices, &indices[0], numIndices);

  std::vector<double> coefs(numIndices);

  src_vector.copyOut(numIndices, &indices[0], &coefs[0]);

  for(size_t j=0; j<coefs.size(); ++j) {
    coefs[j] *= scalar;
  }

  dest_vector.sumIn(numIndices, &indices[0], &coefs[0]);
}

}//namespace linsys
}//namespace stk


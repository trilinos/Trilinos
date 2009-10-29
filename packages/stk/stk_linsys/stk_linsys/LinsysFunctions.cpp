
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

void add_connectivities(stk::linsys::LinearSystem& ls,
                        stk::mesh::EntityType entity_type,
                        stk::mesh::EntityType connected_entity_type,
                        const stk::mesh::FieldBase& field,
                        const stk::mesh::PartVector& part_intersection,
                        const stk::mesh::BulkData& mesh_bulk)
{
  const std::vector<mesh::Bucket*>& mesh_buckets = mesh_bulk.buckets(entity_type);
  stk::mesh::Selector selector(part_intersection);
  std::vector<mesh::Bucket*> part_buckets;
  stk::mesh::get_buckets(selector, mesh_buckets, part_buckets);

  if (part_buckets.empty()) return;

  DofMapper& dof_mapper = ls.get_DofMapper();

  dof_mapper.add_dof_mappings(mesh_bulk, part_intersection, connected_entity_type, field);

  int field_id = dof_mapper.get_field_id(field);

  stk::mesh::Entity& first_entity = *(part_buckets[0]->begin());
  stk::mesh::PairIterRelation rel = first_entity.relations(connected_entity_type);
  int num_connected = rel.second - rel.first;

  fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();

  int pattern_id = matgraph->definePattern(num_connected, connected_entity_type, field_id);

  int num_entities = 0;
  for(size_t i=0; i<part_buckets.size(); ++i) {
    num_entities += part_buckets[i]->size();
  }

  int part_ordinal_sum = 0;
  for(size_t i=0; i<part_intersection.size(); ++i) {
    part_ordinal_sum += part_intersection[i]->mesh_meta_data_ordinal();
  }

  int block_id = 1000*(part_ordinal_sum+1) + entity_type;

  matgraph->initConnectivityBlock(block_id, num_entities, pattern_id);

  std::vector<int> connected_ids(num_connected);

  for(size_t i=0; i<part_buckets.size(); ++i) {
    stk::mesh::Bucket::iterator
      b_iter = part_buckets[i]->begin(),
      b_end  = part_buckets[i]->end();
    for(; b_iter != b_end; ++b_iter) {
      stk::mesh::Entity& entity = *b_iter;
      rel = entity.relations(connected_entity_type);
      for(int j=0; rel.first != rel.second; ++rel.first, ++j) {
        connected_ids[j] = rel.first->entity()->identifier();
      }
      int conn_id = entity.identifier();
      matgraph->initConnectivity(block_id, conn_id, &connected_ids[0]);
    }
  }
}

void dirichlet_bc(stk::linsys::LinearSystem& ls,
                  const stk::mesh::BulkData& mesh,
                  const stk::mesh::Part& bcpart,
                  stk::mesh::EntityType entity_type,
                  const stk::mesh::FieldBase& field,
                  unsigned field_component,
                  double prescribed_value)
{
  std::vector<stk::mesh::Bucket*> buckets;
  stk::mesh::Selector selector(bcpart);
  stk::mesh::get_buckets(selector, mesh.buckets(entity_type), buckets);

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
                                       stk::linsys::impl::entitytype_to_int(entity_type),
                                       field_id, field_component,
                                       &prescribed_values[0]);
}

int fei_solve(int & status, fei::LinearSystem &fei_ls, const Teuchos::ParameterList & params )
{
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

  stk::mesh::EntityType ent_type;
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
      oss << " field (" << field->name() << ") not found on entity with type " << entity.entity_type();
      oss << " and ID " << entity.identifier();
      std::string str = oss.str();
      throw std::runtime_error(str.c_str());
    }

    double * double_data = reinterpret_cast<double *>(data);

    double_data[offset_into_field] = values[i];

  }
}

}//namespace linsys
}//namespace stk


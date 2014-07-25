#ifndef STK_SIERRA_MESH_MESH_CLASS_HPP
#define STK_SIERRA_MESH_MESH_CLASS_HPP

namespace sierra {
namespace mesh {

template <class Mesh>
class mesh_class
{
 public:
  mesh_class(const Mesh& input_mesh) : m_mesh(input_mesh) {}

  ~mesh_class() {}

  inline
  entity_descriptor
  target_entity( const relation_descriptor& relation ) const
  { return sierra::mesh::target_entity(relation, m_mesh); }
  
  inline
  bucket_entity_range
  get_entities(bucket_key bucket) const
  { return sierra::mesh::get_entities(bucket, m_mesh); }
  
  inline
  relation_range
  get_relations(const entity_descriptor& entity,
                const entity_rank& rank) const
  { return sierra::mesh::get_relations(entity, rank, m_mesh); }
  
  inline
  relation_range
  get_node_relations(const entity_descriptor& entity,
                     unsigned num_nodes) const
  { return sierra::mesh::get_node_relations(entity, num_nodes, m_mesh); }

  inline
  bucket_key
  get_bucket(const entity_descriptor& entity) const
  { return sierra::mesh::get_bucket(entity, m_mesh); }
  
  inline
  selected_bucket_range
  get_buckets( const selector & select) const
  { return sierra::mesh::get_buckets(select, m_mesh); }
  
  template<class Field>
  inline
  double*
  get_field_data(Field& field,
                 const entity_descriptor& entity) const
  { return sierra::mesh::get_field_data(field, entity, m_mesh); }
  
  inline
  bool
  is_selected(bucket_key bucket,
              const selector& select) const
  { return sierra::mesh::is_selected(bucket, select, m_mesh); }
  
 private:
  const Mesh& m_mesh;
};


} // namespace mesh
} // namespace sierra

#endif // STK_SIERRA_MESH_MESH_CLASS_HPP

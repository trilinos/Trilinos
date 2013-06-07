
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

namespace stk {

template <unsigned DIM> class STKMesh {
public :
  typedef DIM                         Dimension;
  typedef stk::mesh::Entity           Entity;
  typedef stk::mesh::EntityVec        EntityVec
  typedef stk::mesh::EntityKey        EntityKey;
  typedef std::set<EntityKey>      EntityKeySet;
  typedef std::pair<EntityKey, unsigned>   IdentProc
  typedef std::vector<IdentProc>     IdentProcVec;

  typedef stk::search::box::SphereBoundingBox<IdentProc,float,DIM> BoundingBox;

  STKMesh(const EntityVec                          &ent,
          const stk::mesh::FieldBase             &coord,
          const std::vector<stk::mesh::FieldBase*> &val);
  ~STKMesh();
  stk::ParallelMachine comm() {return m_comm;}
  unsigned      Keys(EntityKeySet &keys) const;

  BoundingBox boundingbox (const EntityKey Id, const double radius);

  void copy_entities(const IdentProcVec    &entities_to_copy,
                     const std::string         &transfer_name);
  
  void update_values();

  // Needed for LinearInterpoate
  const double *coord(const EntityKey k) const;
  const double *value(const EntityKey k, const unsigned i=0) const;
        double *value(const EntityKey k, const unsigned i=0);
  unsigned      value_size(const EntityKey e, const unsigned i=0) const;
  unsigned      num_values() const;

private :
  STKMesh (); 
  STKMesh(const STKMesh &M);
  STKMesh &operator=(const STKMesh&);

  stk::mesh::MetaData                        *m_meta_data;
  stk::mesh::BulkData                        *m_bulk_data;
  const stk::ParallelMachine                       m_comm;
  const EntityKeySet                        m_entity_keys;
  const stk::mesh::FieldBase         &m_coordinates_field;
  const std::vector<stk::mesh::FieldBase*> m_values_field;

  stk::mesh::Ghosting *transfer_entity_ghosting;

  Entity entity(const EntityKey k) const;
  static const EntityKeySet EntityKeys (const stk::mesh::EntityVec &ent);
};

template<unsigned DIM> const EntityKeySet entity_keys (const stk::mesh::BulkData  *bulk_data,
                                                       const stk::mesh::EntityVec &entities);
  EntityKeySet entity_keys;
  for (EntityVec::const_iterator e=entities.begin(); e!=entities.end(); ++e) {
    const stk::mesh::EntityKey k = bulk_data->entity_key(*e);
    entity_keys.insert(k);
  }
  return entity_keys;
}
template<unsigned DIM> STKMesh<DIM>::STKMesh(
          const stk::mesh::EntityVec          &entities,
          const stk::mesh::FieldBase             &coord,
          const std::vector<stk::mesh::FieldBase*> &val) :
    m_meta_data         (&coordinates_field.get_mesh()),
    m_bulk_data         (&stk::mesh::MetaData::get(m_bulk_data)),
    m_comm              (m_bulk_data->parallel()),
    m_entity_keys       (entity_keys(m_bulk_data, entities)), 
    m_coordinates_field (coord), 
    m_values_field      (val) {}

template<unsigned DIM> unsigned STKMesh<DIM>::keys(EntityKeySet &k) const {
  k = m_entity_keys;
  return k.size();
}

template<DIM> BoundingBox STKMesh<DIM>::boundingbox (const EntityKey Id, const double radius) {
  float center[DIM];
  const double *c = coord(Id);
  for (unsigned j=0; j<DIM; ++j) center[j] = c[j];
  const BoundingBox::Key key(Id, stk::parallel_machine_rank(m_comm));
  BoundingBox B(center, radius, key);
  return B;
}

template<unsigned NUM> void STKMesh<NUM>::copy_entities(
                     const IdentProcVec  &keys_to_copy,
                     const std::string         &transfer_name) {

  m_bulk_data->modification_begin();

  stk::mesh::EntityProcVec entities_to_copy(keys_to_copy.size());
  for (size_t i=0; i<keys_to_copy.size(); ++i) {
    stk::mesh::Entity entity = Entity(keys_to_copy[i].first);
    stk::mesh::EntityProc ep( entity, keys_to_copy[i].second);
    entities_to_copy[i] = ep;
  }
  const std::string name = "Transfer "+transfer_name+" Ghosting";
  stk::mesh::Ghosting transfer_entity_ghosting = &m_bulk_data->create_ghosting(name);
  {
    std::vector<stk::mesh::EntityKey> receive;
    transfer_entity_ghosting->receive_list( receive );
    m_bulk_data->change_ghosting(*transfer_entity_ghosting ,
                               entities_to_copy ,
                               receive );
  }
  m_bulk_data->modification_end();

  // Copy coordinates to the newly ghosted nodes
  std::vector<const stk::mesh::FieldBase *> fields;
  fields.push_back(&coordinates_field);
  fields.insert(fields.end(), m_values_field.begin(), m_values_field.end());
  
  stk::mesh::communicate_field_data( *transfer_entity_ghosting , fields);
}

template<DIM> BoundingBox STKMesh<DIM>::update_values () {
  stk::mesh::communicate_field_data( *transfer_entity_ghosting , fields);
}
  
template<DIM> const double *STKMesh<DIM>::coord(const EntityKey k) const {
  const stk::mesh::Entity e = entity(k);
  const double *c = static_cast<const double*>(stk::mesh::field_data(m_coordinates_field, e));
  return  c;
}

template<DIM> unsigned  STKMesh<DIM>::num_values() const {
 const unsigned s = m_values_field.size();
 return s;
}

template<DIM> unsigned  STKMesh<DIM>::value_size(const EntityKey k, const unsigned i) const {
  const stk::mesh::Entity         e = entity(k);
  const stk::mesh::FieldBase &field = *m_values_field[i];
  const stk::mesh::Bucket    &bucket= m_bulk_data->bucket(e);

  const unsigned bytes = m_bulk_data->field_data_size_per_entity(field, bucket);
  const unsigned bytes_per_entry = field.data_traits().size_of;
  const unsigned num_entry = bytes/bytes_per_entry;

  ThrowRequireMsg (bytes == num_entry * bytes_per_entry,    
    __FILE__<<":"<<__LINE__<<" Error:" <<"  bytes:" <<bytes<<"  num_entry:" <<num_entry
         <<"  bytes_per_entry:" <<bytes_per_entry);
  }
  return  num_entry;
}

temlate<DIM> const double *STKMesh<DIM>::value(const EntityKey k, const unsigned i) const {
  const stk::mesh::Entity  e = entity(k);
  stk::mesh::FieldBase *val=m_values_field[i];
  const double *value = static_cast<const double*>(stk::mesh::field_data(*val, e));
  return  value;
}

temlate<DIM> double *STKMesh<DIM>::value(const EntityKey k, const unsigned i) {
  const stk::mesh::Entity e = entity(k);
  stk::mesh::FieldBase *val=values_field[i];
  double *value = static_cast<double*>(stk::mesh::field_data(*val, e));
  return  value;
}

template<DIM> stk::mesh::Entity STKMesh<DIM>::entity(const stk::mesh::EntityKey k) const {
  const stk::mesh::Entity  e = m_bulk_data->get_entity(k);
  return e;
}




}

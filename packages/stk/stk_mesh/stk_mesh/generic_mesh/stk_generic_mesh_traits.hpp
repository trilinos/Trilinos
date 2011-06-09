#ifndef STK_MESH_GENERIC_MESH_TRAITS_HPP
#define STK_MESH_GENERIC_MESH_TRAITS_HPP

#include "stk_mesh/base/Types.hpp"

namespace stk {
namespace mesh {

struct Internal_STK_Relation_Descriptor_Adapter
{
  Internal_STK_Relation_Descriptor_Adapter(Entity * e_from, Entity * e_to, RelationIdentifier rel_id)
    : m_entity_from(e_from), m_entity_to(e_to), m_relation_id(rel_id) { }

  Entity * m_entity_from;
  Entity * m_entity_to;
  RelationIdentifier m_relation_id;

};

struct Internal_STK_Generic_Mesh_Adapter
{

  //****************************************
  // Generic Mesh API BEGIN
  //****************************************

  // persistant through modifications (up to change owner):
  //typedef typename Mesh::entity_global_id                entity_global_id;

  // persistant through modifications (up to change owner):
  typedef typename Entity *     entity_local_id;
  // not-persistant through modifications (up to change owner):
  typedef typename Entity *     entity_descriptor;
  // entity_descriptor_iterator de-references to an entity_descriptor:
  typedef typename Entity *     entity_descriptor_iterator;
  typedef typename Entity       entity_value;
  typedef typename EntityId     entity_size_type;


  typedef typename Part *       part_descriptor;
  // part_descriptor_iterator de-references to a part_descriptor:
  typedef typename Part *       part_descriptor_iterator;
  typedef typename Part         part_value;
  typedef typename PartOridnal  part_size_type; // not used yet

  typedef typename Bucket *     bucket_descriptor;
  // bucket_descriptor_iterator de-references to a bucket_descriptor:
  typedef typename Bucket *     bucket_descriptor_iterator;
  typedef typename Bucket       bucket_value; // not used yet
  typedef typename Ordinal      bucket_size_type; // not used yet
  // potentially heavy bucket_descriptor_iterator:
  typedef typename Mesh::selected_bucket_descriptor_iterator selected_bucket_descriptor_iterator;

  // potentially heavy bucket_descriptor_iterator:
  typedef typename Bucket *     part_bucket_descriptor_iterator;
  typedef typename Ordinal      part_bucket_size_type; // not used yet

  // potentially heavy part_descriptor_iterator:
  typedef typename Part *       bucket_part_descriptor_iterator;
  typedef typename Ordinal      bucket_part_size_type; // not used yet

  // potentially heavy entity_descriptor_iterator:
  typedef typename Entity *     bucket_entity_descriptor_iterator;
  typedef typename Ordinal      bucket_entity_size_type; // not used yet

  typedef typename Internal_STK_Relation_Descriptor_Adapter    relation_descriptor;
  typedef typename Internal_STK_Relation_Descriptor_Adapter *  relation_descriptor_iterator;
  typedef typename Internal_STK_Relation_Descriptor_Adapter *  selected_relation_descriptor_iterator
  typedef typename Relation          relation_value; // not used yet
  typedef typename Ordinal           relation_size_type; // not used yet

  inline part_descriptor             universal_part() { return & m_bulk_data.universal_part(); }

  static inline part_descriptor      null_part() { return NULL; };
  static inline bucket_descriptor    null_bucket() { return NULL; };
  static inline entity_local_id      null_entity() { return NULL; };
  static inline relation_descriptor  null_relation() { return NULL; };


  //****************************************
  // Generic Mesh API END
  //****************************************

  Internal_STK_Generic_Mesh_Adapter(BulkData & bulk_data) : m_bulk_data(bulk_data) {}

  // Data on the struct:
  BulkData & m_bulk_data;
};



// Use the type STKGenericMesh as the mesh type for STK so we can change this in the future.
typedef Internal_STK_Generic_Mesh_Adapter STKGenericMesh;


} // namespace mesh
} // namespace stk


#endif // STK_MESH_GENERIC_MESH_TRAITS_HPP

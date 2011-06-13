#ifndef STK_MESH_GENERIC_MESH_TRAITS_HPP
#define STK_MESH_GENERIC_MESH_TRAITS_HPP

#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/BulkData.hpp"

namespace stk {
namespace mesh {

struct Internal_STK_Relation_Descriptor_Adapter
{
  Internal_STK_Relation_Descriptor_Adapter(Entity * e_from = NULL, Entity * e_to = NULL, RelationIdentifier rel_id = InvalidOrdinal)
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
  typedef Entity *     entity_local_id;
  // not-persistant through modifications (up to change owner):
  typedef Entity *     entity_descriptor;
  // entity_descriptor_iterator de-references to an entity_descriptor:
  typedef Entity *     entity_descriptor_iterator;
  typedef Entity       entity_value;
  typedef EntityId     entity_size_type;


  typedef Part *       part_descriptor;
  // part_descriptor_iterator de-references to a part_descriptor:
  typedef Part *       part_descriptor_iterator;
  typedef Part         part_value;
  typedef PartOrdinal  part_size_type; // not used yet

  typedef Bucket *     bucket_descriptor;
  // bucket_descriptor_iterator de-references to a bucket_descriptor:
  typedef Bucket *     bucket_descriptor_iterator;
  typedef Bucket       bucket_value; // not used yet
  typedef Ordinal      bucket_size_type; // not used yet
  // potentially heavy bucket_descriptor_iterator:
  typedef Bucket *     selected_bucket_descriptor_iterator;

  // potentially heavy bucket_descriptor_iterator:
  typedef Bucket *     part_bucket_descriptor_iterator;
  typedef Ordinal      part_bucket_size_type; // not used yet

  // potentially heavy part_descriptor_iterator:
  typedef Part *       bucket_part_descriptor_iterator;
  typedef Ordinal      bucket_part_size_type; // not used yet

  // potentially heavy entity_descriptor_iterator:
  typedef Entity *     bucket_entity_descriptor_iterator;
  typedef Ordinal      bucket_entity_size_type; // not used yet

  typedef Internal_STK_Relation_Descriptor_Adapter    relation_descriptor;
  typedef Internal_STK_Relation_Descriptor_Adapter *  relation_descriptor_iterator;
  typedef Internal_STK_Relation_Descriptor_Adapter *  selected_relation_descriptor_iterator;
  typedef Relation          relation_value; // not used yet
  typedef Ordinal           relation_size_type; // not used yet

  inline part_descriptor             universal_part() { return & m_bulk_data.mesh_meta_data().universal_part(); }

  static inline part_descriptor      null_part() { return NULL; };
  static inline bucket_descriptor    null_bucket() { return NULL; };
  static inline entity_local_id      null_entity() { return NULL; };
  static inline relation_descriptor  null_relation() { return relation_descriptor(); };


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

#ifndef SIERRA_MESH_STKMESH_MESH_TRAITS_HPP
#define SIERRA_MESH_STKMESH_MESH_TRAITS_HPP

#include <stk_mesh/base/BulkData.hpp>

#include <sierra/mesh/mesh_traits.hpp>

#include <boost/iterator/transform_iterator.hpp>

namespace sierra {
namespace mesh {

// Specialization and adaptors for stkmesh mesh.

struct To_Ordinal : std::unary_function<const stk::mesh::Part*, unsigned>
{
  unsigned operator()(const stk::mesh::Part* part) const
  {
    return part->mesh_meta_data_ordinal();
  }
};

template<>
struct mesh_traits<stk::mesh::BulkData>
{
  typedef stk::mesh::Entity*                                                   entity_key;
  typedef stk::mesh::Entity*                                                   entity_descriptor;
  typedef stk::mesh::EntityId                                                  entity_property;
  typedef stk::mesh::EntityRank                                                entity_rank;
  typedef stk::mesh::BucketVectorEntityIteratorRange                           entity_descriptor_range;
  typedef stk::mesh::BucketPtrIterator                                         bucket_entity_iterator;
  typedef std::pair<bucket_entity_iterator, bucket_entity_iterator>            bucket_entity_range;

  typedef stk::mesh::Relation                                        relation_descriptor;
  typedef stk::mesh::RelationIdentifier                              relation_position;
  typedef stk::mesh::RelationVector::const_iterator           relation_iterator;
  typedef std::pair<relation_iterator, relation_iterator>            relation_range;

  typedef stk::mesh::Bucket*                           bucket_key;
  typedef std::vector<bucket_key>::const_iterator      bucket_iterator;
  typedef std::pair<bucket_iterator, bucket_iterator>  bucket_range;
  typedef unsigned                                     bucket_location;

  typedef stk::mesh::AllSelectedBucketsRange    selected_bucket_range;
  typedef stk::mesh::AllSelectedBucketsIterator selected_bucket_iterator;
  typedef stk::mesh::AllSelectedBucketsRange    part_bucket_range;

  typedef unsigned                                                                              part_key;
  typedef boost::transform_iterator<To_Ordinal, std::vector<stk::mesh::Part*>::const_iterator>  transform_iter;
  typedef std::pair<transform_iter, transform_iter>                                             part_range;
  typedef void*                                                                                 part_property;

  typedef const part_key* bucket_part_iterator;
  typedef std::pair<bucket_part_iterator, bucket_part_iterator>                                           bucket_part_range;

  typedef stk::mesh::Selector            selector;
};

} // namespace mesh
} // namespace sierra

#endif // SIERRA_MESH_STKMESH_MESH_TRAITS_HPP

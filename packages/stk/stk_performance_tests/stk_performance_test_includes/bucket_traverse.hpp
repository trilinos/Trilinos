#ifndef BucketTraverse_hpp
#define BucketTraverse_hpp

#include <sierra/mesh/mesh_traits.hpp>

namespace sierra {
namespace mesh {

template <class Mesh, class EntitySelector, class RelationSelector, class Visitor>
void bucket_traverse(const EntitySelector &entity_selector, const RelationSelector &relation_selector, Visitor &visitor, Mesh &mesh)
{
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_key bucket_key;
  typedef typename sierra::mesh::mesh_traits<Mesh>::selected_bucket_range selected_bucket_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_key bucket_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_entity_range bucket_entity_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_entity_iterator bucket_entity_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::selected_bucket_iterator selected_bucket_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_descriptor entity_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_descriptor_range entity_descriptor_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_range relation_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_descriptor relation_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_iterator relation_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_rank entity_rank;

  visitor.initialize(mesh);

  selected_bucket_range buckets = get_buckets(entity_selector, mesh);
  for (selected_bucket_iterator b_iter = boost::begin(buckets), b_end = boost::end(buckets);  b_iter != b_end; ++b_iter) {
    bucket_descriptor bucket = *b_iter;
    visitor.discover_bucket(bucket, mesh);

    bucket_entity_range entities = get_entities(*b_iter, mesh);
    for (bucket_entity_iterator ent_iter = boost::begin(entities), ent_end = boost::end(entities); ent_iter != ent_end; ++ent_iter) {
      entity_descriptor entity = *ent_iter;
      visitor.discover_entity(entity, mesh);

      relation_range entity_relations = sierra::mesh::get_relations(entity, relation_selector, mesh);
      visitor.examine_relation_range(entity_relations, mesh);

      for (relation_iterator rel_iter = boost::begin(entity_relations), rel_end = boost::end(entity_relations); rel_iter != rel_end; ++rel_iter) {
        relation_descriptor relation = *rel_iter;
        visitor.examine_relation(relation, mesh);
      }

      visitor.finish_entity(entity, mesh);
    }
    visitor.finish_bucket(bucket, mesh);
  }
  visitor.finalize(mesh);
}

template <class Mesh>
class BucketVisitor
{
public:
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_entity_iterator bucket_entity_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_entity_range bucket_entity_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_key bucket_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::bucket_key bucket_key;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_descriptor entity_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_descriptor_range entity_descriptor_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::entity_rank entity_rank;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_descriptor relation_descriptor;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_iterator relation_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::relation_range relation_range;
  typedef typename sierra::mesh::mesh_traits<Mesh>::selected_bucket_iterator selected_bucket_iterator;
  typedef typename sierra::mesh::mesh_traits<Mesh>::selected_bucket_range selected_bucket_range;

protected:
  virtual ~BucketVisitor()
  {}

public:
  void initialize(Mesh &/* mesh */)
  {}

  void finalize(Mesh &/* mesh */)
  {}

  void discover_bucket(bucket_descriptor /* bucket */, Mesh &/* mesh */)
  {}

  void finish_bucket(bucket_descriptor /* bucket */, Mesh &/* mesh */)
  {}


  void discover_entity(entity_descriptor /* entity */, Mesh &/* mesh */)
  {}

  void finish_entity(entity_descriptor /* entity */, Mesh &/* mesh */)
  {}


  void examine_relation_range(const relation_range & /* entity_relations */, Mesh &/* mesh */)
  {}


  void examine_relation(relation_descriptor /* relation */, Mesh &/* mesh */)
  {}
};

} // namespace mesh
} // namespace sierra

#endif // BucketTraverse_hpp

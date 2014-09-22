// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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

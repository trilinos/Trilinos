#ifndef STK_SIERRA_MESH_CSR_MESH_FACTORY_HPP
#define STK_SIERRA_MESH_CSR_MESH_FACTORY_HPP

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>
#include <sierra/mesh/csr/csr_mesh.hpp>

namespace sierra {
namespace mesh {

class csr_mesh_factory {
  public:
    static boost::shared_ptr<csr_mesh> create_from_modifiable(const modifiable_mesh& mmesh)
    {
      boost::shared_ptr<csr_mesh> csrptr(new csr_mesh);
      csr_mesh& csr = *csrptr;

      csr.m_num_entity_ranks = mmesh.m_num_entity_ranks;
      csr.m_entity_key_manager = mmesh.m_entity_key_manager;
      csr.m_bucket_key_manager = mmesh.m_bucket_key_manager;
      csr.m_part_key_manager = mmesh.m_part_key_manager;
      csr.m_bucket_map = mmesh.m_bucket_map;
      csr.m_part_to_buckets_vector = mmesh.m_part_to_buckets_vector;
      csr.m_bucket_to_parts_vector = mmesh.m_bucket_to_parts_vector;
      csr.m_part_property_vector = mmesh.m_part_property_vector;
      csr.m_entity_property_vector = mmesh.m_entity_property_vector;
      csr.m_sharing_procs = mmesh.m_sharing_procs;
      csr.m_offset_into_sharing_procs = mmesh.m_offset_into_sharing_procs;
      csr.m_universal_part = mmesh.m_universal_part;
      csr.m_locally_owned_part = mmesh.m_locally_owned_part;
      csr.m_shared_part = mmesh.m_shared_part;

     size_t n_entities = csr.m_entity_key_manager->size();
     csr.m_entity_descriptor_map.resize(n_entities);
     csr.m_entity_key_map.resize(n_entities);

     csr.m_entity_bucket_location.resize(n_entities);

     csr.m_entity_relations.resize(csr.m_num_entity_ranks);
     for(size_t i=0; i<csr.m_num_entity_ranks; ++i) {
       csr.m_entity_relations[i].resize( n_entities+1);
     }


     csr.m_relation_targets.reserve( mmesh.m_intrusive_relation_container.size());
     csr.m_relation_positions.reserve( mmesh.m_intrusive_relation_container.size());

     csr.m_bucket_begin.resize( mmesh.num_buckets()+1);

     {
       size_t entity_count=0;
       size_t bucket_count=0;
       BOOST_FOREACH( csr_mesh::bucket_key bucket, mmesh.get_buckets() ) {

         csr.m_bucket_begin[bucket_count++] = csr_mesh::bucket_entity_iterator(entity_count, csr_mesh::to_descriptor<csr_mesh::entity_descriptor>());

         BOOST_FOREACH( modifiable_mesh::entity_key key, mmesh.get_entities(bucket) ) {
           csr_mesh::entity_descriptor descriptor(entity_count++);
           csr.m_entity_bucket_location[descriptor] = mmesh.get_bucket_location(key);
           csr.m_entity_descriptor_map[key] = descriptor;
           csr.m_entity_key_map[descriptor] = key;
         }
       }
       //setup the past the end bucket value
       csr.m_bucket_begin[bucket_count] = csr_mesh::bucket_entity_iterator(entity_count, csr_mesh::to_descriptor<csr_mesh::entity_descriptor>());
     }

     //copy the relations
     for(size_t i=0; i<csr.m_num_entity_ranks; ++i)
     {
       const csr_mesh::entity_rank relation_rank(i);

       BOOST_FOREACH( csr_mesh::entity_descriptor descriptor, csr.get_entities()) {
         csr_mesh::entity_key key = csr.get_entity_key(descriptor);

         //csr_mesh::entity_rank current_entity_rank = mmesh[key].m_rank;

         csr_mesh::relation_iterator first(csr.m_relation_targets.size(),csr_mesh::to_descriptor<csr_mesh::relation_descriptor>());

         csr.m_entity_relations[relation_rank][descriptor] = first;

           BOOST_FOREACH( modifiable_mesh::relation_descriptor relation, mmesh.get_out_relations(key)) {
             csr_mesh::entity_key target_key = mmesh.target_entity(relation);
             csr_mesh::entity_rank target_rank = mmesh[target_key].m_rank;
             if ( relation_rank == target_rank) {
               csr.m_relation_targets.push_back(csr.get_entity_descriptor(target_key));
               csr.m_relation_positions.push_back(mmesh.position(relation));
             }
           }
       }

       // setup the past the end relation value
       csr_mesh::relation_iterator last(csr.m_relation_targets.size(),csr_mesh::to_descriptor<csr_mesh::relation_descriptor>());
       csr.m_entity_relations[relation_rank][n_entities] = last;
     }

      return csrptr;
    }

};

} // mesh
} // sierra

#endif // STK_SIERRA_MESH_CSR_MESH_FACTORY_HPP

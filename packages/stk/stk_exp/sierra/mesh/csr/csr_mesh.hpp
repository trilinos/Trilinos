#ifndef STK_SIERRA_MESH_CSR_MESH_HPP
#define STK_SIERRA_MESH_CSR_MESH_HPP

#include <boost/iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/unordered_map.hpp>

#include <sierra/mesh/details/entity_descriptor.hpp>
#include <sierra/mesh/details/entity_key.hpp>
#include <sierra/mesh/details/entity_rank.hpp>
#include <sierra/mesh/details/entity_property.hpp>
#include <sierra/mesh/details/part_property.hpp>
#include <sierra/mesh/details/relation_descriptor.hpp>
#include <sierra/mesh/details/relation_position.hpp>
#include <sierra/mesh/details/bucket_location.hpp>
#include <sierra/mesh/details/selected_buckets.hpp>
#include <sierra/mesh/details/selected_entities.hpp>

#include <algorithm>
#include <iostream>
#include <string>

namespace sierra {
namespace mesh {

class csr_mesh_factory;

class csr_mesh {
  public:
    typedef std::string                   part_identifier;
    typedef details::entity_key           entity_key;
    typedef details::entity_descriptor    entity_descriptor;
    typedef details::entity_property      entity_property;
    typedef details::part_property        part_property;
    typedef details::relation_position    relation_position;
    typedef details::relation_descriptor  relation_descriptor;
    typedef details::bucket_key           bucket_key;
    typedef details::bucket_location      bucket_location;
    typedef details::entity_rank          entity_rank;
    typedef details::part_key             part_key;
    typedef details::bucket_key_manager   bucket_key_manager;
    typedef details::part_key_manager     part_key_manager;
    typedef details::entity_key_manager   entity_key_manager;
    typedef details::selector             selector;

    typedef std::vector<relation_position>                relation_position_container;
    typedef std::vector<entity_descriptor>                entity_descriptor_container;
    typedef std::vector<entity_key>                       entity_key_container;
    typedef std::vector<bucket_location>                  bucket_location_container;

    typedef boost::counting_iterator<size_t>             base_iterator;

    template< class Descriptor>
    struct to_descriptor : std::unary_function<size_t,Descriptor>
    {
      Descriptor operator() (size_t base) const
      { return Descriptor(base); }
    };

    typedef boost::transform_iterator< to_descriptor<entity_descriptor>,
                                       base_iterator
                                     >                           bucket_entity_iterator;
    typedef std::pair<bucket_entity_iterator,bucket_entity_iterator> bucket_entity_range;
    typedef std::vector<bucket_entity_iterator>                    bucket_pool;

    typedef boost::transform_iterator< to_descriptor<relation_descriptor>,
                                       base_iterator
                                     >                     relation_iterator;
    typedef std::pair<relation_iterator,relation_iterator>   relation_range;
    typedef std::vector<std::vector<relation_iterator> >    Entityrelation_pool;

    typedef std::set< bucket_key >                          bucket_set;
    typedef std::set< part_key >                            part_set;

    typedef std::vector< bucket_set >                       part_bucket_vector;
    typedef std::vector< part_set >                         bucket_part_vector;

    typedef boost::unordered_map< part_set, bucket_key>      part_bucket_map;

    typedef std::vector<details::part_property>             part_properties;
    typedef std::vector<details::entity_property>           entity_properties;

    typedef std::vector<int>                                int_vector;
    typedef std::vector<int>::const_iterator                sharing_proc_iterator;
    typedef std::pair<sharing_proc_iterator,
                      sharing_proc_iterator>                sharing_proc_range;

    typedef bucket_key_manager::iterator                     bucket_iterator;
    typedef bucket_key_manager::iterator_range               bucket_range;

    typedef part_key_manager::iterator                       part_iterator;
    typedef part_key_manager::iterator_range                 part_range;

    typedef bucket_entity_range               entity_descriptor_range;

    typedef details::selected_bucket_iterator<csr_mesh> selected_bucket_iterator;
    typedef std::pair< selected_bucket_iterator, selected_bucket_iterator> selected_bucket_range;

    typedef details::selected_entity_iterator<csr_mesh> selected_entity_iterator;
    typedef std::pair< selected_entity_iterator, selected_entity_iterator> selected_entity_range;

    typedef std::pair< bucket_set::const_iterator, bucket_set::const_iterator >                                   part_bucket_range;
    typedef part_set::const_iterator bucket_part_iterator;
    typedef std::pair< bucket_part_iterator, bucket_part_iterator >                                       bucket_part_range;

    //-------------------------------------------------------------------------

    bucket_range get_buckets() const
    { return m_bucket_key_manager->range(); }

    part_range get_parts() const
    { return m_part_key_manager->range(); }

    entity_descriptor_range get_entities() const
    { return std::make_pair(m_bucket_begin.front(), m_bucket_begin.back()); }

    const bucket_entity_range get_entities( const bucket_key & bucket) const
    {
      return std::make_pair(m_bucket_begin[bucket],m_bucket_begin[bucket+1]);
    }

    const bucket_set & get_buckets( const part_key & part ) const {
      return (*m_part_to_buckets_vector)[part];
    }

    const part_set & get_parts( const bucket_key & bucket ) const {
      return (*m_bucket_to_parts_vector)[bucket];
    }

    const bucket_location & get_bucket_location( const entity_key & entity) const {
      return get_bucket_location( get_entity_descriptor(entity));
    }

    const bucket_location & get_bucket_location( const entity_descriptor & entity) const {
      return m_entity_bucket_location[entity];
    }

    const entity_descriptor get_entity_descriptor( const entity_key & entity) const {
      return m_entity_descriptor_map[entity];
    }

    const entity_key get_entity_key( const entity_descriptor & entity) const {
      return m_entity_key_map[entity];
    }

    const relation_range get_relations( const entity_key & entity, const entity_rank & rank ) const {
      return get_relations(get_entity_descriptor(entity),rank);
    }

    const relation_range get_relations( const entity_descriptor & entity, const entity_rank & rank) const {
      return std::make_pair(m_entity_relations[rank][entity],m_entity_relations[rank][1+entity]);
    }

    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    const entity_descriptor  target_entity( const relation_descriptor & rel ) const {
      return m_relation_targets[rel];
    }

#if 0
    const entity_descriptor  source( const relation_descriptor & rel ) const {
      return static_cast<const details::relation*>(rel)->source();
    }

    const relation_position  position( const relation_descriptor & rel ) const {
      return static_cast<const details::relation*>(rel)->position();
    }
#endif

    //-------------------------------------------------------------------------
  public:
    //-------------------------------------------------------------------------
    bool part_contains_bucket( const part_key & part, const bucket_key & bucket) const {
      return (*m_part_to_buckets_vector)[part].count(bucket) > 0;
    }

    bool bucket_contains_part( const part_key & part, const bucket_key & bucket) const {
      return (*m_bucket_to_parts_vector)[bucket].count(part) > 0;
    }
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    size_t num_entities() const
    { return m_entity_key_manager->size(); }

    size_t num_entities(const bucket_key & b) const
    { return m_bucket_begin[b+1] - m_bucket_begin[b]; }

    size_t num_relations( const entity_key & key, const entity_rank & rank) const
    { return num_relations( get_entity_descriptor(key),rank ); }

    size_t num_relations( const entity_descriptor & entity, const entity_rank & rank) const
    { return m_entity_relations[rank][1+entity] - m_entity_relations[rank][entity]; }

    size_t num_buckets() const
    { return m_bucket_key_manager->size(); }

    size_t num_parts() const
    { return m_part_key_manager->size(); }
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    const entity_property & operator[](const entity_key & key) const
    { return (*m_entity_property_vector)[key]; }

    const entity_property & operator[](const entity_descriptor & entity) const
    { return (*this)[get_entity_key(entity)]; }

    const part_property & operator[](const part_key & key) const
    { return (*m_part_property_vector)[key]; }

    const relation_position & operator[](const relation_descriptor rel) const {
      return m_relation_positions[rel];
    }
    //-------------------------------------------------------------------------
    part_key universal_part() const
    { return m_universal_part; }

    part_key locally_owned_part() const
    { return m_locally_owned_part; }

    part_key shared_part() const
    { return m_shared_part; }

    //-------------------------------------------------------------------------
    int get_owner_proc(const entity_key& key) const
    { return (*m_entity_property_vector)[key].m_proc; }

    int get_owner_proc(const entity_descriptor& ent_desc) const
    { return get_owner_proc(get_entity_key(ent_desc)); }

    sharing_proc_range get_sharing_procs(const entity_descriptor& ent_desc) const
    { return get_sharing_procs(get_entity_key(ent_desc)); }

    sharing_proc_range get_sharing_procs(const entity_key& key) const
    {
      sharing_proc_iterator begin = m_sharing_procs->begin()+(*m_offset_into_sharing_procs)[key];
      sharing_proc_iterator end   = m_sharing_procs->begin()+(*m_offset_into_sharing_procs)[key+1];
      return std::make_pair(begin,end);
    }
    //-------------------------------------------------------------------------

    // \TODO fix these! temporary placeholders to get IO working.
    entity_rank node_rank() const {return entity_rank(0);}

    entity_rank element_rank() const {return entity_rank(3);}

    //-------------------------------------------------------------------------
  private:

  public:
  csr_mesh()
      : m_num_entity_ranks()
      , m_entity_key_manager()
      , m_bucket_key_manager()
      , m_part_key_manager()
      , m_bucket_map()
      , m_part_to_buckets_vector()
      , m_bucket_to_parts_vector()
      , m_part_property_vector()
      , m_entity_property_vector()
      , m_sharing_procs()
      , m_offset_into_sharing_procs()
      , m_entity_descriptor_map()
      , m_entity_key_map()
      , m_entity_bucket_location()
      , m_relation_targets()
      , m_relation_positions()
      , m_bucket_begin()
      , m_entity_relations()
      , m_universal_part()
      , m_locally_owned_part()
      , m_shared_part()
   {
   }

  private:
    friend class csr_mesh_factory;

    size_t                                       m_num_entity_ranks;
    boost::shared_ptr<const entity_key_manager>  m_entity_key_manager;
    boost::shared_ptr<const bucket_key_manager>  m_bucket_key_manager;
    boost::shared_ptr<const part_key_manager>    m_part_key_manager;
    boost::shared_ptr<const part_bucket_map>     m_bucket_map;
    boost::shared_ptr<const part_bucket_vector>  m_part_to_buckets_vector;
    boost::shared_ptr<const bucket_part_vector>  m_bucket_to_parts_vector;
    boost::shared_ptr<const part_properties>     m_part_property_vector;
    boost::shared_ptr<const entity_properties>   m_entity_property_vector;
    boost::shared_ptr<const int_vector>          m_sharing_procs;
    boost::shared_ptr<const int_vector>          m_offset_into_sharing_procs;
    entity_descriptor_container                  m_entity_descriptor_map;
    entity_key_container                         m_entity_key_map;
    bucket_location_container                    m_entity_bucket_location;
    entity_descriptor_container                  m_relation_targets;
    relation_position_container                  m_relation_positions;
    bucket_pool                                  m_bucket_begin;
    Entityrelation_pool                          m_entity_relations;
    part_key                                     m_universal_part;
    part_key                                     m_locally_owned_part;
    part_key                                     m_shared_part;
};

} // mesh
} // sierra

#endif // STK_SIERRA_MESH_CSR_MESH_HPP

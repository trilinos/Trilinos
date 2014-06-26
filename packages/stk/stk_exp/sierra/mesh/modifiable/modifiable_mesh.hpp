#ifndef SIERRA_SIERRA_MESH_DETAILS_MODIFIABLE_MESH_HPP
#define SIERRA_SIERRA_MESH_DETAILS_MODIFIABLE_MESH_HPP

#include <sierra/mesh/details/entity_key.hpp>
#include <sierra/mesh/details/bucket_key.hpp>
#include <sierra/mesh/details/part_key.hpp>
#include <sierra/mesh/details/part_property.hpp>
#include <sierra/mesh/details/entity_property.hpp>
#include <sierra/mesh/details/intrusive_entity.hpp>
#include <sierra/mesh/details/intrusive_relation.hpp>
#include <sierra/mesh/details/intrusive_relation_descriptor.hpp>
#include <sierra/mesh/details/selected_buckets.hpp>

#include <boost/unordered_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>

#include <functional>
#include <vector>
#include <deque>
#include <set>

namespace sierra {
namespace mesh {

class csr_mesh_factory;


class modifiable_mesh {
  private:

    struct intrusive_bucket {
      details::intrusive_entityList list;
      intrusive_bucket() : list() {}
      intrusive_bucket(const intrusive_bucket &) : list() {}
      intrusive_bucket & operator=(const intrusive_bucket &) { list.clear(); return *this; }
    };

    typedef std::deque<intrusive_bucket>                    intrusive_bucket_pool;
    typedef std::deque<details::intrusive_entity>           intrusive_entity_pool;
    typedef std::deque<details::intrusive_relation>         intrusive_relation_pool;

    typedef std::vector<details::part_property>             part_properties;
    typedef std::vector<details::entity_property>           entity_properties;
    typedef std::vector<int>                                int_vector;

  public:
    //-------------------------------------------------------------------------
    typedef std::string                   part_identifier;

    typedef details::entity_key           entity_key;
    typedef details::entity_key           entity_descriptor;
    typedef details::entity_property      entity_property;
    typedef details::part_property        part_property;
    typedef details::relation_position    relation_position;
    typedef details::bucket_key           bucket_key;
    typedef details::bucket_location      bucket_location;
    typedef details::entity_rank          entity_rank;
    typedef details::part_key             part_key;
    typedef details::selector             selector;

    typedef details::bucket_key_manager bucket_key_manager;
    typedef details::entity_key_manager entity_key_manager;
    typedef details::part_key_manager   part_key_manager;

    typedef details::intrusive_entity_iterator              bucket_entity_iterator;
    typedef details::intrusive_entityRange                  bucket_entity_range;

    typedef bucket_key_manager::iterator                    bucket_iterator;
    typedef bucket_key_manager::iterator_range              bucket_range;

    typedef details::part_key_manager::iterator             part_iterator;
    typedef details::part_key_manager::iterator_range       part_range;

    typedef details::intrusive_relation_descriptor          relation_descriptor;

    typedef std::set< bucket_key >                          bucket_set;
    typedef std::set< part_key >                            part_set;

    typedef entity_key_manager::iterator_range              entity_descriptor_range;

    typedef details::out_relation_iterator         out_relation_iterator;
    typedef details::in_relation_iterator         in_relation_iterator;
    typedef details::out_relation_range            out_relation_range;
    typedef details::in_relation_range            in_relation_range;

    typedef details::selected_bucket_iterator<modifiable_mesh> selected_bucket_iterator;
    typedef std::pair< selected_bucket_iterator, selected_bucket_iterator > selected_bucket_range;
    typedef std::pair< bucket_set::const_iterator, bucket_set::const_iterator >                                           part_bucket_range;

    typedef part_set::const_iterator bucket_part_iterator;
    typedef std::pair< bucket_part_iterator, bucket_part_iterator >                                               bucket_part_range;

    typedef std::vector<int>::iterator             sharing_proc_iterator;
    typedef std::pair<sharing_proc_iterator,
                      sharing_proc_iterator>       sharing_proc_range;
    //-------------------------------------------------------------------------

 private:
    typedef std::vector< bucket_set >                       part_bucket_vector;
    typedef std::vector< part_set >                         bucket_part_vector;

    typedef boost::unordered_map< part_set, bucket_key>     part_bucket_map;

  public:
    //-------------------------------------------------------------------------
    bucket_range get_buckets() const
    { return m_bucket_key_manager->range(); }

    part_range get_parts() const
    { return m_part_key_manager->range(); }

    const entity_descriptor get_entity_descriptor( const entity_key & ent_key) const {
      return ent_key;
    }

    const entity_key get_entity_key( const entity_descriptor & ent_desc) const {
      return ent_desc;
    }

    entity_descriptor_range get_entities() const
    { return m_entity_key_manager->range(); }

    bucket_entity_range get_entities( const bucket_key & bucket) const
    { return details::range(m_bucket_container[bucket].list); }

    const bucket_set & get_buckets( const part_key & part ) const {
      return (*m_part_to_buckets_vector)[part];
    }

    const part_set & get_parts( const bucket_key & bucket ) const {
      return (*m_bucket_to_parts_vector)[bucket];
    }

    const bucket_location & get_bucket_location( const entity_key & entity) const {
      return m_ientity_container[entity].location();
    }

    const out_relation_range get_out_relations( const entity_key & entity) const {
      return details::range(m_ientity_container[entity].out_relations());
    }

    const in_relation_range get_in_relations( const entity_key & entity) const {
      return details::range(m_ientity_container[entity].in_relations());
    }
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    const entity_key  target_entity( const relation_descriptor & rel ) const {
      return static_cast<const details::intrusive_relation*>(rel)->target();
    }

    const entity_key  source( const relation_descriptor & rel ) const {
      return static_cast<const details::intrusive_relation*>(rel)->source();
    }

    const relation_position  position( const relation_descriptor & rel ) const {
      return static_cast<const details::intrusive_relation*>(rel)->position();
    }
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

  public:
    //-------------------------------------------------------------------------
    entity_key add_entity( const entity_property & prop = entity_property() )
    {
      entity_key key = m_entity_key_manager->create();

      //resize the container if necessary
      if ( m_ientity_container.size() <= static_cast<size_t>(key) ) {
        m_ientity_container.resize( static_cast<size_t>(key) +1);
        (*m_entity_property_vector).resize( static_cast<size_t>(key)+1);
        //m_offset_into_sharing_procs needs to have size num-keys+1:
        (*m_offset_into_sharing_procs).resize( static_cast<size_t>(key)+2);
      }

      m_ientity_container[key].set_key(key);
      (*m_entity_property_vector)[key] = prop;

      //we know that the newly-created key is at the end of the list,
      //here we are setting the last two entries in the m_offset_into_sharing_procs vector:
      (*m_offset_into_sharing_procs)[key] = m_sharing_procs->size();
      (*m_offset_into_sharing_procs)[key+1] = m_sharing_procs->size();

      add_entity_to_bucket( m_initial_bucket, key );

      change_entity_parts(key, &m_universal_part, &m_universal_part+1);

      if (prop.m_proc == m_local_proc) {
        change_entity_parts(key, &m_locally_owned_part, &m_locally_owned_part+1);
      }

      return key;
    }

    bool remove_entity( const entity_key & key )
    {
      bool result = m_entity_key_manager->remove(key);

      if ( result ) {
        //clear out the entities relations
        details::out_relation_set & out_relations = m_ientity_container[key].out_relations();
        details::in_relation_set & in_relations = m_ientity_container[key].in_relations();

        BOOST_FOREACH( details::intrusive_relation & rel, out_relations) {
          rel.reset();
        }
        BOOST_FOREACH( details::intrusive_relation & rel, in_relations) {
          rel.reset();
        }

        out_relations.clear();
        in_relations.clear();

        //set the entity's key and bucket_locations to invalid values
        m_ientity_container[key].set_key(entity_key());
        m_ientity_container[key].set_location(bucket_location());

        //clear property
        (*m_entity_property_vector)[key] = entity_property();
      }
      return result;
    }
    //-------------------------------------------------------------------------

  public:
    //-------------------------------------------------------------------------
    relation_descriptor add_relation(const entity_key & from,
                                    const entity_key & to,
                                    const relation_position & position)
    {
      details::to_intrusive_relation_descriptor convert;

      details::intrusive_relation_value rel(from,position,to);

      details::out_relation_set & out_relations = m_ientity_container[from].out_relations();

      // \TODO see if the relation already exists

      m_intrusive_relation_container.push_back(rel);

      details::intrusive_relation & intrusive_relation = m_intrusive_relation_container.back();

      // \TODO can do much better than push_front and sort can be combined with search for existance
      out_relations.push_front(intrusive_relation);
      out_relations.sort(details::intrusive_out_relation_order());

      details::in_relation_set & in_relations = m_ientity_container[to].in_relations();

      // \TODO can do much better than push_front and sort
      in_relations.push_front(intrusive_relation);
      in_relations.sort(details::intrusive_in_relation_order());

      return convert(intrusive_relation);
    }

    bool remove_relation( const relation_descriptor & rel_descriptor )
    {
      details::intrusive_relation & rel = *const_cast<details::intrusive_relation*>(static_cast<const details::intrusive_relation*>(rel_descriptor));

      if (rel.source() != entity_key::invalid()) {
        details::out_relation_set & out_relations = m_ientity_container[rel.source()].out_relations();
        details::in_relation_set & in_relations = m_ientity_container[rel.target()].in_relations();

        out_relations.erase(out_relations.iterator_to(rel));
        in_relations.erase(in_relations.iterator_to(rel));

        rel.reset();

        return true;
      }

      return false;
    }
    //-------------------------------------------------------------------------

  private:
    //-------------------------------------------------------------------------
    bucket_key add_bucket()
    {
      bucket_key bucket = m_bucket_key_manager->create();

      //resize the container if necessary
      if ( m_bucket_container.size() <= static_cast<size_t>(bucket) ) {
        m_bucket_container.resize( static_cast<size_t>(bucket) +1);
        (*m_bucket_to_parts_vector).resize( static_cast<size_t>(bucket)+1);
      }

      return bucket;
    }

    bool remove_bucket( const bucket_key & bucket )
    {
      const bool result = m_bucket_container[bucket].list.empty();

      if ( result ) {
        m_bucket_key_manager->remove(bucket);
        m_bucket_container[bucket].list.clear();
        part_set tmp;
        (*m_bucket_to_parts_vector)[bucket].swap(tmp);
      }
      return result;
    }
    //-------------------------------------------------------------------------

  public:
    //-------------------------------------------------------------------------
    part_key declare_part( const part_identifier & part_name )
    {
      part_properties::const_iterator i =  (*m_part_property_vector).begin();
      const part_properties::const_iterator end =  (*m_part_property_vector).end();

      for( ; i != end; ++i) {
        if( i->name() == part_name ) return i->key();
      }

      // create a new part with the given name
      part_key key = m_part_key_manager->create();

      //resize the container if necessary
      if ( (*m_part_to_buckets_vector).size() <= static_cast<size_t>(key) ) {
        (*m_part_to_buckets_vector).resize( static_cast<size_t>(key)+1);
        (*m_part_property_vector).resize( static_cast<size_t>(key)+1);
      }

      (*m_part_property_vector)[key] = part_property(part_name,key);

      return key;
    }

    part_key find_part( const part_identifier & part_name ) const
    {
      part_properties::const_iterator i =  (*m_part_property_vector).begin();
      const part_properties::const_iterator end =  (*m_part_property_vector).end();

      for( ; i != end; ++i) {
        if( i->name() == part_name ) return i->key();
      }

      //didn't find it, return invalid (default-constructed) key
      return part_key();
    }

    template<class T>
    part_key declare_part( const part_identifier & part_name, const T & property)
    {
      part_key key = declare_part(part_name);
      (*m_part_property_vector)[key].add_property(property);
      return key;
    }

    //-------------------------------------------------------------------------

  public:
    //-------------------------------------------------------------------------
    size_t num_entities() const
    { return m_entity_key_manager->size(); }

    size_t num_entities(const bucket_key & b) const
    { return m_bucket_container[b].list.size(); }

    size_t num_out_relations( const entity_key & key) const
    { return m_ientity_container[key].out_relations().size(); }

    size_t num_in_relations( const entity_key & key) const
    { return m_ientity_container[key].in_relations().size(); }

    size_t num_buckets() const
    { return m_bucket_key_manager->size(); }

    size_t num_parts() const
    { return m_part_key_manager->size(); }

    size_t get_num_entity_ranks() const
    { return m_num_entity_ranks; }

    void set_num_entity_ranks(size_t num_ranks)
    { m_num_entity_ranks = num_ranks; }

    //-------------------------------------------------------------------------
    part_key universal_part() const
    { return m_universal_part; }

    part_key locally_owned_part() const
    { return m_locally_owned_part; }

    part_key shared_part() const
    { return m_shared_part; }

    //-------------------------------------------------------------------------
    void set_owner_proc(const entity_key& key, int proc)
    { (*m_entity_property_vector)[key].m_proc = proc; }

    int get_owner_proc(const entity_key& key) const
    { return (*m_entity_property_vector)[key].m_proc; }

    template<class procInputIterator>
    void add_sharing_procs(const entity_key& key,
                           procInputIterator proc_begin,
                           procInputIterator proc_end)
    {
      size_t num_inserted = 0;
      for(procInputIterator input_it=proc_begin; input_it!=proc_end; ++input_it) {
        sharing_proc_range sharing_procs = get_sharing_procs(key);
        bool already_present = false;
        sharing_proc_iterator insert_pos = sharing_procs.first;
        for(sharing_proc_iterator it=sharing_procs.first; it!=sharing_procs.second; ++it) {
          if (*it == *input_it) {
            already_present = true;
            break;
          }
         insert_pos = it;
          if (*it > *input_it) break;
        }

        if (!already_present) {
          (*m_sharing_procs).insert(insert_pos, *input_it);
          ++num_inserted;
        }
      }

      //now adjust offsets to account for the sharing-procs we just inserted:
      for(size_t i=key+1; i<m_offset_into_sharing_procs->size(); ++i) {
        (*m_offset_into_sharing_procs)[i] += num_inserted;
      }
    }

    sharing_proc_range get_sharing_procs(const entity_key& key)
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
    bucket_entity_range operator[](const bucket_key & b) const
    { return details::range(m_bucket_container[b].list); }

    const entity_property & operator[](const entity_key & key) const
    { return (*m_entity_property_vector)[key]; }

    const part_property & operator[](const part_key & key) const
    { return (*m_part_property_vector)[key]; }

    part_property& get(const part_key& key)
    { return (*m_part_property_vector)[key]; }

    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    template <class part_keyInputIterator>
    void change_entity_parts( const entity_key & key,
                         part_keyInputIterator add_first,
                         part_keyInputIterator add_last
                        )
    {
      change_entity_parts(key, add_first, add_last, add_last, add_last);
    }

    template <class Addpart_keyInputIterator, class Removepart_keyInputIterator>
    void change_entity_parts( const entity_key & entity,
                         Addpart_keyInputIterator add_first,
                         Addpart_keyInputIterator add_last,
                         Removepart_keyInputIterator remove_first,
                         Removepart_keyInputIterator remove_last
                        )
    {
      details::intrusive_entity & ientity = m_ientity_container[entity];
      part_set entity_parts = (*m_bucket_to_parts_vector)[ientity.location().bucket()];

      entity_parts.insert(add_first,add_last);
      for(Removepart_keyInputIterator it=remove_first, end=remove_last; it!=end; ++it) {
        entity_parts.erase(*it);
      }

      const part_set & prev_entity_parts = (*m_bucket_to_parts_vector)[ientity.location().bucket()];
      if (entity_parts.size() == prev_entity_parts.size()) {
        if (std::equal(entity_parts.begin(), entity_parts.end(), prev_entity_parts.begin())) {
          //entity_parts is the same as prev_entity_parts, so nothing to do.
          return;
        }
      }

      part_bucket_map::iterator itr = (*m_bucket_map).find(entity_parts);

      bucket_key bucket;
      if ( itr == (*m_bucket_map).end() ) { // bucket not found
        bucket = add_bucket();
        //add to the bucket map
        (*m_bucket_map)[entity_parts] = bucket;

        //store the entity_parts on the new bucket
        (*m_bucket_to_parts_vector)[bucket] = entity_parts;

        //add the bucket to each part
        BOOST_FOREACH( part_key part, entity_parts) {
          (*m_part_to_buckets_vector)[part].insert(bucket);
        }
      }
      else {
        bucket = itr->second;
      }

      move_entity_to_bucket(bucket,entity);
    }

    //-------------------------------------------------------------------------

  private:
    //-------------------------------------------------------------------------
    const bucket_location
    move_entity_to_bucket (
        const bucket_key  & bucket,
        const entity_key  & entity
          )
      {
        bucket_key current_bucket_key = get_bucket_location(entity).bucket();

        remove_entity_from_bucket( current_bucket_key, entity );

        return add_entity_to_bucket( bucket, entity );
      }

    const bucket_location
      add_entity_to_bucket (
          const bucket_key  & bucket_id,
          const entity_key  & entity
          )
      {
        details::intrusive_entityList & bucket = m_bucket_container[bucket_id].list;

        details::intrusive_entity & ientity = m_ientity_container[entity];

        ientity.set_location(bucket_location(bucket_id, details::bucket_ordinal(bucket.size())));

        bucket.push_back(ientity);

        return ientity.location();
      }

    void remove_entity_from_bucket (
        const bucket_key & bucket_id,
        const entity_key & entity
        )
    {
      details::intrusive_entityList & bucket = m_bucket_container[bucket_id].list;
      details::intrusive_entity & ientity = m_ientity_container[entity];

      // move the back of the bucket to the current location
      details::intrusive_entity & back = bucket.back();
      back.set_location(ientity.location());
      back.swap_nodes(ientity);

      bucket.pop_back();
    }
    //-------------------------------------------------------------------------

  public:
    //-------------------------------------------------------------------------
    modifiable_mesh(
        int local_proc = 0,
        size_t num_entity_ranks = 4,
        const boost::shared_ptr<entity_key_manager>    & arg_entity_key_manager =
        boost::shared_ptr<entity_key_manager>(new entity_key_manager),
        const boost::shared_ptr<bucket_key_manager>    & arg_bucket_key_manager =
        boost::shared_ptr<bucket_key_manager>(new bucket_key_manager),
        const boost::shared_ptr<part_key_manager>      & arg_part_key_manager =
        boost::shared_ptr<part_key_manager>(new part_key_manager),
        const boost::shared_ptr<part_bucket_map>           & arg_bucket_map =
        boost::shared_ptr<part_bucket_map>(new part_bucket_map),
        const boost::shared_ptr<part_bucket_vector>    & arg_part_to_buckets_vector =
        boost::shared_ptr<part_bucket_vector>(new part_bucket_vector),
        const boost::shared_ptr<bucket_part_vector>    & arg_bucket_to_parts_vector =
        boost::shared_ptr<bucket_part_vector>(new bucket_part_vector),
        const boost::shared_ptr<part_properties>  & arg_part_property_vector =
        boost::shared_ptr<part_properties>(new part_properties),
        const boost::shared_ptr<entity_properties>  & arg_entity_property_vector =
        boost::shared_ptr<entity_properties>(new entity_properties)
        )
      : m_local_proc(local_proc)
      , m_num_entity_ranks(num_entity_ranks)
      , m_entity_key_manager(arg_entity_key_manager)
      , m_bucket_key_manager(arg_bucket_key_manager)
      , m_part_key_manager(arg_part_key_manager)
      , m_bucket_map(arg_bucket_map)
      , m_part_to_buckets_vector(arg_part_to_buckets_vector)
      , m_bucket_to_parts_vector(arg_bucket_to_parts_vector)
      , m_part_property_vector(arg_part_property_vector)
      , m_entity_property_vector(arg_entity_property_vector)
      , m_sharing_procs(new int_vector)
      , m_offset_into_sharing_procs(new int_vector(1,0))
      , m_ientity_container()
      , m_intrusive_relation_container()
      , m_bucket_container()
      , m_universal_part()
      , m_locally_owned_part()
      , m_shared_part()
      , m_initial_bucket()
    {
      part_set empty_key;
      m_initial_bucket = add_bucket();
      (*m_bucket_map)[empty_key] = m_initial_bucket;
      (*m_bucket_to_parts_vector).push_back(empty_key);

      m_universal_part =     declare_part("SIERRA::MESH::UNIVERSAL_PART");
      m_locally_owned_part = declare_part("SIERRA::MESH::LOCALLY_OWNED_PART");
      m_shared_part =        declare_part("SIERRA::MESH::SHARED_PART");
    }

    ~modifiable_mesh()
    {}
    //-------------------------------------------------------------------------

  private:
    int                                       m_local_proc;
    size_t                                    m_num_entity_ranks;
    boost::shared_ptr<entity_key_manager>       m_entity_key_manager;
    boost::shared_ptr<bucket_key_manager>       m_bucket_key_manager;
    boost::shared_ptr<part_key_manager>         m_part_key_manager;
    boost::shared_ptr<part_bucket_map>              m_bucket_map;
    boost::shared_ptr<part_bucket_vector>       m_part_to_buckets_vector;
    boost::shared_ptr<bucket_part_vector>       m_bucket_to_parts_vector;
    boost::shared_ptr<part_properties>     m_part_property_vector;
    boost::shared_ptr<entity_properties>   m_entity_property_vector;
    boost::shared_ptr<int_vector>          m_sharing_procs;
    boost::shared_ptr<int_vector>          m_offset_into_sharing_procs;
    intrusive_entity_pool                  m_ientity_container;
    intrusive_relation_pool                m_intrusive_relation_container;
    intrusive_bucket_pool                      m_bucket_container;
    part_key                                   m_universal_part;
    part_key                                   m_locally_owned_part;
    part_key                                   m_shared_part;
    bucket_key                                 m_initial_bucket;

    friend class sierra::mesh::csr_mesh_factory;

    //disable copy constructor and assignment operator
    modifiable_mesh( const modifiable_mesh &);
    modifiable_mesh & operator=(const modifiable_mesh);
};

} // mesh
} // sierra

#endif // SIERRA_SIERRA_MESH_DETAILS_MODIFIABLE_MESH_HPP

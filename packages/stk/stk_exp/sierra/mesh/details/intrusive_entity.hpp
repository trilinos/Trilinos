#ifndef SIERRA_SIERRA_MESH_DETAILS_INTRUSIVE_ENTITY_HPP
#define SIERRA_SIERRA_MESH_DETAILS_INTRUSIVE_ENTITY_HPP


#include <sierra/mesh/details/entity_key.hpp>
#include <sierra/mesh/details/bucket_location.hpp>
#include <sierra/mesh/details/intrusive_relation.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/intrusive/list.hpp>

namespace sierra {
namespace mesh {
namespace details {

namespace bi = boost::intrusive;

//-------------------------------------------------------------------------
typedef bi::list_base_hook< bi::link_mode<bi::normal_link> > ListHook;

class intrusive_entity : public ListHook
{
  public:

    explicit
      intrusive_entity(const entity_key & key = entity_key(),
                      const bucket_location & des = bucket_location()
                     )
      : m_key(key)
      , m_location(des)
      , m_out_relations()
      , m_in_relations()
      {}

    intrusive_entity( const intrusive_entity & rhs)
      : m_key(rhs.m_key)
      , m_location(rhs.m_location)
      , m_out_relations()
      , m_in_relations()
    {}

    intrusive_entity & operator=(const intrusive_entity & rhs)
    {
      if (this != &rhs) {
        m_key = rhs.m_key;
        m_location = rhs.m_location;
        m_out_relations.clear();
        m_in_relations.clear();
      }
      return *this;
    }

    entity_key key() const { return m_key; }
    const bucket_location & location() const { return m_location; }

    void set_key( const entity_key & key ) {
      m_key = key;
    }

    void set_location( const bucket_location & des ) {
      m_location = des;
    }

          out_relation_set & out_relations()       { return m_out_relations; }
    const out_relation_set & out_relations() const { return m_out_relations; }

          in_relation_set & in_relations()       { return m_in_relations; }
    const in_relation_set & in_relations() const { return m_in_relations; }

  private:
    entity_key         m_key;
    bucket_location    m_location;
    out_relation_set    m_out_relations;
    in_relation_set     m_in_relations;


};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
struct to_entity_key : std::unary_function<intrusive_entity,entity_key>
{
  entity_key operator() (const intrusive_entity & entity) const
  { return entity.key(); }
};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
typedef boost::intrusive::list< intrusive_entity,
                                bi::base_hook<ListHook>
                              > intrusive_entityList;


typedef boost::transform_iterator< to_entity_key,
                                   intrusive_entityList::const_iterator
                                 > intrusive_entity_iterator;


typedef std::pair<intrusive_entity_iterator,intrusive_entity_iterator>
          intrusive_entityRange;


inline
const intrusive_entityRange range(const intrusive_entityList & list)
{
  return std::make_pair(intrusive_entity_iterator(list.begin(), to_entity_key()) ,
                        intrusive_entity_iterator(list.end(),   to_entity_key()) );
}
//-------------------------------------------------------------------------

} // details
} // mesh
} // sierra

#endif //SIERRA_SIERRA_MESH_DETAILS_INTRUSIVE_ENTITY_HPP

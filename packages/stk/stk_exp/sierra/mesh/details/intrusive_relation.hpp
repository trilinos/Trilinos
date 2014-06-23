#ifndef SIERRA_SIERRA_MESH_DETAILS_INTRUSIVE_RELATION_HPP
#define SIERRA_SIERRA_MESH_DETAILS_INTRUSIVE_RELATION_HPP

// \TODO Template relation on position concept
#include <sierra/mesh/details/entity_key.hpp>
#include <sierra/mesh/details/relation_position.hpp>
#include <sierra/mesh/details/intrusive_relation_descriptor.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/intrusive/slist.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <functional>

namespace sierra {
namespace mesh {
namespace details {

namespace bi = boost::intrusive;
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
struct out_relations_tag{};
struct in_relations_tag{};

typedef bi::slist_base_hook <
                             bi::tag<out_relations_tag>,
                             bi::link_mode<bi::normal_link>
                          > out_hook;
typedef bi::slist_base_hook <
                             bi::tag<in_relations_tag>,
                             bi::link_mode<bi::normal_link>
                          > in_hook;

typedef boost::tuple<entity_key,relation_position,entity_key> intrusive_relation_value;
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
class intrusive_relation
  : public out_hook, public in_hook
{
  public:

    intrusive_relation( const entity_key & arg_source,
                       const entity_key & arg_target,
                       const relation_position & arg_position )
      : m_value(arg_source,arg_position,arg_target)
    {}

    intrusive_relation( const intrusive_relation_value & arg_value)
      : m_value(arg_value)
    {}

    const entity_key & source() const { return m_value.get<0>(); }
    const relation_position & position() const { return m_value.get<1>(); }
    const entity_key & target() const { return m_value.get<2>(); }

    operator intrusive_relation_value() const { return m_value; }

    void reset() {
      intrusive_relation_value temp;
      m_value = temp;
    }


  private:
    intrusive_relation_value m_value;

};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
struct intrusive_out_relation_order: std::binary_function<intrusive_relation_value,intrusive_relation_value,bool>
{
  bool operator()(const intrusive_relation_value & a, const intrusive_relation_value & b) const
  { return a < b; }
};

struct intrusive_in_relation_order: std::binary_function<intrusive_relation_value,intrusive_relation_value,bool>
{
  bool operator()(const intrusive_relation_value & a, const intrusive_relation_value & b) const
  {
    if(a.get<2>() < b.get<2>()) return true;
    if(b.get<2>() < a.get<2>()) return false;

    if(a.get<1>() < b.get<1>()) return true;
    if(b.get<1>() < a.get<1>()) return false;

    return a.get<0>() < b.get<0>();
  }
};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
typedef bi::slist<
                intrusive_relation,
                bi::base_hook<out_hook>,
                bi::constant_time_size<false>
               > out_relation_set;

typedef bi::slist<
                intrusive_relation,
                bi::base_hook<in_hook>,
                bi::constant_time_size<false>
               > in_relation_set;
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
struct to_intrusive_relation_descriptor
  : std::unary_function<intrusive_relation,intrusive_relation_descriptor>
{
  intrusive_relation_descriptor operator() (const intrusive_relation & rel) const
  { return intrusive_relation_descriptor(&rel); }
};
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
typedef boost::transform_iterator< to_intrusive_relation_descriptor,
                                   out_relation_set::const_iterator
                                 > out_relation_iterator;


typedef std::pair<out_relation_iterator,out_relation_iterator>
          out_relation_range;


inline
const out_relation_range range(const out_relation_set & set)
{
  return std::make_pair(out_relation_iterator(set.begin(), to_intrusive_relation_descriptor()) ,
                        out_relation_iterator(set.end(),   to_intrusive_relation_descriptor()) );
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
typedef boost::transform_iterator< to_intrusive_relation_descriptor,
                                   in_relation_set::const_iterator
                                 > in_relation_iterator;


typedef std::pair<in_relation_iterator,in_relation_iterator>
          in_relation_range;


inline
const in_relation_range range(const in_relation_set & set)
{
  return std::make_pair(in_relation_iterator(set.begin(), to_intrusive_relation_descriptor()) ,
                        in_relation_iterator(set.end(),   to_intrusive_relation_descriptor()) );
}
//-------------------------------------------------------------------------

} // details
} // mesh
} // sierra

#endif //SIERRA_SIERRA_MESH_INTRUSIVE_RELATION_HPP`

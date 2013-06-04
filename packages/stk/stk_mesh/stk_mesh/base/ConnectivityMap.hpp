#ifndef STK_STK_CONNECTIVITY_MAP_HPP
#define STK_STK_CONNECTIVITY_MAP_HPP

#include <stk_mesh/base/Types.hpp>

#include <boost/array.hpp>

namespace stk  {
namespace mesh {

// TODO - enforce invariants (such as only makes sense to have fixed connectivity to an element from other elements)

struct ConnectivityMap
{
  typedef boost::array<boost::array<ConnectivityType, 4>, 4> map_type;

  static ConnectivityMap const& none()
  {
    static const map_type map =
    {{
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& none_2d()
  {
    static const map_type map =
    {{
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& all_dynamic()
  {
    static const map_type map =
    {{
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& all_dynamic_2d()
  {
    static const map_type map =
    {{
      {{dynamic(),dynamic(),invalid(),dynamic()}},
      {{dynamic(),dynamic(),invalid(),dynamic()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{dynamic(),dynamic(),invalid(),dynamic()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& default_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  dynamic(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  dynamic(), dynamic()}},
      /*face*/  {{ fixed()  , dynamic(),  invalid(), dynamic()}},
      /*elem*/  {{ fixed()  , dynamic(),  dynamic(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& minimal_back_relations_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), invalid(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), invalid()}},
      /*face*/  {{ fixed()  , dynamic(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , dynamic(),  dynamic(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& default_map_2d()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), dynamic()}},
      /*face*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , dynamic(),  invalid(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& fixed_edges_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  dynamic(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  dynamic(), dynamic()}},
      /*face*/  {{ fixed()  , fixed()  ,  invalid(), dynamic()}},
      /*elem*/  {{ fixed()  , fixed()  ,  dynamic(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& fixed_edges_map_2d()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), dynamic()}},
      /*face*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , fixed()  ,  invalid(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& fixed_downward_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), dynamic(),  dynamic(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  dynamic(), dynamic()}},
      /*face*/  {{ fixed()  , fixed()  ,  invalid(), dynamic()}},
      /*elem*/  {{ fixed()  , fixed()  ,  fixed()  , invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityMap const& element_node()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), invalid(),  invalid(), dynamic()}},
      /*edge*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*face*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , invalid(),  invalid(), invalid()}}
    }};
    static ConnectivityMap r = {map};
    return r;
  }

  static ConnectivityType invalid() { return INVALID_CONNECTIVITY_TYPE; }
  static ConnectivityType fixed()   { return FIXED_CONNECTIVITY;   }
  static ConnectivityType dynamic() { return DYNAMIC_CONNECTIVITY; }

  ConnectivityType & operator()(EntityRank from, EntityRank to)
  { return m_map[from][to]; }

  ConnectivityType operator()(EntityRank from, EntityRank to) const
  { return m_map[from][to]; }

  map_type m_map;
};

} //namespace mesh
} //namespace stk

#endif //STK_STK_CONNECTIVITY_MAP_HPP

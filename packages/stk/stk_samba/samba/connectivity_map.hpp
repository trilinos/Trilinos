#ifndef SAMBA_SAMBA_CONNECTIVITY_MAP_HPP
#define SAMBA_SAMBA_CONNECTIVITY_MAP_HPP

#include <samba/entity_rank.hpp>

#include <samba/connectivity_kind.hpp>
#include <samba/spatial_dimension.hpp>

#include <boost/array.hpp>

namespace samba {

// TODO - enforce invariants (such as only makes sense to have fixed connectivity to an element from other elements)

struct connectivity_map
{
  typedef boost::array<boost::array<connectivity_kind,4>,4> map_type;

  static connectivity_map const& none()
  {
    static const map_type map =
    {{
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}}
    }};
    static samba::spatial_dimension sd = samba::spatial_dimension::create(3);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_map const& none_2d()
  {
    static const map_type map =
    {{
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{invalid(),invalid(),invalid(),invalid()}}
    }};
    static samba::spatial_dimension sd = samba::spatial_dimension::create(2);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_map const& all_dynamic()
  {
    static const map_type map =
    {{
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}},
      {{dynamic(),dynamic(),dynamic(),dynamic()}}
    }};
    static samba::spatial_dimension sd = samba::spatial_dimension::create(3);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_map const& all_dynamic_2d()
  {
    static const map_type map =
    {{
      {{dynamic(),dynamic(),invalid(),dynamic()}},
      {{dynamic(),dynamic(),invalid(),dynamic()}},
      {{invalid(),invalid(),invalid(),invalid()}},
      {{dynamic(),dynamic(),invalid(),dynamic()}}
    }};
    static samba::spatial_dimension sd = samba::spatial_dimension::create(2);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_map const& default_map()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), invalid(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), invalid()}},
      /*face*/  {{ fixed()  , dynamic(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , invalid(),  dynamic(), invalid()}}
    }};
    static samba::spatial_dimension sd = samba::spatial_dimension::create(3);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_map const& default_map_2d()
  {
    static const map_type map =
    {{
      //         TO
      //FROM       node       edge        face       element
      /*node*/  {{ invalid(), invalid(),  invalid(), dynamic()}},
      /*edge*/  {{ fixed()  , invalid(),  invalid(), invalid()}},
      /*face*/  {{ invalid(), invalid(),  invalid(), invalid()}},
      /*elem*/  {{ fixed()  , dynamic(),  invalid(), invalid()}}
    }};
    static samba::spatial_dimension sd = samba::spatial_dimension::create(2);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_map const& element_node()
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
    static samba::spatial_dimension sd = samba::spatial_dimension::create(3);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_map const& element_node_2d()
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
    static samba::spatial_dimension sd = samba::spatial_dimension::create(2);
    static connectivity_map r = {map,sd};
    return r;
  }

  static connectivity_kind invalid() { return connectivity_kind::invalid(); }
  static connectivity_kind fixed()   { return connectivity_kind::fixed();   }
  static connectivity_kind dynamic() { return connectivity_kind::dynamic(); }

  connectivity_kind & operator()(entity_rank from, entity_rank to)
  { return m_map[from()][to()]; }

  connectivity_kind operator()(entity_rank from, entity_rank to) const
  { return m_map[from()][to()]; }

  samba::spatial_dimension & spatial_dimension()
  { return m_spatial_dimension; }

  samba::spatial_dimension spatial_dimension() const
  { return m_spatial_dimension; }

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version)
  {
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
         ar & m_map[i][j];
    ar & m_spatial_dimension;
  }

  map_type m_map;
  samba::spatial_dimension m_spatial_dimension;
};

} //namespace samba

#endif //SAMBA_SAMBA_CONNECTIVITY_MAP_HPP

#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_PENTAGON_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_PENTAGON_HPP

namespace samba {

//entity_rank::pentagon_5_type: dimension = 2, vertices = 5, nodes = 5
inline std::ostream& operator<<(std::ostream& out, entity_topology::pentagon_5_type)
{ return out << "pentagon_5"; }

template<> struct entity_topology::is_valid_topology<entity_topology::pentagon_5_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::pentagon_5_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::pentagon_5_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_nodes<entity_topology::pentagon_5_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_edges<entity_topology::pentagon_5_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_sides<entity_topology::pentagon_5_type>
  : public boost::mpl::int_<5>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::pentagon_5_type, SideId>
{ typedef entity_topology::line_2_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::pentagon_5_type /*t*/, unsigned edge_id)
{
  switch (edge_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{2}, {3}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{3}, {4}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{4}, {0}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_PENTAGON_HPP

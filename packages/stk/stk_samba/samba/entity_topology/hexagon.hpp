#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_HEXAGON_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_HEXAGON_HPP

namespace samba {

//entity_rank::hexagon_6_type: dimension = 2, vertices = 6, nodes = 6
inline std::ostream& operator<<(std::ostream& out, entity_topology::hexagon_6_type)
{ return out << "hexagon_6"; }

template<> struct entity_topology::is_valid_topology<entity_topology::hexagon_6_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::hexagon_6_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::hexagon_6_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_nodes<entity_topology::hexagon_6_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_edges<entity_topology::hexagon_6_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_sides<entity_topology::hexagon_6_type>
  : public boost::mpl::int_<6>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::hexagon_6_type, SideId>
{ typedef entity_topology::line_2_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::hexagon_6_type /*t*/, unsigned edge_id)
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
      static const connectivity_ordinal nodes[] = {{4}, {5}};
      return nodes;
    }
  case 5:
    {
      static const connectivity_ordinal nodes[] = {{5}, {0}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_HEXAGON_HPP

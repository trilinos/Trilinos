#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_LINE_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_LINE_HPP

namespace samba {

//entity_rank::shell_line_2_type: dimension = 2, vertices = 2, nodes = 2
inline std::ostream& operator<<(std::ostream& out, entity_topology::shell_line_2_type)
{ return out << "shell_line_2"; }

template<> struct entity_topology::is_valid_topology<entity_topology::shell_line_2_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::shell_line_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::shell_line_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_nodes<entity_topology::shell_line_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_edges<entity_topology::shell_line_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_sides<entity_topology::shell_line_2_type>
  : public boost::mpl::int_<2>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::shell_line_2_type, SideId>
{ typedef entity_topology::line_2_type type; };

//entity_rank::shell_line_3_type: dimension = 2, vertices = 2, nodes = 3
inline std::ostream& operator<<(std::ostream& out, entity_topology::shell_line_3_type)
{ return out << "shell_line_3"; }

template<> struct entity_topology::is_valid_topology<entity_topology::shell_line_3_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::shell_line_3_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::shell_line_3_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_nodes<entity_topology::shell_line_3_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_edges<entity_topology::shell_line_3_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_sides<entity_topology::shell_line_3_type>
  : public boost::mpl::int_<2>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::shell_line_3_type, SideId>
{ typedef entity_topology::line_3_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::shell_line_3_type /*t*/, unsigned edge_id)
{
  switch (edge_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {2}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {0}, {2}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology::shell_line_2_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::shell_line_3_type(), edge_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_LINE_HPP

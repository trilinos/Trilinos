#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_TRIANGLE_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_TRIANGLE_HPP

#include <samba/entity_topology/triangle.hpp>

namespace samba {

//entity_rank::shell_triangle_3_type: dimension = 3, vertices = 3, nodes = 3
inline std::ostream& operator<<(std::ostream& out, entity_topology::shell_triangle_3_type)
{ return out << "shell_triangle_3"; }

template<> struct entity_topology::is_valid_topology<entity_topology::shell_triangle_3_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::shell_triangle_3_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::shell_triangle_3_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_nodes<entity_topology::shell_triangle_3_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_edges<entity_topology::shell_triangle_3_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_sides<entity_topology::shell_triangle_3_type>
  : public boost::mpl::int_<3>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::shell_triangle_3_type, SideId>
{ typedef entity_topology::line_2_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::shell_triangle_3_type, SideId>
{ typedef entity_topology::triangle_3_type type; };

//entity_rank::shell_triangle_6_type: dimension = 3, vertices = 3, nodes = 6
inline std::ostream& operator<<(std::ostream& out, entity_topology::shell_triangle_6_type)
{ return out << "shell_triangle_6"; }

template<> struct entity_topology::is_valid_topology<entity_topology::shell_triangle_6_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::shell_triangle_6_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::shell_triangle_6_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_nodes<entity_topology::shell_triangle_6_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_edges<entity_topology::shell_triangle_6_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_sides<entity_topology::shell_triangle_6_type>
  : public boost::mpl::int_<3>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::shell_triangle_6_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::shell_triangle_6_type, SideId>
{ typedef entity_topology::triangle_6_type type; };

inline
connectivity_ordinal const* const face_nodes(entity_topology::shell_triangle_6_type /*t*/, unsigned face_id)
{
  switch (face_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {2}, {3}, {4}, {5}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{0}, {2}, {1}, {5}, {4}, {3}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const face_nodes(entity_topology::shell_triangle_3_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::shell_triangle_6_type(), face_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::shell_triangle_6_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::triangle_6_type(), edge_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::shell_triangle_3_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::shell_triangle_6_type(), edge_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_TRIANGLE_HPP

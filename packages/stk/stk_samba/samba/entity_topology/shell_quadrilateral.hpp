
#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_QUAD_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_QUAD_HPP

#include <samba/entity_topology/quadrilateral.hpp>

namespace samba {

//entity_rank::shell_quadrilateral_4_type: dimension = 3, vertices = 4, nodes = 4
inline std::ostream& operator<<(std::ostream& out, entity_topology::shell_quadrilateral_4_type)
{ return out << "shell_quadrilateral_4"; }

template<> struct entity_topology::is_valid_topology<entity_topology::shell_quadrilateral_4_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::shell_quadrilateral_4_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::shell_quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::shell_quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_edges<entity_topology::shell_quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::shell_quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::shell_quadrilateral_4_type, SideId>
{ typedef entity_topology::line_2_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::shell_quadrilateral_4_type, SideId>
{ typedef entity_topology::quadrilateral_4_type type; };

//entity_rank::shell_quadrilateral_8_type: dimension = 3, vertices = 4, nodes = 8
inline std::ostream& operator<<(std::ostream& out, entity_topology::shell_quadrilateral_8_type)
{ return out << "shell_quadrilateral_8"; }

template<> struct entity_topology::is_valid_topology<entity_topology::shell_quadrilateral_8_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::shell_quadrilateral_8_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::shell_quadrilateral_8_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::shell_quadrilateral_8_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_edges<entity_topology::shell_quadrilateral_8_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::shell_quadrilateral_8_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::shell_quadrilateral_8_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::shell_quadrilateral_8_type, SideId>
{ typedef entity_topology::quadrilateral_8_type type; };

//entity_rank::shell_quadrilateral_9_type: dimension = 3, vertices = 4, nodes = 9
inline std::ostream& operator<<(std::ostream& out, entity_topology::shell_quadrilateral_9_type)
{ return out << "shell_quadrilateral_9"; }

template<> struct entity_topology::is_valid_topology<entity_topology::shell_quadrilateral_9_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::shell_quadrilateral_9_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::shell_quadrilateral_9_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::shell_quadrilateral_9_type>
  : public boost::mpl::int_<9>  {};

template<> struct entity_topology::num_edges<entity_topology::shell_quadrilateral_9_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::shell_quadrilateral_9_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::shell_quadrilateral_9_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::shell_quadrilateral_9_type, SideId>
{ typedef entity_topology::quadrilateral_9_type type; };

inline
connectivity_ordinal const* const face_nodes(entity_topology::shell_quadrilateral_9_type /*t*/, unsigned face_id)
{
  switch (face_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{0}, {3}, {2}, {1}, {7}, {6}, {5}, {4}, {8}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const face_nodes(entity_topology::shell_quadrilateral_8_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::shell_quadrilateral_9_type(), face_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::shell_quadrilateral_4_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::shell_quadrilateral_9_type(), face_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::shell_quadrilateral_9_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::quadrilateral_9_type(), edge_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::shell_quadrilateral_4_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::shell_quadrilateral_9_type(), edge_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::shell_quadrilateral_8_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::shell_quadrilateral_9_type(), edge_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_SHELL_QUAD_HPP

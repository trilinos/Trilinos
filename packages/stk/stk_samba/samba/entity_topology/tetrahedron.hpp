#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_TET_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_TET_HPP

namespace samba {

//entity_rank::tetrahedron_4_type: dimension = 3, vertices = 4, nodes = 4
inline std::ostream& operator<<(std::ostream& out, entity_topology::tetrahedron_4_type)
{ return out << "tetrahedron_4"; }

template<> struct entity_topology::is_valid_topology<entity_topology::tetrahedron_4_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::tetrahedron_4_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::tetrahedron_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::tetrahedron_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_edges<entity_topology::tetrahedron_4_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_faces<entity_topology::tetrahedron_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::tetrahedron_4_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::tetrahedron_4_type, SideId>
{ typedef entity_topology::line_2_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::tetrahedron_4_type, SideId>
{ typedef entity_topology::triangle_3_type type; };

//entity_rank::tetrahedron_8_type: dimension = 3, vertices = 4, nodes = 8
inline std::ostream& operator<<(std::ostream& out, entity_topology::tetrahedron_8_type)
{ return out << "tetrahedron_8"; }

template<> struct entity_topology::is_valid_topology<entity_topology::tetrahedron_8_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::tetrahedron_8_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::tetrahedron_8_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::tetrahedron_8_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_edges<entity_topology::tetrahedron_8_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_faces<entity_topology::tetrahedron_8_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::tetrahedron_8_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::tetrahedron_8_type, SideId>
{ typedef entity_topology::line_2_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::tetrahedron_8_type, SideId>
{ typedef entity_topology::triangle_4_type type; };

//entity_rank::tetrahedron_10_type: dimension = 3, vertices = 4, nodes = 10
inline std::ostream& operator<<(std::ostream& out, entity_topology::tetrahedron_10_type)
{ return out << "tetrahedron_10"; }

template<> struct entity_topology::is_valid_topology<entity_topology::tetrahedron_10_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::tetrahedron_10_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::tetrahedron_10_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::tetrahedron_10_type>
  : public boost::mpl::int_<10>  {};

template<> struct entity_topology::num_edges<entity_topology::tetrahedron_10_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_faces<entity_topology::tetrahedron_10_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::tetrahedron_10_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::tetrahedron_10_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::tetrahedron_10_type, SideId>
{ typedef entity_topology::triangle_6_type type; };

//entity_rank::tetrahedron_11_type: dimension = 3, vertices = 4, nodes = 11
inline std::ostream& operator<<(std::ostream& out, entity_topology::tetrahedron_11_type)
{ return out << "tetrahedron_11"; }

template<> struct entity_topology::is_valid_topology<entity_topology::tetrahedron_11_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::tetrahedron_11_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::tetrahedron_11_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::tetrahedron_11_type>
  : public boost::mpl::int_<11>  {};

template<> struct entity_topology::num_edges<entity_topology::tetrahedron_11_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_faces<entity_topology::tetrahedron_11_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::tetrahedron_11_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::tetrahedron_11_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::tetrahedron_11_type, SideId>
{ typedef entity_topology::triangle_6_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::tetrahedron_11_type /*t*/, unsigned edge_id)
{
  switch (edge_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {4}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {5}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{2}, {0}, {6}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{0}, {3}, {7}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{1}, {3}, {8}};
      return nodes;
    }
  case 5:
    {
      static const connectivity_ordinal nodes[] = {{2}, {3}, {9}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const face_nodes(entity_topology::tetrahedron_11_type /*t*/, unsigned face_id)
{
  switch (face_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {3}, {4}, {8}, {7}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {3}, {5}, {9}, {8}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{0}, {3}, {2}, {7}, {9}, {6}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{0}, {2}, {1}, {6}, {5}, {4}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology::tetrahedron_4_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::tetrahedron_11_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::tetrahedron_4_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::tetrahedron_11_type(), face_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::tetrahedron_8_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::tetrahedron_11_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::tetrahedron_8_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::tetrahedron_11_type(), face_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::tetrahedron_10_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::tetrahedron_11_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::tetrahedron_10_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::tetrahedron_11_type(), face_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_TET_HPP

#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_HEX_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_HEX_HPP

namespace samba {

//entity_rank::hexahedron_8_type: dimension = 3, vertices = 6, nodes = 8
inline std::ostream& operator<<(std::ostream& out, entity_topology::hexahedron_8_type)
{ return out << "hexahedron_8"; }

template<> struct entity_topology::is_valid_topology<entity_topology::hexahedron_8_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::hexahedron_8_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::hexahedron_8_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_nodes<entity_topology::hexahedron_8_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_edges<entity_topology::hexahedron_8_type>
  : public boost::mpl::int_<12>  {};

template<> struct entity_topology::num_faces<entity_topology::hexahedron_8_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_sides<entity_topology::hexahedron_8_type>
  : public boost::mpl::int_<6>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::hexahedron_8_type, SideId>
{ typedef entity_topology::line_2_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::hexahedron_8_type, SideId>
{ typedef entity_topology::quadrilateral_4_type type; };

//entity_rank::hexahedron_20_type: dimension = 3, vertices = 8, nodes = 20
inline std::ostream& operator<<(std::ostream& out, entity_topology::hexahedron_20_type)
{ return out << "hexahedron_20"; }

template<> struct entity_topology::is_valid_topology<entity_topology::hexahedron_20_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::hexahedron_20_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::hexahedron_20_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_nodes<entity_topology::hexahedron_20_type>
  : public boost::mpl::int_<20>  {};

template<> struct entity_topology::num_edges<entity_topology::hexahedron_20_type>
  : public boost::mpl::int_<12>  {};

template<> struct entity_topology::num_faces<entity_topology::hexahedron_20_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_sides<entity_topology::hexahedron_20_type>
  : public boost::mpl::int_<6>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::hexahedron_20_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::hexahedron_20_type, SideId>
{ typedef entity_topology::quadrilateral_8_type type; };

//entity_rank::hexahedron_27_type: dimension = 3, vertices = 8, nodes = 27
inline std::ostream& operator<<(std::ostream& out, entity_topology::hexahedron_27_type)
{ return out << "hexahedron_27"; }

template<> struct entity_topology::is_valid_topology<entity_topology::hexahedron_27_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::hexahedron_27_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::hexahedron_27_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_nodes<entity_topology::hexahedron_27_type>
  : public boost::mpl::int_<27>  {};

template<> struct entity_topology::num_edges<entity_topology::hexahedron_27_type>
  : public boost::mpl::int_<12>  {};

template<> struct entity_topology::num_faces<entity_topology::hexahedron_27_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_sides<entity_topology::hexahedron_27_type>
  : public boost::mpl::int_<6>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::hexahedron_27_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::hexahedron_27_type, SideId>
{ typedef entity_topology::quadrilateral_9_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::hexahedron_27_type /*t*/, unsigned edge_id)
{
  switch (edge_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {8}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {9}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{2}, {3}, {10}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{3}, {0}, {11}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{4}, {5}, {16}};
      return nodes;
    }
  case 5:
    {
      static const connectivity_ordinal nodes[] = {{5}, {6}, {17}};
      return nodes;
    }
  case 6:
    {
      static const connectivity_ordinal nodes[] = {{6}, {7}, {18}};
      return nodes;
    }
  case 7:
    {
      static const connectivity_ordinal nodes[] = {{7}, {4}, {19}};
      return nodes;
    }
  case 8:
    {
      static const connectivity_ordinal nodes[] = {{0}, {4}, {12}};
      return nodes;
    }
  case 9:
    {
      static const connectivity_ordinal nodes[] = {{1}, {5}, {13}};
      return nodes;
    }
  case 10:
    {
      static const connectivity_ordinal nodes[] = {{2}, {6}, {14}};
      return nodes;
    }
  case 11:
    {
      static const connectivity_ordinal nodes[] = {{3}, {7}, {15}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const face_nodes(entity_topology::hexahedron_27_type /*t*/, unsigned face_id)
{
  switch (face_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {5}, {4}, {8}, {13}, {16}, {12}, {25}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {6}, {5}, {9}, {14}, {17}, {13}, {24}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{2}, {3}, {7}, {6}, {10}, {15}, {18}, {14}, {26}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{0}, {4}, {7}, {3}, {12}, {19}, {15}, {11}, {23}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{0}, {3}, {2}, {1}, {11}, {10}, {9}, {8}, {21}};
      return nodes;
    }
  case 5:
    {
      static const connectivity_ordinal nodes[] = {{4}, {5}, {6}, {7}, {16}, {17}, {18}, {19}, {22}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology::hexahedron_8_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::hexahedron_27_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::hexahedron_8_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::hexahedron_27_type(), face_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::hexahedron_20_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::hexahedron_27_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::hexahedron_20_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::hexahedron_27_type(), face_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_HEX_HPP

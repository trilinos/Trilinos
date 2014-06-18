#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_WEDGE_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_WEDGE_HPP

namespace samba {

//entity_rank::wedge_6_type: dimension = 3, vertices = 6, nodes = 6
inline std::ostream& operator<<(std::ostream& out, entity_topology::wedge_6_type)
{ return out << "wedge_6"; }

template<> struct entity_topology::is_valid_topology<entity_topology::wedge_6_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::has_homogeneous_sides<entity_topology::wedge_6_type>
  : public boost::mpl::false_ {};

template<> struct entity_topology::dimension<entity_topology::wedge_6_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::wedge_6_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_nodes<entity_topology::wedge_6_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_edges<entity_topology::wedge_6_type>
  : public boost::mpl::int_<9>  {};

template<> struct entity_topology::num_faces<entity_topology::wedge_6_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_sides<entity_topology::wedge_6_type>
  : public boost::mpl::int_<5>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::wedge_6_type, SideId>
{ typedef entity_topology::line_2_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::wedge_6_type, SideId>
{ typedef entity_topology::quadrilateral_4_type type; };

template<> struct entity_topology::face_topology<entity_topology::wedge_6_type, 3>
{ typedef entity_topology::triangle_3_type type; };

template<> struct entity_topology::face_topology<entity_topology::wedge_6_type, 4>
{ typedef entity_topology::triangle_3_type type; };

//entity_rank::wedge_15_type: dimension = 3, vertices = 6, nodes = 15
inline std::ostream& operator<<(std::ostream& out, entity_topology::wedge_15_type)
{ return out << "wedge_15"; }

template<> struct entity_topology::is_valid_topology<entity_topology::wedge_15_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::has_homogeneous_sides<entity_topology::wedge_15_type>
  : public boost::mpl::false_ {};

template<> struct entity_topology::dimension<entity_topology::wedge_15_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::wedge_15_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_nodes<entity_topology::wedge_15_type>
  : public boost::mpl::int_<15>  {};

template<> struct entity_topology::num_edges<entity_topology::wedge_15_type>
  : public boost::mpl::int_<9>  {};

template<> struct entity_topology::num_faces<entity_topology::wedge_15_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_sides<entity_topology::wedge_15_type>
  : public boost::mpl::int_<5>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::wedge_15_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::wedge_15_type, SideId>
{ typedef entity_topology::quadrilateral_8_type type; };

template<> struct entity_topology::face_topology<entity_topology::wedge_15_type, 3>
{ typedef entity_topology::triangle_6_type type; };

template<> struct entity_topology::face_topology<entity_topology::wedge_15_type, 4>
{ typedef entity_topology::triangle_6_type type; };

//entity_rank::wedge_18_type: dimension = 3, vertices = 6, nodes = 18
inline std::ostream& operator<<(std::ostream& out, entity_topology::wedge_18_type)
{ return out << "wedge_18"; }

template<> struct entity_topology::is_valid_topology<entity_topology::wedge_18_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::has_homogeneous_sides<entity_topology::wedge_18_type>
  : public boost::mpl::false_ {};

template<> struct entity_topology::dimension<entity_topology::wedge_18_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::wedge_18_type>
  : public boost::mpl::int_<6>  {};

template<> struct entity_topology::num_nodes<entity_topology::wedge_18_type>
  : public boost::mpl::int_<18>  {};

template<> struct entity_topology::num_edges<entity_topology::wedge_18_type>
  : public boost::mpl::int_<9>  {};

template<> struct entity_topology::num_faces<entity_topology::wedge_18_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_sides<entity_topology::wedge_18_type>
  : public boost::mpl::int_<5>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::wedge_18_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::wedge_18_type, SideId>
{ typedef entity_topology::quadrilateral_9_type type; };

template<> struct entity_topology::face_topology<entity_topology::wedge_18_type, 3>
{ typedef entity_topology::triangle_6_type type; };

template<> struct entity_topology::face_topology<entity_topology::wedge_18_type, 4>
{ typedef entity_topology::triangle_6_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::wedge_18_type /*t*/, unsigned edge_id)
{
  switch (edge_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {6}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {7}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{2}, {0}, {8}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{3}, {4}, {12}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{4}, {5}, {13}};
      return nodes;
    }
  case 5:
    {
      static const connectivity_ordinal nodes[] = {{5}, {3}, {14}};
      return nodes;
    }
  case 6:
    {
      static const connectivity_ordinal nodes[] = {{0}, {3}, {9}};
      return nodes;
    }
  case 7:
    {
      static const connectivity_ordinal nodes[] = {{1}, {4}, {10}};
      return nodes;
    }
  case 8:
    {
      static const connectivity_ordinal nodes[] = {{2}, {5}, {11}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const face_nodes(entity_topology::wedge_18_type /*t*/, unsigned face_id)
{
  switch (face_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {4}, {3}, {6}, {10}, {12}, {9}, {15}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {5}, {4}, {7}, {11}, {13}, {10}, {16}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{0}, {3}, {5}, {2}, {9}, {14}, {11}, {8}, {17}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{0}, {2}, {1}, {8}, {7}, {6}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{3}, {4}, {5}, {12}, {13}, {14}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology::wedge_6_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::wedge_18_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::wedge_6_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::wedge_18_type(), face_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::wedge_15_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::wedge_18_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::wedge_15_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::wedge_18_type(), face_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_WEDGE_HPP

#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_PYRAMID_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_PYRAMID_HPP

namespace samba {

//entity_rank::pyramid_5_type: dimension = 3, vertices = 5, nodes = 5
inline std::ostream& operator<<(std::ostream& out, entity_topology::pyramid_5_type)
{ return out << "pyramid_5"; }

template<> struct entity_topology::is_valid_topology<entity_topology::pyramid_5_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::has_homogeneous_sides<entity_topology::pyramid_5_type>
  : public boost::mpl::false_ {};

template<> struct entity_topology::dimension<entity_topology::pyramid_5_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::pyramid_5_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_nodes<entity_topology::pyramid_5_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_edges<entity_topology::pyramid_5_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_faces<entity_topology::pyramid_5_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_sides<entity_topology::pyramid_5_type>
  : public boost::mpl::int_<5>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::pyramid_5_type, SideId>
{ typedef entity_topology::line_2_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::pyramid_5_type, SideId>
{ typedef entity_topology::triangle_3_type type; };

template<> struct entity_topology::face_topology<entity_topology::pyramid_5_type, 4>
{ typedef entity_topology::quadrilateral_4_type type; };

//entity_rank::pyramid_13_type: dimension = 3, vertices = 5, nodes = 13
inline std::ostream& operator<<(std::ostream& out, entity_topology::pyramid_13_type)
{ return out << "pyramid_13"; }

template<> struct entity_topology::is_valid_topology<entity_topology::pyramid_13_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::has_homogeneous_sides<entity_topology::pyramid_13_type>
  : public boost::mpl::false_ {};

template<> struct entity_topology::dimension<entity_topology::pyramid_13_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::pyramid_13_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_nodes<entity_topology::pyramid_13_type>
  : public boost::mpl::int_<13>  {};

template<> struct entity_topology::num_edges<entity_topology::pyramid_13_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_faces<entity_topology::pyramid_13_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_sides<entity_topology::pyramid_13_type>
  : public boost::mpl::int_<5>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::pyramid_13_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::pyramid_13_type, SideId>
{ typedef entity_topology::triangle_6_type type; };

template<> struct entity_topology::face_topology<entity_topology::pyramid_13_type, 4>
{ typedef entity_topology::quadrilateral_8_type type; };

//entity_rank::pyramid_14_type: dimension = 3, vertices = 5, nodes = 14
inline std::ostream& operator<<(std::ostream& out, entity_topology::pyramid_14_type)
{ return out << "pyramid_14"; }

template<> struct entity_topology::is_valid_topology<entity_topology::pyramid_14_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::has_homogeneous_sides<entity_topology::pyramid_14_type>
  : public boost::mpl::false_ {};

template<> struct entity_topology::dimension<entity_topology::pyramid_14_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_vertices<entity_topology::pyramid_14_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_nodes<entity_topology::pyramid_14_type>
  : public boost::mpl::int_<14>  {};

template<> struct entity_topology::num_edges<entity_topology::pyramid_14_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_faces<entity_topology::pyramid_14_type>
  : public boost::mpl::int_<5>  {};

template<> struct entity_topology::num_sides<entity_topology::pyramid_14_type>
  : public boost::mpl::int_<5>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::pyramid_14_type, SideId>
{ typedef entity_topology::line_3_type type; };

template<unsigned SideId> struct entity_topology::face_topology<entity_topology::pyramid_14_type, SideId>
{ typedef entity_topology::triangle_6_type type; };

template<> struct entity_topology::face_topology<entity_topology::pyramid_14_type, 4>
{ typedef entity_topology::quadrilateral_9_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::pyramid_14_type /*t*/, unsigned edge_id)
{
  switch (edge_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {5}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {6}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{2}, {3}, {7}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{3}, {0}, {8}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{0}, {4}, {9}};
      return nodes;
    }
  case 5:
    {
      static const connectivity_ordinal nodes[] = {{1}, {4}, {10}};
      return nodes;
    }
  case 6:
    {
      static const connectivity_ordinal nodes[] = {{2}, {4}, {11}};
      return nodes;
    }
  case 7:
    {
      static const connectivity_ordinal nodes[] = {{3}, {4}, {12}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const face_nodes(entity_topology::pyramid_14_type /*t*/, unsigned face_id)
{
  switch (face_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {4}, {5}, {10}, {9}};
      return nodes;
    }
  case 1:
    {
      static const connectivity_ordinal nodes[] = {{1}, {2}, {4}, {6}, {11}, {10}};
      return nodes;
    }
  case 2:
    {
      static const connectivity_ordinal nodes[] = {{2}, {3}, {4}, {7}, {12}, {11}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{3}, {0}, {4}, {8}, {9}, {12}};
      return nodes;
    }
  case 4:
    {
      static const connectivity_ordinal nodes[] = {{0}, {3}, {2}, {1}, {8}, {7}, {6}, {5}, {13}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology::pyramid_5_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::pyramid_14_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::pyramid_5_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::pyramid_14_type(), face_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::pyramid_13_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::pyramid_14_type(), edge_id); }

inline
connectivity_ordinal const* const face_nodes(entity_topology::pyramid_13_type /*t*/, unsigned face_id)
{ return face_nodes(entity_topology::pyramid_14_type(), face_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_PYRAMID_HPP

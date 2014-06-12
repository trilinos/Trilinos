#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_QUAD_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_QUAD_HPP

namespace samba {

//entity_rank::quadrilateral_4_type: dimension = 2, vertices = 4, nodes = 4
inline std::ostream& operator<<(std::ostream& out, entity_topology::quadrilateral_4_type)
{ return out << "quadrilateral_4"; }

template<> struct entity_topology::is_valid_topology<entity_topology::quadrilateral_4_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::quadrilateral_4_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_edges<entity_topology::quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::quadrilateral_4_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::quadrilateral_4_type, SideId>
{ typedef entity_topology::line_2_type type; };

//entity_rank::quadrilateral_8_type: dimension = 2, vertices = 4, nodes = 8
inline std::ostream& operator<<(std::ostream& out, entity_topology::quadrilateral_8_type)
{ return out << "quadrilateral_8"; }

template<> struct entity_topology::is_valid_topology<entity_topology::quadrilateral_8_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::quadrilateral_8_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::quadrilateral_8_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::quadrilateral_8_type>
  : public boost::mpl::int_<8>  {};

template<> struct entity_topology::num_edges<entity_topology::quadrilateral_8_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::quadrilateral_8_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::quadrilateral_8_type, SideId>
{ typedef entity_topology::line_3_type type; };

//entity_rank::quadrilateral_9_type: dimension = 2, vertices = 4, nodes = 9
inline std::ostream& operator<<(std::ostream& out, entity_topology::quadrilateral_9_type)
{ return out << "quadrilateral_9"; }

template<> struct entity_topology::is_valid_topology<entity_topology::quadrilateral_9_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::quadrilateral_9_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::quadrilateral_9_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_nodes<entity_topology::quadrilateral_9_type>
  : public boost::mpl::int_<9>  {};

template<> struct entity_topology::num_edges<entity_topology::quadrilateral_9_type>
  : public boost::mpl::int_<4>  {};

template<> struct entity_topology::num_sides<entity_topology::quadrilateral_9_type>
  : public boost::mpl::int_<4>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::quadrilateral_9_type, SideId>
{ typedef entity_topology::line_3_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::quadrilateral_9_type /*t*/, unsigned edge_id)
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
      static const connectivity_ordinal nodes[] = {{2}, {3}, {6}};
      return nodes;
    }
  case 3:
    {
      static const connectivity_ordinal nodes[] = {{3}, {0}, {7}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology::quadrilateral_4_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::quadrilateral_9_type(), edge_id); }

inline
connectivity_ordinal const* const edge_nodes(entity_topology::quadrilateral_8_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::quadrilateral_9_type(), edge_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_QUAD_HPP

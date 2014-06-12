#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_BEAM_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_BEAM_HPP

namespace samba {

//entity_rank::beam_2_type: dimension = 2, vertices = 2, nodes = 2
inline std::ostream& operator<<(std::ostream& out, entity_topology::beam_2_type)
{ return out << "beam_2"; }

template<> struct entity_topology::is_valid_topology<entity_topology::beam_2_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::beam_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::beam_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_nodes<entity_topology::beam_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_edges<entity_topology::beam_2_type>
  : public boost::mpl::int_<1>  {};

template<> struct entity_topology::num_sides<entity_topology::beam_2_type>
  : public boost::mpl::int_<1>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::beam_2_type, SideId>
{ typedef entity_topology::line_2_type type; };


//entity_rank::beam_3_type: dimension = 2, vertices = 2, nodes = 3
inline std::ostream& operator<<(std::ostream& out, entity_topology::beam_3_type)
{ return out << "beam_3"; }

template<> struct entity_topology::is_valid_topology<entity_topology::beam_3_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::beam_3_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_vertices<entity_topology::beam_3_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_nodes<entity_topology::beam_3_type>
  : public boost::mpl::int_<3>  {};

template<> struct entity_topology::num_edges<entity_topology::beam_3_type>
  : public boost::mpl::int_<1>  {};

template<> struct entity_topology::num_sides<entity_topology::beam_3_type>
  : public boost::mpl::int_<1>  {};

template<unsigned SideId> struct entity_topology::edge_topology<entity_topology::beam_3_type, SideId>
{ typedef entity_topology::line_3_type type; };

inline
connectivity_ordinal const* const edge_nodes(entity_topology::beam_3_type /*t*/, unsigned edge_id)
{
  switch (edge_id) {
  case 0:
    {
      static const connectivity_ordinal nodes[] = {{0}, {1}, {2}};
      return nodes;
    }
  default:
    {}
  }

  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology::beam_2_type /*t*/, unsigned edge_id)
{ return edge_nodes(entity_topology::beam_3_type(), edge_id); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_BEAM_HPP

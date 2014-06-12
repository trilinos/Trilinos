#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_LINE_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_LINE_HPP

namespace samba {

//entity_rank::line_2_type: dimension = 1, vertices = 2, nodes = 2
inline std::ostream& operator<<(std::ostream& out, entity_topology::line_2_type)
{ return out << "line_2"; }

template<> struct entity_topology::is_valid_topology<entity_topology::line_2_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::line_2_type>
  : public boost::mpl::int_<1>  {};

template<> struct entity_topology::num_vertices<entity_topology::line_2_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_nodes<entity_topology::line_2_type>
  : public boost::mpl::int_<2>  {};


//entity_rank::line_3_type: dimension = 1, vertices = 2, nodes = 3
inline std::ostream& operator<<(std::ostream& out, entity_topology::line_3_type)
{ return out << "line_3"; }

template<> struct entity_topology::is_valid_topology<entity_topology::line_3_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::line_3_type>
  : public boost::mpl::int_<1>  {};

template<> struct entity_topology::num_vertices<entity_topology::line_3_type>
  : public boost::mpl::int_<2>  {};

template<> struct entity_topology::num_nodes<entity_topology::line_3_type>
  : public boost::mpl::int_<3>  {};

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_LINE_HPP

#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_NODE_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_NODE_HPP

namespace samba {

//entity_rank::node_type: dimension = 0, vertices = 0, nodes = 0
inline std::ostream& operator<<(std::ostream& out, entity_topology::node_type)
{ return out << "node"; }

template<> struct entity_topology::is_valid_topology<entity_topology::node_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::node_type>
  : public boost::mpl::int_<0>  {};

template<> struct entity_topology::num_vertices<entity_topology::node_type>
  : public boost::mpl::int_<0>  {};

template<> struct entity_topology::num_nodes<entity_topology::node_type>
  : public boost::mpl::int_<0>  {};

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_NODE_HPP

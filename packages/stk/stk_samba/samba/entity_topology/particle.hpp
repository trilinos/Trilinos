#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_PARTICLE_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_PARTICLE_HPP

namespace samba {

//entity_rank::particle_type: dimension = 1, vertices = 1, nodes = 1
inline std::ostream& operator<<(std::ostream& out, entity_topology::particle_type)
{ return out << "particle"; }

template<> struct entity_topology::is_valid_topology<entity_topology::particle_type>
  : public boost::mpl::true_ {};

template<> struct entity_topology::dimension<entity_topology::particle_type>
  : public boost::mpl::int_<1>  {};

template<> struct entity_topology::num_vertices<entity_topology::particle_type>
  : public boost::mpl::int_<1>  {};

template<> struct entity_topology::num_nodes<entity_topology::particle_type>
  : public boost::mpl::int_<1>  {};

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_PARTICLE_HPP

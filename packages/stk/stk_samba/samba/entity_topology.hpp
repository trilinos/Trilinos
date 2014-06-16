#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_HPP
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_HPP

#include <samba/utility.hpp>
#include <samba/entity_rank.hpp>
#include <samba/spatial_dimension.hpp>
#include <samba/connectivity_ordinal.hpp>

#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>

namespace samba {

struct entity_topology
{
  struct tag
  {
    typedef int value_type;
    static const int value = 1;
    friend inline std::ostream& operator<<(std::ostream& out,tag)
    { return out << "entity_topology"; }
  };

  typedef uint8_t value_type;

  static const int num_bits = 8;

  template <entity_topology::value_type Topology>
  struct topology_type
  {
    typedef entity_topology::value_type value_type;
    static const value_type value = Topology;
  };

  typedef topology_type<0 > node_type;
  typedef topology_type<1 > particle_type;
  typedef topology_type<2 > line_2_type;
  typedef topology_type<3 > line_3_type;
  typedef topology_type<4 > beam_2_type;
  typedef topology_type<5 > beam_3_type;
  typedef topology_type<6 > shell_line_2_type;
  typedef topology_type<7 > shell_line_3_type;
  typedef topology_type<8 > triangle_3_type;
  typedef topology_type<9 > triangle_4_type;
  typedef topology_type<10> triangle_6_type;
  typedef topology_type<11> quadrilateral_4_type;
  typedef topology_type<12> quadrilateral_8_type;
  typedef topology_type<13> quadrilateral_9_type;
  typedef topology_type<14> pentagon_5_type;
  typedef topology_type<15> hexagon_6_type;
  typedef topology_type<16> shell_triangle_3_type;
  typedef topology_type<17> shell_triangle_6_type;
  typedef topology_type<18> shell_quadrilateral_4_type;
  typedef topology_type<19> shell_quadrilateral_8_type;
  typedef topology_type<20> shell_quadrilateral_9_type;
  typedef topology_type<21> tetrahedron_4_type;
  typedef topology_type<22> tetrahedron_8_type;
  typedef topology_type<23> tetrahedron_10_type;
  typedef topology_type<24> tetrahedron_11_type;
  typedef topology_type<25> pyramid_5_type;
  typedef topology_type<26> pyramid_13_type;
  typedef topology_type<27> pyramid_14_type;
  typedef topology_type<28> wedge_6_type;
  typedef topology_type<29> wedge_15_type;
  typedef topology_type<30> wedge_18_type;
  typedef topology_type<31> hexahedron_8_type;
  typedef topology_type<32> hexahedron_20_type;
  typedef topology_type<33> hexahedron_27_type;
  typedef topology_type<34> invalid_type;

  static const value_type max_possible_topology = invalid_type::value;

  static const entity_topology node()
  { static entity_topology d = {node_type::value}; return d; }

  static const entity_topology particle()
  { static entity_topology d = {particle_type::value}; return d; }

  static const entity_topology line_2()
  { static entity_topology d = {line_2_type::value}; return d; }

  static const entity_topology line_3()
  { static entity_topology d = {line_3_type::value}; return d; }

  static const entity_topology beam_2()
  { static entity_topology d = {beam_2_type::value}; return d; }

  static const entity_topology beam_3()
  { static entity_topology d = {beam_3_type::value}; return d; }

  static const entity_topology shell_line_2()
  { static entity_topology d = {shell_line_2_type::value}; return d; }

  static const entity_topology shell_line_3()
  { static entity_topology d = {shell_line_3_type::value}; return d; }

  static const entity_topology triangle_3()
  { static entity_topology d = {triangle_3_type::value}; return d; }

  static const entity_topology triangle_4()
  { static entity_topology d = {triangle_4_type::value}; return d; }

  static const entity_topology triangle_6()
  { static entity_topology d = {triangle_6_type::value}; return d; }

  static const entity_topology quadrilateral_4()
  { static entity_topology d = {quadrilateral_4_type::value}; return d; }

  static const entity_topology quadrilateral_8()
  { static entity_topology d = {quadrilateral_8_type::value}; return d; }

  static const entity_topology quadrilateral_9()
  { static entity_topology d = {quadrilateral_9_type::value}; return d; }

  static const entity_topology pentagon_5()
  { static entity_topology d = {pentagon_5_type::value}; return d; }

  static const entity_topology hexagon_6()
  { static entity_topology d = {hexagon_6_type::value}; return d; }

  static const entity_topology shell_triangle_3()
  { static entity_topology d = {shell_triangle_3_type::value}; return d; }

  static const entity_topology shell_triangle_6()
  { static entity_topology d = {shell_triangle_6_type::value}; return d; }

  static const entity_topology shell_quadrilateral_4()
  { static entity_topology d = {shell_quadrilateral_4_type::value}; return d; }

  static const entity_topology shell_quadrilateral_8()
  { static entity_topology d = {shell_quadrilateral_8_type::value}; return d; }

  static const entity_topology shell_quadrilateral_9()
  { static entity_topology d = {shell_quadrilateral_9_type::value}; return d; }

  static const entity_topology tetrahedron_4()
  { static entity_topology d = {tetrahedron_4_type::value}; return d; }

  static const entity_topology tetrahedron_8()
  { static entity_topology d = {tetrahedron_8_type::value}; return d; }

  static const entity_topology tetrahedron_10()
  { static entity_topology d = {tetrahedron_10_type::value}; return d; }

  static const entity_topology tetrahedron_11()
  { static entity_topology d = {tetrahedron_11_type::value}; return d; }

  static const entity_topology pyramid_5()
  { static entity_topology d = {pyramid_5_type::value}; return d; }

  static const entity_topology pyramid_13()
  { static entity_topology d = {pyramid_13_type::value}; return d; }

  static const entity_topology pyramid_14()
  { static entity_topology d = {pyramid_14_type::value}; return d; }

  static const entity_topology wedge_6()
  { static entity_topology d = {wedge_6_type::value}; return d; }

  static const entity_topology wedge_15()
  { static entity_topology d = {wedge_15_type::value}; return d; }

  static const entity_topology wedge_18()
  { static entity_topology d = {wedge_18_type::value}; return d; }

  static const entity_topology hexahedron_8()
  { static entity_topology d = {hexahedron_8_type::value}; return d; }

  static const entity_topology hexahedron_20()
  { static entity_topology d = {hexahedron_20_type::value}; return d; }

  static const entity_topology hexahedron_27()
  { static entity_topology d = {hexahedron_27_type::value}; return d; }

  static const entity_topology invalid()
  { static entity_topology d = {invalid_type::value}; return d; }

  //****************************************************
  //Abbreviations
  //****************************************************

  static const entity_topology tri_3() { return triangle_3(); }
  static const entity_topology tri_4() { return triangle_4(); }
  static const entity_topology tri_6() { return triangle_6(); }

  static const entity_topology shell_tri_3() { return shell_triangle_3(); }
  static const entity_topology shell_tri_6() { return shell_triangle_6(); }

  static const entity_topology quad_4() { return quadrilateral_4(); }
  static const entity_topology quad_8() { return quadrilateral_8(); }
  static const entity_topology quad_9() { return quadrilateral_9(); }

  static const entity_topology shell_quad_4() { return shell_quadrilateral_4(); }
  static const entity_topology shell_quad_8() { return shell_quadrilateral_8(); }
  static const entity_topology shell_quad_9() { return shell_quadrilateral_9(); }

  static const entity_topology tet_4() { return tetrahedron_4(); }
  static const entity_topology tet_8() { return tetrahedron_8(); }
  static const entity_topology tet_10() { return tetrahedron_10(); }
  static const entity_topology tet_11() { return tetrahedron_11(); }

  static const entity_topology hex_8() { return hexahedron_8(); }
  static const entity_topology hex_20() { return hexahedron_20(); }
  static const entity_topology hex_27() { return hexahedron_27(); }

  //****************************************************

  static const entity_topology create(value_type v)
  { entity_topology d = {v}; return d; }

  value_type operator()() const { return m_value; }

  SAMBA_PRIMITIVE_COMPARISONS(entity_topology,value_type)
  SAMBA_ARITHMETIC_OPERATORS(entity_topology,value_type)

  entity_topology & operator=(entity_topology::value_type v)
  { m_value = v; return *this; }

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version)
  { ar & m_value; }

  //***************************************************************************
  //meta_functions
  //***************************************************************************

  template <typename T> struct is_valid_topology     : public boost::mpl::false_ {};
  template <typename T> struct has_homogeneous_sides : public boost::mpl::true_ {};
  template <typename T> struct dimension             : public boost::mpl::int_<0> {};
  template <typename T> struct side_dimension        : public boost::mpl::int_<dimension<T>::value-1> {};
  template <typename T> struct num_nodes             : public boost::mpl::int_<0> {};
  template <typename T> struct num_vertices          : public boost::mpl::int_<0> {};
  template <typename T> struct num_edges             : public boost::mpl::int_<0> {};
  template <typename T> struct num_faces             : public boost::mpl::int_<0> {};
  template <typename T> struct num_sides             : public boost::mpl::int_<0> {};

  template <typename T, unsigned SideId=0> struct face_topology
  {
    typedef invalid_type type;
  };

  template <typename T, unsigned SideId=0> struct edge_topology
  {
    typedef invalid_type type;
  };

  value_type m_value;
};

bool is_valid_topology(entity_topology t);
bool has_homogeneous_sides(entity_topology t);
unsigned dimension(entity_topology t);
unsigned side_dimension(entity_topology t);
unsigned num_nodes(entity_topology t);
unsigned num_vertices(entity_topology t);
unsigned num_edges(entity_topology t);
unsigned num_faces(entity_topology t);
unsigned num_sides(entity_topology t);
void print_topology_detail(std::ostream& out, entity_topology t);
entity_topology side_topology(entity_topology t, unsigned side_id);
entity_topology face_topology(entity_topology t, unsigned side_id);
entity_topology edge_topology(entity_topology t, unsigned side_id);
entity_topology connectivity_topology(entity_topology from, entity_rank to, connectivity_ordinal to_ordinal);

entity_rank topology_rank(entity_topology t,spatial_dimension sp);
entity_rank side_rank(entity_topology t,spatial_dimension sp);

connectivity_ordinal const* const side_nodes(entity_topology t, unsigned side_id); // side = topo rank - 1
connectivity_ordinal const* const face_nodes(entity_topology t, unsigned face_id);
connectivity_ordinal const* const edge_nodes(entity_topology t, unsigned edge_id);

} //namespace samba

SAMBA_IS_PRIMITIVE(samba::entity_topology)

#include <samba/entity_topology/entity_topology.tcc>

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_HPP

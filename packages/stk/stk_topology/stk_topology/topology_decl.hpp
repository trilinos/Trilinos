#ifndef TOPOLOGY_DECL_HPP
#define TOPOLOGY_DECL_HPP

// IWYU pragma: private, include "stk_topology/topology.hpp"
#include "stk_util/stk_config.h"
#include <string>
#include <cstdint>
#include <limits>

namespace stk {

struct EquivalentPermutation
{
  STK_FUNCTION
  EquivalentPermutation()
    : is_equivalent(false),
      permutation_number(0) {}

  STK_FUNCTION
  EquivalentPermutation(bool _is_equivalent, unsigned _permutation_number)
    : is_equivalent(_is_equivalent),
      permutation_number(_permutation_number) {}

  bool     is_equivalent;
  unsigned permutation_number;
};

struct topology
{
  enum rank_t : int8_t
  {
    BEGIN_RANK,
    NODE_RANK = BEGIN_RANK,
    EDGE_RANK,
    FACE_RANK,
    ELEM_RANK, ELEMENT_RANK = ELEM_RANK,
    CONSTRAINT_RANK,
    END_RANK,
    NUM_RANKS = END_RANK,
    INVALID_RANK = std::numeric_limits<int8_t>::max()
  };

  //To add new topologies consult the toolkit team
  enum topology_t
  {
    INVALID_TOPOLOGY,
    BEGIN_TOPOLOGY,

    // NODE_RANK
    NODE = BEGIN_TOPOLOGY,

    // EDGE_RANK
    LINE_2, BEGIN_EDGE_RANK = LINE_2,
    LINE_3,

    // FACE_RANK
    TRI_3, TRIANGLE_3 = TRI_3, END_EDGE_RANK = TRI_3, BEGIN_FACE_RANK = TRI_3,
    TRI_4, TRIANGLE_4 = TRI_4,
    TRI_6, TRIANGLE_6 = TRI_6,
    QUAD_4, QUADRILATERAL_4 = QUAD_4,
    QUAD_6, QUADRILATERAL_6 = QUAD_6,
    QUAD_8, QUADRILATERAL_8 = QUAD_8,
    QUAD_9, QUADRILATERAL_9 = QUAD_9,
    SHELL_SIDE_BEAM_2,
    SHELL_SIDE_BEAM_3,

    // ELEMENT_RANK
    PARTICLE, END_FACE_RANK = PARTICLE, BEGIN_ELEMENT_RANK = PARTICLE,
    LINE_2_1D,
    LINE_3_1D,
    BEAM_2,
    BEAM_3,
    SHELL_LINE_2,
    SHELL_LINE_3,
    SPRING_2,
    SPRING_3,
    TRI_3_2D, TRIANGLE_3_2D = TRI_3_2D,
    TRI_4_2D, TRIANGLE_4_2D = TRI_4_2D,
    TRI_6_2D, TRIANGLE_6_2D = TRI_6_2D,
    QUAD_4_2D, QUADRILATERAL_4_2D = QUAD_4_2D,
    QUAD_8_2D, QUADRILATERAL_8_2D = QUAD_8_2D,
    QUAD_9_2D, QUADRILATERAL_9_2D = QUAD_9_2D,
    SHELL_TRI_3, SHELL_TRIANGLE_3 = SHELL_TRI_3,
    SHELL_TRI_4, SHELL_TRIANGLE_4 = SHELL_TRI_4,
    SHELL_TRI_6, SHELL_TRIANGLE_6 = SHELL_TRI_6,
    SHELL_TRI_3_ALL_FACE_SIDES, SHELL_TRIANGLE_3_ALL_FACE_SIDES = SHELL_TRI_3_ALL_FACE_SIDES,
    SHELL_TRI_4_ALL_FACE_SIDES, SHELL_TRIANGLE_4_ALL_FACE_SIDES = SHELL_TRI_4_ALL_FACE_SIDES,
    SHELL_TRI_6_ALL_FACE_SIDES, SHELL_TRIANGLE_6_ALL_FACE_SIDES = SHELL_TRI_6_ALL_FACE_SIDES,
    SHELL_QUAD_4, SHELL_QUADRILATERAL_4 = SHELL_QUAD_4,
    SHELL_QUAD_8, SHELL_QUADRILATERAL_8 = SHELL_QUAD_8,
    SHELL_QUAD_9, SHELL_QUADRILATERAL_9 = SHELL_QUAD_9,
    SHELL_QUAD_4_ALL_FACE_SIDES, SHELL_QUADRILATERAL_4_ALL_FACE_SIDES = SHELL_QUAD_4_ALL_FACE_SIDES,
    SHELL_QUAD_8_ALL_FACE_SIDES, SHELL_QUADRILATERAL_8_ALL_FACE_SIDES = SHELL_QUAD_8_ALL_FACE_SIDES,
    SHELL_QUAD_9_ALL_FACE_SIDES, SHELL_QUADRILATERAL_9_ALL_FACE_SIDES = SHELL_QUAD_9_ALL_FACE_SIDES,
    TET_4,  TETRAHEDRON_4  = TET_4,
    TET_8,  TETRAHEDRON_8  = TET_8,
    TET_10, TETRAHEDRON_10 = TET_10,
    TET_11, TETRAHEDRON_11 = TET_11,
    PYRAMID_5,
    PYRAMID_13,
    PYRAMID_14,
    WEDGE_6,
    WEDGE_12,
    WEDGE_15,
    WEDGE_18,
    HEX_8,  HEXAHEDRON_8  = HEX_8,
    HEX_20, HEXAHEDRON_20 = HEX_20,
    HEX_27, HEXAHEDRON_27 = HEX_27,

    END_TOPOLOGY, END_ELEMENT_RANK = END_TOPOLOGY,
    NUM_TOPOLOGIES = END_TOPOLOGY - BEGIN_TOPOLOGY,

    SUPEREDGE_START = END_TOPOLOGY+1,
    SUPEREDGE_END = SUPEREDGE_START + 1000,
    SUPERFACE_START = SUPEREDGE_END+1,
    SUPERFACE_END = SUPERFACE_START + 1000,
    SUPERELEMENT_START = SUPERFACE_END+1,
    FORCE_TOPOLOGY_TO_UNSIGNED = ~0U // max unsigned int
  };

  //***************************************************************************
  //member functions
  //***************************************************************************

  /// is this topology valid
  STK_INLINE_FUNCTION
  bool is_valid() const { return m_value != INVALID_TOPOLOGY; }

  /// get the name of this topology
  std::string name() const;

  /// get the name of this topology
  STK_FUNCTION
  const char * char_name() const;

  /// does this topology have homogeneous faces
  STK_INLINE_FUNCTION
  bool has_homogeneous_faces() const;

  /// is this topology a shell topology (i.e. an element with only two sides)
  STK_INLINE_FUNCTION
  bool is_shell() const;

  STK_INLINE_FUNCTION
  bool is_shell_side_ordinal(unsigned ord) const;

  STK_INLINE_FUNCTION
  bool is_shell_with_face_sides() const;

  /// what is the rank of this topology
  STK_INLINE_FUNCTION
  rank_t rank() const;

  /// what is the side rank of this topology
  STK_INLINE_FUNCTION
  rank_t side_rank() const;

  /// what is the topological dimension of this topology
  STK_INLINE_FUNCTION
  unsigned dimension() const;

  /// how many nodes define this topology
  STK_INLINE_FUNCTION
  unsigned num_nodes() const;

  /// how many nodes are vertices
  STK_INLINE_FUNCTION
  unsigned num_vertices() const;

  /// how many edges does this topology have
  STK_INLINE_FUNCTION
  unsigned num_edges() const;

  /// how many faces does this topology have
  STK_INLINE_FUNCTION
  unsigned num_faces() const;

  /// how many different node permutations does this topology have
  STK_INLINE_FUNCTION
  unsigned num_permutations() const;

  /// how many different positive node permutations does this topology have
  STK_INLINE_FUNCTION
  unsigned num_positive_permutations() const;

  STK_INLINE_FUNCTION
  bool is_positive_polarity(unsigned permutation_ordinal) const;

  /// is this topology defined on the given spatial dimension
  STK_INLINE_FUNCTION
  bool defined_on_spatial_dimension(unsigned spatial_dimension) const;

  /// what is the base topology (i.e. topology where num_nodes == num_vertices)
  STK_INLINE_FUNCTION
  topology base() const;

  /// what is the topology of the given edge
  STK_INLINE_FUNCTION
  topology edge_topology() const;

  /// what is the topology of the given edge
  STK_INLINE_FUNCTION
  topology edge_topology(unsigned edge_orginal) const;

  /// what is the topology of the given face
  STK_INLINE_FUNCTION
  topology face_topology(unsigned face_ordinal) const;

  STK_INLINE_FUNCTION
  topology shell_side_topology(unsigned shell_side_ordinal = 0) const;

  /// fill the output ordinals with the ordinals that make up the given edge
  template <typename OrdinalOutputIterator>
  STK_INLINE_FUNCTION
  void edge_node_ordinals(unsigned edge_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output ordinals with the ordinals that make up the given face
  template <typename OrdinalOutputIterator>
  STK_INLINE_FUNCTION
  void face_node_ordinals(unsigned face_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output ordinals with the ordinals that make up the given permutation
  template <typename OrdinalOutputIterator>
  STK_INLINE_FUNCTION
  void permutation_node_ordinals(unsigned permutation_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output nodes with the nodes that make up the given edge
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  STK_INLINE_FUNCTION
  void edge_nodes(const NodeArray & nodes, unsigned edge_ordinal, NodeOutputIterator output_nodes) const;

  /// fill the output nodes with the nodes that make up the given face
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  STK_INLINE_FUNCTION
  void face_nodes(const NodeArray & nodes, unsigned face_ordinal, NodeOutputIterator output_nodes) const;

  /// fill the output nodes with the nodes that make up the given permutation
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  STK_INLINE_FUNCTION
  void permutation_nodes(const NodeArray & nodes, unsigned permutation_ordinal, NodeOutputIterator output_nodes) const;

  /// do the two arrays define equivalent entities (same nodes, but maybe a different permutation)
  /// return a struct containing a bool and permutation number from a to b
  template <typename NodeArrayA, typename NodeArrayB>
  STK_INLINE_FUNCTION
  EquivalentPermutation is_equivalent(const NodeArrayA & a, const NodeArrayB & possible_permutation_of_a) const;

  /// return the permutation index which gives the lowest lexicographical ordering of the nodes
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray>
  STK_INLINE_FUNCTION
  unsigned lexicographical_smallest_permutation(const NodeArray &nodes, bool only_positive_permutations = false) const;

  /// return the permutation index which gives the lowest lexicographical ordering of the nodes that preserves polarity
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray>
  STK_INLINE_FUNCTION
  unsigned lexicographical_smallest_permutation_preserve_polarity(const NodeArray &nodes, const NodeArray &element_nodes) const;

  /// fill the output ordinals with the ordinals that make up the given sub topology
  template <typename OrdinalOutputIterator>
  STK_INLINE_FUNCTION
  void sub_topology_node_ordinals(unsigned sub_rank, unsigned sub_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output nodes with the nodes that make up the given sub topology
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  STK_INLINE_FUNCTION
  void sub_topology_nodes(const NodeArray & nodes, unsigned sub_rank, unsigned sub_ordinal, NodeOutputIterator output_nodes) const;

  /// how many 'sub topologies' does this topology have
  STK_INLINE_FUNCTION
  unsigned num_sub_topology(unsigned sub_rank) const;

  /// what is the topology of the given sub topology
  STK_INLINE_FUNCTION
  topology sub_topology(unsigned sub_rank, unsigned sub_ordinal = 0) const;

  /// fill the output ordinals with the ordinals that make up the given side topology
  template <typename OrdinalOutputIterator>
  STK_INLINE_FUNCTION
  void side_node_ordinals(unsigned side_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output nodes with the nodes that make up the given side topology
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  STK_INLINE_FUNCTION
  void side_nodes(const NodeArray & nodes, unsigned side_ordinal, NodeOutputIterator output_nodes) const;

  /// how many 'side topologies' does this topology have
  STK_INLINE_FUNCTION
  unsigned num_sides() const;


  /// what is the topology of the given side topology
  STK_INLINE_FUNCTION
  topology side_topology(unsigned side_ordinal = 0) const;

  STK_INLINE_FUNCTION
  bool is_superelement() const;

  STK_INLINE_FUNCTION
  bool is_superface() const;

  STK_INLINE_FUNCTION
  bool is_superedge() const;

  STK_INLINE_FUNCTION
  bool is_super_topology() const;

  //***************************************************************************
  //cast to integer type
  //***************************************************************************

  /// implicit cast to topology_t enum type
  STK_INLINE_FUNCTION
  operator topology_t() const
  { return m_value; }

  /// return topology_t enum type
  STK_INLINE_FUNCTION
  topology_t operator()() const
  { return m_value; }

  /// return topology_t enum type
  STK_INLINE_FUNCTION
  topology_t value() const
  { return m_value; }

  //***************************************************************************
  //constructors
  //***************************************************************************
  /// default construct to invalid
  STK_INLINE_FUNCTION
  topology()
    : m_value(INVALID_TOPOLOGY)
  {}

  /// implicit construct from a topology_t
  STK_INLINE_FUNCTION
  topology(topology_t topo)
    : m_value(topo)
  {}

  //copy constructor needs to be defined for GPU
  STK_INLINE_FUNCTION
  topology(const topology& topo)
    : m_value(topo.m_value)
  {}

  STK_INLINE_FUNCTION
  stk::topology& operator=(const topology& rhs)
  {
    if (&rhs != this) {
      m_value = rhs.m_value;
    }
    return *this;
  }

  //***************************************************************************
  // comparison operator
  //***************************************************************************
  STK_INLINE_FUNCTION
  bool operator==(const topology& rhs) const
  {
    return m_value == rhs.m_value;
  }

  STK_INLINE_FUNCTION
  bool operator==(const topology_t& rhs) const
  {
    return m_value == rhs;
  }

  //***************************************************************************
  // compile time for topology
  // used with apply_functor
  //***************************************************************************
  template <topology_t Topology>
  struct topology_type;

  struct types;

  //***************************************************************************
  //used to convert from a runtime to a compile time topology
  //  Functor is of the form
  //  struct Functor {
  //    typedef ... result_type;
  //
  //    template <topology_t Topology>
  //    result_type operator(topology_type<Topology> t)
  //    { ... }
  //  };
  //***************************************************************************
  template <typename Functor>
  struct apply_host_functor;
  template <typename Functor>
  struct apply_functor;

  //***************************************************************************
  // topology_type conversion constructors
  //***************************************************************************
  template <topology_t Topology>
  STK_INLINE_FUNCTION
  topology(topology_type<Topology> /* t */ )
    : m_value(Topology)
  {}

  //***************************************************************************
  //data member
  //***************************************************************************
  topology_t m_value;
};


//***************************************************************************
//increment and decrement rank_t
//***************************************************************************
STK_INLINE_FUNCTION
topology::rank_t operator++(stk::topology::rank_t &r)
{
  r = static_cast<topology::rank_t>(r+1);
  return r;
}

STK_INLINE_FUNCTION
topology::rank_t operator++(stk::topology::rank_t &r,int)
{
  topology::rank_t tmp = r;
  r = static_cast<topology::rank_t>(r+1);
  return tmp;
}

STK_INLINE_FUNCTION
topology::rank_t operator--(stk::topology::rank_t &r)
{
  r = static_cast<topology::rank_t>(r-1);
  return r;
}

STK_INLINE_FUNCTION
topology::rank_t operator--(stk::topology::rank_t &r,int)
{
  topology::rank_t tmp = r;
  r = static_cast<topology::rank_t>(r-1);
  return tmp;
}

//***************************************************************************
//increment and decrement topology_t
//***************************************************************************
STK_INLINE_FUNCTION
topology::topology_t operator++(stk::topology::topology_t &t)
{
  t = static_cast<topology::topology_t>(t+1);
  return t;
}

STK_INLINE_FUNCTION
topology::topology_t operator++(stk::topology::topology_t &t,int)
{
  topology::topology_t tmp = t;
  t = static_cast<topology::topology_t>(t+1);
  return tmp;
}

STK_INLINE_FUNCTION
topology::topology_t operator--(stk::topology::topology_t &t)
{
  t = static_cast<topology::topology_t>(t-1);
  return t;
}

STK_INLINE_FUNCTION
topology::topology_t operator--(stk::topology::topology_t &t,int)
{
  topology::topology_t tmp = t;
  t = static_cast<topology::topology_t>(t-1);
  return tmp;
}

//***************************************************************************
//increment and decrement topology
//***************************************************************************
STK_INLINE_FUNCTION
topology operator++(topology &t)
{
  ++t.m_value;
  return t;
}

STK_INLINE_FUNCTION
topology operator++(topology &t,int)
{
  topology tmp = t;
  ++t.m_value;
  return tmp;
}

STK_INLINE_FUNCTION
topology operator--(topology &t)
{
  --t.m_value;
  return t;
}

STK_INLINE_FUNCTION
topology operator--(topology &t,int)
{
  topology tmp = t;
  --t.m_value;
  return tmp;
}

//***************************************************************************
//create superelement
//***************************************************************************
STK_INLINE_FUNCTION
topology create_superedge_topology(unsigned num_nodes)
{
  if ( num_nodes < 1u ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPEREDGE_START);
}

STK_INLINE_FUNCTION
topology create_superedge_topology(int num_nodes)
{
  if ( num_nodes < 1 ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPEREDGE_START);
}

STK_INLINE_FUNCTION
topology create_superface_topology(unsigned num_nodes)
{
  if ( num_nodes < 1u ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPERFACE_START);
}

STK_INLINE_FUNCTION
topology create_superface_topology(int num_nodes)
{
  if ( num_nodes < 1 ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPERFACE_START);
}

STK_INLINE_FUNCTION
topology create_superelement_topology(unsigned num_nodes)
{
  if ( num_nodes < 1u ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPERELEMENT_START);
}

STK_INLINE_FUNCTION
topology create_superelement_topology(int num_nodes)
{
  if ( num_nodes < 1 ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPERELEMENT_START);
}

//***************************************************************************
//output operators and verbose topology printing
//***************************************************************************
std::ostream & operator<<(std::ostream &out, topology::rank_t r);
std::ostream & operator<<(std::ostream &out, topology t);
void verbose_print_topology(std::ostream &out, topology t);

bool isTriangleElement (topology topo);
bool isQuadrilateralElement (topology topo);
bool isTetrahedronElement (topology topo);
bool isHexahedronElement (topology topo);
bool is_solid_element(stk::topology t);

bool is_quad_side(topology topo);
bool is_tri_side(topology topo);

} //namespace stk

#endif // TOPOLOGY_DECL_HPP

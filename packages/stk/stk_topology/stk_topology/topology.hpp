#ifndef STKTOPOLOGY_TOPOLOGY_HPP
#define STKTOPOLOGY_TOPOLOGY_HPP

#include <stk_topology/topology_detail/macros.hpp>

#include <utility>
#include <iosfwd>

namespace stk {


struct topology
{

  enum rank_t
  {
      BEGIN_RANK
    , NODE_RANK = BEGIN_RANK
    , EDGE_RANK
    , FACE_RANK
    , ELEM_RANK, ELEMENT_RANK = ELEM_RANK
    , CONSTRAINT_RANK
    , END_RANK
    , NUM_RANKS = END_RANK
    , INVALID_RANK = ~0
  };

  //To add new topologies consult the toolkit team
  enum topology_t
  {
      BEGIN_TOPOLOGY
    //NODE_RANK
    , NODE = BEGIN_TOPOLOGY
    //EDGE_RANK
    , LINE_2
    , LINE_3
    //FACE_RANK
    , TRI_3, TRIANGLE_3 = TRI_3
    , TRI_4, TRIANGLE_4 = TRI_4
    , TRI_6, TRIANGLE_6 = TRI_6
    , QUAD_4, QUADRILATERAL_4 = QUAD_4
    , QUAD_8, QUADRILATERAL_8 = QUAD_8
    , QUAD_9, QUADRILATERAL_9 = QUAD_9
    //ELEMENT_RANK
    , PARTICLE
    , LINE_2_1D
    , LINE_3_1D
    , BEAM_2
    , BEAM_3
    , SHELL_LINE_2
    , SHELL_LINE_3
    , TRI_3_2D, TRIANGLE_3_2D = TRI_3_2D
    , TRI_4_2D, TRIANGLE_4_2D = TRI_4_2D
    , TRI_6_2D, TRIANGLE_6_2D = TRI_6_2D
    , QUAD_4_2D, QUADRILATERAL_4_2D = QUAD_4_2D
    , QUAD_8_2D, QUADRILATERAL_8_2D = QUAD_8_2D
    , QUAD_9_2D, QUADRILATERAL_9_2D = QUAD_9_2D
    , SHELL_TRI_3, SHELL_TRIANGLE_3 = SHELL_TRI_3
    , SHELL_TRI_4, SHELL_TRIANGLE_4 = SHELL_TRI_4
    , SHELL_TRI_6, SHELL_TRIANGLE_6 = SHELL_TRI_6
    , SHELL_QUAD_4, SHELL_QUADRILATERAL_4 = SHELL_QUAD_4
    , SHELL_QUAD_8, SHELL_QUADRILATERAL_8 = SHELL_QUAD_8
    , SHELL_QUAD_9, SHELL_QUADRILATERAL_9 = SHELL_QUAD_9
    , TET_4,  TETRAHEDRON_4  = TET_4
    , TET_8,  TETRAHEDRON_8  = TET_8
    , TET_10, TETRAHEDRON_10 = TET_10
    , TET_11, TETRAHEDRON_11 = TET_11
    , PYRAMID_5
    , PYRAMID_13
    , PYRAMID_14
    , WEDGE_6
    , WEDGE_15
    , WEDGE_18
    , HEX_8,  HEXAHEDRON_8  = HEX_8
    , HEX_20, HEXAHEDRON_20 = HEX_20
    , HEX_27, HEXAHEDRON_27 = HEX_27
    , END_TOPOLOGY
    , INVALID_TOPOLOGY = END_TOPOLOGY
    , NUM_TOPOLOGIES = END_TOPOLOGY
  };

  static const char * rank_names[];     // indexed by rank_t
  static const char * topology_names[]; // indexed by topology_t

  //***************************************************************************
  //member functions
  //***************************************************************************

  /// is this topology valid
  bool is_valid() const { return m_value != INVALID_TOPOLOGY; }

  /// get the name of this topology
  const char * name() const { return topology_names[m_value]; }

  /// does this topology have homogeneous edges
  bool has_homogeneous_edges() const;

  /// does this topology have homogeneous faces
  bool has_homogeneous_faces() const;

  /// does this topology have homogeneous faces
  bool has_homogeneous_sides() const;

  /// is this topology a shell topology (i.e. an element with only two sides)
  bool is_shell() const;

  /// what is the rank of this topology
  rank_t rank() const;

  /// what is the side rank of this topology
  rank_t side_rank() const;

  /// what is the topological dimension of this topology
  int dimension() const;

  /// how many nodes define this topology
  int num_nodes() const;

  /// how many nodes are vertices
  int num_vertices() const;

  /// how many edges does this topology have
  int num_edges() const;

  /// how many faces does this topology have
  int num_faces() const;

  /// how many sides does this topology have
  int num_sides() const;

  /// how many different node permutations does this topology have
  int num_permutations() const;

  /// how many different positive node permutations does this topology have
  int num_positive_permutations() const;

  /// is this topology defined on the given spatial dimension
  bool defined_on_spatial_dimension(int spatial_dimension) const;

  /// what is the base topology (i.e. topology where num_nodes == num_vertices)
  topology base() const;

  /// what is the topology of the given edge
  topology edge_topology(int edge_ordinal = 0) const;

  /// what is the topology of the given face
  topology face_topology(int face_ordinal = 0) const;

  /// what is the topology of the given side
  topology side_topology(int side_ordinal = 0) const;

  /// fill the output ordinals with the ordinals that make up the given edge
  template <typename OrdinalOutputIterator>
  void edge_node_ordinals(int edge_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output ordinals with the ordinals that make up the given face
  template <typename OrdinalOutputIterator>
  void face_node_ordinals(int face_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output ordinals with the ordinals that make up the given side
  template <typename OrdinalOutputIterator>
  void side_node_ordinals(int side_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output ordinals with the ordinals that make up the given permutation
  template <typename OrdinalOutputIterator>
  void permutation_node_ordinals(int permutation_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output nodes with the nodes that make up the given edge
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  void edge_nodes(const NodeArray & nodes, int edge_ordinal, NodeOutputIterator output_nodes) const;

  /// fill the output nodes with the nodes that make up the given face
  /// input 'nodes' is expected to be of length num_nodes.
 template <typename NodeArray, typename NodeOutputIterator>
  void face_nodes(const NodeArray & nodes, int face_ordinal, NodeOutputIterator output_nodes) const;

  /// fill the output nodes with the nodes that make up the given side
 /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  void side_nodes(const NodeArray & nodes, int side_ordinal, NodeOutputIterator output_nodes) const;

  /// fill the output nodes with the nodes that make up the given permutation
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  void permutation_nodes(const NodeArray & nodes, int permutation_ordinal, NodeOutputIterator output_nodes) const;

  /// do the two arrays define equivalent entities (same nodes, but maybe a different permutation)
  /// return a pair<bool, permutation_number> bool and permutation number from a to b
  template <typename NodeArrayA, typename NodeArrayB>
  std::pair<bool,int> equivalent(const NodeArrayA & a, const NodeArrayB & b) const;


  /// return the permutation index which gives the lowest lexicographical ordering of the nodes
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray>
  int lexicographical_smallest_permutation(const NodeArray &nodes, bool only_positive_permutations = false) const;

  /// fill the output ordinals with the ordinals that make up the given sub topology
  template <typename OrdinalOutputIterator>
  STKTOPOLOGY_INLINE_FUNCTION
  void sub_topology_node_ordinals(int sub_rank, int sub_ordinal, OrdinalOutputIterator output_ordinals) const
  {
    switch(sub_rank)
    {
    case NODE_RANK: *output_ordinals = sub_ordinal;                    break;
    case EDGE_RANK: edge_node_ordinals(sub_ordinal, output_ordinals);  break;
    case FACE_RANK: face_node_ordinals(sub_ordinal, output_ordinals);  break;
    default: break;
    }
  }

  /// fill the output nodes with the nodes that make up the given sub topology
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  STKTOPOLOGY_INLINE_FUNCTION
  void sub_topology_nodes(const NodeArray & nodes, int sub_rank, int sub_ordinal, NodeOutputIterator output_nodes) const
  {
    switch(sub_rank)
    {
    case NODE_RANK: *output_nodes = nodes[sub_ordinal];                    break;
    case EDGE_RANK: edge_nodes(nodes, sub_ordinal, output_nodes);  break;
    case FACE_RANK: face_nodes(nodes, sub_ordinal, output_nodes);  break;
    default: break;
    }
  }

  /// how many 'sub topologies' does this topology have
  STKTOPOLOGY_INLINE_FUNCTION
  int num_sub_topology(int sub_rank) const
  {
    switch(sub_rank)
    {
    case NODE_RANK: return num_nodes();
    case EDGE_RANK: return num_edges();
    case FACE_RANK: return num_faces();
    default: break;
    }
    return 0;

  }


  /// what is the topology of the given sub topology
  STKTOPOLOGY_INLINE_FUNCTION
  topology sub_topology(int sub_rank, int sub_ordinal = 0) const
  {
    switch(sub_rank)
    {
    case NODE_RANK: return NODE;
    case EDGE_RANK: return edge_topology(sub_ordinal);
    case FACE_RANK: return face_topology(sub_ordinal);
    default: break;
    }
    return INVALID_TOPOLOGY;
  }

  //***************************************************************************
  //cast to integer type
  //***************************************************************************

  /// implicit cast to topology_t enum type
  STKTOPOLOGY_INLINE_FUNCTION
  operator topology_t() const
  { return m_value; }

  /// return topology_t enum type
  STKTOPOLOGY_INLINE_FUNCTION
  topology_t operator()() const
  { return m_value; }

  /// return topology_t enum type
  STKTOPOLOGY_INLINE_FUNCTION
  topology_t value() const
  { return m_value; }

  //***************************************************************************
  //constructors
  //***************************************************************************
  /// default construct to invalid
  topology()
    : m_value(INVALID_TOPOLOGY)
  {}

  /// implicit construct from a topology_t
  topology(topology_t topo)
    : m_value(topo)
  {}

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
  struct apply_functor;

  //***************************************************************************
  // topology_type conversion constructors
  //***************************************************************************
  template <topology_t Topology>
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
STKTOPOLOGY_INLINE_FUNCTION
topology::rank_t operator++(topology::rank_t &r)
{
  r = static_cast<topology::rank_t>(r+1);
  return r;
}

STKTOPOLOGY_INLINE_FUNCTION
topology::rank_t operator++(topology::rank_t &r,int)
{
  topology::rank_t tmp = r;
  r = static_cast<topology::rank_t>(r+1);
  return tmp;
}

STKTOPOLOGY_INLINE_FUNCTION
topology::rank_t operator--(topology::rank_t &r)
{
  r = static_cast<topology::rank_t>(r-1);
  return r;
}

STKTOPOLOGY_INLINE_FUNCTION
topology::rank_t operator--(topology::rank_t &r,int)
{
  topology::rank_t tmp = r;
  r = static_cast<topology::rank_t>(r-1);
  return tmp;
}

//***************************************************************************
//increment and decrement topology_t
//***************************************************************************
STKTOPOLOGY_INLINE_FUNCTION
topology::topology_t operator++(topology::topology_t &t)
{
  t = static_cast<topology::topology_t>(t+1);
  return t;
}

STKTOPOLOGY_INLINE_FUNCTION
topology::topology_t operator++(topology::topology_t &t,int)
{
  topology::topology_t tmp = t;
  t = static_cast<topology::topology_t>(t+1);
  return tmp;
}

STKTOPOLOGY_INLINE_FUNCTION
topology::topology_t operator--(topology::topology_t &t)
{
  t = static_cast<topology::topology_t>(t-1);
  return t;
}

STKTOPOLOGY_INLINE_FUNCTION
topology::topology_t operator--(topology::topology_t &t,int)
{
  topology::topology_t tmp = t;
  t = static_cast<topology::topology_t>(t-1);
  return tmp;
}

//***************************************************************************
//increment and decrement topology
//***************************************************************************
STKTOPOLOGY_INLINE_FUNCTION
topology operator++(topology &t)
{
  ++t.m_value;
  return t;
}

STKTOPOLOGY_INLINE_FUNCTION
topology operator++(topology &t,int)
{
  topology tmp = t;
  ++t.m_value;
  return tmp;
}

STKTOPOLOGY_INLINE_FUNCTION
topology operator--(topology &t)
{
  --t.m_value;
  return t;
}

STKTOPOLOGY_INLINE_FUNCTION
topology operator--(topology &t,int)
{
  topology tmp = t;
  --t.m_value;
  return tmp;
}


//***************************************************************************
//output operators and verbose topology printing
//***************************************************************************
std::ostream & operator<<(std::ostream &out, topology::rank_t r);
std::ostream & operator<<(std::ostream &out, topology t);

void verbose_print_topology(std::ostream &out, topology t);

} //namespace stk

#include <stk_topology/topology_type.tcc>
#include <stk_topology/types.tcc>

#include <stk_topology/apply_functor.tcc>
#include <stk_topology/topology.tcc>

#endif //STKTOPOLOGY_TOPOLOGY_HPP


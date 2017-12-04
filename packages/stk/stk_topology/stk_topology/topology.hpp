// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STKTOPOLOGY_TOPOLOGY_HPP
#define STKTOPOLOGY_TOPOLOGY_HPP

#include <stk_util/stk_config.h>

#include <string>

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
    , INVALID_RANK = 256
    , FORCE_RANK_TO_UNSIGNED = 256
  };

  //To add new topologies consult the toolkit team
  enum topology_t
  {
      INVALID_TOPOLOGY
    , BEGIN_TOPOLOGY
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
    , PARTICLE, BEGIN_ELEMENT_RANK = PARTICLE
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
    , NUM_TOPOLOGIES = END_TOPOLOGY - BEGIN_TOPOLOGY
    , SUPEREDGE_START = END_TOPOLOGY+1
    , SUPEREDGE_END = SUPEREDGE_START + 1000
    , SUPERFACE_START = SUPEREDGE_END+1
    , SUPERFACE_END = SUPERFACE_START + 1000
    , SUPERELEMENT_START = SUPERFACE_END+1
    , FORCE_TOPOLOGY_TO_UNSIGNED = ~0U // max unsigned int
  };

  //***************************************************************************
  //member functions
  //***************************************************************************

  /// is this topology valid
  bool is_valid() const { return m_value != INVALID_TOPOLOGY; }

  /// get the name of this topology
  std::string name() const;

  /// does this topology have homogeneous faces
  bool has_homogeneous_faces() const;

  /// is this topology a shell topology (i.e. an element with only two sides)
  bool is_shell() const;

  /// what is the rank of this topology
  rank_t rank() const;

  /// what is the side rank of this topology
  rank_t side_rank() const;

  /// what is the topological dimension of this topology
  unsigned dimension() const;

  /// how many nodes define this topology
  STK_FUNCTION
  unsigned num_nodes() const;

  /// how many nodes are vertices
  unsigned num_vertices() const;

  /// how many edges does this topology have
  unsigned num_edges() const;

  /// how many faces does this topology have
  unsigned num_faces() const;

  /// how many different node permutations does this topology have
  unsigned num_permutations() const;

  /// how many different positive node permutations does this topology have
  unsigned num_positive_permutations() const;

  bool is_positive_polarity(unsigned permutation_ordinal) const {
      return (permutation_ordinal < num_positive_permutations());
  }

  /// is this topology defined on the given spatial dimension
  bool defined_on_spatial_dimension(unsigned spatial_dimension) const;

  /// what is the base topology (i.e. topology where num_nodes == num_vertices)
  topology base() const;

  /// what is the topology of the given edge
  topology edge_topology() const;

  /// what is the topology of the given face
  topology face_topology(unsigned face_ordinal) const;

  /// fill the output ordinals with the ordinals that make up the given edge
  template <typename OrdinalOutputIterator>
  void edge_node_ordinals(unsigned edge_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output ordinals with the ordinals that make up the given face
  template <typename OrdinalOutputIterator>
  void face_node_ordinals(unsigned face_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output ordinals with the ordinals that make up the given permutation
  template <typename OrdinalOutputIterator>
  void permutation_node_ordinals(unsigned permutation_ordinal, OrdinalOutputIterator output_ordinals) const;

  /// fill the output nodes with the nodes that make up the given edge
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  void edge_nodes(const NodeArray & nodes, unsigned edge_ordinal, NodeOutputIterator output_nodes) const;

  /// fill the output nodes with the nodes that make up the given face
  /// input 'nodes' is expected to be of length num_nodes.
 template <typename NodeArray, typename NodeOutputIterator>
  void face_nodes(const NodeArray & nodes, unsigned face_ordinal, NodeOutputIterator output_nodes) const;

  /// fill the output nodes with the nodes that make up the given permutation
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  void permutation_nodes(const NodeArray & nodes, unsigned permutation_ordinal, NodeOutputIterator output_nodes) const;

  /// do the two arrays define equivalent entities (same nodes, but maybe a different permutation)
  /// return a pair<bool, permutation_number> bool and permutation number from a to b
  template <typename NodeArrayA, typename NodeArrayB>
  std::pair<bool,unsigned> equivalent(const NodeArrayA & a, const NodeArrayB & possible_permutation_of_a) const;

  /// return the permutation index which gives the lowest lexicographical ordering of the nodes
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray>
  unsigned lexicographical_smallest_permutation(const NodeArray &nodes, bool only_positive_permutations = false) const;

  /// return the permutation index which gives the lowest lexicographical ordering of the nodes that preserves polarity
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray>
  unsigned lexicographical_smallest_permutation_preserve_polarity(const NodeArray &nodes, const NodeArray &element_nodes) const;

  /// fill the output ordinals with the ordinals that make up the given sub topology
#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
  template <typename OrdinalOutputIterator>
  void sub_topology_node_ordinals(unsigned sub_rank, unsigned sub_ordinal, OrdinalOutputIterator output_ordinals) const
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
#ifdef __CUDACC__
#pragma hd_warning_disable
#endif
  template <typename NodeArray, typename NodeOutputIterator>
  void sub_topology_nodes(const NodeArray & nodes, unsigned sub_rank, unsigned sub_ordinal, NodeOutputIterator output_nodes) const
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
  unsigned num_sub_topology(unsigned sub_rank) const
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
  topology sub_topology(unsigned sub_rank, unsigned sub_ordinal = 0) const
  {
    switch(sub_rank)
    {
    case NODE_RANK: return NODE;
    case EDGE_RANK: return edge_topology();
    case FACE_RANK: return face_topology(sub_ordinal);
    default: break;
    }
    return INVALID_TOPOLOGY;
  }



  /// fill the output ordinals with the ordinals that make up the given side topology
  template <typename OrdinalOutputIterator>
  void side_node_ordinals(unsigned side_ordinal, OrdinalOutputIterator output_ordinals) const
  {
    sub_topology_node_ordinals( side_rank(), side_ordinal, output_ordinals);
  }

  /// fill the output nodes with the nodes that make up the given side topology
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  void side_nodes(const NodeArray & nodes, unsigned side_ordinal, NodeOutputIterator output_nodes) const
  {
    sub_topology_nodes( nodes, side_rank(), side_ordinal, output_nodes);
  }

  /// how many 'side topologies' does this topology have
  unsigned num_sides() const
  {
    unsigned num_sides_out = 0u;
    if (side_rank() != INVALID_RANK) {
      num_sides_out = side_rank() > NODE_RANK? num_sub_topology(side_rank()) : num_vertices();
    }
    return num_sides_out;
  }


  /// what is the topology of the given side topology
  topology side_topology(unsigned side_ordinal = 0) const
  {
    return sub_topology(side_rank(), side_ordinal);
  }

  bool is_superelement() const
  {
    return m_value > SUPERELEMENT_START;
  }

  STK_FUNCTION
  bool is_superface() const
  {
    return m_value > SUPERFACE_START && m_value < SUPERFACE_END;
  }

  STK_FUNCTION
  bool is_superedge() const
  {
    return m_value > SUPEREDGE_START && m_value < SUPEREDGE_END;
  }

  bool is_super_topology() const
  {
    return is_superelement() || is_superface() || is_superedge();
  }

  //***************************************************************************
  //cast to integer type
  //***************************************************************************

  /// implicit cast to topology_t enum type
  operator topology_t() const
  { return m_value; }

  /// return topology_t enum type
  topology_t operator()() const
  { return m_value; }

  /// return topology_t enum type
  topology_t value() const
  { return m_value; }

  //***************************************************************************
  //constructors
  //***************************************************************************
  /// default construct to invalid
  STK_FUNCTION
  topology()
    : m_value(INVALID_TOPOLOGY)
  {}

  /// implicit construct from a topology_t
  STK_FUNCTION
  topology(topology_t topo)
    : m_value(topo)
  {}

  //copy constructor needs to be defined for GPU
  STK_FUNCTION
  topology(const topology& topo)
    : m_value(topo.m_value)
  {}

  //***************************************************************************
  // comparison operator
  //***************************************************************************
  STK_FUNCTION
  bool operator==(const topology& rhs) const
  {
    return m_value == rhs.m_value;
  }

  STK_FUNCTION
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
inline
topology::rank_t operator++(stk::topology::rank_t &r)
{
  r = static_cast<topology::rank_t>(r+1);
  return r;
}

inline
topology::rank_t operator++(stk::topology::rank_t &r,int)
{
  topology::rank_t tmp = r;
  r = static_cast<topology::rank_t>(r+1);
  return tmp;
}

inline
topology::rank_t operator--(stk::topology::rank_t &r)
{
  r = static_cast<topology::rank_t>(r-1);
  return r;
}

inline
topology::rank_t operator--(stk::topology::rank_t &r,int)
{
  topology::rank_t tmp = r;
  r = static_cast<topology::rank_t>(r-1);
  return tmp;
}

//***************************************************************************
//increment and decrement topology_t
//***************************************************************************
inline
topology::topology_t operator++(stk::topology::topology_t &t)
{
  t = static_cast<topology::topology_t>(t+1);
  return t;
}

inline
topology::topology_t operator++(stk::topology::topology_t &t,int)
{
  topology::topology_t tmp = t;
  t = static_cast<topology::topology_t>(t+1);
  return tmp;
}

inline
topology::topology_t operator--(stk::topology::topology_t &t)
{
  t = static_cast<topology::topology_t>(t-1);
  return t;
}

inline
topology::topology_t operator--(stk::topology::topology_t &t,int)
{
  topology::topology_t tmp = t;
  t = static_cast<topology::topology_t>(t-1);
  return tmp;
}

//***************************************************************************
//increment and decrement topology
//***************************************************************************
inline
topology operator++(topology &t)
{
  ++t.m_value;
  return t;
}

inline
topology operator++(topology &t,int)
{
  topology tmp = t;
  ++t.m_value;
  return tmp;
}

inline
topology operator--(topology &t)
{
  --t.m_value;
  return t;
}

inline
topology operator--(topology &t,int)
{
  topology tmp = t;
  --t.m_value;
  return tmp;
}

//***************************************************************************
//create superelement
//***************************************************************************
inline
topology create_superedge_topology(unsigned num_nodes)
{
  if ( num_nodes < 1u ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPEREDGE_START);
}

inline
topology create_superedge_topology(int num_nodes)
{
  if ( num_nodes < 1 ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPEREDGE_START);
}

inline
topology create_superface_topology(unsigned num_nodes)
{
  if ( num_nodes < 1u ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPERFACE_START);
}

inline
topology create_superface_topology(int num_nodes)
{
  if ( num_nodes < 1 ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPERFACE_START);
}

inline
topology create_superelement_topology(unsigned num_nodes)
{
  if ( num_nodes < 1u ) return topology::INVALID_TOPOLOGY;
  return static_cast<topology::topology_t>(num_nodes + topology::SUPERELEMENT_START);
}

inline
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
bool isTriangle (topology topo);
bool isQuadrilateral (topology topo);
bool isTetrahedron (topology topo);
bool isHexahedron (topology topo);
bool is_solid_element(stk::topology t);

} //namespace stk

#include <stk_topology/topology_type.tcc>
#include <stk_topology/types.tcc>

#include <stk_topology/apply_functor.tcc>
#include <stk_topology/topology.tcc>

#endif //STKTOPOLOGY_TOPOLOGY_HPP

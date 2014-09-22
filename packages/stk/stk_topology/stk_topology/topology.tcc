# Copyright (c) 2013, Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Governement retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of Sandia Corporation nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#ifndef STKTOPOLOGY_TOPOLOGY_TCC
#define STKTOPOLOGY_TOPOLOGY_TCC

#define STKTOPOLOGY_ORDINAL_NODES_MEMBER(name)                                     \
  namespace stk { namespace topology_detail {                                      \
  template <typename OrdinalOutputIterator>                                        \
  struct name##_impl {                                                             \
    typedef void result_type;                                                      \
    BOOST_GPU_ENABLED                                                              \
    name##_impl(unsigned ordinal, OrdinalOutputIterator output_ordinals)           \
      : m_ordinal(ordinal)                                                         \
      , m_output_ordinals(output_ordinals)                                         \
    {}                                                                             \
    template <typename Topology>                                                   \
    BOOST_GPU_ENABLED                                                              \
    void operator()(Topology) const                                         \
    { Topology::name(m_ordinal,m_output_ordinals); }                               \
    unsigned                   m_ordinal;                                          \
    OrdinalOutputIterator m_output_ordinals;                                       \
  };                                                                               \
  }} /*namespace stk::topology_detail*/                                            \
  namespace stk {                                                                  \
  template <typename OrdinalOutputIterator>                                        \
  BOOST_GPU_ENABLED inline                                                         \
  void topology::name( unsigned ordinal, OrdinalOutputIterator output_ordinals) const   \
  { typedef topology_detail::name##_impl<OrdinalOutputIterator> functor;           \
    functor f(ordinal,output_ordinals);                                            \
    topology::apply_functor< functor > apply( f );                                 \
    apply(m_value);                                                                \
  }                                                                                \
  } /*namespace stk*/

#define STKTOPOLOGY_NODES_MEMBER(name)                                             \
  namespace stk { namespace topology_detail {                                      \
  template <typename NodeArray, typename NodeOutputIterator>                       \
  struct name##_impl {                                                             \
    typedef void result_type;                                                      \
    BOOST_GPU_ENABLED                                                              \
    name##_impl(   const NodeArray &nodes                                          \
                 , unsigned ordinal                                                \
                 , NodeOutputIterator output_ordinals)                             \
      : m_nodes(nodes)                                                             \
      , m_ordinal(ordinal)                                                         \
      , m_output_ordinals(output_ordinals)                                         \
    {}                                                                             \
    template <typename Topology>                                                   \
    BOOST_GPU_ENABLED                                                              \
    void operator()(Topology) const                                         \
    { Topology::name(m_nodes,m_ordinal,m_output_ordinals); }                \
    const NodeArray    & m_nodes;                                                  \
    unsigned                  m_ordinal;                                           \
    NodeOutputIterator   m_output_ordinals;                                        \
  };                                                                               \
  }} /*namespace stk::topology_detail*/                                            \
  namespace stk {                                                                  \
  template <typename NodeArray, typename NodeOutputIterator>                       \
  BOOST_GPU_ENABLED inline                                                         \
  void topology::name(   const NodeArray & nodes                                   \
                       , unsigned ordinal                                          \
                       , NodeOutputIterator output_ordinals) const                 \
  { typedef topology_detail::name##_impl<NodeArray,NodeOutputIterator> functor;    \
    functor f(nodes,ordinal,output_ordinals);                                      \
    topology::apply_functor< functor > apply( f );                                 \
    apply(m_value);                                                                \
  }                                                                                \
  } /*namespace stk*/

#define STKTOPOLOGY_SIMPLE_MEMBER(name,result)         \
  namespace stk { namespace topology_detail {          \
  struct name##_impl {                                 \
    typedef result result_type;                        \
    template <typename Topology>                       \
    BOOST_GPU_ENABLED                                  \
    result_type operator()(Topology) const             \
    { return Topology::name; }                         \
  };                                                   \
  }} /*namespace stk::topology_detail*/                \
  namespace stk {                                      \
  BOOST_GPU_ENABLED inline                             \
  result topology::name() const                        \
  { typedef topology_detail::name##_impl functor;      \
    topology::apply_functor< functor > apply;          \
    return apply(m_value);                             \
  }                                                    \
  } /*namespace stk*/

#define STKTOPOLOGY_ORDINAL_MEMBER(name,result)        \
  namespace stk { namespace topology_detail {          \
  struct name##_impl {                                 \
    typedef result result_type;                        \
    BOOST_GPU_ENABLED                                  \
    name##_impl(unsigned ordinal)                      \
      : m_ordinal(ordinal)                             \
    {}                                                 \
    template <typename Topology>                       \
    BOOST_GPU_ENABLED                                  \
    result_type operator()(Topology) const             \
    { return Topology::name(m_ordinal); }              \
    unsigned m_ordinal;                                \
  };                                                   \
  }} /*namespace stk::topology_detail*/                \
  namespace stk {                                      \
  BOOST_GPU_ENABLED inline                             \
  result topology::name(unsigned ordinal) const        \
  { typedef topology_detail::name##_impl functor;      \
    functor f(ordinal);                                \
    topology::apply_functor< functor > apply( f );     \
    return apply(m_value);                             \
  }                                                    \
  } /*namespace stk*/

STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_faces,bool)
STKTOPOLOGY_SIMPLE_MEMBER(is_shell,bool)
STKTOPOLOGY_SIMPLE_MEMBER(side_rank,stk::topology::rank_t)
STKTOPOLOGY_SIMPLE_MEMBER(dimension,unsigned)
STKTOPOLOGY_SIMPLE_MEMBER(num_vertices,unsigned)
STKTOPOLOGY_SIMPLE_MEMBER(num_edges,unsigned)
STKTOPOLOGY_SIMPLE_MEMBER(num_faces,unsigned)
STKTOPOLOGY_SIMPLE_MEMBER(num_permutations,unsigned)
STKTOPOLOGY_SIMPLE_MEMBER(num_positive_permutations,unsigned)
STKTOPOLOGY_SIMPLE_MEMBER(base,stk::topology)
STKTOPOLOGY_SIMPLE_MEMBER(edge_topology,stk::topology)

STKTOPOLOGY_ORDINAL_MEMBER(defined_on_spatial_dimension,bool)
STKTOPOLOGY_ORDINAL_MEMBER(face_topology,stk::topology)

STKTOPOLOGY_ORDINAL_NODES_MEMBER(edge_node_ordinals)
STKTOPOLOGY_ORDINAL_NODES_MEMBER(face_node_ordinals)
STKTOPOLOGY_ORDINAL_NODES_MEMBER(permutation_node_ordinals)

STKTOPOLOGY_NODES_MEMBER(edge_nodes)
STKTOPOLOGY_NODES_MEMBER(face_nodes)
STKTOPOLOGY_NODES_MEMBER(permutation_nodes)

#undef STKTOPOLOGY_SIMPLE_MEMBER
#undef STKTOPOLOGY_ORDINAL_MEMBER

#undef STKTOPOLOGY_ORDINAL_NODES_MEMBER
#undef STKTOPOLOGY_NODES_MEMBER

namespace stk { namespace topology_detail {

struct num_nodes_impl {
  typedef unsigned result_type;
  template <typename Topology>
  BOOST_GPU_ENABLED
  result_type operator()(Topology) const
  { return Topology::num_nodes; }
};

struct rank_impl {
  typedef topology::rank_t result_type;
  template <typename Topology>
  BOOST_GPU_ENABLED
  result_type operator()(Topology) const
  { return Topology::rank; }
};

template <typename NodeArrayA, typename NodeArrayB>
struct equivalent_impl {
  typedef std::pair<bool,unsigned> result_type;

  BOOST_GPU_ENABLED
  equivalent_impl( const NodeArrayA &a , const NodeArrayB &b )
    : m_a(a), m_b(b)
  {}

  template <typename Topology>
  BOOST_GPU_ENABLED
  result_type operator()(Topology) const
  { return Topology::equivalent(m_a, m_b); }

  const NodeArrayA & m_a;
  const NodeArrayB & m_b;
};

template <typename NodeArray>
struct lexicographical_smallest_permutation_impl {
  typedef unsigned result_type;

  BOOST_GPU_ENABLED
  lexicographical_smallest_permutation_impl( const NodeArray &nodes , bool only_positive_permutations )
    : m_nodes(nodes), m_only_positive_permutations(only_positive_permutations)
  {}

  template <typename Topology>
  BOOST_GPU_ENABLED
  result_type operator()(Topology) const
  { return Topology::lexicographical_smallest_permutation(m_nodes, m_only_positive_permutations); }

  const NodeArray & m_nodes;
  bool              m_only_positive_permutations;
};

}} /*namespace stk::topology_detail*/

namespace stk {

BOOST_GPU_ENABLED inline
unsigned topology::num_nodes() const
{
  typedef topology_detail::num_nodes_impl functor;
  topology::apply_functor< functor > apply;
  return m_value < SUPERELEMENT_START ? apply(m_value) : m_value - SUPERELEMENT_START;
}

BOOST_GPU_ENABLED inline
topology::rank_t topology::rank() const
{
  typedef topology_detail::rank_impl functor;
  topology::apply_functor< functor > apply;
  return m_value < SUPERELEMENT_START ? apply(m_value) : ( m_value > SUPERELEMENT_START ? topology::ELEMENT_RANK : topology::INVALID_RANK);
}

template <typename NodeArrayA, typename NodeArrayB>
BOOST_GPU_ENABLED inline
std::pair<bool,unsigned> topology::equivalent( const NodeArrayA &a, const NodeArrayB &b) const
{ typedef topology_detail::equivalent_impl<NodeArrayA,NodeArrayB> functor;
  functor f(a,b);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
BOOST_GPU_ENABLED inline
unsigned topology::lexicographical_smallest_permutation( const NodeArray &nodes, bool only_positive_permutations) const
{
  typedef topology_detail::lexicographical_smallest_permutation_impl< NodeArray > functor;
  functor f(nodes, only_positive_permutations);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

} /*namespace stk*/

#endif //STKTOPOLOGY_TOPOLOGY_TCC


#ifndef STKTOPOLOGY_TOPOLOGY_TCC
#define STKTOPOLOGY_TOPOLOGY_TCC

#define STKTOPOLOGY_SIMPLE_MEMBER(name,result)                                     \
  namespace stk { namespace detail {                                               \
  struct name##_impl {                                                             \
    typedef result result_type;                                                    \
    template <typename Topology>                                                   \
    STKTOPOLOGY_INLINE_FUNCTION                                                    \
    result_type operator()(Topology) const                                         \
    { return Topology::name; }                                                     \
  };                                                                               \
  }} /*namespace stk::detail*/                                                     \
  namespace stk {                                                                  \
  STKTOPOLOGY_INLINE_FUNCTION                                                      \
  result topology::name() const                                                    \
  { topology::apply_functor< detail::name##_impl > apply;                          \
    return apply(m_value);                                                         \
  }                                                                                \
  } /*namespace stk*/

#define STKTOPOLOGY_ORDINAL_MEMBER(name,result)                                    \
  namespace stk { namespace detail {                                               \
  struct name##_impl {                                                             \
    typedef result result_type;                                                    \
    STKTOPOLOGY_INLINE_FUNCTION                                                    \
    name##_impl(int ordinal)                                                       \
      : m_ordinal(ordinal)                                                         \
    {}                                                                             \
    template <typename Topology>                                                   \
    STKTOPOLOGY_INLINE_FUNCTION                                                    \
    result_type operator()(Topology) const                                         \
    { return Topology::name(m_ordinal); }                                          \
    int m_ordinal;                                                                 \
  };                                                                               \
  }} /*namespace stk::detail*/                                                     \
  namespace stk {                                                                  \
  STKTOPOLOGY_INLINE_FUNCTION                                                      \
  result topology::name(int ordinal) const                                         \
  { detail::name##_impl f(ordinal);                                                \
    topology::apply_functor< detail::name##_impl > apply( f );                     \
    return apply(m_value);                                                         \
  }                                                                                \
  } /*namespace stk*/

#define STKTOPOLOGY_ORDINAL_NODES_MEMBER(name)                                     \
  namespace stk { namespace detail {                                               \
  template <typename OrdinalOutputIterator>                                        \
  struct name##_impl {                                                             \
    typedef void result_type;                                                      \
    STKTOPOLOGY_INLINE_FUNCTION                                                    \
    name##_impl(int ordinal, OrdinalOutputIterator output_ordinals)                \
      : m_ordinal(ordinal)                                                         \
      , m_output_ordinals(output_ordinals)                                         \
    {}                                                                             \
    template <typename Topology>                                                   \
    STKTOPOLOGY_INLINE_FUNCTION                                                    \
    result_type operator()(Topology) const                                         \
    { return Topology::name(m_ordinal,m_output_ordinals); }                        \
    int                   m_ordinal;                                               \
    OrdinalOutputIterator m_output_ordinals;                                       \
  };                                                                               \
  }} /*namespace stk::detail*/                                                     \
  namespace stk {                                                                  \
  template <typename OrdinalOutputIterator>                                        \
  STKTOPOLOGY_INLINE_FUNCTION                                                      \
  void topology::name( int ordinal, OrdinalOutputIterator output_ordinals) const   \
  { typedef detail::name##_impl<OrdinalOutputIterator> functor;                    \
    functor f(ordinal,output_ordinals);                                            \
    topology::apply_functor< functor > apply( f );                                 \
    apply(m_value);                                                                \
  }                                                                                \
  } /*namespace stk*/

#define STKTOPOLOGY_NODES_MEMBER(name)                                             \
  namespace stk { namespace detail {                                               \
  template <typename NodeArray, typename NodeOutputIterator>                       \
  struct name##_impl {                                                             \
    typedef void result_type;                                                      \
    STKTOPOLOGY_INLINE_FUNCTION                                                    \
    name##_impl(   const NodeArray &nodes                                          \
                 , int ordinal                                                     \
                 , NodeOutputIterator output_ordinals)                             \
      : m_nodes(nodes)                                                             \
      , m_ordinal(ordinal)                                                         \
      , m_output_ordinals(output_ordinals)                                         \
    {}                                                                             \
    template <typename Topology>                                                   \
    STKTOPOLOGY_INLINE_FUNCTION                                                    \
    result_type operator()(Topology) const                                         \
    { return Topology::name(m_nodes,m_ordinal,m_output_ordinals); }                \
    const NodeArray    & m_nodes;                                                  \
    int                  m_ordinal;                                                \
    NodeOutputIterator   m_output_ordinals;                                        \
  };                                                                               \
  }} /*namespace stk::detail*/                                                     \
  namespace stk {                                                                  \
  template <typename NodeArray, typename NodeOutputIterator>                       \
  STKTOPOLOGY_INLINE_FUNCTION                                                      \
  void topology::name(   const NodeArray & nodes                                   \
                       , int ordinal                                               \
                       , NodeOutputIterator output_ordinals) const                 \
  { typedef detail::name##_impl<NodeArray,NodeOutputIterator> functor;             \
    functor f(nodes,ordinal,output_ordinals);                                      \
    topology::apply_functor< functor > apply( f );                                 \
    apply(m_value);                                                                \
  }                                                                                \
  } /*namespace stk*/

STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_edges,bool)
STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_faces,bool)
STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_sides,bool)
STKTOPOLOGY_SIMPLE_MEMBER(is_shell,bool)
STKTOPOLOGY_SIMPLE_MEMBER(rank,topology::rank_t)
STKTOPOLOGY_SIMPLE_MEMBER(side_rank,topology::rank_t)
STKTOPOLOGY_SIMPLE_MEMBER(dimension,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_nodes,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_vertices,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_edges,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_faces,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_sides,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_permutations,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_positive_permutations,int)
STKTOPOLOGY_SIMPLE_MEMBER(base,topology)

STKTOPOLOGY_ORDINAL_MEMBER(defined_on_spatial_dimension,bool)
STKTOPOLOGY_ORDINAL_MEMBER(edge_topology,topology)
STKTOPOLOGY_ORDINAL_MEMBER(face_topology,topology)
STKTOPOLOGY_ORDINAL_MEMBER(side_topology,topology)

STKTOPOLOGY_ORDINAL_NODES_MEMBER(edge_node_ordinals)
STKTOPOLOGY_ORDINAL_NODES_MEMBER(face_node_ordinals)
STKTOPOLOGY_ORDINAL_NODES_MEMBER(side_node_ordinals)
STKTOPOLOGY_ORDINAL_NODES_MEMBER(permutation_node_ordinals)

STKTOPOLOGY_NODES_MEMBER(edge_nodes)
STKTOPOLOGY_NODES_MEMBER(face_nodes)
STKTOPOLOGY_NODES_MEMBER(side_nodes)
STKTOPOLOGY_NODES_MEMBER(permutation_nodes)

#undef STKTOPOLOGY_SIMPLE_MEMBER
#undef STKTOPOLOGY_ORDINAL_MEMBER
#undef STKTOPOLOGY_ORDINAL_NODES_MEMBER
#undef STKTOPOLOGY_NODES_MEMBER

#endif //STKTOPOLOGY_TOPOLOGY_TCC


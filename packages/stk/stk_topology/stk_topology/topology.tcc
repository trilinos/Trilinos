#ifndef STKTOPOLOGY_TOPOLOGY_TCC
#define STKTOPOLOGY_TOPOLOGY_TCC

#define STKTOPOLOGY_ORDINAL_NODES_MEMBER(name)                                     \
  namespace stk { namespace topology_detail {                                      \
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
  }} /*namespace stk::topology_detail*/                                            \
  namespace stk {                                                                  \
  template <typename OrdinalOutputIterator>                                        \
  STKTOPOLOGY_INLINE_FUNCTION                                                      \
  void topology::name( int ordinal, OrdinalOutputIterator output_ordinals) const   \
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
  }} /*namespace stk::topology_detail*/                                            \
  namespace stk {                                                                  \
  template <typename NodeArray, typename NodeOutputIterator>                       \
  STKTOPOLOGY_INLINE_FUNCTION                                                      \
  void topology::name(   const NodeArray & nodes                                   \
                       , int ordinal                                               \
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
    STKTOPOLOGY_INLINE_FUNCTION                        \
    result_type operator()(Topology) const             \
    { return Topology::name; }                         \
  };                                                   \
  }} /*namespace stk::topology_detail*/                \
  namespace stk {                                      \
  STKTOPOLOGY_INLINE_FUNCTION                          \
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
    STKTOPOLOGY_INLINE_FUNCTION                        \
    name##_impl(int ordinal)                           \
      : m_ordinal(ordinal)                             \
    {}                                                 \
    template <typename Topology>                       \
    STKTOPOLOGY_INLINE_FUNCTION                        \
    result_type operator()(Topology) const             \
    { return Topology::name(m_ordinal); }              \
    int m_ordinal;                                     \
  };                                                   \
  }} /*namespace stk::topology_detail*/                \
  namespace stk {                                      \
  STKTOPOLOGY_INLINE_FUNCTION                          \
  result topology::name(int ordinal) const             \
  { typedef topology_detail::name##_impl functor;      \
    functor f(ordinal);                                \
    topology::apply_functor< functor > apply( f );     \
    return apply(m_value);                             \
  }                                                    \
  } /*namespace stk*/

STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_faces,bool)
STKTOPOLOGY_SIMPLE_MEMBER(is_shell,bool)
STKTOPOLOGY_SIMPLE_MEMBER(side_rank,stk::topology::rank_t)
STKTOPOLOGY_SIMPLE_MEMBER(dimension,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_vertices,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_edges,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_faces,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_permutations,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_positive_permutations,int)
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
  typedef int result_type;
  template <typename Topology>
  STKTOPOLOGY_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::num_nodes; }
};

struct rank_impl {
  typedef topology::rank_t result_type;
  template <typename Topology>
  STKTOPOLOGY_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::rank; }
};

template <typename NodeArrayA, typename NodeArrayB>
struct equivalent_impl {
  typedef std::pair<bool,int> result_type;

  STKTOPOLOGY_INLINE_FUNCTION
  equivalent_impl( const NodeArrayA &a , const NodeArrayB &b )
    : m_a(a), m_b(b)
  {}

  template <typename Topology>
  STKTOPOLOGY_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::equivalent(m_a, m_b); }

  const NodeArrayA & m_a;
  const NodeArrayB & m_b;
};

template <typename NodeArray>
struct lexicographical_smallest_permutation_impl {
  typedef int result_type;

  STKTOPOLOGY_INLINE_FUNCTION
  lexicographical_smallest_permutation_impl( const NodeArray &nodes , bool only_positive_permutations )
    : m_nodes(nodes), m_only_positive_permutations(only_positive_permutations)
  {}

  template <typename Topology>
  STKTOPOLOGY_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::lexicographical_smallest_permutation(m_nodes, m_only_positive_permutations); }

  const NodeArray & m_nodes;
  bool              m_only_positive_permutations;
};

}} /*namespace stk::topology_detail*/

namespace stk {

STKTOPOLOGY_INLINE_FUNCTION
int topology::num_nodes() const
{
  typedef topology_detail::num_nodes_impl functor;
  topology::apply_functor< functor > apply;
  return m_value <= END_TOPOLOGY ? apply(m_value) : m_value - SUPERELEMENT_START;
}

STKTOPOLOGY_INLINE_FUNCTION
topology::rank_t topology::rank() const
{
  typedef topology_detail::rank_impl functor;
  topology::apply_functor< functor > apply;
  return m_value <= END_TOPOLOGY ? apply(m_value) : topology::ELEMENT_RANK;
}

template <typename NodeArrayA, typename NodeArrayB>
STKTOPOLOGY_INLINE_FUNCTION
std::pair<bool,int> topology::equivalent( const NodeArrayA &a, const NodeArrayB &b) const
{ typedef topology_detail::equivalent_impl<NodeArrayA,NodeArrayB> functor;
  functor f(a,b);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
STKTOPOLOGY_INLINE_FUNCTION
int topology::lexicographical_smallest_permutation( const NodeArray &nodes, bool only_positive_permutations) const
{
  typedef topology_detail::lexicographical_smallest_permutation_impl< NodeArray > functor;
  functor f(nodes, only_positive_permutations);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

} /*namespace stk*/

#endif //STKTOPOLOGY_TOPOLOGY_TCC


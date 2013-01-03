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


STKTOPOLOGY_ORDINAL_NODES_MEMBER(edge_node_ordinals)
STKTOPOLOGY_ORDINAL_NODES_MEMBER(face_node_ordinals)
STKTOPOLOGY_ORDINAL_NODES_MEMBER(permutation_node_ordinals)

STKTOPOLOGY_NODES_MEMBER(edge_nodes)
STKTOPOLOGY_NODES_MEMBER(face_nodes)
STKTOPOLOGY_NODES_MEMBER(permutation_nodes)

#undef STKTOPOLOGY_ORDINAL_NODES_MEMBER
#undef STKTOPOLOGY_NODES_MEMBER


namespace stk { namespace topology_detail {

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


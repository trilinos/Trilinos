#include <stk_topology/topology.hpp>


#define STKTOPOLOGY_SIMPLE_MEMBER(name,result)         \
  namespace {                                          \
  struct name##_impl {                                 \
    typedef result result_type;                        \
    template <typename Topology>                       \
    result_type operator()(Topology) const             \
    { return Topology::name; }                         \
  };                                                   \
  }  /*namespace */                                    \
  namespace stk {                                      \
  result topology::name() const                        \
  { topology::apply_functor< name##_impl > apply;      \
    return apply(m_value);                             \
  }                                                    \
  } /*namespace stk*/

#define STKTOPOLOGY_ORDINAL_MEMBER(name,result)        \
  namespace {                                          \
  struct name##_impl {                                 \
    typedef result result_type;                        \
    name##_impl(int ordinal)                           \
      : m_ordinal(ordinal)                             \
    {}                                                 \
    template <typename Topology>                       \
    result_type operator()(Topology) const             \
    { return Topology::name(m_ordinal); }              \
    int m_ordinal;                                     \
  };                                                   \
  }  /*namespace */                                    \
  namespace stk {                                      \
  result topology::name(int ordinal) const             \
  { name##_impl f(ordinal);                            \
    topology::apply_functor< name##_impl > apply( f ); \
    return apply(m_value);                             \
  }                                                    \
  } /*namespace stk*/

STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_edges,bool)
STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_faces,bool)
STKTOPOLOGY_SIMPLE_MEMBER(has_homogeneous_sides,bool)
STKTOPOLOGY_SIMPLE_MEMBER(is_shell,bool)
STKTOPOLOGY_SIMPLE_MEMBER(rank,stk::topology::rank_t)
STKTOPOLOGY_SIMPLE_MEMBER(side_rank,stk::topology::rank_t)
STKTOPOLOGY_SIMPLE_MEMBER(dimension,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_nodes,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_vertices,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_edges,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_faces,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_sides,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_permutations,int)
STKTOPOLOGY_SIMPLE_MEMBER(num_positive_permutations,int)
STKTOPOLOGY_SIMPLE_MEMBER(base,stk::topology)

STKTOPOLOGY_ORDINAL_MEMBER(defined_on_spatial_dimension,bool)
STKTOPOLOGY_ORDINAL_MEMBER(edge_topology,stk::topology)
STKTOPOLOGY_ORDINAL_MEMBER(face_topology,stk::topology)
STKTOPOLOGY_ORDINAL_MEMBER(side_topology,stk::topology)

#undef STKTOPOLOGY_SIMPLE_MEMBER
#undef STKTOPOLOGY_ORDINAL_MEMBER

namespace stk {

const char * topology::rank_names[] =
{
    "NODE_RANK"
  , "EDGE_RANK"
  , "FACE_RANK"
  , "ELEMENT_RANK"
  , "CONSTRAINT_RANK"
  , "INVALID_RANK"
};

const char * topology::topology_names[] =
{
    "NODE"
  , "LINE_2"
  , "LINE_3"
  , "TRIANGLE_3"
  , "TRIANGLE_4"
  , "TRIANGLE_6"
  , "QUADRILATERAL_4"
  , "QUADRILATERAL_8"
  , "QUADRILATERAL_9"
  , "PARTICLE"
  , "LINE_2_1D"
  , "LINE_3_1D"
  , "BEAM_2"
  , "BEAM_3"
  , "SHELL_LINE_2"
  , "SHELL_LINE_3"
  , "TRIANGLE_3_2D"
  , "TRIANGLE_4_2D"
  , "TRIANGLE_6_2D"
  , "QUADRILATERAL_4_2D"
  , "QUADRILATERAL_8_2D"
  , "QUADRILATERAL_9_2D"
  , "SHELL_TRIANGLE_3"
  , "SHELL_TRIANGLE_4"
  , "SHELL_TRIANGLE_6"
  , "SHELL_QUADRILATERAL_4"
  , "SHELL_QUADRILATERAL_8"
  , "SHELL_QUADRILATERAL_9"
  , "TETRAHEDRON_4"
  , "TETRAHEDRON_8"
  , "TETRAHEDRON_10"
  , "TETRAHEDRON_11"
  , "PYRAMID_5"
  , "PYRAMID_13"
  , "PYRAMID_14"
  , "WEDGE_6"
  , "WEDGE_15"
  , "WEDGE_18"
  , "HEXAHEDRON_8"
  , "HEXAHEDRON_20"
  , "HEXAHEDRON_27"
  , "INVALID_TOPOLOGY"
};

} //namespace stk



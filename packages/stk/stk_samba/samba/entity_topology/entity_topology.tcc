#ifndef SAMBA_SAMBA_ENTITY_TOPOLOGY_ENTITY_TOPOLOGY_TCC
#define SAMBA_SAMBA_ENTITY_TOPOLOGY_ENTITY_TOPOLOGY_TCC

#include <samba/detail/apply_uint_8_functor.hpp>

#include <samba/connectivity_ordinal.hpp>

#include <samba/entity_topology/beam.hpp>
#include <samba/entity_topology/hexahedron.hpp>
#include <samba/entity_topology/hexagon.hpp>
#include <samba/entity_topology/line.hpp>
#include <samba/entity_topology/node.hpp>
#include <samba/entity_topology/particle.hpp>
#include <samba/entity_topology/pentagon.hpp>
#include <samba/entity_topology/pyramid.hpp>
#include <samba/entity_topology/quadrilateral.hpp>
#include <samba/entity_topology/shell_line.hpp>
#include <samba/entity_topology/shell_quadrilateral.hpp>
#include <samba/entity_topology/shell_triangle.hpp>
#include <samba/entity_topology/tetrahedron.hpp>
#include <samba/entity_topology/triangle.hpp>
#include <samba/entity_topology/wedge.hpp>

namespace samba {

namespace detail {
struct entity_topology_output_impl
{
  std::ostream& m_out;

  entity_topology_output_impl(std::ostream& out)
    : m_out(out)
  {}

  typedef void result_type;

  template <entity_topology::value_type Topology>
  void operator()(entity_topology::topology_type<Topology> t) const
  { m_out << t; }
};
} //namespace detail

template <entity_topology::value_type Topology>
inline std::ostream& operator<<(std::ostream& out, entity_topology::topology_type<Topology>)
{ return out << static_cast<int>(Topology); }

inline std::ostream& operator<<(std::ostream& out, entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::entity_topology_output_impl
                       > apply;

  detail::entity_topology_output_impl f(out);
  out << "{" << entity_topology::tag() << ":";

  apply(f, t());

  out << "}";
  return out;
}

//dimension
namespace detail {
struct dimension_impl
{
  typedef unsigned result_type;
  template <typename T>
  unsigned operator()(T t) const
  { return entity_topology::dimension<T>(); }
};
} //namespace detail

inline unsigned dimension(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::dimension_impl
                       > apply;
  return apply(detail::dimension_impl(),t());
}

//side_dimension
namespace detail {
struct side_dimension_impl
{
  typedef unsigned result_type;
  template <typename T>
  unsigned operator()(T t) const
  { return entity_topology::side_dimension<T>(); }
};
} //namespace detail

inline unsigned side_dimension(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::side_dimension_impl
                       > apply;
  return apply(detail::side_dimension_impl(),t());
}

//num_nodes
namespace detail {
struct num_nodes_impl
{
  typedef unsigned result_type;
  template <typename T>
  unsigned operator()(T t) const
  { return entity_topology::num_nodes<T>(); }
};
} //namespace detail

inline unsigned num_nodes(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::num_nodes_impl
                       > apply;
  return apply(detail::num_nodes_impl(),t());
}

//num_vertices
namespace detail {
struct num_vertices_impl
{
  typedef unsigned result_type;
  template <typename T>
  unsigned operator()(T t) const
  { return entity_topology::num_vertices<T>(); }
};
} //namespace detail

inline unsigned num_vertices(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::num_vertices_impl
                       > apply;
  return apply(detail::num_vertices_impl(),t());
}

//num_edges
namespace detail {
struct num_edges_impl
{
  typedef unsigned result_type;
  template <typename T>
  unsigned operator()(T t) const
  { return entity_topology::num_edges<T>(); }
};
} //namespace detail

inline unsigned num_edges(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::num_edges_impl
                       > apply;
  return apply(detail::num_edges_impl(),t());
}

//num_faces
namespace detail {
struct num_faces_impl
{
  typedef unsigned result_type;
  template <typename T>
  unsigned operator()(T t) const
  { return entity_topology::num_faces<T>(); }
};
} //namespace detail

inline unsigned num_faces(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::num_faces_impl
                       > apply;
  return apply(detail::num_faces_impl(),t());
}

//num_sides
namespace detail {
struct num_sides_impl
{
  typedef unsigned result_type;
  template <typename T>
  unsigned operator()(T t) const
  { return entity_topology::num_sides<T>(); }
};
} //namespace detail

inline unsigned num_sides(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::num_sides_impl
                       > apply;
  return apply(detail::num_sides_impl(),t());
}

//has_homogeneous_sides
namespace detail {
struct has_homogeneous_sides_impl
{
  typedef bool result_type;
  template <typename T>
  bool operator()(T t) const
  { return entity_topology::has_homogeneous_sides<T>(); }
};
} //namespace detail

inline bool has_homogeneous_sides(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::has_homogeneous_sides_impl
                       > apply;
  return apply(detail::has_homogeneous_sides_impl(),t());
}

inline void print_topology_detail(std::ostream& out, entity_topology t)
{
  out << t << ":" << std::endl;
  out << "           dimension: " << dimension(t)    << std::endl;
  out << "        num_vertices: " << num_vertices(t) << std::endl;
  out << "           num_nodes: " << num_nodes(t)    << std::endl;
  out << "           num_edges: " << num_edges(t)    << std::endl;
  out << "           num_faces: " << num_faces(t)    << std::endl;
  out << "           num_sides: " << num_sides(t)    << std::endl;
  out << "   homogeneous_sides: " << has_homogeneous_sides(t) << std::endl;
  out << std::endl;

}

inline entity_rank topology_rank(entity_topology t, spatial_dimension sp)
{
  entity_rank rank = entity_rank::invalid();

  if( sp()==2) {
    if (dimension(t)==2) {
      rank = entity_rank::element();
    }
    else {
      rank = entity_rank::create(dimension(t));
    }
  }
  else if(sp() == 3) {
    rank = entity_rank::create(dimension(t));
  }

  return rank;
};

inline entity_rank side_rank(entity_topology t, spatial_dimension sp)
{ return entity_rank::create(side_dimension(t)); }

//
// is_valid
//

namespace detail {
struct entity_topology_is_valid_topology_impl
{
  typedef bool result_type;

  template <typename T>
  result_type operator()(T /*t*/) const
  { return entity_topology::is_valid_topology<T>(); }
};
} //namespace detail

inline bool is_valid_topology(entity_topology t)
{
  detail::apply_functor< entity_topology::topology_type
                        ,detail::entity_topology_is_valid_topology_impl
                       > apply;
  return apply(detail::entity_topology_is_valid_topology_impl(), t());
}

//
// getting nodes for sides/faces/edges
//

template <typename UndefinedTopology>
inline
connectivity_ordinal const* const face_nodes(UndefinedTopology t, unsigned face_id)
{
  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

template <typename UndefinedTopology>
inline
connectivity_ordinal const* const edge_nodes(UndefinedTopology t, unsigned edge_id)
{
  static const connectivity_ordinal* nodes = NULL;
  return nodes;
}

namespace detail {

struct entity_topology_face_nodes_impl
{
  unsigned m_face_id;

  entity_topology_face_nodes_impl(unsigned face_id) : m_face_id(face_id) {}

  typedef connectivity_ordinal const* const result_type;

  template <typename T>
  result_type operator()(T t) const
  { return face_nodes(t, m_face_id); }
};

struct entity_topology_edge_nodes_impl
{
  unsigned m_edge_id;

  entity_topology_edge_nodes_impl(unsigned edge_id) : m_edge_id(edge_id) {}

  typedef connectivity_ordinal const* const result_type;

  template <typename T>
  result_type operator()(T t) const
  { return edge_nodes(t, m_edge_id); }
};

} //namespace detail

inline
connectivity_ordinal const* const side_nodes(entity_topology t, unsigned side_id)
{
  entity_rank rank_of_side = entity_rank::create(dimension(t) - 1);
  if (rank_of_side == entity_rank::face()) {
    return face_nodes(t, side_id);
  }
  else if (rank_of_side == entity_rank::edge()) {
    return edge_nodes(t, side_id);
  }
  else {
    static const connectivity_ordinal* nodes = NULL;
    return nodes;
  }
}

inline
connectivity_ordinal const* const face_nodes(entity_topology t, unsigned face_id)
{
  detail::apply_functor< entity_topology::topology_type
                         ,detail::entity_topology_face_nodes_impl
                         > apply;

  detail::entity_topology_face_nodes_impl f(face_id);

  return apply(f, t());
}

inline
connectivity_ordinal const* const edge_nodes(entity_topology t, unsigned edge_id)
{
  detail::apply_functor< entity_topology::topology_type
                         ,detail::entity_topology_edge_nodes_impl
                         > apply;

  detail::entity_topology_edge_nodes_impl f(edge_id);

  return apply(f, t());
}

//
// getting topologies of sides
//

namespace detail {

struct entity_topology_edge_topology_impl
{
  unsigned m_edge_id;

  entity_topology_edge_topology_impl(unsigned edge_id) : m_edge_id(edge_id) {}

  typedef entity_topology result_type;

  template <typename T>
  result_type operator()(T t) const
  {
    if (m_edge_id > static_cast<unsigned>(entity_topology::num_edges<T>())) {
      return entity_topology::invalid();
    }

    switch(m_edge_id) {

    case 0:  return entity_topology::create(entity_topology::edge_topology<T, 0 >::type::value);

    case 1:  return entity_topology::create(entity_topology::edge_topology<T, 1 >::type::value);
    case 2:  return entity_topology::create(entity_topology::edge_topology<T, 2 >::type::value);
    case 3:  return entity_topology::create(entity_topology::edge_topology<T, 3 >::type::value);
    case 4:  return entity_topology::create(entity_topology::edge_topology<T, 4 >::type::value);
    case 5:  return entity_topology::create(entity_topology::edge_topology<T, 5 >::type::value);
    case 6:  return entity_topology::create(entity_topology::edge_topology<T, 6 >::type::value);
    case 7:  return entity_topology::create(entity_topology::edge_topology<T, 7 >::type::value);
    case 8:  return entity_topology::create(entity_topology::edge_topology<T, 8 >::type::value);
    case 9:  return entity_topology::create(entity_topology::edge_topology<T, 9 >::type::value);
    case 10: return entity_topology::create(entity_topology::edge_topology<T, 10>::type::value);
    case 11: return entity_topology::create(entity_topology::edge_topology<T, 11>::type::value);
    case 12: return entity_topology::create(entity_topology::edge_topology<T, 12>::type::value);
    case 13: return entity_topology::create(entity_topology::edge_topology<T, 13>::type::value);
    case 14: return entity_topology::create(entity_topology::edge_topology<T, 14>::type::value);
    case 15: return entity_topology::create(entity_topology::edge_topology<T, 15>::type::value);
    case 16: return entity_topology::create(entity_topology::edge_topology<T, 16>::type::value);
    case 17: return entity_topology::create(entity_topology::edge_topology<T, 17>::type::value);
    case 18: return entity_topology::create(entity_topology::edge_topology<T, 18>::type::value);
    case 19: return entity_topology::create(entity_topology::edge_topology<T, 19>::type::value);
    case 20: return entity_topology::create(entity_topology::edge_topology<T, 20>::type::value);
    default:
      BOOST_ASSERT_MSG(false, "Need more case statements");
    }
    return entity_topology::invalid();
  }
};

struct entity_topology_face_topology_impl
{
  unsigned m_face_id;

  entity_topology_face_topology_impl(unsigned face_id) : m_face_id(face_id) {}

  typedef entity_topology result_type;

  template <typename T>
  result_type operator()(T t) const
  {
    if (m_face_id > static_cast<unsigned>(entity_topology::num_faces<T>())) {
      return entity_topology::invalid();
    }

    switch(m_face_id) {

    case 0:  return entity_topology::create(entity_topology::face_topology<T, 0 >::type::value);

    case 1:  return entity_topology::create(entity_topology::face_topology<T, 1 >::type::value);
    case 2:  return entity_topology::create(entity_topology::face_topology<T, 2 >::type::value);
    case 3:  return entity_topology::create(entity_topology::face_topology<T, 3 >::type::value);
    case 4:  return entity_topology::create(entity_topology::face_topology<T, 4 >::type::value);
    case 5:  return entity_topology::create(entity_topology::face_topology<T, 5 >::type::value);
    case 6:  return entity_topology::create(entity_topology::face_topology<T, 6 >::type::value);
    case 7:  return entity_topology::create(entity_topology::face_topology<T, 7 >::type::value);
    case 8:  return entity_topology::create(entity_topology::face_topology<T, 8 >::type::value);
    case 9:  return entity_topology::create(entity_topology::face_topology<T, 9 >::type::value);
    case 10: return entity_topology::create(entity_topology::face_topology<T, 10>::type::value);
    case 11: return entity_topology::create(entity_topology::face_topology<T, 11>::type::value);
    case 12: return entity_topology::create(entity_topology::face_topology<T, 12>::type::value);
    case 13: return entity_topology::create(entity_topology::face_topology<T, 13>::type::value);
    case 14: return entity_topology::create(entity_topology::face_topology<T, 14>::type::value);
    case 15: return entity_topology::create(entity_topology::face_topology<T, 15>::type::value);
    case 16: return entity_topology::create(entity_topology::face_topology<T, 16>::type::value);
    case 17: return entity_topology::create(entity_topology::face_topology<T, 17>::type::value);
    case 18: return entity_topology::create(entity_topology::face_topology<T, 18>::type::value);
    case 19: return entity_topology::create(entity_topology::face_topology<T, 19>::type::value);
    case 20: return entity_topology::create(entity_topology::face_topology<T, 20>::type::value);
    default:
      BOOST_ASSERT_MSG(false, "Need more case statements");
    }
    return entity_topology::invalid();
  }
};

} // namespace detail

inline
entity_topology side_topology(entity_topology t, unsigned side_id)
{
  entity_rank rank_of_side = entity_rank::create(dimension(t) - 1);
  if (rank_of_side == entity_rank::face()) {
    return face_topology(t, side_id);
  }
  else if (rank_of_side == entity_rank::edge()) {
    return edge_topology(t, side_id);
  }
  else {
    return entity_topology::invalid();
  }
}

inline
entity_topology face_topology(entity_topology t, unsigned face_id)
{
  detail::apply_functor< entity_topology::topology_type
                         ,detail::entity_topology_face_topology_impl
                         > apply;

  detail::entity_topology_face_topology_impl f(face_id);

  return apply(f, t());
}

inline
entity_topology edge_topology(entity_topology t, unsigned edge_id)
{
  detail::apply_functor< entity_topology::topology_type
                         ,detail::entity_topology_edge_topology_impl
                         > apply;

  detail::entity_topology_edge_topology_impl f(edge_id);

  return apply(f, t());
}

inline
entity_topology connectivity_topology(entity_topology from, entity_rank to, connectivity_ordinal to_ordinal)
{
  // TODO
  return entity_topology::invalid();
}

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_TOPOLOGY_ENTITY_TOPOLOGY_TCC

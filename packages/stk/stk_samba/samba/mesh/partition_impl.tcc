#ifndef SAMBA_SAMBA_MESH_PARTITION_IMPL_TCC
#define SAMBA_SAMBA_MESH_PARTITION_IMPL_TCC

#include <samba/types.hpp>

#include <boost/iterator/iterator_traits.hpp>

#include <boost/iterator/zip_iterator.hpp>



#define SAMBA_PARTITION_IMPL_DEFAULT_RANGE(rank_name)                                                               \
inline partition_index_range partition_impl::rank_name##s(partition_offset offset) const                            \
  {                                                                                                                 \
    BOOST_ASSERT(!m_in_modification_cycle);                                                                         \
                                                                                                                    \
    BOOST_ASSERT_MSG( m_##rank_name##_kind != connectivity_kind::invalid(),                                         \
        (debug_message() << "Invalid connectivity from " << rank() << " to " << entity_rank::rank_name()));         \
                                                                                                                    \
    return connectivity_kind::fixed() == m_##rank_name##_kind                                                       \
      ? m_fixed_##rank_name##_connectivity.range<partition_index>(offset)                                           \
      : m_dynamic_##rank_name##_connectivity.range<partition_index>(offset);                                        \
  }

#define SAMBA_PARTITION_IMPL_DEFAULT_BEGIN(rank_name)                                                               \
inline partition_index_iterator partition_impl::begin_##rank_name##s(partition_offset offset) const                 \
  {                                                                                                                 \
    BOOST_ASSERT(!m_in_modification_cycle);                                                                         \
                                                                                                                    \
    BOOST_ASSERT_MSG( m_##rank_name##_kind != connectivity_kind::invalid(),                                         \
        (debug_message() << "Invalid connectivity from " << rank() << " to " << entity_rank::rank_name()));         \
                                                                                                                    \
    return connectivity_kind::fixed() == m_##rank_name##_kind                                                       \
      ? m_fixed_##rank_name##_connectivity.begin<partition_index>(offset)                                           \
      : m_dynamic_##rank_name##_connectivity.begin<partition_index>(offset);                                        \
  }

#define SAMBA_PARTITION_IMPL_DEFAULT_END(rank_name)                                                               \
inline partition_index_iterator partition_impl::end_##rank_name##s(partition_offset offset) const                 \
  {                                                                                                               \
    BOOST_ASSERT(!m_in_modification_cycle);                                                                       \
                                                                                                                  \
    BOOST_ASSERT_MSG( m_##rank_name##_kind != connectivity_kind::invalid(),                                       \
        (debug_message() << "Invalid connectivity from " << rank() << " to " << entity_rank::rank_name()));       \
                                                                                                                  \
    return connectivity_kind::fixed() == m_##rank_name##_kind                                                     \
      ? m_fixed_##rank_name##_connectivity.end<partition_index>(offset)                                           \
      : m_dynamic_##rank_name##_connectivity.end<partition_index>(offset);                                        \
  }

#define SAMBA_PARTITION_IMPL_NUM_CONNECTIVITY(rank_name)                                                          \
inline size_t partition_impl::num_##rank_name##s(partition_offset offset) const                                   \
  {                                                                                                               \
    BOOST_ASSERT_MSG( m_##rank_name##_kind != connectivity_kind::invalid(),                                       \
        (debug_message() << "Invalid connectivity from " << rank() << " to " << entity_rank::rank_name()));       \
                                                                                                                  \
    return connectivity_kind::fixed() == m_##rank_name##_kind                                                     \
      ? m_fixed_##rank_name##_connectivity.num_connectivity(offset)                                               \
      : m_dynamic_##rank_name##_connectivity.num_connectivity(offset);                                            \
  }

#define SAMBA_PARTITION_IMPL_RANGE(rank_name)                                                                       \
template <typename T>                                                                                               \
inline std::pair<const T*, const T*> partition_impl::rank_name##s(partition_offset offset) const                    \
  {                                                                                                                 \
    BOOST_ASSERT_MSG( m_##rank_name##_kind != connectivity_kind::invalid(),                                         \
        (debug_message() << "Invalid connectivity from " << rank() << " to " << entity_rank::rank_name()));         \
                                                                                                                    \
    return connectivity_kind::fixed() == m_##rank_name##_kind                                                       \
      ? m_fixed_##rank_name##_connectivity.range<T>(offset)                                                         \
      : m_dynamic_##rank_name##_connectivity.range<T>(offset);                                                      \
  }

#define SAMBA_PARTITION_IMPL_BEGIN(rank_name)                                                                       \
template <typename T>                                                                                               \
inline const T * partition_impl::begin_##rank_name##s(partition_offset offset) const                                \
  {                                                                                                                 \
    BOOST_ASSERT_MSG( m_##rank_name##_kind != connectivity_kind::invalid(),                                         \
        (debug_message() << "Invalid connectivity from " << rank() << " to " << entity_rank::rank_name()));         \
                                                                                                                    \
    return connectivity_kind::fixed() == m_##rank_name##_kind                                                       \
      ? m_fixed_##rank_name##_connectivity.begin<T>(offset)                                                         \
      : m_dynamic_##rank_name##_connectivity.begin<T>(offset);                                                      \
  }

#define SAMBA_PARTITION_IMPL_END(rank_name)                                                                       \
template <typename T>                                                                                             \
inline const T * partition_impl::end_##rank_name##s(partition_offset offset) const                                \
  {                                                                                                               \
    BOOST_ASSERT_MSG( m_##rank_name##_kind != connectivity_kind::invalid(),                                       \
        (debug_message() << "Invalid connectivity from " << rank() << " to " << entity_rank::rank_name()));       \
                                                                                                                  \
    return connectivity_kind::fixed() == m_##rank_name##_kind                                                     \
      ? m_fixed_##rank_name##_connectivity.end<T>(offset)                                                         \
      : m_dynamic_##rank_name##_connectivity.end<T>(offset);                                                      \
  }



namespace samba { namespace detail {

//defined here to break cyclic dependency

inline
partition_impl::partition_impl( mesh_impl * arg_mesh
                                ,partition_id arg_partition
                                ,entity_part_vector const& arg_parts
                                )
  : m_mesh(arg_mesh)
  , m_partition(arg_partition)
  , m_parts(arg_parts)
  , m_num_entities(0)
  , m_starting_rank_count(0)
  , m_in_modification_cycle(true)
  , m_node_kind(m_mesh->connectivity_map()(rank(),entity_rank::node()))
  , m_edge_kind(m_mesh->connectivity_map()(rank(),entity_rank::edge()))
  , m_face_kind(m_mesh->connectivity_map()(rank(),entity_rank::face()))
  , m_element_kind(m_mesh->connectivity_map()(rank(),entity_rank::element()))
  , m_fixed_node_connectivity(topology(),m_mesh->spatial_dimension())
  , m_fixed_edge_connectivity(topology(),m_mesh->spatial_dimension())
  , m_fixed_face_connectivity(topology(),m_mesh->spatial_dimension())
  , m_fixed_element_connectivity(topology(),m_mesh->spatial_dimension())
  , m_dynamic_node_connectivity(topology(),m_mesh->spatial_dimension())
  , m_dynamic_edge_connectivity(topology(),m_mesh->spatial_dimension())
  , m_dynamic_face_connectivity(topology(),m_mesh->spatial_dimension())
  , m_dynamic_element_connectivity(topology(),m_mesh->spatial_dimension())
{
  BOOST_ASSERT_MSG( m_parts.size() >= 2, "rank and topology must be specified" );
  BOOST_ASSERT_MSG( m_parts[0].which() == entity_part::Rank, "first part must be the rank");
  BOOST_ASSERT_MSG( m_parts[1].which() == entity_part::Topology, "second part must be the topology");
}

//*************************************************************************
//connectivity accessors
//*************************************************************************

// descriptor

inline
size_t partition_impl::num_connectivity(entity_rank to_rank, partition_offset offset) const
{
  switch(to_rank()) {
  case entity_rank::node_type::value:    return num_nodes(offset);
  case entity_rank::edge_type::value:    return num_edges(offset);
  case entity_rank::face_type::value:    return num_faces(offset);
  case entity_rank::element_type::value: return num_elements(offset);
  default:
    BOOST_ASSERT_MSG(false, "Invalid rank");
    return 0;
  }
}

template <typename T>
inline
std::pair<const T*, const T*> partition_impl::connectivity(entity_rank rank, partition_offset offset) const
{
  switch(rank()) {
  case entity_rank::node_type::value:    return nodes<T>(offset);
  case entity_rank::edge_type::value:    return edges<T>(offset);
  case entity_rank::face_type::value:    return faces<T>(offset);
  case entity_rank::element_type::value: return elements<T>(offset);
  default:
    BOOST_ASSERT_MSG(false, "Invalid rank");
    return std::make_pair<const T*,const T*>(NULL,NULL);
  }
  return std::make_pair<const T*,const T*>(NULL,NULL);
}

template <typename T>
inline
const T* partition_impl::begin_connectivity(entity_rank to_rank, partition_offset offset) const
{
  switch(to_rank()) {
  case entity_rank::node_type::value:    return begin_nodes<T>(offset);
  case entity_rank::edge_type::value:    return begin_edges<T>(offset);
  case entity_rank::face_type::value:    return begin_faces<T>(offset);
  case entity_rank::element_type::value: return begin_elements<T>(offset);
  default:
    BOOST_ASSERT_MSG(false, "Invalid rank");
    return NULL;
  }
}

template <typename T>
inline
const T* partition_impl::end_connectivity(entity_rank to_rank, partition_offset offset) const
{
  switch(to_rank()) {
  case entity_rank::node_type::value:    return end_nodes<T>(offset);
  case entity_rank::edge_type::value:    return end_edges<T>(offset);
  case entity_rank::face_type::value:    return end_faces<T>(offset);
  case entity_rank::element_type::value: return end_elements<T>(offset);
  default:
    BOOST_ASSERT_MSG(false, "Invalid rank");
    return NULL;
  }
}

// partition_index range
SAMBA_PARTITION_IMPL_DEFAULT_RANGE(node)
SAMBA_PARTITION_IMPL_DEFAULT_RANGE(edge)
SAMBA_PARTITION_IMPL_DEFAULT_RANGE(face)
SAMBA_PARTITION_IMPL_DEFAULT_RANGE(element)

SAMBA_PARTITION_IMPL_DEFAULT_BEGIN(node)
SAMBA_PARTITION_IMPL_DEFAULT_BEGIN(edge)
SAMBA_PARTITION_IMPL_DEFAULT_BEGIN(face)
SAMBA_PARTITION_IMPL_DEFAULT_BEGIN(element)

SAMBA_PARTITION_IMPL_DEFAULT_END(node)
SAMBA_PARTITION_IMPL_DEFAULT_END(edge)
SAMBA_PARTITION_IMPL_DEFAULT_END(face)
SAMBA_PARTITION_IMPL_DEFAULT_END(element)

SAMBA_PARTITION_IMPL_RANGE(node)
SAMBA_PARTITION_IMPL_RANGE(edge)
SAMBA_PARTITION_IMPL_RANGE(face)
SAMBA_PARTITION_IMPL_RANGE(element)

SAMBA_PARTITION_IMPL_BEGIN(node)
SAMBA_PARTITION_IMPL_BEGIN(edge)
SAMBA_PARTITION_IMPL_BEGIN(face)
SAMBA_PARTITION_IMPL_BEGIN(element)

SAMBA_PARTITION_IMPL_END(node)
SAMBA_PARTITION_IMPL_END(edge)
SAMBA_PARTITION_IMPL_END(face)
SAMBA_PARTITION_IMPL_END(element)

SAMBA_PARTITION_IMPL_NUM_CONNECTIVITY(node)
SAMBA_PARTITION_IMPL_NUM_CONNECTIVITY(edge)
SAMBA_PARTITION_IMPL_NUM_CONNECTIVITY(face)
SAMBA_PARTITION_IMPL_NUM_CONNECTIVITY(element)


//*************************************************************************
//partition querys
//*************************************************************************

inline samba::spatial_dimension partition_impl::spatial_dimension() const
{ return m_mesh->spatial_dimension(); }

inline partition_id partition_impl::partition() const
{ return m_partition; }

inline entity_topology partition_impl::topology() const
{ return m_parts[1].topology(); }

inline entity_rank partition_impl::rank() const
{ return m_parts[0].rank(); }

inline entity_part_vector const& partition_impl::parts() const
{ return m_parts; }

inline size_t partition_impl::size() const
{ return m_num_entities; }

inline bool partition_impl::empty() const
{ return size() == 0; }

inline convertable_to_entity_proxy partition_impl::first() const
{
  BOOST_ASSERT_MSG( !empty(), "out of range" );
  return convertable_to_entity_proxy(this, partition_offset::create(0));
}

inline convertable_to_entity_proxy partition_impl::last() const
{
  BOOST_ASSERT_MSG( !empty(), "out of range" );
  return convertable_to_entity_proxy(this, partition_offset::create(size()-1));
}

inline
samba::connectivity_kind partition_impl::connectivity_kind(entity_rank r) const
{
  switch(r()) {
  case entity_rank::node_type::value:    return m_node_kind;
  case entity_rank::edge_type::value:    return m_edge_kind;
  case entity_rank::face_type::value:    return m_face_kind;
  case entity_rank::element_type::value: return m_element_kind;
  default:
    BOOST_ASSERT_MSG(false, "Invalid rank");
  }

  return connectivity_kind::invalid();
}

template <typename EntityRank>
inline
samba::connectivity_kind partition_impl::connectivity_kind() const
{ return connectivity_kind(entity_rank::create(EntityRank::value)); }

//*************************************************************************
//entity querys
//*************************************************************************

inline convertable_to_entity_proxy partition_impl::operator[](partition_index d) const
{ return (*this)[d.offset()]; }

inline convertable_to_entity_proxy partition_impl::operator[](partition_offset o) const
{
  BOOST_ASSERT_MSG( o < size()
                    , (debug_message() << o << "out of range")
                   );
  return convertable_to_entity_proxy(this, o);
}

inline
convertable_to_entity_proxy partition_impl::operator[](unsigned i) const
{ return (*this)[partition_offset::create(i)]; }

inline entity_key partition_impl::key(partition_offset offset) const
{ return m_mesh->convert(descriptor(offset)); }

inline partition_index partition_impl::descriptor(partition_offset offset) const
{
  BOOST_ASSERT_MSG( offset < size()
                    , (debug_message() << offset << "out of range")
                   );
  return partition_index::create(rank(),m_partition,offset);
}

template <typename RankIndex>
inline
RankIndex partition_impl::rank_index(partition_offset offset) const
{
  BOOST_ASSERT_MSG(rank()() == index_type_to_entity_rank<RankIndex>::type::value,
                   "Error: Cannot convert offset to a rank index of that rank");
  BOOST_ASSERT_MSG(!m_in_modification_cycle,
                   "Cannot get rank index during modification");

  return RankIndex::create(m_starting_rank_count + offset());
}


//*************************************************************************
//modification api
//*************************************************************************

inline void partition_impl::add_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal, connectivity_orientation orientation)
{
  entity_rank to_rank  = topology_rank(to.topology(), m_mesh->spatial_dimension());

  switch(to_rank()) {
  case entity_rank::node_type::value:
    switch(m_node_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  case entity_rank::edge_type::value:
    switch(m_edge_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  case entity_rank::face_type::value:
    switch(m_face_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  case entity_rank::element_type::value:
    switch(m_element_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.add_connectivity(from,to,ordinal,orientation); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  default:
    BOOST_ASSERT_MSG(false, "Invalid rank"); // TODO:  clean up error messages
  }

  check_invariant_helper();
}

inline void partition_impl::remove_connectivity(partition_offset from, entity_key to, connectivity_ordinal ordinal)
{
  entity_rank to_rank  = topology_rank(to.topology(), m_mesh->spatial_dimension());

  switch(to_rank()) {
  case entity_rank::node_type::value:
    switch(m_node_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.remove_connectivity(from,to,ordinal); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.remove_connectivity(from,to,ordinal); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  case entity_rank::edge_type::value:
    switch(m_edge_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.remove_connectivity(from,to,ordinal); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.remove_connectivity(from,to,ordinal); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  case entity_rank::face_type::value:
    switch(m_face_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.remove_connectivity(from,to,ordinal); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.remove_connectivity(from,to,ordinal); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  case entity_rank::element_type::value:
    switch(m_element_kind()) {
    case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.remove_connectivity(from,to,ordinal); break;
    case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.remove_connectivity(from,to,ordinal); break;
    default:
      BOOST_ASSERT_MSG(false, "Invalid connetivity"); // TODO:  clean up error messages
    } break;
  default:
    BOOST_ASSERT_MSG(false, "Invalid rank"); // TODO:  clean up error messages
  }

  check_invariant_helper();
}

inline
void partition_impl::swap(partition_offset first, partition_offset second)
{
  switch(m_node_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.swap(first,second); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.swap(first,second); break;
  default:
    break;
  }

  switch(m_edge_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.swap(first,second); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.swap(first,second); break;
  default:
    break;
  }

  switch(m_face_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.swap(first,second); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.swap(first,second); break;
  default:
    break;
  }

  switch(m_element_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.swap(first,second); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.swap(first,second); break;
  default:
    break;
  }

  check_invariant_helper();
}

inline
void partition_impl::begin_modification_impl()
{
  BOOST_ASSERT_MSG(!m_in_modification_cycle, "Already in modification cycle");

  m_in_modification_cycle = true;

  switch(m_node_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.begin_modification(); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.begin_modification(); break;
  default:
    break;
  }

  switch(m_edge_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.begin_modification(); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.begin_modification(); break;
  default:
    break;
  }

  switch(m_face_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.begin_modification(); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.begin_modification(); break;
  default:
    break;
  }

  switch(m_element_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.begin_modification(); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.begin_modification(); break;
  default:
    break;
  }

  check_invariant_helper();
}

inline
void partition_impl::end_modification_impl()
{
  BOOST_ASSERT_MSG(m_in_modification_cycle, "Not in modification cycle");

  m_in_modification_cycle = false;

  switch(m_node_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.end_modification(m_mesh); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.end_modification(m_mesh); break;
  default:
    break;
  }

  switch(m_edge_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.end_modification(m_mesh); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.end_modification(m_mesh); break;
  default:
    break;
  }

  switch(m_face_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.end_modification(m_mesh); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.end_modification(m_mesh); break;
  default:
    break;
  }

  switch(m_element_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.end_modification(m_mesh); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.end_modification(m_mesh); break;
  default:
    break;
  }

  check_invariant_helper();
}

inline
void partition_impl::add_entities_impl(size_t how_many)
{
  m_num_entities += how_many;

  switch(m_node_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.add_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.add_entities(how_many); break;
  default:
    break;
  }

  switch(m_edge_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.add_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.add_entities(how_many); break;
  default:
    break;
  }

  switch(m_face_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.add_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.add_entities(how_many); break;
  default:
    break;
  }

  switch(m_element_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.add_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.add_entities(how_many); break;
  default:
    break;
  }

  check_invariant_helper();
}

inline
void partition_impl::remove_entities_impl(size_t how_many)
{
  m_num_entities -= how_many;

  switch(m_node_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.remove_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.remove_entities(how_many); break;
  default:
    break;
  }

  switch(m_edge_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.remove_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.remove_entities(how_many); break;
  default:
    break;
  }

  switch(m_face_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.remove_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.remove_entities(how_many); break;
  default:
    break;
  }

  switch(m_element_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.remove_entities(how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.remove_entities(how_many); break;
  default:
    break;
  }

  check_invariant_helper();
}

inline
void partition_impl::move_entities_impl(partition_impl& to, size_t how_many)
{
  m_num_entities    -= how_many;
  to.m_num_entities += how_many;

  switch(m_node_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_node_connectivity.move_entities(to.m_fixed_node_connectivity,how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_node_connectivity.move_entities(to.m_dynamic_node_connectivity,how_many); break;
  default:
    break;
  }

  switch(m_edge_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_edge_connectivity.move_entities(to.m_fixed_edge_connectivity,how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_edge_connectivity.move_entities(to.m_dynamic_edge_connectivity,how_many); break;
  default:
    break;
  }

  switch(m_face_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_face_connectivity.move_entities(to.m_fixed_face_connectivity,how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_face_connectivity.move_entities(to.m_dynamic_face_connectivity,how_many); break;
  default:
    break;
  }

  switch(m_element_kind()) {
  case connectivity_kind::fixed_type::value:   m_fixed_element_connectivity.move_entities(to.m_fixed_element_connectivity,how_many); break;
  case connectivity_kind::dynamic_type::value: m_dynamic_element_connectivity.move_entities(to.m_dynamic_element_connectivity,how_many); break;
  default:
    break;
  }

  check_invariant_helper();
  to.check_invariant_helper();
}

inline
void partition_impl::check_invariant_helper() const
{
  switch(m_node_kind()) {
  case connectivity_kind::fixed_type::value:
    BOOST_ASSERT_MSG(m_fixed_node_connectivity.size() == size(),
                     (debug_message() << "Size of fixed-node connectivity " << m_fixed_node_connectivity.size() << " does not match partition size " << size()));
    break;
  case connectivity_kind::dynamic_type::value:
    BOOST_ASSERT_MSG(m_dynamic_node_connectivity.size() == size(),
                     (debug_message() << "Size of dynamic-node connectivity " << m_dynamic_node_connectivity.size() << " does not match partition size " << size()));
    break;
  default:
    break;
  }

  switch(m_edge_kind()) {
  case connectivity_kind::fixed_type::value:
    BOOST_ASSERT_MSG(m_fixed_edge_connectivity.size() == size(),
                     (debug_message() << "Size of fixed-edge connectivity " << m_fixed_edge_connectivity.size() << " does not match partition size " << size()));
    break;
  case connectivity_kind::dynamic_type::value:
    BOOST_ASSERT_MSG(m_dynamic_edge_connectivity.size() == size(),
                     (debug_message() << "Size of dynamic-edge connectivity " << m_dynamic_edge_connectivity.size() << " does not match partition size " << size()));
    break;
  default:
    break;
  }

  switch(m_face_kind()) {
  case connectivity_kind::fixed_type::value:
    BOOST_ASSERT_MSG(m_fixed_face_connectivity.size() == size(),
                     (debug_message() << "Size of fixed-face connectivity " << m_fixed_face_connectivity.size() << " does not match partition size " << size()));
    break;
  case connectivity_kind::dynamic_type::value:
    BOOST_ASSERT_MSG(m_dynamic_face_connectivity.size() == size(),
                     (debug_message() << "Size of dynamic-face connectivity " << m_dynamic_face_connectivity.size() << " does not match partition size " << size()));
    break;
  default:
    break;
  }

  switch(m_element_kind()) {
  case connectivity_kind::fixed_type::value:
    BOOST_ASSERT_MSG(m_fixed_element_connectivity.size() == size(),
                     (debug_message() << "Size of fixed-element connectivity " << m_fixed_element_connectivity.size() << " does not match partition size " << size()));
    break;
  case connectivity_kind::dynamic_type::value:
    BOOST_ASSERT_MSG(m_dynamic_element_connectivity.size() == size(),
                     (debug_message() << "Size of dynamic-element connectivity " << m_dynamic_element_connectivity.size() << " does not match partition size " << size()));
    break;
  default:
    break;
  }
}

}} //namespace samba::detail

#undef SAMBA_PARTITION_IMPL_DEFAULT_RANGE
#undef SAMBA_PARTITION_IMPL_DEFAULT_BEGIN
#undef SAMBA_PARTITION_IMPL_DEFAULT_END
#undef SAMBA_PARTITION_IMPL_RANGE
#undef SAMBA_PARTITION_IMPL_BEGIN
#undef SAMBA_PARTITION_IMPL_END
#undef SAMBA_PARTITION_IMPL_NUM_CONNECTIVITY

#endif //SAMBA_SAMBA_MESH_PARTITION_IMPL_TCC

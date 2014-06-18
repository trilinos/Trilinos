#ifndef SAMBA_SAMBA_ENTITY_ENTITY_PROXY_TCC
#define SAMBA_SAMBA_ENTITY_ENTITY_PROXY_TCC

#include <samba/types.hpp>

namespace samba {

//defined here to break cyclic dependency

//*************************************************************************
//connectivity accessors
//*************************************************************************

inline partition_index_range entity_proxy::nodes() const
{ return m_partition->nodes(offset()); }
inline partition_index_range entity_proxy::edges() const
{ return m_partition->edges(offset()); }
inline partition_index_range entity_proxy::faces() const
{ return m_partition->faces(offset()); }
inline partition_index_range entity_proxy::elements() const
{ return m_partition->elements(offset()); }

inline partition_index_iterator entity_proxy::begin_nodes() const
{ return m_partition->begin_nodes(offset()); }
inline partition_index_iterator entity_proxy::begin_edges() const
{ return m_partition->begin_edges(offset()); }
inline partition_index_iterator entity_proxy::begin_faces() const
{ return m_partition->begin_faces(offset()); }
inline partition_index_iterator entity_proxy::begin_elements() const
{ return m_partition->begin_elements(offset()); }

inline partition_index_iterator entity_proxy::end_nodes() const
{ return m_partition->end_nodes(offset()); }
inline partition_index_iterator entity_proxy::end_edges() const
{ return m_partition->end_edges(offset()); }
inline partition_index_iterator entity_proxy::end_faces() const
{ return m_partition->end_faces(offset()); }
inline partition_index_iterator entity_proxy::end_elements() const
{ return m_partition->end_elements(offset()); }

inline size_t entity_proxy::num_connectivity(entity_rank rank) const
{ return m_partition->num_connectivity(rank, offset()); }
inline size_t entity_proxy::num_nodes() const
{ return m_partition->num_nodes(offset()); }
inline size_t entity_proxy::num_edges() const
{ return m_partition->num_edges(offset()); }
inline size_t entity_proxy::num_faces() const
{ return m_partition->num_faces(offset()); }
inline size_t entity_proxy::num_elements() const
{ return m_partition->num_elements(offset()); }

template <typename T>
inline std::pair<const T*, const T*> entity_proxy::connectivity(entity_rank rank) const
{ return m_partition->connectivity<T>(rank, offset()); }
template <typename T>
inline std::pair<const T*, const T*> entity_proxy::nodes() const
{ return m_partition->nodes<T>(offset()); }
template <typename T>
inline std::pair<const T*, const T*> entity_proxy::edges() const
{ return m_partition->edges<T>(offset()); }
template <typename T>
inline std::pair<const T*, const T*> entity_proxy::faces() const
{ return m_partition->faces<T>(offset()); }
template <typename T>
inline std::pair<const T*, const T*> entity_proxy::elements() const
{ return m_partition->elements<T>(offset()); }

template <typename T>
inline const T* entity_proxy::begin_connectivity(entity_rank rank) const
{ return m_partition->begin_connectivity<T>(rank, offset()); }
template <typename T>
inline const T* entity_proxy::begin_nodes() const
{ return m_partition->begin_nodes<T>(offset()); }
template <typename T>
inline const T* entity_proxy::begin_edges() const
{ return m_partition->begin_edges<T>(offset()); }
template <typename T>
inline const T* entity_proxy::begin_faces() const
{ return m_partition->begin_faces<T>(offset()); }
template <typename T>
inline const T* entity_proxy::begin_elements() const
{ return m_partition->begin_elements<T>(offset()); }

template <typename T>
inline const T* entity_proxy::end_connectivity(entity_rank rank) const
{ return m_partition->end_connectivity<T>(rank, offset()); }
template <typename T>
inline const T* entity_proxy::end_nodes() const
{ return m_partition->end_nodes<T>(offset()); }
template <typename T>
inline const T* entity_proxy::end_edges() const
{ return m_partition->end_edges<T>(offset()); }
template <typename T>
inline const T* entity_proxy::end_faces() const
{ return m_partition->end_faces<T>(offset()); }
template <typename T>
inline const T* entity_proxy::end_elements() const
{ return m_partition->end_elements<T>(offset()); }

//*************************************************************************
//entity querys
//*************************************************************************
inline entity_part_range entity_proxy::parts() const
{ return m_partition->parts(); }

inline entity_key entity_proxy::key() const
{ return m_partition->key(m_offset); }

inline partition_index entity_proxy::descriptor() const
{ return m_partition->descriptor(m_offset); }

template <>
inline
entity_key entity_proxy::get<entity_key>() const
{ return key(); }

template <>
inline
partition_index entity_proxy::get<partition_index>() const
{ return descriptor(); }

template <typename RankIndex>
inline
RankIndex entity_proxy::get() const
{ return m_partition->rank_index<RankIndex>(m_offset); }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_ENTITY_PROXY_TCC

#ifndef SAMBA_SAMBA_MESH_PARTITION_PROXY_TCC
#define SAMBA_SAMBA_MESH_PARTITION_PROXY_TCC

#include <samba/types.hpp>

namespace samba {

//defined here to break cyclic dependency

//*************************************************************************
//connectivity accessors
//*************************************************************************

inline partition_index_range partition_proxy::nodes(partition_offset offset) const
{ return m_partition->nodes(offset); }
inline partition_index_range partition_proxy::edges(partition_offset offset) const
{ return m_partition->edges(offset); }
inline partition_index_range partition_proxy::faces(partition_offset offset) const
{ return m_partition->faces(offset); }
inline partition_index_range partition_proxy::elements(partition_offset offset) const
{ return m_partition->elements(offset); }

inline partition_index_iterator partition_proxy::begin_nodes(partition_offset offset) const
{ return m_partition->begin_nodes(offset); }
inline partition_index_iterator partition_proxy::begin_edges(partition_offset offset) const
{ return m_partition->begin_edges(offset); }
inline partition_index_iterator partition_proxy::begin_faces(partition_offset offset) const
{ return m_partition->begin_faces(offset); }
inline partition_index_iterator partition_proxy::begin_elements(partition_offset offset) const
{ return m_partition->begin_elements(offset); }

inline partition_index_iterator partition_proxy::end_nodes(partition_offset offset) const
{ return m_partition->end_nodes(offset); }
inline partition_index_iterator partition_proxy::end_edges(partition_offset offset) const
{ return m_partition->end_edges(offset); }
inline partition_index_iterator partition_proxy::end_faces(partition_offset offset) const
{ return m_partition->end_faces(offset); }
inline partition_index_iterator partition_proxy::end_elements(partition_offset offset) const
{ return m_partition->end_elements(offset); }

inline size_t partition_proxy::num_connectivity(entity_rank rank, partition_offset offset) const
{ return m_partition->num_connectivity(rank, offset); }
inline size_t partition_proxy::num_nodes(partition_offset offset) const
{ return m_partition->num_nodes(offset); }
inline size_t partition_proxy::num_edges(partition_offset offset) const
{ return m_partition->num_edges(offset); }
inline size_t partition_proxy::num_faces(partition_offset offset) const
{ return m_partition->num_faces(offset); }
inline size_t partition_proxy::num_elements(partition_offset offset) const
{ return m_partition->num_elements(offset); }


template <typename T>
inline std::pair<const T*, const T*> partition_proxy::connectivity(entity_rank rank, partition_offset offset) const
{ return m_partition->connectivity<T>(rank, offset); }
template <typename T>
inline std::pair<const T*, const T*> partition_proxy::nodes(partition_offset offset) const
{ return m_partition->nodes<T>(offset); }
template <typename T>
inline std::pair<const T*, const T*> partition_proxy::edges(partition_offset offset) const
{ return m_partition->edges<T>(offset); }
template <typename T>
inline std::pair<const T*, const T*> partition_proxy::faces(partition_offset offset) const
{ return m_partition->faces<T>(offset); }
template <typename T>
inline std::pair<const T*, const T*> partition_proxy::elements(partition_offset offset) const
{ return m_partition->elements<T>(offset); }

template <typename T>
inline const T* partition_proxy::begin_connectivity(entity_rank rank, partition_offset offset) const
{ return m_partition->begin_connectivity<T>(rank, offset); }
template <typename T>
inline const T* partition_proxy::begin_nodes(partition_offset offset) const
{ return m_partition->begin_nodes<T>(offset); }
template <typename T>
inline const T* partition_proxy::begin_edges(partition_offset offset) const
{ return m_partition->begin_edges<T>(offset); }
template <typename T>
inline const T* partition_proxy::begin_faces(partition_offset offset) const
{ return m_partition->begin_faces<T>(offset); }
template <typename T>
inline const T* partition_proxy::begin_elements(partition_offset offset) const
{ return m_partition->begin_elements<T>(offset); }

template <typename T>
inline const T* partition_proxy::end_connectivity(entity_rank rank, partition_offset offset) const
{ return m_partition->end_connectivity<T>(rank, offset); }
template <typename T>
inline const T* partition_proxy::end_nodes(partition_offset offset) const
{ return m_partition->end_nodes<T>(offset); }
template <typename T>
inline const T* partition_proxy::end_edges(partition_offset offset) const
{ return m_partition->end_edges<T>(offset); }
template <typename T>
inline const T* partition_proxy::end_faces(partition_offset offset) const
{ return m_partition->end_faces<T>(offset); }
template <typename T>
inline const T* partition_proxy::end_elements(partition_offset offset) const
{ return m_partition->end_elements<T>(offset); }


//*************************************************************************
//partition querys
//*************************************************************************
inline samba::spatial_dimension partition_proxy::spatial_dimension() const
{ return m_partition->spatial_dimension(); }

inline samba::partition_id partition_proxy::partition_id() const
{ return m_partition->partition(); }

inline entity_topology partition_proxy::topology() const
{ return m_partition->topology(); }

inline entity_rank partition_proxy::rank() const
{ return m_partition->rank(); }

inline entity_part_range partition_proxy::parts() const
{ return m_partition->parts(); }

inline size_t partition_proxy::size() const
{ return m_partition->size(); }

inline bool partition_proxy::empty() const
{ return m_partition->empty(); }

inline entity_proxy partition_proxy::first() const
{ return m_partition->first(); }

inline entity_proxy partition_proxy::last() const
{ return m_partition->last(); }

inline samba::connectivity_kind partition_proxy::connectivity_kind(entity_rank r) const
{ return m_partition->connectivity_kind(r); }

//*************************************************************************
//entity querys
//*************************************************************************

inline entity_proxy partition_proxy::operator[](partition_index d) const
{ return (*m_partition)[d]; }

inline entity_proxy partition_proxy::operator[](partition_offset o) const
{ return (*m_partition)[o]; }

inline entity_proxy partition_proxy::operator[](unsigned i) const
{ return (*m_partition)[i]; }

inline entity_key partition_proxy::key(partition_offset offset) const
{ return m_partition->key(offset); }

inline partition_index partition_proxy::descriptor(partition_offset offset) const
{ return m_partition->descriptor(offset); }

template <typename RankIndex>
inline
RankIndex partition_proxy::rank_index(partition_offset offset) const
{ return m_partition->rank_index<RankIndex>(offset); }

} //namespace samba

#endif //SAMBA_SAMBA_MESH_PARTITION_PROXY_TCC

/*
 * Partition_inl.hpp
 *
 * Not really template cc file, but needs to be included like one.
 *
 */
#ifndef STK_MESH_IMPL_PARTITION_TCC_
#define STK_MESH_IMPL_PARTITION_TCC_

#include  <stk_mesh/baseImpl/BucketRepository.hpp>

namespace stk {
namespace mesh {
namespace impl {

inline std::vector<Bucket *>::iterator Partition::begin()
{
    if (m_modifyingBucketSet)
    {
        return m_buckets.begin();
    }

    std::vector<Bucket *> &buckets = m_repository->m_buckets[m_rank];
    if (buckets.empty())
    {
        return buckets.end();
    }
    return buckets.begin() + m_beginBucketIndex;
}


inline std::vector<Bucket *>::iterator Partition::end()
{
    if (m_modifyingBucketSet)
    {
        return m_buckets.end();
    }

    std::vector<Bucket *> &buckets = m_repository->m_buckets[m_rank];
    if (buckets.empty())
    {
        return buckets.end();
    }
    return buckets.begin() + m_endBucketIndex;
}


inline std::vector<Bucket *>::const_iterator Partition::begin() const
{
    if (m_modifyingBucketSet)
    {
        return m_buckets.begin();
    }

    const std::vector<Bucket *> &buckets = m_repository->m_buckets[m_rank];
    if (buckets.empty())
    {
        return buckets.end();
    }
    return buckets.begin() + m_beginBucketIndex;
}


inline std::vector<Bucket *>::const_iterator Partition::end() const
{
    if (m_modifyingBucketSet)
    {
        return m_buckets.end();
    }

    const std::vector<Bucket *> &buckets = m_repository->m_buckets[m_rank];
    if (buckets.empty())
    {
        return buckets.end();
    }
    return buckets.begin() + m_endBucketIndex;
}


} // impl
} // mesh
} // stk

#endif

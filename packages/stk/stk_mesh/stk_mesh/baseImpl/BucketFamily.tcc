/*
 * BucketFamilyImpl.tcc
 *
 * Not really template cc file, but needs to be included like one.
 *
 */
#ifndef STK_MESH_IMPL_BUCKETFAMILY_TCC_
#define STK_MESH_IMPL_BUCKETFAMILY_TCC_

#include  <stk_mesh/baseImpl/BucketRepository.hpp>

namespace stk {
namespace mesh {
namespace impl {

inline std::vector<Bucket *>::iterator BucketFamily::begin() const
{
	std::vector<Bucket *> &buckets = m_repository->m_buckets[m_rank];
	if (buckets.empty())
	{
		return buckets.end();
	}
	return buckets.begin() + m_beginBucketIndex;
}


inline std::vector<Bucket *>::iterator BucketFamily::end() const
{
	std::vector<Bucket *> &buckets = m_repository->m_buckets[m_rank];
	if (buckets.empty())
	{
		return buckets.end();
	}
	return buckets.begin() + m_endBucketIndex;
}

inline bool BucketFamily::belongs(Bucket *bkt) const
{
	return bkt->getBucketFamily() == this;
}

} // impl
} // mesh
} // stk

#endif

/*
 * BucketFamily.hpp
 *
 */

#ifndef STK_MESH_IMPL_BUCKETFAMILY_HPP_
#define STK_MESH_IMPL_BUCKETFAMILY_HPP_

#include <stk_mesh/base/Types.hpp>


namespace stk {
namespace mesh {

namespace impl {

class BucketFamily
{
    friend class BucketRepository;
public:
    BucketFamily(EntityRank rank = 0)
    : m_rank(rank), m_beginBucketIndex(0), m_endBucketIndex(0) { }

    virtual ~BucketFamily();

private:
    EntityRank m_rank;
    std::vector<PartOrdinal> m_stkPartition;
    unsigned m_beginBucketIndex;
    unsigned m_endBucketIndex;
};

} // impl
} // mesh
} // stk

#endif /* BUCKETFAMILY_HPP_ */

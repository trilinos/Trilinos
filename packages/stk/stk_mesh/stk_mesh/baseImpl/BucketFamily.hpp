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

class BucketRepository;

class BucketFamily
{
    friend class BucketRepository;
public:
    BucketFamily(BucketRepository *repo, EntityRank rank)
    : m_repository(repo), m_rank(rank), m_beginBucketIndex(0), m_endBucketIndex(0) { }

    virtual ~BucketFamily();

    std::ostream &streamit(std::ostream &os) const;

    const EntityRank get_rank() const { return m_rank; }

    inline std::vector<Bucket *>::iterator begin() const;
    inline std::vector<Bucket *>::iterator end() const;

    const std::vector<PartOrdinal> &get_legacy_partition_id() const
    {
        return m_stkPartition;
    }

    inline bool belongs(Bucket *bkt) const;

    // Just for unit testing.  Remove after refactor.
    static BucketRepository &getRepository(stk::mesh::BulkData &mesh);

private:

    BucketRepository *m_repository;
    EntityRank m_rank;
    std::vector<PartOrdinal> m_stkPartition;
    unsigned m_beginBucketIndex;
    unsigned m_endBucketIndex;
};

std::ostream &operator<<(std::ostream &, const stk::mesh::impl::BucketFamily &);

} // impl
} // mesh
} // stk



#endif /* BUCKETFAMILY_HPP_ */

/*
 * BucketFamily.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: pgxavie
 */
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/baseImpl/BucketFamily.hpp>

#include <iostream>

namespace stk {
namespace mesh {
namespace impl {


std::ostream &operator<<(std::ostream &os, const stk::mesh::impl::BucketFamily &bf)
{
    return bf.streamit(os);
}

} // impl
} // mesh
} // stk

using namespace stk::mesh::impl;

std::ostream &BucketFamily::streamit(std::ostream &os) const
{
    os << "{BucketFamily " << std::endl
       << "  m_repository = " << m_repository << "  m_rank = " << m_rank << std::endl
       << "  (" << m_beginBucketIndex << ", " << m_endBucketIndex << ")" << std::endl;

    os << "  legacy partition id : {";
    const std::vector<unsigned> &family_key = get_legacy_partition_id();
    for (size_t i = 0; i < family_key.size(); ++i)
    {
        os << " " << family_key[i];
    }
    os << " }"  << std::endl << "}";

    return os;
}

BucketFamily::~BucketFamily()
{
    // TODO Auto-generated destructor stub
}

BucketRepository &BucketFamily::getRepository(stk::mesh::BulkData &mesh)
{
    return mesh.m_bucket_repository;
}


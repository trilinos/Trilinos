/*
 * Partition.cpp
 *
 *  Created on: Oct 22, 2012
 */
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>

#include <iostream>

namespace stk {
namespace mesh {
namespace impl {


std::ostream &operator<<(std::ostream &os, const stk::mesh::impl::Partition &bf)
{
    return bf.streamit(os);
}

} // impl
} // mesh
} // stk

using namespace stk::mesh::impl;

std::ostream &Partition::streamit(std::ostream &os) const
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

Partition::~Partition()
{
    // TODO Auto-generated destructor stub
}

BucketRepository &Partition::getRepository(stk::mesh::BulkData &mesh)
{
    return mesh.m_bucket_repository;
}

void Partition::sort()
{
    // Gather up the entities in this partition.

    // Sort the entities.

    // More to come...
}

/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stk_ngp/Ngp.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/baseImpl/ForEachEntityLoopAbstractions.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/StkVector.hpp>
#include "stk_mesh/base/FieldParallel.hpp"

#include <limits>

class NgpMeshTest : public stk::unit_test_util::MeshFixture
{
public:
    void run_get_nodes_using_FastMeshIndex_test()
    {
        setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

        stk::Vector<double> numNodesVec("numNodes", 1);

        ngp::Mesh ngpMesh(get_bulk());
        Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
        {
            ngp::Mesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, stk::mesh::FastMeshIndex{0,0});
            numNodesVec.device_get(0) = nodes.size();
        });

        numNodesVec.copy_device_to_host();
        ASSERT_EQ(8u, numNodesVec[0]);
    }
};
TEST_F(NgpMeshTest, get_nodes_using_FastMeshIndex)
{
    run_get_nodes_using_FastMeshIndex_test();
}

class EntityIndexSpace : public stk::unit_test_util::MeshFixture {};
TEST_F(EntityIndexSpace, accessingLocalData_useLocalOffset)
{
    setup_mesh("generated:1x1x1", stk::mesh::BulkData::NO_AUTO_AURA);
    std::vector<unsigned> entityToLocalOffset(get_bulk().get_size_of_entity_index_space(), 0);

    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<get_meta().entity_rank_count(); ++rank)
    {
        unsigned localOffset = 0;
        const stk::mesh::BucketVector &buckets = get_bulk().buckets(stk::topology::NODE_RANK);
        for(const stk::mesh::Bucket *bucket : buckets)
        {
            for(stk::mesh::Entity entity : *bucket)
            {
                entityToLocalOffset[entity.local_offset()] = localOffset;
                localOffset++;
            }
        }
    }

    std::vector<unsigned> gold {0,0,1,2,3,4,5,6,7,0};
    ASSERT_EQ(gold.size(), entityToLocalOffset.size());
    for(size_t i=0; i<gold.size(); i++)
    {
        EXPECT_EQ(gold[i], entityToLocalOffset[i]);
    }
}

//void fails_build() { Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) { std::max(1, 2); }); }
//TEST_F(NodalVolumeCalculator, fails)  { fails_build(); }

//KOKKOS_INLINE_FUNCTION double my_max(double a, double b) { return std::max(a, b); }
//void builds_ok() { Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& i) { my_max(1, 2); }); }
//TEST_F(NodalVolumeCalculator, builds) { builds_ok(); }

void run_vector_gpu_test()
{
    size_t n = 10;
    stk::Vector<double> vec("vec", n);
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int i)
    {
        vec.device_get(i) = i;
    });
    vec.copy_device_to_host();
    for(size_t i=0; i<n; i++)
        EXPECT_EQ(i, vec[i]);
}
TEST(StkVectorGpuTest, gpu_runs)
{
    run_vector_gpu_test();
}



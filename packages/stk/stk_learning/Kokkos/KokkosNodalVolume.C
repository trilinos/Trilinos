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

#include "mtk_kokkos.h"
#include "StaticMesh.h"
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/stk_config.h>
#include "../../stk_mesh/stk_mesh/base/FieldParallel.hpp"

#include <limits>

namespace ngp {



template<typename T>
class StaticField {
public:
    typedef Kokkos::View<T*> FieldDataType;

    STK_FUNCTION
    StaticField(stk::mesh::EntityRank rank, const T& initialValue, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    : ngpBulk(bulk), deviceData(), numPerEntity(1)
    {
        unsigned alloc_size = compute_alloc_size(rank, bulk, selector);
        create_device_field_data(alloc_size, initialValue);
    }

    StaticField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &field)
    : ngpBulk(bulk), deviceData(), numPerEntity(0)
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);

        if(!buckets.empty())
        {
            numPerEntity = field.max_size(field.entity_rank());
            size_t allocationSize = buckets.size()*buckets[0]->capacity()*numPerEntity;
            deviceData = FieldDataType("deviceData"+std::to_string(numPerEntity), allocationSize);

            hostData = Kokkos::create_mirror_view(deviceData);

            for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            {
                const stk::mesh::Bucket &bucket = *buckets[iBucket];
                const T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
                for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
                {
                    for(unsigned j=0; j<numPerEntity; j++)
                    {
                        unsigned index = iBucket*bucket.capacity()*numPerEntity + iEntity*numPerEntity + j;
                        hostData(index) = data[iEntity*numPerEntity + j];
                    }
                }
            }

            Kokkos::deep_copy(deviceData, hostData);
        }
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
        stk::mesh::Selector selector = stk::mesh::selectField(field);
        const stk::mesh::BucketVector& buckets = bulk.get_buckets(field.entity_rank(), selector);

        if(!buckets.empty())
        {
            Kokkos::deep_copy(hostData, deviceData);

            for(size_t iBucket=0; iBucket<buckets.size(); iBucket++)
            {
                const stk::mesh::Bucket &bucket = *buckets[iBucket];
                T* data = static_cast<T*>(stk::mesh::field_data(field, bucket));
                for(size_t iEntity=0; iEntity<bucket.size(); iEntity++)
                {
                    for(unsigned j=0; j<numPerEntity; j++)
                    {
                        unsigned index = iBucket*bucket.capacity()*numPerEntity + iEntity*numPerEntity + j;
                        data[iEntity*numPerEntity + j] = hostData(index);
                    }
                }
            }
        }
    }

    void clear()
    {
        ngpBulk.clear();
    }

    STK_FUNCTION
    ~StaticField(){}

    STK_FUNCTION
    unsigned get_index(stk::mesh::Entity entity) const
    {
        const StaticMesh& localRef = ngpBulk;
        unsigned bkt_id = localRef.mesh_index(entity).bucket_id;
        unsigned bkt_ord = localRef.mesh_index(entity).bucket_ord;
        unsigned idx = bkt_id*512*numPerEntity + bkt_ord*numPerEntity;
        return idx;
    }

    STK_FUNCTION
    T* operator[](stk::mesh::Entity entity) const
    {
        return &deviceData(get_index(entity));
    }

    STK_FUNCTION
    T& get(stk::mesh::Entity entity, int component) const
    {
        unsigned idx = get_index(entity);
        return deviceData(idx+component);
    }

private:
    unsigned compute_alloc_size(stk::mesh::EntityRank rank, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector)
    {
        const stk::mesh::BucketVector& bkts = bulk.get_buckets(rank, selector);
        unsigned alloc_size = 0;
        for(const stk::mesh::Bucket* bktptr : bkts) {
            const stk::mesh::Bucket& bkt = *bktptr;
            alloc_size += bkt.capacity();
        }
        return alloc_size*numPerEntity;
    }

    void create_device_field_data(unsigned allocSize, const T& initialValue)
    {
        deviceData = FieldDataType("deviceData", allocSize);
        hostData = Kokkos::create_mirror_view(deviceData);

        for(size_t i=0; i<allocSize; ++i)
            hostData(i) = initialValue;
        Kokkos::deep_copy(deviceData, hostData);
    }

    StaticMesh ngpBulk;
    typename FieldDataType::HostMirror hostData;
    FieldDataType deviceData;
    unsigned numPerEntity;
};

}

STK_FUNCTION
double calculate_element_volume(ngp::ConnectedNodesType elemNodes,
                                unsigned numElemNodes,
                                const ngp::StaticField<double> &staticCoords)
{
    double min[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    double max[3] = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    for(unsigned i=0; i<numElemNodes; ++i) {
        for(int j=0; j<3; ++j) {
            double val = staticCoords.get(elemNodes(i),j);
            if (val > max[j]) {
                max[j] = val;
            }
                
            if (val < min[j]) {
                min[j] = val;
            }
        }
    }
    return (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
}

typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;
typedef Kokkos::TeamPolicy<ExecSpace, ScheduleType> TeamPolicyType;
typedef TeamPolicyType::member_type TeamHandleType;

void calculate_nodal_volume(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
    ngp::StaticField<double> staticCoords(mesh, coords);
    ngp::StaticField<double> staticNodalVolume(mesh, nodalVolumeField);
    ngp::StaticMesh staticMesh(mesh);
    unsigned numBuckets = staticMesh.num_buckets(stk::topology::ELEM_RANK);

    Kokkos::parallel_for(Kokkos::TeamPolicy< ExecSpace >(numBuckets, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamHandleType& team) {
        const int elementBucketIndex = team.league_rank();
        const ngp::StaticBucket &bucket = staticMesh.get_bucket(stk::topology::ELEM_RANK, elementBucketIndex);
        const unsigned numNodesPerElem = bucket.get_num_nodes_per_entity();
        unsigned numElements = bucket.size();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i) {
            ngp::ConnectedNodesType elemNodes = bucket.get_nodes(i);
            double elemVolumePerNode = calculate_element_volume(elemNodes, numNodesPerElem, staticCoords)/numNodesPerElem;
            for(unsigned j=0; j<numNodesPerElem; ++j) {
                double* nodalVolume = staticNodalVolume[elemNodes(j)];
                Kokkos::atomic_add(nodalVolume, elemVolumePerNode);
            }
        });
    });
    Kokkos::fence();
    staticNodalVolume.copy_device_to_host(mesh, nodalVolumeField);
    staticMesh.clear();
    staticCoords.clear();
    staticNodalVolume.clear();
    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

void expect_nodal_volume(const stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, stk::mesh::Field<int>& numElemsPerNodeField)
{
    stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part();
    const stk::mesh::BucketVector& nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, selector);
    for(const stk::mesh::Bucket* bucket : nodeBuckets) {
        const double* nodalVolume = stk::mesh::field_data(nodalVolumeField, *bucket);
        const int *numElemsPerNode = stk::mesh::field_data(numElemsPerNodeField, *bucket);
        for(size_t i=0; i<bucket->size(); ++i) {
            double expectedVolume = 1.0;
            if (numElemsPerNode[i] == 8) {
                expectedVolume = 1.0;
            }
            else if (numElemsPerNode[i] == 4) {
                expectedVolume = 0.5;
            }
            else if (numElemsPerNode[i] == 2) {
                expectedVolume = 0.25;
            }
            else if (numElemsPerNode[i] == 1) {
                expectedVolume = 0.125;
            }
            else {
                ASSERT_TRUE(false) << "numElemsPerNode = " << numElemsPerNode[i];
            }

            EXPECT_NEAR(expectedVolume, nodalVolume[i], 1.e-9);
        }
    }
}

void count_num_elems_per_node(const stk::mesh::BulkData& mesh, stk::mesh::Field<int>& numElemsPerNodeField)
{
    stk::mesh::Selector ownedOrShared = mesh.mesh_meta_data().locally_owned_part() | mesh.mesh_meta_data().globally_shared_part();
    const stk::mesh::BucketVector& nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, ownedOrShared);
    for(const stk::mesh::Bucket* bucket : nodeBuckets) {
        int* numElemsPerNode = stk::mesh::field_data(numElemsPerNodeField, *bucket);
        for(size_t i=0; i<bucket->size(); ++i) {
            unsigned numElems = bucket->num_elements(i);
            const stk::mesh::Entity *elems = bucket->begin_elements(i);
            for(size_t j=0; j<numElems; ++j)
            {
                if(mesh.bucket(elems[j]).owned())
                {
                    numElemsPerNode[i]++;
                }
            }
        }
    }
    stk::mesh::parallel_sum(mesh, {&numElemsPerNodeField});
}

class NodalVolumeCalculator : public stk::unit_test_util::MeshFixture
{
protected:
    void test_nodal_volume(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        auto& nodalVolumeField = get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "nodal_volume");
        stk::mesh::put_field(nodalVolumeField, get_meta().universal_part());
        auto& numElemsPerNodeField = get_meta().declare_field<stk::mesh::Field<int> >(stk::topology::NODE_RANK, "numElemsPerNode");
        stk::mesh::put_field(numElemsPerNodeField, get_meta().universal_part());
        setup_mesh("generated:20x20x20", auraOption);
        count_num_elems_per_node(get_bulk(), numElemsPerNodeField);

        calculate_nodal_volume(get_bulk(), nodalVolumeField);

        expect_nodal_volume(get_bulk(), nodalVolumeField, numElemsPerNodeField);
    }
};
TEST_F(NodalVolumeCalculator, nodalVolumeWithAura) {
    test_nodal_volume(stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeWithoutAura) {
    test_nodal_volume(stk::mesh::BulkData::NO_AUTO_AURA);
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

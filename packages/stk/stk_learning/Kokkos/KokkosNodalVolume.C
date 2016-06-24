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
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/environment/WallTime.hpp>
#include "../../stk_mesh/stk_mesh/base/FieldParallel.hpp"

#include <limits>

STK_FUNCTION
double calculate_element_volume(const ngp::StaticMesh& ngpMesh, ngp::ConnectedNodesType elemNodes,
                                unsigned numElemNodes,
                                const ngp::StaticField<double> &staticCoords)
{
    double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double max[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
    for(unsigned i=0; i<numElemNodes; ++i) {
        for(int j=0; j<3; ++j) {
            double val = staticCoords.const_get(ngpMesh, elemNodes(i),j);
            if (val > max[j])
                max[j] = val;
            if (val < min[j])
                min[j] = val;
        }
    }
    return (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
}

double calculate_element_volume(const stk::mesh::Entity* elemNodes,
                                unsigned numElemNodes,
                                const stk::mesh::FieldBase& coords)
{
    double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double max[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
    for(unsigned i=0; i<numElemNodes; ++i) {
        const double* val = static_cast<const double*>(stk::mesh::field_data(coords, elemNodes[i]));
        for(int j=0; j<3; ++j) {
            if (val[j] > max[j])
                max[j] = val[j];
            if (val[j] < min[j])
                min[j] = val[j];
        }
    }
    return (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
}

STK_FUNCTION
void calculate_nodal_volume_given_elem_nodes(const ngp::StaticMesh &staticMesh,
                                             const ngp::ConnectedNodesType& elemNodes,
                                             const ngp::StaticField<double>& staticCoords,
                                             const ngp::StaticField<double>& staticNodalVolume)
{
    const unsigned numNodesPerElem = elemNodes.size();
    double elemVolumePerNode = calculate_element_volume(staticMesh, elemNodes, numNodesPerElem, staticCoords) / numNodesPerElem;
    for(unsigned j = 0; j < numNodesPerElem; ++j)
    {
        double* nodalVolume = &staticNodalVolume.get(staticMesh, elemNodes(j), 0);
        Kokkos::atomic_add(nodalVolume, elemVolumePerNode);
    }
}

void calculate_nodal_volume_given_elem_nodes_stkmesh(const stk::mesh::Entity* elemNodes,
                                                     const unsigned numNodesPerElem,
                                                     const stk::mesh::FieldBase& coords,
                                                     stk::mesh::Field<double>& nodalVolumeField)
{
    double elemVolumePerNode = calculate_element_volume(elemNodes, numNodesPerElem, coords) / numNodesPerElem;
    for(unsigned j = 0; j < numNodesPerElem; ++j)
    {
        double* nodalVolume = stk::mesh::field_data(nodalVolumeField, elemNodes[j]);
        Kokkos::atomic_add(nodalVolume, elemVolumePerNode);
    }
}

void calculate_nodal_volume_staticmesh(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, int numRepeat)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
    ngp::StaticField<double> staticCoords(mesh, coords);
    ngp::StaticField<double> staticNodalVolume(mesh, nodalVolumeField);
    ngp::StaticMesh staticMesh(mesh);

    unsigned numBuckets = staticMesh.num_buckets(stk::topology::ELEM_RANK);
    for(int n=0; n<numRepeat; ++n)
    {
        Kokkos::parallel_for(Kokkos::TeamPolicy< ExecSpace >(numBuckets, Kokkos::AUTO), KOKKOS_LAMBDA (const ngp::TeamHandleType& team)
        {
            const int elementBucketIndex = team.league_rank();
            const ngp::StaticBucket &bucket = staticMesh.get_bucket(stk::topology::ELEM_RANK, elementBucketIndex);
            unsigned numElements = bucket.size();
    
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i)
            {
                ngp::ConnectedNodesType elemNodes = bucket.get_nodes(i);
                calculate_nodal_volume_given_elem_nodes(staticMesh, elemNodes, staticCoords, staticNodalVolume);
            });
        });
    }

    staticNodalVolume.copy_device_to_host(mesh, nodalVolumeField);
    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

void calculate_nodal_volume_static_mesh_entity_loop(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, int numRepeat)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
    ngp::StaticField<double> staticCoords(mesh, coords);
    ngp::StaticField<double> staticNodalVolume(mesh, nodalVolumeField);
    ngp::StaticMesh staticMesh(mesh);

    for(int n=0; n<numRepeat; ++n)
    {
        ngp::for_each_entity_run(staticMesh, stk::topology::ELEM_RANK, KOKKOS_LAMBDA(stk::mesh::FastMeshIndex elem)
        {
            ngp::ConnectedNodesType elemNodes = staticMesh.get_nodes(elem);
            calculate_nodal_volume_given_elem_nodes(staticMesh, elemNodes, staticCoords, staticNodalVolume);
        });
    }

    staticNodalVolume.copy_device_to_host(mesh, nodalVolumeField);
    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

void calculate_nodal_volume_stkmesh(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, int numRepeat)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();

    const stk::mesh::BucketVector& elemBuckets = mesh.get_buckets(stk::topology::ELEM_RANK, mesh.mesh_meta_data().locally_owned_part());
    unsigned numBuckets = elemBuckets.size();
    for(int n=0; n<numRepeat; ++n)
    {
        Kokkos::parallel_for(Kokkos::TeamPolicy< HostExecSpace >(numBuckets, Kokkos::AUTO), [&] (const ngp::HostTeamHandleType& team)
        {
            const int elementBucketIndex = team.league_rank();
            const stk::mesh::Bucket &bucket = *elemBuckets[elementBucketIndex];
            const unsigned numNodesPerElem = bucket.topology().num_nodes();
            unsigned numElements = bucket.size();
    
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i)
            {
                const stk::mesh::Entity* elemNodes = bucket.begin_nodes(i);
                calculate_nodal_volume_given_elem_nodes_stkmesh(elemNodes, numNodesPerElem, coords, nodalVolumeField);
            });
        });
    }

    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

void calculate_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, int numRepeat)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();

    for(int n=0; n<numRepeat; ++n)
    {
        ngp::for_each_entity_run_stkmesh(mesh, stk::topology::ELEM_RANK, [&](stk::mesh::FastMeshIndex elem)
        {
            const stk::mesh::Bucket *bucket = mesh.buckets(stk::topology::ELEM_RANK)[elem.bucket_id];
            const stk::mesh::Entity* elemNodes = bucket->begin_nodes(elem.bucket_ord);
            calculate_nodal_volume_given_elem_nodes_stkmesh(elemNodes, bucket->topology().num_nodes(), coords, nodalVolumeField);
        });
    }

    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

void expect_nodal_volume(const stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, stk::mesh::Field<int>& numElemsPerNodeField, int numRepeat)
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

            EXPECT_NEAR(numRepeat*expectedVolume, nodalVolume[i], 1.e-9);
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
    void create_mesh_and_fields(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        int dim = unitTestUtils::get_command_line_option<int>("-dim", "20");
        nodalVolumeField = &get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "nodal_volume");
        stk::mesh::put_field(*nodalVolumeField, get_meta().universal_part());
        numElemsPerNodeField = &get_meta().declare_field<stk::mesh::Field<int> >(stk::topology::NODE_RANK, "numElemsPerNode");
        stk::mesh::put_field(*numElemsPerNodeField, get_meta().universal_part());
        std::string sDim = std::to_string(dim);
        std::string meshName("generated:");
        meshName += sDim+"x"+sDim+"x"+sDim;
        setup_mesh(meshName, auraOption);
        std::cout << "Using mesh: "<<meshName<<", numRepeat: "<<get_num_repeat()<<std::endl;
        count_num_elems_per_node(get_bulk(), *numElemsPerNodeField);
    }

    int get_num_repeat()
    {
        return unitTestUtils::get_command_line_option<int>("-n", "20");
    }

    void test_nodal_volume_staticmesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        calculate_nodal_volume_staticmesh(get_bulk(), *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"calculate_nodal_volume_staticmesh: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    void test_nodal_volume_staticmesh_entity_loop(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        calculate_nodal_volume_static_mesh_entity_loop(get_bulk(), *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"calculate_nodal_volume_staticmesh_entity_loop: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    void test_nodal_volume_stkmesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        calculate_nodal_volume_stkmesh(get_bulk(), *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"calculate_nodal_volume_stkmesh: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    void test_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        calculate_nodal_volume_stkmesh_entity_loop(get_bulk(), *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"calculate_nodal_volume_stkmesh_entity_loop: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    stk::mesh::Field<double> *nodalVolumeField;
    stk::mesh::Field<int> *numElemsPerNodeField;
};
TEST_F(NodalVolumeCalculator, nodalVolumeWithAuraStaticMesh) {
    test_nodal_volume_staticmesh(stk::mesh::BulkData::AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeWithoutAuraStaticMesh) {
    test_nodal_volume_staticmesh(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeWithoutAuraStaticMeshForEachEntityLoop) {
    test_nodal_volume_staticmesh_entity_loop(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeWithoutAuraStkMesh) {
    test_nodal_volume_stkmesh(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeWithoutAuraStkMeshForEachEntityLoop) {
    test_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData::NO_AUTO_AURA);
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


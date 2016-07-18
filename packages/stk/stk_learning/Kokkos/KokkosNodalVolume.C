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
#include <stk_util/util/StkVector.hpp>
#include "../../stk_mesh/stk_mesh/base/FieldParallel.hpp"

#include <limits>



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

typedef Kokkos::TeamPolicy<HostExecSpace, ngp::ScheduleType> HostTeamPolicyType;
typedef HostTeamPolicyType::member_type HostTeamHandleType;

void calculate_nodal_volume_stkmesh(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, int numRepeat)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();

    const stk::mesh::BucketVector& elemBuckets = mesh.get_buckets(stk::topology::ELEM_RANK,
                                                                  mesh.mesh_meta_data().locally_owned_part());
    unsigned numBuckets = elemBuckets.size();
    for(int n=0; n<numRepeat; ++n)
    {
        Kokkos::parallel_for(Kokkos::TeamPolicy< HostExecSpace >(numBuckets, Kokkos::AUTO),
        [&] (const HostTeamHandleType& team)
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

void calculate_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData& mesh,
                                                stk::mesh::Field<double>& nodalVolumeField,
                                                int numRepeat)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
    for(int n=0; n<numRepeat; ++n)
    {
        ngp::for_each_entity_run_stkmesh(mesh, stk::topology::ELEM_RANK, [&](stk::mesh::FastMeshIndex elem)
        {
            const stk::mesh::Bucket *bucket = mesh.buckets(stk::topology::ELEM_RANK)[elem.bucket_id];
            const stk::mesh::Entity* elemNodes = bucket->begin_nodes(elem.bucket_ord);
            const unsigned numNodesPerElem = bucket->topology().num_nodes();
            calculate_nodal_volume_given_elem_nodes_stkmesh(elemNodes, numNodesPerElem, coords, nodalVolumeField);
        });
    }
    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}





template <typename Mesh, typename Field> STK_FUNCTION
void calculate_nodal_volume_given_elem_volume(const Mesh &ngpMesh,
                                             typename Mesh::ConnectedNodes nodes,
                                             double elemVolumePerNode,
                                             const Field& nodalVolume)
{
    const unsigned numNodesPerElem = nodes.size();
    for(unsigned j = 0; j < numNodesPerElem; ++j)
    {
        double* v = &nodalVolume.get(ngpMesh, nodes[j], 0);
        Kokkos::atomic_add(v, elemVolumePerNode);
    }
}

template <typename Mesh, typename Field> STK_FUNCTION
double calculate_element_volume(const Mesh &ngpMesh,
                                typename Mesh::ConnectedNodes nodes,
                                const Field &coords)
{
    double min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double max[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
    const unsigned numElemNodes = nodes.size();
    for(unsigned i=0; i<numElemNodes; ++i) {
        stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[i]);
        for(int j=0; j<3; ++j) {
            double val = coords.const_get(nodeIndex, j);
            if (val > max[j])
                max[j] = val;
            if (val < min[j])
                min[j] = val;
        }
    }
    return (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
}

template <typename Mesh, typename Field> STK_FUNCTION
void calculate_nodal_volume_device(const Mesh &ngpMesh,
                                   typename Mesh::MeshIndex elem,
                                   const Field &coords,
                                   const Field &nodalVolume)
{
    typename Mesh::ConnectedNodes nodes = ngpMesh.get_nodes(elem);
    const unsigned numNodesPerElem = nodes.size();
    double elemVolumePerNode = calculate_element_volume(ngpMesh, nodes, coords) / numNodesPerElem;
    calculate_nodal_volume_given_elem_volume(ngpMesh, nodes, elemVolumePerNode, nodalVolume);
}

template <typename Mesh, typename Field>
void calculate_nodal_volume(Mesh &ngpMesh, const Field &coords, Field &nodalVolume)
{
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, KOKKOS_LAMBDA(typename Mesh::MeshIndex elem)
    {
        calculate_nodal_volume_device(ngpMesh, elem, coords, nodalVolume);
    });
}

void calculate_nodal_volume_entity_loop(stk::mesh::BulkData& mesh, stk::mesh::Field<double>& nodalVolumeField, int numRepeat)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
    ngp::StkNgpField ngpCoords(mesh, coords);
    ngp::StkNgpField ngpNodalVolume(mesh, nodalVolumeField);
    ngp::StkNgpMesh ngpMesh(mesh);

    for(int n=0; n<numRepeat; ++n)
        calculate_nodal_volume(ngpMesh, ngpCoords, ngpNodalVolume);

    ngpNodalVolume.copy_device_to_host(mesh, nodalVolumeField);
    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}





template <typename Algorithm>
void repeat_for_each_entity_loop_for_algorithm(stk::mesh::BulkData& mesh,
                                               stk::mesh::Field<double>& elemVolumePerNodeField,
                                               stk::mesh::Field<double>& nodalVolumeField,
                                               int numRepeat,
                                               const Algorithm &algorithm)
{
    ngp::StkNgpField staticNodalVolume(mesh, nodalVolumeField);
    ngp::StkNgpField staticElemVolume(mesh, elemVolumePerNodeField);
    ngp::StkNgpMesh ngpMesh(mesh);

    for(int n=0; n<numRepeat; ++n)
        algorithm(ngpMesh, staticElemVolume, staticNodalVolume);

    staticNodalVolume.copy_device_to_host(mesh, nodalVolumeField);
    stk::mesh::parallel_sum(mesh, {&nodalVolumeField});
}

template <typename Mesh, typename Field> STK_FUNCTION
void assemble_nodal_volume_from_elem_volume_device(const Mesh &ngpMesh,
                                   typename Mesh::MeshIndex elem,
                                   const double elemVolumePerNode,
                                   const Field &staticNodalVolume)
{
    typename Mesh::ConnectedNodes nodes = ngpMesh.get_nodes(elem);
    calculate_nodal_volume_given_elem_volume(ngpMesh, nodes, elemVolumePerNode, staticNodalVolume);
}

template <typename Mesh, typename Field>
void assemble_nodal_volume_from_elem_volume_bucket_field_access(Mesh &ngpMesh,
                                                                const Field &elemVolumePerNode,
                                                                Field &staticNodalVolume)
{
    typedef typename Kokkos::TeamPolicy<typename Mesh::MeshExecSpace, ngp::ScheduleType>::member_type TeamHandleType;
    unsigned numBuckets = ngpMesh.num_buckets(stk::topology::ELEM_RANK);
    Kokkos::parallel_for(Kokkos::TeamPolicy<typename Mesh::MeshExecSpace>(numBuckets, Kokkos::AUTO), KOKKOS_LAMBDA (const TeamHandleType& team)
    {
        const int bucketIndex = team.league_rank();
        const typename Mesh::BucketType &bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
        unsigned numElements = bucket.size();
        double *volumePerNodePtr = &elemVolumePerNode.get(stk::mesh::FastMeshIndex{bucket.bucket_id(), 0}, 0);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& i)
        {
            typename Mesh::MeshIndex elem{&bucket, static_cast<unsigned>(i)};
            assemble_nodal_volume_from_elem_volume_device(ngpMesh, elem, volumePerNodePtr[i], staticNodalVolume);
        });
    });
}

void assemble_nodal_volume_bucket_field_access(stk::mesh::BulkData& mesh,
                                               stk::mesh::Field<double>& elemVolumePerNodeField,
                                               stk::mesh::Field<double>& nodalVolumeField,
                                               int numRepeat)
{
    repeat_for_each_entity_loop_for_algorithm(mesh, elemVolumePerNodeField, nodalVolumeField, numRepeat,
        [=](ngp::StkNgpMesh &ngpMesh, const ngp::StkNgpField &elemVolumePerNode, ngp::StkNgpField &staticNodalVolume)
        {
            assemble_nodal_volume_from_elem_volume_bucket_field_access(ngpMesh, elemVolumePerNode, staticNodalVolume);
        });
}

template <typename Mesh, typename Field>
void assemble_nodal_volume_from_elem_volume_entity_field_access(Mesh &ngpMesh,
                                                                const Field &elemVolumePerNode,
                                                                Field &staticNodalVolume)
{
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, KOKKOS_LAMBDA(typename Mesh::MeshIndex elem)
    {
        double volumePerNode = elemVolumePerNode.get(elem, 0);
        assemble_nodal_volume_from_elem_volume_device(ngpMesh, elem, volumePerNode, staticNodalVolume);
    });
}

void assemble_nodal_volume_entity_field_access(stk::mesh::BulkData& mesh,
                                               stk::mesh::Field<double>& elemVolumePerNodeField,
                                               stk::mesh::Field<double>& nodalVolumeField,
                                               int numRepeat)
{
    repeat_for_each_entity_loop_for_algorithm(mesh, elemVolumePerNodeField, nodalVolumeField, numRepeat,
        [=](ngp::StkNgpMesh &ngpMesh, const ngp::StkNgpField &elemVolumePerNode, ngp::StkNgpField &staticNodalVolume)
        {
            assemble_nodal_volume_from_elem_volume_entity_field_access(ngpMesh, elemVolumePerNode, staticNodalVolume);
        });
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

void fill_elem_volume_per_node(const stk::mesh::BulkData& mesh, stk::mesh::Field<double>& elemVolumePerNodeField)
{
    const stk::mesh::FieldBase& coords = *mesh.mesh_meta_data().coordinate_field();
    for(const stk::mesh::Bucket *bucket : mesh.buckets(stk::topology::ELEM_RANK))
    {
        unsigned numNodesPerElem = bucket->topology().num_nodes();
        double *volume = stk::mesh::field_data(elemVolumePerNodeField, *bucket);
        for(size_t j=0; j<bucket->size(); j++)
        {
            const stk::mesh::Entity *nodes = bucket->begin_nodes(j);
            volume[j] = calculate_element_volume(nodes, bucket->topology().num_nodes(), coords) / numNodesPerElem;
        }
    }
}

class NodalVolumeCalculator : public stk::unit_test_util::MeshFixture
{
protected:
    void create_mesh_and_fields(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        elemVolumePerNodeField = &get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::ELEM_RANK, "elemVolumePerNodeField");
        stk::mesh::put_field(*elemVolumePerNodeField, get_meta().universal_part());
        nodalVolumeField = &get_meta().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "nodal_volume");
        stk::mesh::put_field(*nodalVolumeField, get_meta().universal_part());
        numElemsPerNodeField = &get_meta().declare_field<stk::mesh::Field<int> >(stk::topology::NODE_RANK, "numElemsPerNode");
        stk::mesh::put_field(*numElemsPerNodeField, get_meta().universal_part());
        int dim = unitTestUtils::get_command_line_option<int>("-dim", "20");
        std::string sDim = std::to_string(dim);
        std::string meshName("generated:");
        meshName += sDim+"x"+sDim+"x"+sDim;
        setup_mesh(meshName, auraOption);
        std::cout << "Using mesh: "<<meshName<<", numRepeat: "<<get_num_repeat()<<std::endl;
        count_num_elems_per_node(get_bulk(), *numElemsPerNodeField);
        fill_elem_volume_per_node(get_bulk(), *elemVolumePerNodeField);
    }

    int get_num_repeat()
    {
        return unitTestUtils::get_command_line_option<int>("-n", "20");
    }

    void test_nodal_volume_stkmesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        calculate_nodal_volume_stkmesh(get_bulk(), *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"\t\tstkmesh: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    void test_nodal_volume_bucket_field_access(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        assemble_nodal_volume_bucket_field_access(get_bulk(), *elemVolumePerNodeField, *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"\t\tbucket field access: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    void test_nodal_volume_entity_field_access(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        assemble_nodal_volume_entity_field_access(get_bulk(), *elemVolumePerNodeField, *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"\t\tentity field access: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    void test_nodal_volume_entity_loop(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        calculate_nodal_volume_entity_loop(get_bulk(), *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"\t\tentity_loop: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    void test_nodal_volume_stkmesh_entity_loop(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        create_mesh_and_fields(auraOption);
        int numRepeat = get_num_repeat();
        double startTime = stk::wall_time();
        calculate_nodal_volume_stkmesh_entity_loop(get_bulk(), *nodalVolumeField, numRepeat);
        double elapsedTime = stk::wall_time() - startTime;
        std::cout<<"\t\tstkmesh_entity_loop: "<<elapsedTime<<" seconds"<<std::endl;
        expect_nodal_volume(get_bulk(), *nodalVolumeField, *numElemsPerNodeField, numRepeat);
    }

    stk::mesh::Field<double> *elemVolumePerNodeField;
    stk::mesh::Field<double> *nodalVolumeField;
    stk::mesh::Field<int> *numElemsPerNodeField;
};
TEST_F(NodalVolumeCalculator, primerTest) {
    test_nodal_volume_stkmesh(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, stkMeshNodalVolume) {
    test_nodal_volume_stkmesh(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeBucketFieldAccess) {
    test_nodal_volume_bucket_field_access(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeEntityFieldAccess) {
    test_nodal_volume_entity_field_access(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, nodalVolumeForEachEntityLoop) {
    test_nodal_volume_entity_loop(stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(NodalVolumeCalculator, stkMeshNodalVolumeForEachEntityLoop) {
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



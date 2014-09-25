// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>

#include <stk_util/stk_config.h>

// restrict this file to only build if KokkosCore is enabled and if STKIO is enabled
#if defined(HAVE_STK_KokkosCore) && defined(HAVE_STKIO)

#include <Kokkos_Threads.hpp>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_View.hpp>

#include <iostream>

namespace {

struct ForceFunctor
{
    typedef typename Kokkos::Threads device_type;

    ForceFunctor(const double* disp, const double* vel, const double* acc, double *force) :
            mDisp(disp), mVel(vel), mAcc(acc), mForce(force), alpha(-1.4), beta(0.3333333), gamma(3.14159)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t rowIndex) const
    {
        mForce[rowIndex] = alpha * mDisp[rowIndex] + beta * mVel[rowIndex] + gamma * mAcc[rowIndex];
    }

private:
    const double *mDisp;
    const double *mVel;
    const double *mAcc;
    double *mForce;
    const double alpha;
    const double beta;
    const double gamma;
};

} // unnamed namespace

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldDataManager.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <mpi.h>                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <vector>
#include <string>
#include <optionParsing/getOption.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

double initial_value1[3] = {-1, 2, -0.3};

void createNodalVectorField(stk::mesh::MetaData& meshMetaData, const std::string &field_name)
{
    stk::mesh::Field<double, stk::mesh::Cartesian3d> &field1 = meshMetaData.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, field_name);
    stk::mesh::put_field_on_entire_mesh_with_initial_value(field1, initial_value1);
}

void createNodalVectorFields(stk::mesh::MetaData& meshMetaData)
{
    createNodalVectorField(meshMetaData, "disp");
    createNodalVectorField(meshMetaData, "vel");
    createNodalVectorField(meshMetaData, "acc");
    createNodalVectorField(meshMetaData, "force");
    meshMetaData.commit();
}

void createMetaAndBulkData(stk::io::StkMeshIoBroker &exodusFileReader, stk::mesh::FieldDataManager *fieldDataManager)
{
    std::string exodusFileName = unitTestUtils::getOption("-i", "NO_FILE_SPECIFIED");
    ASSERT_NE(exodusFileName, "NO_FILE_SPECIFIED");

    exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

    std::cerr << "Starting To Read Mesh: " << exodusFileName << std::endl;

    exodusFileReader.create_input_mesh();
    stk::mesh::MetaData &stkMeshMetaData = exodusFileReader.meta_data();
    createNodalVectorFields(stkMeshMetaData);

    Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data(new stk::mesh::BulkData(stkMeshMetaData, MPI_COMM_WORLD, false, NULL, fieldDataManager));
    exodusFileReader.set_bulk_data(arg_bulk_data);
    stk::mesh::BulkData& stkMeshBulkData = *arg_bulk_data;

    bool delay_field_data_allocation = true;
    exodusFileReader.populate_mesh(delay_field_data_allocation);
    exodusFileReader.populate_field_data();

    std::cerr << "Finished Reading Mesh: " << exodusFileName << std::endl;

    stk::mesh::Selector allEntities = stkMeshMetaData.universal_part();
    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(allEntities, stkMeshBulkData, entityCounts);
    size_t numElements = entityCounts[stk::topology::ELEMENT_RANK];
    size_t numNodes = entityCounts[stk::topology::NODE_RANK];

    std::cerr << "Number of elements in " << exodusFileName << " = " << numElements << std::endl;
    std::cerr << "Number of nodes in " << exodusFileName << " = " << numNodes << std::endl;
}

//void timeFieldOperations(stk::mesh::MetaData &stkMeshMetaData, stk::mesh::BulkData &stkMeshBulkData, double alpha, double beta, double gamma)
//{
//    std::string numIterationsString = getOption("-numIter", "1");
//    const int numIterations = std::atoi(numIterationsString.c_str());
//    stk::mesh::Field<double, stk::mesh::Cartesian3d> &disp_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "disp");
//    stk::mesh::Field<double, stk::mesh::Cartesian3d> &vel_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "vel");
//    stk::mesh::Field<double, stk::mesh::Cartesian3d> &acc_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "acc");
//    stk::mesh::Field<double, stk::mesh::Cartesian3d> &force_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "force");
//
//    const stk::mesh::BucketVector& allNodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
//    const size_t numNodeBuckets = allNodeBuckets.size();
//
//    std::cerr << "Starting timer for " << numIterations << " iterations of nodal addition (vec4 = vec1 + vec2 + vec3) test." << std::endl;
//    double startTime = stk::cpu_time();
//    double startWallTime = stk::wall_time();
//
//    for(int i = 0; i < numIterations; i++)
//    {
//#if defined(_OPENMP)
//#pragma omp parallel for schedule(dynamic, 5)
//#endif
//        for(size_t bucketIndex = 0; bucketIndex < numNodeBuckets; ++bucketIndex)
//        {
//            const stk::mesh::Bucket &bucket = *allNodeBuckets[bucketIndex];
//            double *disp = stk::mesh::field_data(disp_field, bucket);
//            double *vel = stk::mesh::field_data(vel_field, bucket);
//            double *acc = stk::mesh::field_data(acc_field, bucket);
//            double *force = stk::mesh::field_data(force_field, bucket);
//            size_t bucketLoopEnd = 3 * bucket.size();
//            for(size_t nodeIndex = 0; nodeIndex < bucketLoopEnd; nodeIndex ++)
//            {
//                force[nodeIndex] = alpha * disp[nodeIndex] + beta * vel[nodeIndex] + gamma * acc[nodeIndex];
//            }
//        }
//    }
//
//    double elapsedTime = stk::cpu_time() - startTime;
//    double elapsedWallTime = stk::wall_time() - startWallTime;
//
//    std::cerr << "That took CPU time of " << elapsedTime << " seconds." << std::endl;
//    std::cerr << "That took wall time of " << elapsedWallTime << " seconds." << std::endl;
//}

void testGoldValues(stk::mesh::MetaData &stkMeshMetaData, stk::mesh::BulkData &stkMeshBulkData, double alpha, double beta, double gamma)
{
    double goldX = initial_value1[0] * (alpha + beta + gamma);
    double goldY = initial_value1[1] * (alpha + beta + gamma);
    double goldZ = initial_value1[2] * (alpha + beta + gamma);

    stk::mesh::Field<double, stk::mesh::Cartesian3d> &force_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "force");
    double tol = 1.0e-4;
    const stk::mesh::BucketVector& buckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    for(size_t bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex)
    {
        const stk::mesh::Bucket &bucket = *buckets[bucketIndex];
        double *force = stk::mesh::field_data(force_field, bucket);
        for(size_t nodeIndex = 0; nodeIndex < 3 * bucket.size(); nodeIndex += 3)
        {
            EXPECT_NEAR(goldX, force[nodeIndex+0], tol);
            EXPECT_NEAR(goldY, force[nodeIndex+1], tol);
            EXPECT_NEAR(goldZ, force[nodeIndex+2], tol);
        }
    }
}

void test1ToNSumOfNodalFields(stk::mesh::ContiguousFieldDataManager *fieldDataManager)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    createMetaAndBulkData(exodusFileReader, fieldDataManager);

    stk::mesh::MetaData &stkMeshMetaData = exodusFileReader.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = exodusFileReader.bulk_data();

    const stk::mesh::BucketVector& buckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
    std::cerr << "Number of node buckets: " << buckets.size() << std::endl;

    double alpha = -1.4;
    double beta = 0.3333333;
    double gamma = 3.14159;

    std::string numIterationsString = unitTestUtils::getOption("-numIter", "1");
    const int numIterations = std::atoi(numIterationsString.c_str());

    stk::mesh::Field<double, stk::mesh::Cartesian3d> &disp_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "disp");
    stk::mesh::Field<double, stk::mesh::Cartesian3d> &vel_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "vel");
    stk::mesh::Field<double, stk::mesh::Cartesian3d> &acc_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "acc");
    stk::mesh::Field<double, stk::mesh::Cartesian3d> &force_field = *stkMeshMetaData.get_field<stk::mesh::Field<double, stk::mesh::Cartesian3d> >(stk::topology::NODE_RANK, "force");

    const stk::mesh::BucketVector& allNodeBuckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);

    size_t numBytesForField = fieldDataManager->get_num_bytes_allocated_on_field(disp_field.mesh_meta_data_ordinal());
    stk::mesh::FieldMetaDataVector& field_meta_data_vector = const_cast<stk::mesh::FieldMetaDataVector&>(disp_field.get_meta_data_for_field());
            int numBytesPerEntity = field_meta_data_vector[0].m_bytes_per_entity;
    size_t numNodes = numBytesForField / numBytesPerEntity;
    size_t nodeLoopEnd = 3 * numNodes;

    std::cerr << "Starting timer for " << numIterations << " iterations of nodal addition (vec4 = vec1 + vec2 + vec3) test." << std::endl;
    double startTime = stk::cpu_time();
    double startWallTime = stk::wall_time();

    size_t numThreads = 8;
//    const size_t chunkSize = 118720;

    for(int i = 0; i < numIterations; i++)
    {
#if defined(_OPENMP)
#pragma omp parallel firstprivate(alpha, beta, gamma)
#endif
        {
            const stk::mesh::Bucket &bucket = *allNodeBuckets[0];
            double *disp = stk::mesh::field_data(disp_field, bucket);
            double *vel = stk::mesh::field_data(vel_field, bucket);
            double *acc = stk::mesh::field_data(acc_field, bucket);
            double *force = stk::mesh::field_data(force_field, bucket);
#if defined(_OPENMP)
#pragma omp for
#endif
            Kokkos::Threads::initialize( numThreads );
            ForceFunctor forceF(disp, vel, acc, force);
            Kokkos::parallel_for(nodeLoopEnd, forceF);
            Kokkos::Threads::finalize();
#if defined(_OPENMP)
            numThreads = omp_get_num_threads();
#endif
        }
    }
    double elapsedTime = stk::cpu_time() - startTime;
    double elapsedWallTime = stk::wall_time() - startWallTime;

    std::cerr << "That took CPU time of " << elapsedTime << " seconds." << std::endl;
    std::cerr << "That took wall time of " << elapsedWallTime << " seconds." << std::endl;
    std::cerr << "Using " << numThreads << " threads" << std::endl;

    testGoldValues(stkMeshMetaData, stkMeshBulkData, alpha, beta, gamma);
}

//void testSumOfNodalFields(stk::mesh::FieldDataManager *fieldDataManager)
//{
//    MPI_Comm communicator = MPI_COMM_WORLD;
//    stk::io::StkMeshIoBroker exodusFileReader(communicator);
//
//    createMetaAndBulkData(exodusFileReader, fieldDataManager);
//
//    stk::mesh::MetaData &stkMeshMetaData = exodusFileReader.meta_data();
//    stk::mesh::BulkData &stkMeshBulkData = exodusFileReader.bulk_data();
//
//    const stk::mesh::BucketVector& buckets = stkMeshBulkData.buckets(stk::topology::NODE_RANK);
//    std::cerr << "Number of node buckets: " << buckets.size() << std::endl;
//
//    double alpha = -1.4;
//    double beta = 0.3333333;
//    double gamma = 3.14159;
//
//    timeFieldOperations(stkMeshMetaData, stkMeshBulkData, alpha, beta, gamma);
//
//    testGoldValues(stkMeshMetaData, stkMeshBulkData, alpha, beta, gamma);
//}

//TEST(NodalFieldPerformance, addNodalFieldsDefaultFieldManager)
//{
//    const int weKnowThereAreFiveRanks = 5;
//    stk::mesh::DefaultFieldDataManager fieldDataManager(weKnowThereAreFiveRanks);
//    testSumOfNodalFields(&fieldDataManager);
//}
//
//TEST(NodalFieldPerformance, addNodalFieldsContiguousFieldManager)
//{
//    stk::mesh::ContiguousFieldDataManager fieldDataManager;
//    testSumOfNodalFields(&fieldDataManager);
//}

TEST(NodalFieldPerformance, addNodalFields1ToNOrderingContiguousFieldManager)
{
    stk::mesh::ContiguousFieldDataManager fieldDataManager;
    test1ToNSumOfNodalFields(&fieldDataManager);
}

#endif


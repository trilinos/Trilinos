// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <stk_mesh/base/FieldBLAS.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <string>
#include <limits>

#if defined(_OPENMP) && !defined(__INTEL_COMPILER)
#define OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
// there seems to be an issue with OpenMP combined with GoogleTest macros with Intel compilers
// example of error:
//    openMP.C(206): internal error: assertion failed at: "shared/cfe/edgcpfe/checkdir.c", line 5531
#include <omp.h>
#endif

//#ifndef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//TEST(OMP,single_failing) {
//    std::cout << "!defined(_OPENMP) || !defined(__INTEL_COMPILER)" <<std::endl;
//    EXPECT_TRUE(false);
//}
//
//TEST(OMP,all_failing) {
//    std::cout << "!defined(_OPENMP) || !defined(__INTEL_COMPILER)" <<std::endl;
//    EXPECT_TRUE(false);
//}
//#endif
//
//template<class A>
//struct BLASFixture {
//
//    A initial_value1;
//    A initial_value2;
//    int numEntities;
//    int numBuckets;
//    stk::io::StkMeshIoBroker * stkMeshIoBroker;
//    stk::mesh::BulkData * stkMeshBulkData;
//    stk::mesh::Field<A> * field1;
//    stk::mesh::FieldBase * fieldBase1;
//    stk::mesh::Field<A> * field2;
//    stk::mesh::FieldBase * fieldBase2;
//
//    BLASFixture(const A init1 ,const A init2 = A(), const int MeshSize = 64);
//    ~BLASFixture();
//
//};
//
//template<class A>
//BLASFixture<A>::BLASFixture(const A init1,const A init2, const int MeshSize)
//{
//    MPI_Comm communicator = MPI_COMM_WORLD;
//
//    initial_value1 = init1;
//    initial_value2 = init2;
//
//    stkMeshIoBroker = new stk::io::StkMeshIoBroker(communicator);
//    stk::io::StkMeshIoBroker & io = *stkMeshIoBroker;
//    char generatedMeshSpecification [25];
//    sprintf(generatedMeshSpecification,"generated:%dx%dx%d",MeshSize,MeshSize,MeshSize);
//    io.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
//    io.create_input_mesh();
//    stk::mesh::MetaData &meta_data = io.meta_data();
//
//    field1 = &meta_data.declare_field<stk::mesh::Field<A> >(stk::topology::NODE_RANK, "field1",3);
//    stk::mesh::put_field(*field1,field1->mesh_meta_data().universal_part(),&initial_value1);
//    fieldBase1 = dynamic_cast<stk::mesh::FieldBase*>(field1);
//
//    field2 = &meta_data.declare_field<stk::mesh::Field<A> >(stk::topology::NODE_RANK, "field2",3);
//    stk::mesh::put_field(*field2,field2->mesh_meta_data().universal_part(),&initial_value2);
//    fieldBase2 = dynamic_cast<stk::mesh::FieldBase*>(field2);
//
//    io.populate_bulk_data(); // THIS IS THE SLOW LINE
//    stkMeshBulkData = &io.bulk_data();
//
//    const stk::mesh::Selector selector = meta_data.universal_part() & stk::mesh::selectField(*field1);
//    const stk::mesh::BucketVector & buckets = stkMeshBulkData->get_buckets(field1->entity_rank(),selector);
//
//    numBuckets = buckets.size();
//    numEntities = 0;
//    for (int i=0; i<numBuckets;++i){
//        numEntities += (buckets[i]->size());
//    }
//}
//
//template<class A>BLASFixture<A>::~BLASFixture() {
//    delete stkMeshIoBroker;
//}
//
//template<class Scalar>
//void testFieldValidation(BLASFixture<Scalar> & Fixture,Scalar val1,Scalar val2,double tol=1.0e-3) {
//    /*const stk::mesh::BucketVector& buckets = Fixture.stkMeshBulkData->buckets(stk::topology::NODE_RANK);
//    for(size_t j = 0; j < buckets.size(); j++)
//    {
//        const stk::mesh::Bucket& bucket = *buckets[j];
//        for(size_t i=0; i<bucket.size(); ++i)*/
//    const stk::mesh::BucketVector& buckets = Fixture.stkMeshBulkData->buckets(stk::topology::NODE_RANK);
//    for(size_t j = 0; j < 1u; j++)
//    {
//        const stk::mesh::Bucket& bucket = *buckets[j];
//        for(size_t i=0; i < 1u; i++)
//        {
//            Scalar* field_value1 = reinterpret_cast<Scalar *>(stk::mesh::field_data(*Fixture.field1, bucket[i]));
//            EXPECT_NEAR(val1,*field_value1,tol);
//
//            Scalar* field_value2 = reinterpret_cast<Scalar *>(stk::mesh::field_data(*Fixture.field2, bucket[i]));
//            EXPECT_NEAR(val2,*field_value2,tol);
//        }
//
//    }
//}
//
//void axpyCompare(BLASFixture<double> * Fixture,const int numThreads,const int repTestMAX)
//{
//    const double alpha = -1.15;
//    ThrowAssert(alpha*repTestMAX*Fixture.initial_value1+Fixture.initial_value2<std::numeric_limits<double>::max());
//
//    //test init
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//
//    wallTime_field -= stk::wall_time();
//    for (int i=0;i<repTestMAX;i++)
//    {
//        stk::mesh::field_axpy(alpha,*Fixture->field1,*Fixture->field2);
//    }
//    wallTime_field += stk::wall_time();
//
//    testFieldValidation(*Fixture,Fixture->initial_value1,alpha*repTestMAX*Fixture->initial_value1+Fixture->initial_value2);
//
//    //return to init
//    stk::mesh::field_axpy(-repTestMAX*alpha,*Fixture->field1,*Fixture->field2);
//
//    wallTime_fieldBase -= stk::wall_time();
//    for (int i=0;i<repTestMAX;i++)
//    {
//        stk::mesh::field_axpy(alpha,*Fixture->fieldBase1,*Fixture->fieldBase2);
//    }
//    wallTime_fieldBase += stk::wall_time();
//
//    //test post
//    testFieldValidation(*Fixture,Fixture->initial_value1,alpha*repTestMAX*Fixture->initial_value1+Fixture->initial_value2);
//
//    //return to init
//    stk::mesh::field_axpy(-repTestMAX*alpha,*Fixture->field1,*Fixture->field2);
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_axpyCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            axpyCompare(&Fixture,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void scaleCompare(BLASFixture<double> * Fixture,const int numThreads,const int repTestMAX)
//{
//    const double alpha1 = 1.1;
//    const double alpha2 = 0.92;
//
//    const int inner_repeat = 5;
//
//    //test init
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//
//    for (int j=0;j<repTestMAX;j++)
//    {
//        wallTime_field -= stk::wall_time();
//        for (int i=0;i<inner_repeat;i++)
//        {
//            stk::mesh::field_scale(alpha1,*Fixture->field1);
//            stk::mesh::field_scale(alpha2,*Fixture->field1);
//        }
//        wallTime_field += stk::wall_time();
//
//        testFieldValidation(*Fixture,pow(alpha1*alpha2,inner_repeat)*double(Fixture->initial_value1),Fixture->initial_value2);
//
//        stk::mesh::field_scale(pow(alpha1*alpha2,-inner_repeat),*Fixture->field1);
//
//        wallTime_fieldBase -= stk::wall_time();
//        for (int rep = 0; rep<inner_repeat;rep++)
//        {
//            stk::mesh::field_scale(alpha1,*Fixture->fieldBase1);
//            stk::mesh::field_scale(alpha2,*Fixture->fieldBase1);
//        }
//        wallTime_fieldBase += stk::wall_time();
//
//        testFieldValidation(*Fixture,pow(alpha1*alpha2,inner_repeat)*double(Fixture->initial_value1),Fixture->initial_value2);
//
//        stk::mesh::field_scale(pow(alpha1*alpha2,-inner_repeat),*Fixture->field1);
//    }
//
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_scaleCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            scaleCompare(&Fixture,numThread_list[j],repTestMAX/3);
//        }
//    }
//}
//
//void copyCompare(BLASFixture<double> * Fixture1,BLASFixture<double> * Fixture2,const int numThreads,const int repTestMAX)
//{
//    //test init
//    testFieldValidation(*Fixture1,Fixture1->initial_value1,Fixture1->initial_value2);
//    testFieldValidation(*Fixture2,Fixture2->initial_value1,Fixture2->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        stk::mesh::field_copy(*Fixture1->field1,*Fixture1->field2);
//        stk::mesh::field_copy(*Fixture2->field2,*Fixture2->field1);
//        wallTime_field += stk::wall_time();
//
//        testFieldValidation(*Fixture1,Fixture1->initial_value1,Fixture1->initial_value1);
//        testFieldValidation(*Fixture2,Fixture2->initial_value2,Fixture2->initial_value2);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_copy(*Fixture1->fieldBase1,*Fixture2->fieldBase1);
//        stk::mesh::field_copy(*Fixture2->fieldBase2,*Fixture1->fieldBase2);
//        wallTime_fieldBase += stk::wall_time();
//
//        testFieldValidation(*Fixture1,Fixture1->initial_value1,Fixture1->initial_value2);
//        testFieldValidation(*Fixture2,Fixture2->initial_value1,Fixture2->initial_value2);
//    }
//
//    //cerr
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture1->numEntities,Fixture1->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture1->numEntities,Fixture1->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_copyCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture1 (initial_value1,initial_value2,meshSize_list[i]);
//        BLASFixture<double> Fixture2 (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            copyCompare(&Fixture1,&Fixture2,numThread_list[j],repTestMAX/3);
//        }
//    }
//}
//
//void dotCompare(BLASFixture<double> * Fixture,const int numThreads,const int repTestMAX)
//{
//    ThrowAssert(std::abs(initial_value1*initial_value2)*Fixture.numEntities<std::numeric_limits<double>::max());
//
//    //test init
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    const double tol = 1.0e-3;
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//    double result_field = 0.0;
//    double result_fieldBase = 0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        result_field = stk::mesh::field_dot(*Fixture->field1,*Fixture->field2);
//        wallTime_field += stk::wall_time();
//
//        EXPECT_NEAR(Fixture->initial_value1*Fixture->initial_value2*Fixture->numEntities,result_field,tol);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_dot(result_fieldBase,*Fixture->fieldBase1,*Fixture->fieldBase2);
//        wallTime_fieldBase += stk::wall_time();
//
//        EXPECT_NEAR(Fixture->initial_value1*Fixture->initial_value2*Fixture->numEntities,result_fieldBase,tol);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_dotCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            dotCompare(&Fixture,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void fillCompare(BLASFixture<double> * Fixture,const int numThreads,const int repTestMAX)
//{
//    const double other_value1 = pow(Fixture->initial_value1,2.0)+1;
//    const double other_value2 = -pow(Fixture->initial_value2,4.0)-1;
//
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        stk::mesh::field_fill(other_value1,*Fixture->field1);
//        stk::mesh::field_fill(other_value2,*Fixture->field2);
//        wallTime_field += stk::wall_time();
//
//        testFieldValidation(*Fixture,other_value1,other_value2);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_fill(Fixture->initial_value1,*Fixture->fieldBase1);
//        stk::mesh::field_fill(Fixture->initial_value2,*Fixture->fieldBase2);
//        wallTime_fieldBase += stk::wall_time();
//
//        testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_fillCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            fillCompare(&Fixture,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void swapCompare(BLASFixture<double> * Fixture,const int numThreads,const int repTestMAX)
//{
//    //test init
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        stk::mesh::field_swap(*Fixture->field1,*Fixture->field2);
//        wallTime_field += stk::wall_time();
//
//        testFieldValidation(*Fixture,Fixture->initial_value2,Fixture->initial_value1);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_swap(*Fixture->fieldBase1,*Fixture->fieldBase2);
//        wallTime_fieldBase += stk::wall_time();
//
//        testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_swapCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            swapCompare(&Fixture,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void nrm2Compare(BLASFixture<double> * Fixture,const int numThreads,const int repTestMAX)
//{
//    EXPECT_LT(std::abs(Fixture->initial_value1*Fixture->initial_value2)*sqrt(double(Fixture->numEntities)),std::numeric_limits<double>::max());
//
//    //test init
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    const double tol = 1.0e-3;
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//    double result_field;
//    double result_fieldBase;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        result_field = stk::mesh::field_nrm2(*Fixture->field1);
//        wallTime_field += stk::wall_time();
//
//        EXPECT_NEAR(std::abs(Fixture->initial_value1)*sqrt(double(Fixture->numEntities)),result_field,tol);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_nrm2(result_fieldBase,*Fixture->fieldBase1);
//        wallTime_fieldBase += stk::wall_time();
//
//        EXPECT_NEAR(std::abs(Fixture->initial_value1)*sqrt(double(Fixture->numEntities)),result_fieldBase,tol);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_nrm2Compare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            nrm2Compare(&Fixture,numThread_list[j],repTestMAX*3);
//        }
//    }
//}
//
//void asumCompare(BLASFixture<double> * Fixture,const int numThreads,const int repTestMAX)
//{
//    ThrowAssert(std::abs(initial_value1*initial_value2)*sqrt(double(Fixture.numEntities))<std::numeric_limits<double>::max());
//
//    //test init
//    testFieldValidation(*Fixture,Fixture->initial_value1,Fixture->initial_value2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    const double tol = 1.0e-3;
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//    double result_field = 0.0;
//    double result_fieldBase = 0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        result_field=stk::mesh::field_asum(*Fixture->field1);
//        wallTime_field += stk::wall_time();
//
//        EXPECT_NEAR(std::abs(Fixture->initial_value1)*Fixture->numEntities,result_field,tol);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_asum(result_fieldBase,*Fixture->fieldBase1);
//        wallTime_fieldBase += stk::wall_time();
//
//        EXPECT_NEAR(std::abs(Fixture->initial_value1)*Fixture->numEntities,result_fieldBase,tol);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//
//TEST(OMP,single_asumCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double initial_value1 = 2.24;
//    const double initial_value2 = 4.4;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            asumCompare(&Fixture,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void eamaxCompare(BLASFixture<double> * Fixture,const double low_value,const double high_value,const int numThreads,const int repTestMAX)
//{
//    EXPECT_TRUE(std::abs(low_value)<std::abs(high_value));
//    EXPECT_TRUE(high_value<low_value);
//
//    //reset
//    stk::mesh::field_fill(low_value,*Fixture->field1);
//    stk::mesh::field_fill(low_value,*Fixture->field2);
//
//    const stk::mesh::MetaData &metaData = Fixture->stkMeshBulkData->mesh_meta_data();
//    const stk::mesh::Selector selector = metaData.universal_part() & stk::mesh::selectField(*Fixture->field1);
//    const stk::mesh::BucketVector buckets = Fixture->stkMeshBulkData->get_buckets(Fixture->field1->entity_rank(),selector);
//    int iter_inc = int(floor((buckets.size()-1.0)/10.0));
//    if (iter_inc<1) iter_inc=1;
//    for (int buckets_i=buckets.size()-1; buckets_i>0; buckets_i-=iter_inc)
//    {
//        stk::mesh::Bucket & b = *buckets[buckets_i];
//        double * x = (double*)stk::mesh::field_data(*Fixture->field1, b);
//        x[int(floor((b.size()-1)*(buckets_i)/double(buckets.size()-1)))]=low_value+(high_value-low_value)*(buckets_i)/double(buckets.size()-1);
//    }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//    stk::mesh::Entity entityOfamax_field;
//    stk::mesh::Entity entityOfamax_fieldBase;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        entityOfamax_field = stk::mesh::field_eamax(*Fixture->field1);
//        wallTime_field += stk::wall_time();
//
//        EXPECT_EQ(*stk::mesh::field_data(*Fixture->field1,entityOfamax_field),high_value);
//
//        wallTime_fieldBase -= stk::wall_time();
//        entityOfamax_fieldBase = stk::mesh::field_eamax(*Fixture->fieldBase1);
//        wallTime_fieldBase += stk::wall_time();
//
//        EXPECT_EQ(*stk::mesh::field_data(*Fixture->field1,entityOfamax_fieldBase),high_value);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//
//    //reset
//    stk::mesh::field_fill(Fixture->initial_value1,*Fixture->field1);
//    stk::mesh::field_fill(Fixture->initial_value2,*Fixture->field2);
//}
//
//
//TEST(OMP,single_eamaxCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double low_value = -5.2;
//    const double high_value = -10.7;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (low_value,low_value,meshSize_list[i]);
//
//        for (int j=0;j<num_numThreads;j++)
//        {
//            eamaxCompare(&Fixture,low_value,high_value,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void eaminCompare(BLASFixture<double> * Fixture,const double low_value,const double high_value,const int numThreads,const int repTestMAX)
//{
//    EXPECT_TRUE(std::abs(low_value)<std::abs(high_value));
//    EXPECT_TRUE(high_value<low_value);
//
//    //reset
//    stk::mesh::field_fill(high_value,*Fixture->field1);
//    stk::mesh::field_fill(high_value,*Fixture->field2);
//
//    const stk::mesh::MetaData &metaData = Fixture->stkMeshBulkData->mesh_meta_data();
//    const stk::mesh::Selector selector = metaData.universal_part() & stk::mesh::selectField(*Fixture->field1);
//    const stk::mesh::BucketVector buckets = Fixture->stkMeshBulkData->get_buckets(Fixture->field1->entity_rank(),selector);
//    int iter_inc = int(floor((buckets.size()-1.0)/10.0));
//    if (iter_inc<1) iter_inc=1;
//    for (int buckets_i=buckets.size()-1; buckets_i>0; buckets_i-=iter_inc)
//    {
//        stk::mesh::Bucket & b = *buckets[buckets_i];
//        double * x = (double*)stk::mesh::field_data(*Fixture->field1, b);
//        x[int(floor((b.size()-1)*(buckets_i)/double(buckets.size()-1)))]=high_value+(low_value-high_value)*(buckets_i)/double(buckets.size()-1);
//    }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//    stk::mesh::Entity entityOfamax_field;
//    stk::mesh::Entity entityOfamax_fieldBase;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        entityOfamax_field = stk::mesh::field_eamin(*Fixture->field1);
//        wallTime_field += stk::wall_time();
//
//        EXPECT_EQ(*stk::mesh::field_data(*Fixture->field1,entityOfamax_field),low_value);
//
//        wallTime_fieldBase -= stk::wall_time();
//        entityOfamax_fieldBase = stk::mesh::field_eamin(*Fixture->fieldBase1);
//        wallTime_fieldBase += stk::wall_time();
//
//        EXPECT_EQ(*stk::mesh::field_data(*Fixture->field1,entityOfamax_fieldBase),low_value);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//
//    //reset
//    stk::mesh::field_fill(Fixture->initial_value1,*Fixture->field1);
//    stk::mesh::field_fill(Fixture->initial_value2,*Fixture->field2);
//}
//
//
//TEST(OMP,single_eaminCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double low_value = -5.2;
//    const double high_value = -10.7;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (high_value,high_value,meshSize_list[i]);
//
//        for (int j=0;j<num_numThreads;j++)
//        {
//            eaminCompare(&Fixture,low_value,high_value,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void productCompare(BLASFixture<double> * Fixture1,BLASFixture<double> * Fixture2,const double alpha1, const double alpha2,const int numThreads,const int repTestMAX)
//{
//    EXPECT_NEAR(alpha1,1.0,0.15);
//    EXPECT_NEAR(alpha2,1.0,0.15);
//    EXPECT_LT(std::abs(pow(alpha1*alpha2,double(repTestMAX)))*std::abs(double(Fixture2->initial_value1)),std::numeric_limits<double>::max());
//
//    //test init
//    stk::mesh::field_fill(alpha1,*Fixture1->field1);
//    stk::mesh::field_fill(alpha2,*Fixture1->field2);
//    stk::mesh::field_fill(Fixture2->initial_value1,*Fixture2->field1);
//    stk::mesh::field_fill(Fixture2->initial_value2,*Fixture2->field2);
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        stk::mesh::field_product(*Fixture1->field1,*Fixture2->field1,*Fixture2->field1);
//        stk::mesh::field_product(*Fixture1->field2,*Fixture2->field1,*Fixture2->field1);
//        wallTime_field += stk::wall_time();
//    }
//
//    testFieldValidation(*Fixture1,alpha1,alpha2);
//    testFieldValidation(*Fixture2,pow(alpha1*alpha2,double(repTestMAX))*Fixture2->initial_value1,Fixture2->initial_value2);
//
//    stk::mesh::field_scale(pow(alpha1*alpha2,-double(repTestMAX)),*Fixture2->field1);
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_product(*Fixture1->field1,*Fixture2->field2,*Fixture2->field2);
//        stk::mesh::field_product(*Fixture1->field2,*Fixture2->field2,*Fixture2->field2);
//        wallTime_fieldBase += stk::wall_time();
//    }
//
//    testFieldValidation(*Fixture1,alpha1,alpha2);
//    testFieldValidation(*Fixture2,Fixture2->initial_value1,pow(alpha1*alpha2,double(repTestMAX))*Fixture2->initial_value2);
//
//    stk::mesh::field_scale(pow(alpha1*alpha2,-repTestMAX),*Fixture2->field2);
//
//    //cerr
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture1->numEntities,Fixture1->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture1->numEntities,Fixture1->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//}
//
//TEST(OMP,single_productCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double alpha1 = 1.1;
//    const double alpha2 = 0.92;
//    const double initial_value1 = 13.13;
//    const double initial_value2 = 4.21;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture1 (initial_value1,initial_value2,meshSize_list[i]);
//        BLASFixture<double> Fixture2 (initial_value1,initial_value2,meshSize_list[i]);
//        for (int j=0;j<num_numThreads;j++)
//        {
//            productCompare(&Fixture1,&Fixture2,alpha1,alpha2,numThread_list[j],repTestMAX/2);
//        }
//    }
//}
//
//void amaxCompare(BLASFixture<double> * Fixture,const double low_value,const double high_value,const int numThreads,const int repTestMAX)
//{
//    EXPECT_TRUE(std::abs(low_value)<std::abs(high_value));
//    EXPECT_TRUE(high_value<low_value);
//
//    //reset
//    stk::mesh::field_fill(low_value,*Fixture->field1);
//    stk::mesh::field_fill(low_value,*Fixture->field2);
//
//    const stk::mesh::MetaData &metaData = Fixture->stkMeshBulkData->mesh_meta_data();
//    const stk::mesh::Selector selector = metaData.universal_part() & stk::mesh::selectField(*Fixture->field1);
//    const stk::mesh::BucketVector buckets = Fixture->stkMeshBulkData->get_buckets(Fixture->field1->entity_rank(),selector);
//    int iter_inc = int(floor((buckets.size()-1.0)/10.0));
//    if (iter_inc<1) iter_inc=1;
//    for (int buckets_i=buckets.size()-1; buckets_i>0; buckets_i-=iter_inc)
//    {
//        stk::mesh::Bucket & b = *buckets[buckets_i];
//        double * x = (double*)stk::mesh::field_data(*Fixture->field1, b);
//        x[int(floor((b.size()-1)*(buckets_i)/double(buckets.size()-1)))]=low_value+(high_value-low_value)*(buckets_i)/double(buckets.size()-1);
//    }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//    const double tol = 1.0e-3;
//    double amax_field=0.0;
//    double amax_fieldBase=0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        amax_field=stk::mesh::field_amax(*Fixture->field1);
//        wallTime_field += stk::wall_time();
//
//        EXPECT_NEAR(amax_field,std::abs(high_value),tol);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_amax(amax_fieldBase,*Fixture->fieldBase1);
//        wallTime_fieldBase += stk::wall_time();
//
//        EXPECT_NEAR(amax_fieldBase,std::abs(high_value),tol);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//
//    //reset
//    stk::mesh::field_fill(Fixture->initial_value1,*Fixture->field1);
//    stk::mesh::field_fill(Fixture->initial_value2,*Fixture->field2);
//}
//
//
//TEST(OMP,single_amaxCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double low_value = -5.2;
//    const double high_value = -10.7;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (low_value,low_value,meshSize_list[i]);
//
//        for (int j=0;j<num_numThreads;j++)
//        {
//            amaxCompare(&Fixture,low_value,high_value,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void aminCompare(BLASFixture<double> * Fixture,const double low_value,const double high_value,const int numThreads,const int repTestMAX)
//{
//    EXPECT_TRUE(std::abs(low_value)<std::abs(high_value));
//    EXPECT_TRUE(high_value<low_value);
//
//    //reset
//    stk::mesh::field_fill(high_value,*Fixture->field1);
//    stk::mesh::field_fill(high_value,*Fixture->field2);
//
//    const stk::mesh::MetaData &metaData = Fixture->stkMeshBulkData->mesh_meta_data();
//    const stk::mesh::Selector selector = metaData.universal_part() & stk::mesh::selectField(*Fixture->field1);
//    const stk::mesh::BucketVector buckets = Fixture->stkMeshBulkData->get_buckets(Fixture->field1->entity_rank(),selector);
//    int iter_inc = int(floor((buckets.size()-1.0)/10.0));
//    if (iter_inc<1) iter_inc=1;
//    for (int buckets_i=buckets.size()-1; buckets_i>0; buckets_i-=iter_inc)
//    {
//        stk::mesh::Bucket & b = *buckets[buckets_i];
//        double * x = (double*)stk::mesh::field_data(*Fixture->field1, b);
//        x[int(floor((b.size()-1)*(buckets_i)/double(buckets.size()-1)))]=high_value+(low_value-high_value)*(buckets_i)/double(buckets.size()-1);
//    }
//
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    omp_set_num_threads(numThreads);
//#endif
//
//    double wallTime_field = 0.0;
//    double wallTime_fieldBase = 0.0;
//    const double tol = 1.0e-3;
//    double amin_field=0.0;
//    double amin_fieldBase=0.0;
//
//    for (int i=0;i<repTestMAX;i++)
//    {
//        wallTime_field -= stk::wall_time();
//        amin_field=stk::mesh::field_amin(*Fixture->field1);
//        wallTime_field += stk::wall_time();
//
//        EXPECT_NEAR(amin_field,std::abs(low_value),tol);
//
//        wallTime_fieldBase -= stk::wall_time();
//        stk::mesh::field_amin(amin_fieldBase,*Fixture->fieldBase1);
//        wallTime_fieldBase += stk::wall_time();
//
//        EXPECT_NEAR(amin_fieldBase,std::abs(low_value),tol);
//    }
//
//    char coutBuffer [120];
//#ifdef OPEN_MP_ACTIVE_PERFOPENMPMKL_CPP
//    sprintf(coutBuffer,"%20d , %20d , %20d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,numThreads,wallTime_field,wallTime_fieldBase);
//#else
//    sprintf(coutBuffer,"%20d , %20d , %8d of %8d , %20f , %20f\n",Fixture->numEntities,Fixture->numBuckets,0,numThreads,wallTime_field,wallTime_fieldBase);
//#endif
//    std::cout << coutBuffer;
//
//    //reset
//    stk::mesh::field_fill(Fixture->initial_value1,*Fixture->field1);
//    stk::mesh::field_fill(Fixture->initial_value2,*Fixture->field2);
//}
//
//
//TEST(OMP,single_aminCompare)
//{
//    char coutBuffer [120];
//    sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//    std::cout << coutBuffer;
//
//    const double low_value = -5.2;
//    const double high_value = -10.7;
//    const int repTestMAX = 2000;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,8,16,24,32,40};
//
//    const int num_meshSize_list = 4;
//    int meshSize_list [num_meshSize_list] = {32,48,64,128};
//    for (int i=0;i<num_meshSize_list;i++)
//    {
//        BLASFixture<double> Fixture (high_value,high_value,meshSize_list[i]);
//
//        for (int j=0;j<num_numThreads;j++)
//        {
//            aminCompare(&Fixture,low_value,high_value,numThread_list[j],repTestMAX);
//        }
//    }
//}
//
//void INTERNAL_cout_describe_func(bool begin, std::string name)
//{
//    if (begin) { std::cout<<"begin "; }
//    else { std::cout<<"end "; }
//    std::cout << name << std::endl;
//    if (begin)
//    {
//        char coutBuffer [120];
//        sprintf(coutBuffer,"%20s , %20s , %20s , %20s , %20s\n","Num Entities","Num Buckets","Num Threads","Field Time","FieldBase Time");
//        std::cout << coutBuffer;
//    }
//}
//
//TEST(OMP,all) {
//    const double initial_value1 = 2.24;
//    const double initial_value2 = -4.4;
//    const double low_value = -5.2;
//    const double high_value = -10.7;
//    const double alpha1 = 1.1;
//    const double alpha2 = 0.92;
//
//    const int repTestMAX = 2000;
//
//    const int num_meshSizes = 4;
//    int meshSize_list [num_meshSizes] = {12,14,16,18};//{32,48,64,128};
//    BLASFixture<double>*fixture1s [4];
//    BLASFixture<double>*fixture2s [4];
//    BLASFixture<double> fixture1_0 (initial_value1,initial_value2,meshSize_list[0]);fixture1s[0]=&fixture1_0;
//    BLASFixture<double> fixture2_0 (initial_value1,initial_value2,meshSize_list[0]);fixture2s[0]=&fixture2_0;
//    BLASFixture<double> fixture1_1 (initial_value1,initial_value2,meshSize_list[1]);fixture1s[1]=&fixture1_1;
//    BLASFixture<double> fixture2_1 (initial_value1,initial_value2,meshSize_list[1]);fixture2s[1]=&fixture2_1;
//    BLASFixture<double> fixture1_2 (initial_value1,initial_value2,meshSize_list[2]);fixture1s[2]=&fixture1_2;
//    BLASFixture<double> fixture2_2 (initial_value1,initial_value2,meshSize_list[2]);fixture2s[2]=&fixture2_2;
//    BLASFixture<double> fixture1_3 (initial_value1,initial_value2,meshSize_list[3]);fixture1s[3]=&fixture1_3;
//    BLASFixture<double> fixture2_3 (initial_value1,initial_value2,meshSize_list[3]);fixture2s[3]=&fixture2_3;
//
//    const int num_numThreads = 6;
//    int numThread_list [num_numThreads] = {1,10,20,30,40,50};
//
//    INTERNAL_cout_describe_func(true,"axpy");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            axpyCompare(fixture1s[iter_meshSize],numThread_list[iter_numThreads],repTestMAX);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"axpy");
//
//    INTERNAL_cout_describe_func(true,"scale");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            scaleCompare(fixture1s[iter_meshSize],numThread_list[iter_numThreads],repTestMAX/2);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"scale");
//
//    INTERNAL_cout_describe_func(true,"copy");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            copyCompare(fixture1s[iter_meshSize],fixture2s[iter_meshSize],numThread_list[iter_numThreads],2*repTestMAX/3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"copy");
//
//    INTERNAL_cout_describe_func(true,"dot");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            dotCompare(fixture1s[iter_meshSize],numThread_list[iter_numThreads],repTestMAX);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"dot");
//
//    INTERNAL_cout_describe_func(true,"fill");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            fillCompare(fixture1s[iter_meshSize],numThread_list[iter_numThreads],4*repTestMAX/3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"fill");
//
//    INTERNAL_cout_describe_func(true,"swap");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            swapCompare(fixture1s[iter_meshSize],numThread_list[iter_numThreads],repTestMAX);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"swap");
//
//    INTERNAL_cout_describe_func(true,"nrm2");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            nrm2Compare(fixture1s[iter_meshSize],numThread_list[iter_numThreads],repTestMAX*3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"nrm2");
//
//    INTERNAL_cout_describe_func(true,"asum");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            asumCompare(fixture1s[iter_meshSize],numThread_list[iter_numThreads],repTestMAX*3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"asum");
//
//    INTERNAL_cout_describe_func(true,"eamax");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            eamaxCompare(fixture1s[iter_meshSize],low_value,high_value,numThread_list[iter_numThreads],repTestMAX*3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"eamax");
//
//    INTERNAL_cout_describe_func(true,"eamin");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            eaminCompare(fixture1s[iter_meshSize],low_value,high_value,numThread_list[iter_numThreads],repTestMAX*3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"eamin");
//
//    INTERNAL_cout_describe_func(true,"product");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            productCompare(fixture1s[iter_meshSize],fixture2s[iter_meshSize],alpha1,alpha2,numThread_list[iter_numThreads],repTestMAX/2);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"product");
//
//    INTERNAL_cout_describe_func(true,"amax");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            amaxCompare(fixture1s[iter_meshSize],low_value,high_value,numThread_list[iter_numThreads],repTestMAX*3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"amax");
//
//    INTERNAL_cout_describe_func(true,"amin");
//    for (int iter_meshSize=0;iter_meshSize<num_meshSizes;iter_meshSize++)
//    {
//        for (int iter_numThreads=0;iter_numThreads<num_numThreads;iter_numThreads++)
//        {
//            aminCompare(fixture1s[iter_meshSize],low_value,high_value,numThread_list[iter_numThreads],repTestMAX*3);
//        }
//    }
//    INTERNAL_cout_describe_func(false,"amin");
//}

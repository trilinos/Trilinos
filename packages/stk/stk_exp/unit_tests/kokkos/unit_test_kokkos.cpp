#include <gtest/gtest.h>

// restrict this file to only build on platforms that defined
// STK_EXP_BUILD_KOKKOS
#ifdef STK_EXP_BUILD_KOKKOS
#include <Kokkos_Threads.hpp>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_View.hpp>

#include <iostream>

namespace {

// setup and tear down for the KokkosThreads unit tests
class KokkosThreads : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    unsigned num_threads = 8;

    // if hwloc is present we will get better thread placement and
    // numa aware allocation.
    // Currently sierra does not have the hwloc TPL enabled
    if (Kokkos::hwloc::available()) {
      num_threads = Kokkos::hwloc::get_available_numa_count()
                    * Kokkos::hwloc::get_available_cores_per_numa()
                 // * Kokkos::hwloc::get_available_threads_per_core()
                    ;

    }

    std::cout << "Threads: " << num_threads << std::endl;

    Kokkos::Threads::initialize( num_threads );
  }

  static void TearDownTestCase()
  {
    Kokkos::Threads::finalize();
  }
};

const size_t RUN_TIME_DIMENSION = 4000000;
const size_t COMPILE_TIME_DIMENSION = 3;

} // unnamed namespace


TEST_F( KokkosThreads, SerialInitialize)
{
  // allocate a rank 2 array witn that is RUN_TIME_DIMENSION x COMPILE_TIME_DIMENSION

  // View will default initialize all the values unless it is explicitly disabled, ie,
  // Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> a("node views", RUN_TIME_DIMENSION);
  // zero fills the array, but
  // Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> a( Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);
  // will allocate without initializing the array

  Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> a( Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);

  for (size_t i=0; i < a.dimension_0(); ++i) {
    for (size_t x=0; x < a.dimension_1(); ++x) {
      a(i,x) = i;
    }
  }

  // get a const view to the same array
  // this view shares the same memory as a, but cannot modify the values
  Kokkos::View<const unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> b = a;

  for (size_t i=0; i < b.dimension_0(); ++i) {
    for (size_t x=0; x < b.dimension_1(); ++x) {
      EXPECT_EQ(i, b(i,x));
    }
  }
}

// Not available until c++11 support is enable in Sierra and Kokkos
#if defined (KOKKOS_HAVE_C_PLUS_PLUS_11_LAMBDA)
TEST_F( KokkosThreads, LambdaInitialize)
{
  Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> a( Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);

  Kokkos::parallel_for<Kokkos::Threads>(
    a.dimension_0() ,
    [=](size_t i) {
      for (size_t x=0; x < a.dimension_1(); ++x) {
        a(i,x) = i;
      }
    }
  );

  Kokkos::View<const unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> b = a;

  int num_error = 0;
  // Cannot portably call a GTEST macro in parallel
  // count the errors and test that they are equal to zero
  Kokkos::parallel_reduce<Kokkos::Threads, int /*reduction value type */>(
    b.dimension_0() ,
    [](int & local_errors)                                    // init lambda
    { local_errors = 0; } ,
    [=](size_t i, int & local_errors) {                       // operator() lambda
      for (size_t x=0; x < b.dimension_1(); ++x)
        local_errors += i == b(i,x) ? 0 : 1;
    } ,
    [](volatile int & dst_err, volatile int const& src_err)   // join lambda
    { dst_err += src_err; } ,
    num_errors                                                // where to store the result
  );
  EXPECT_EQ( 0, num_errors);

}
#endif


namespace {

// Functors need to initialize and check the view with out c++11 support

template <typename View>
struct InitializeView
{
  // need a device_type typedef for all parallel functors
  typedef typename View::device_type device_type;

  View a;

  // get a view to the a
  template <typename RhsView>
  InitializeView( RhsView const& arg_a )
    : a(arg_a)
  {}

  void apply()
  {
    // call parallel_for on this functor
    Kokkos::parallel_for( a.dimension_0(), *this);
  }

  // initialize the a
  KOKKOS_INLINE_FUNCTION
  void operator()(size_t i) const
  {
    for (size_t x=0; x < a.dimension_1(); ++x) {
      a(i,x) = i;
    }
  }
};

template <typename View>
struct CheckView
{
  // need a device_type typedef for all parallel functors
  typedef typename View::device_type device_type;

  // need a value_type typedef for the reduction
  typedef int value_type;

  View a;

  // get a view to the a
  template <typename RhsView>
  CheckView( RhsView const& arg_a )
    : a(arg_a)
  {}

  // return the number of errors found
  value_type apply()
  {
    int num_errors = 0;
    // call a parallel_reduce to count the errors
    Kokkos::parallel_reduce( a.dimension_0(), *this, num_errors);
    return num_errors;
  }

  // initialize the reduction type
  KOKKOS_INLINE_FUNCTION
  void init(value_type & v) const
  { v = 0; }

  // this threads contribution to the reduction type
  // check that the value is equal to the expected value
  // otherwise increment the error count
  KOKKOS_INLINE_FUNCTION
  void operator()(size_t i, value_type & error) const
  {
    for (size_t x=0; x < a.dimension_1(); ++x) {
      error += i == a(i,x) ? 0 : 1;
    }
  }

  // join two threads together
  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, volatile value_type const& src) const
  { dst += src; }
};

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

} // unnameed namespace

TEST_F( KokkosThreads, ParallelInitialize)
{
  typedef Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> view_type;
  typedef Kokkos::View<const unsigned*[COMPILE_TIME_DIMENSION], Kokkos::Threads> const_view_type;

  view_type a(Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);

  // call the InitializeView functor
  {
    InitializeView<view_type> f(a);
    f.apply();
  }

  // call the CheckView functor
  // and expect no errors
  {
    const_view_type b = a;
    CheckView<view_type> f(a);
    const int num_errors = f.apply();
    EXPECT_EQ( 0, num_errors);
  }
}

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

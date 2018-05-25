#include "KokkosCentroidCalculation.h"

STK_FUNCTION void get_value(DeviceViewVectorType data, size_t index, my_double& value)
{
    value = data(index);
}

my_double get_average(const std::vector<my_double>& data, int num_repeat)
{
    HostViewVectorType host_data("host_data", data.size());
    DeviceViewVectorType device_data("device_data", data.size());

    for(size_t i=0;i<data.size();++i)
    {
        host_data(i) = data[i];
    }

    Kokkos::deep_copy(device_data, host_data);

    my_double sum = 0;
    for(int j=0;j<num_repeat;++j)
    {
        sum = 0;
        Kokkos::parallel_reduce( data.size(), KOKKOS_LAMBDA ( size_t i, my_double &local_sum )
        {
            my_double value = 0;
            get_value(device_data, i, value);
            local_sum += value;
        },
        sum );
    }

    my_double num = data.size();
    return sum/num;
}

TEST_F(MTK_Kokkos, calculate_average)
{
    size_t size_data = stk::unit_test_util::get_command_line_option<size_t>("-s", 2);
    std::vector<my_double> data(size_data);
    for(size_t i=0;i<data.size();++i)
        data[i]=i;
    int num_repeat = stk::unit_test_util::get_command_line_option<size_t>("-n", 1);

    struct timeval begin,end;
    gettimeofday(&begin,NULL);

    my_double average = get_average(data, num_repeat);

    gettimeofday(&end,NULL);

    double time = 1.0*(end.tv_sec-begin.tv_sec) +
                  1.0e-6*(end.tv_usec-begin.tv_usec);
    double Gbytes = 1.0e-9 * double(sizeof(my_double) * ( data.size() )) ;
    printf("bandwidth( %g GB/s )\n", Gbytes * num_repeat / time );

    my_double num_items = data.size();
    my_double sum = (num_items*(num_items-1))/2.0;
    my_double gold_ave = sum/data.size();
    EXPECT_NEAR(gold_ave, average, 0.0000001);
}

#if defined(KOKKOS_ENABLE_CUDA)
void run_it_uvm()
{
    typedef Kokkos::View<int*, Layout, MemSpace> IntViewType;
    typedef Kokkos::View<int**, Layout, Kokkos::CudaUVMSpace> BucketConnectivityView;
    typedef Kokkos::View<BucketConnectivityView*, Layout, Kokkos::CudaUVMSpace> ViewOfBuckets;


    const unsigned numBuckets = 2;
    unsigned numElemsPerBucket[numBuckets] = {2, 3};
    unsigned numNodesPerElem[numBuckets] = {4, 8};

    IntViewType elemsPerBucket("D_elemsPerBucket", numBuckets);
    IntViewType::HostMirror hostElemsPerBucket = Kokkos::create_mirror_view(elemsPerBucket);

    IntViewType nodesPerElem("D_nodesPerElem", numBuckets);
    IntViewType::HostMirror hostNodesPerElem = Kokkos::create_mirror_view(nodesPerElem);

    for(unsigned i=0; i<numBuckets; ++i) {
        hostElemsPerBucket(i) = numElemsPerBucket[i];
        hostNodesPerElem(i) = numNodesPerElem[i];
    }

    Kokkos::deep_copy(elemsPerBucket, hostElemsPerBucket);
    Kokkos::deep_copy(nodesPerElem, hostNodesPerElem);

    IntViewType expectedBucketOffsets("D_expectedBucketOffsets", numBuckets);
    IntViewType::HostMirror hostExpectedBucketOffsets = Kokkos::create_mirror_view(expectedBucketOffsets);

    ViewOfBuckets viewOfBuckets("viewOfBuckets", numBuckets);
    int offset = 0;
    for(unsigned i=0; i<numBuckets; ++i) {
        viewOfBuckets(i) = BucketConnectivityView("Bucket"+std::to_string(i), numElemsPerBucket[i], numNodesPerElem[i]);
        hostExpectedBucketOffsets(i) = offset;

        for(unsigned elemIndex=0; elemIndex<numElemsPerBucket[i]; ++elemIndex) {
            for(unsigned nodeIndex=0; nodeIndex<numNodesPerElem[i]; ++nodeIndex) {
                viewOfBuckets(i)(elemIndex, nodeIndex) = offset++;
            }
        }
    }

    Kokkos::deep_copy(expectedBucketOffsets, hostExpectedBucketOffsets);

    double errorCheck = 0.0;
    Kokkos::parallel_reduce(numBuckets, KOKKOS_LAMBDA(int i, double& update) {
        BucketConnectivityView bucket = viewOfBuckets(i);
        int expectedOffset = expectedBucketOffsets(i);
        for(unsigned elemIndex=0; elemIndex<elemsPerBucket(i); ++elemIndex) {
            for(unsigned nodeIndex=0; nodeIndex<nodesPerElem(i); ++nodeIndex) {
               int offset = bucket(elemIndex, nodeIndex);
               printf("i=%d, elemIndex=%d, nodeIndex=%d, offset=%d, expectedOffset=%d\n",i, elemIndex,nodeIndex,offset,expectedOffset);
               if (offset != expectedOffset) {
                   update += 1.0;
               }
               ++expectedOffset;
            }
        }
    }, errorCheck);

    EXPECT_EQ(0.0, errorCheck);
}

TEST_F(MTK_Kokkos, view_of_views_uvm)
{
    run_it_uvm();
}
#endif

void run_it_flat()
{
    typedef Kokkos::View<int*, Layout, MemSpace> IntViewType;
    typedef Kokkos::View<int*, Layout, MemSpace> ViewOfBuckets;
    typedef Kokkos::View<int*, Layout, Kokkos::HostSpace> HostViewOfBuckets;

    const unsigned numBuckets = 2;
    unsigned numElemsPerBucket[numBuckets] = {2, 3};
    unsigned numNodesPerElem[numBuckets] = {4, 8};

    IntViewType elemsPerBucket("D_elemsPerBucket", numBuckets);
    IntViewType::HostMirror hostElemsPerBucket = Kokkos::create_mirror_view(elemsPerBucket);

    IntViewType nodesPerElem("D_nodesPerElem", numBuckets);
    IntViewType::HostMirror hostNodesPerElem = Kokkos::create_mirror_view(nodesPerElem);

    IntViewType bucketOffsets("D_bucketOffsets", numBuckets);
    IntViewType::HostMirror hostBucketOffsets = Kokkos::create_mirror_view(bucketOffsets);

    unsigned numConnectivities = 0;
    
    for(unsigned i=0; i<numBuckets; ++i) {
        hostBucketOffsets(i) = numConnectivities;
      
        hostElemsPerBucket(i) = numElemsPerBucket[i];
        hostNodesPerElem(i) = numNodesPerElem[i];

        numConnectivities += numElemsPerBucket[i]*numNodesPerElem[i]; 
    }

    Kokkos::deep_copy(elemsPerBucket, hostElemsPerBucket);
    Kokkos::deep_copy(nodesPerElem, hostNodesPerElem);
    Kokkos::deep_copy(bucketOffsets, hostBucketOffsets);

    IntViewType expectedBucketOffsets("D_expectedBucketOffsets", numBuckets);
    IntViewType::HostMirror hostExpectedBucketOffsets = Kokkos::create_mirror_view(expectedBucketOffsets);

    ViewOfBuckets viewOfBuckets("D_viewOfBuckets", numConnectivities);
    ViewOfBuckets::HostMirror hostViewOfBuckets = Kokkos::create_mirror_view(viewOfBuckets);

    unsigned offset = 0;
    for(unsigned i=0; i<numBuckets; ++i) {
        hostExpectedBucketOffsets(i) = offset;

        unsigned bucketOffset = hostBucketOffsets(i);
  
        for(unsigned elemIndex=0; elemIndex<numElemsPerBucket[i]; ++elemIndex) {
            unsigned elemOffset = elemIndex*numNodesPerElem[i] + bucketOffset;
    
            for(unsigned nodeIndex=0; nodeIndex<numNodesPerElem[i]; ++nodeIndex) {
                unsigned nodeOffset = elemOffset + nodeIndex;
                EXPECT_EQ(nodeOffset, offset);
    
                hostViewOfBuckets(nodeOffset) = offset++;
            }
        }
    }

    Kokkos::deep_copy(expectedBucketOffsets, hostExpectedBucketOffsets);
    Kokkos::deep_copy(viewOfBuckets, hostViewOfBuckets);

    double errorCheck = 0.0;
    Kokkos::parallel_reduce(numBuckets, KOKKOS_LAMBDA(int i, double& update) {
        int expectedConnId = expectedBucketOffsets(i);
        unsigned bucketOffset = bucketOffsets(i);

        for(int elemIndex=0; elemIndex<elemsPerBucket(i); ++elemIndex) {
            int elemOffset = elemIndex*nodesPerElem(i) + bucketOffset;
    
            for(int nodeIndex=0; nodeIndex<nodesPerElem(i); ++nodeIndex) {
               int nodeOffset = elemOffset + nodeIndex;
               int connId = viewOfBuckets(nodeOffset);
               if (connId != expectedConnId) {
                   update += 1.0;
               }
               ++expectedConnId;
            }
        }
    }, errorCheck);

    EXPECT_EQ(0.0, errorCheck);
}

TEST_F(MTK_Kokkos, view_of_views_flat)
{
    run_it_flat();
}

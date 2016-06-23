#include "KokkosBulkDataCentroidCalculation.h"
#include "StaticMesh.h"


struct GpuGatherBucketScratchData
{
    GpuGatherBucketScratchData(const stk::mesh::BulkData &bulk)
    : combineAllBuckets(false), ngpMesh(bulk)
    {
    }

    unsigned getNumParallelItems() const
    {
        return ngpMesh.num_buckets(stk::topology::ELEM_RANK);
    }

    void setup_mesh_indices(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
    {
        const stk::mesh::BucketVector& nodeBuckets = bulk.buckets(stk::topology::NODE_RANK);

        meshIndices = DeviceViewMeshIndicesType("DMeshIndices", bulk.get_size_of_entity_index_space());
        constMeshIndices = meshIndices;
        hostMeshIndices = Kokkos::create_mirror_view(meshIndices);

        for(unsigned bktIndex = 0; bktIndex < nodeBuckets.size(); ++bktIndex)
        {
            const stk::mesh::Bucket& bucket = *nodeBuckets[bktIndex];
            for(unsigned nodeIndex = 0; nodeIndex < bucket.size(); ++nodeIndex)
            {
                stk::mesh::Entity node = bucket[nodeIndex];
                stk::mesh::FastMeshIndex meshIndex(bucket.bucket_id(), nodeIndex);
                hostMeshIndices(node.local_offset()) = meshIndex;
            }
        }

        const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

        for(unsigned bktIndex = 0; bktIndex < elemBuckets.size(); ++bktIndex)
        {
            const stk::mesh::Bucket& bucket = *elemBuckets[bktIndex];
            for(unsigned elemIndex = 0; elemIndex < bucket.size(); ++elemIndex)
            {
                stk::mesh::Entity elem = bucket[elemIndex];
                stk::mesh::FastMeshIndex meshIndex(bucket.bucket_id(), elemIndex);
                hostMeshIndices(elem.local_offset()) = meshIndex;
            }
        }

        Kokkos::deep_copy(meshIndices, hostMeshIndices);
    }

    void setup_node_coords(const stk::mesh::BulkData& bulk, const CoordFieldType& coords, const stk::mesh::Selector& selector)
    {
        const stk::mesh::BucketVector& nodeBuckets = bulk.buckets(stk::topology::NODE_RANK);
        bucketCapacity = nodeBuckets[0]->capacity();
        unsigned nodeAllocSize = nodeBuckets.size()*bucketCapacity;

        nodeCoords = DeviceViewMatrixType("DNodeCoords", nodeAllocSize, bulk.mesh_meta_data().spatial_dimension());
        constNodeCoords = nodeCoords;
        hostNodeCoords =  Kokkos::create_mirror_view(nodeCoords);

        for(unsigned bktIndex = 0; bktIndex < nodeBuckets.size(); ++bktIndex)
        {
            const stk::mesh::Bucket& bucket = *nodeBuckets[bktIndex];
            for(unsigned nodeIndex = 0; nodeIndex < bucket.size(); ++nodeIndex)
            {
                stk::mesh::Entity node = bucket[nodeIndex];
                double *node_coords = stk::mesh::field_data(coords, node);
                unsigned nodeCoordIndex = host_get_index(node);
                for(unsigned k=0; k<bulk.mesh_meta_data().spatial_dimension(); ++k)
                {
                    hostNodeCoords(nodeCoordIndex, k) = node_coords[k];
                }
            }
        }

        Kokkos::deep_copy(nodeCoords, hostNodeCoords);
    }

    void setup_element_centroids(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
    {
        const stk::mesh::BucketVector& elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);
        bucketCapacity = elemBuckets[0]->capacity();
        unsigned elemAllocSize = elemBuckets.size()*bucketCapacity;

        elementCentroids = DeviceViewMatrixType("DElemCentroids", elemAllocSize, bulk.mesh_meta_data().spatial_dimension());
        hostElementCentroids =  Kokkos::create_mirror_view(elementCentroids);
    }

    void initialize(const stk::mesh::BulkData& bulk, const CoordFieldType& coords, CoordFieldType& centroid, const stk::mesh::Selector& selector)
    {
        setup_mesh_indices(bulk, selector);
        setup_node_coords(bulk, coords, selector);
        setup_element_centroids(bulk, selector);
    }

    unsigned host_get_index(stk::mesh::Entity entity) const
    {
        const stk::mesh::FastMeshIndex& meshIndex = hostMeshIndices(entity.local_offset());
        return meshIndex.bucket_id * bucketCapacity + meshIndex.bucket_ord;
    }

    STK_FUNCTION unsigned get_index(stk::mesh::Entity entity) const
    {
        const stk::mesh::FastMeshIndex& meshIndex = constMeshIndices(entity.local_offset());
        return meshIndex.bucket_id * bucketCapacity + meshIndex.bucket_ord;
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, solo, compact), const int elementBucketIndex) const
    {
        const ngp::StaticBucket& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, elementBucketIndex);
        const unsigned numElements = bucket.size();
        const unsigned nodesPerElem = bucket.get_num_nodes_per_entity();
        const unsigned dim = elementCentroids.extent(1);
        double temp[3] = {0.0, 0.0, 0.0};
        for(unsigned elementIndex=0; elementIndex<numElements; ++elementIndex) {
            const unsigned elemFieldIndex = get_index(bucket[elementIndex]);
            for(unsigned k=0;k<dim;++k) // loop over x y z coordinates
            {
                temp[k] = 0.0;
            }
            ngp::ConnectedNodesType nodesView = bucket.get_nodes(elementIndex);
            for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
            {
                unsigned idx = get_index(nodesView(nodeIndex));
                for(unsigned k=0; k<dim; ++k) {
                    temp[k] += nodeCoords(idx, k);
                }
            }
            for(unsigned k=0; k<dim; ++k) {
                elementCentroids(elemFieldIndex, k) = temp[k] * 0.125;
            }
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, solo, unroll), const int elementBucketIndex) const
    {
        const ngp::StaticBucket& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, elementBucketIndex);
        const int numElements = bucket.size();
        const int nodesPerElem = bucket.get_num_nodes_per_entity();
        double tempx = 0.0, tempy = 0.0, tempz = 0.0;
        for(int elementIndex=0; elementIndex<numElements; ++elementIndex) {
            const int elemFieldIndex = get_index(bucket[elementIndex]);
            tempx = 0.0;
            tempy = 0.0;
            tempz = 0.0;
            ngp::ConnectedNodesType nodesView = bucket.get_nodes(elementIndex);
            for(int nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
            {
                int idx = get_index(nodesView(nodeIndex));
                tempx += constNodeCoords(idx, 0);
                tempy += constNodeCoords(idx, 1);
                tempz += constNodeCoords(idx, 2);
            }
            elementCentroids(elemFieldIndex, 0) = tempx * 0.125;
            elementCentroids(elemFieldIndex, 1) = tempy * 0.125;
            elementCentroids(elemFieldIndex, 2) = tempz * 0.125;
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, team, compact), const TeamHandleType& team) const
    {
        const int elementBucketIndex = team.league_rank();
        const ngp::StaticBucket& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, elementBucketIndex);
        const unsigned numElements = bucket.size();
        const unsigned nodesPerElem = bucket.get_num_nodes_per_entity();
        const unsigned dim = elementCentroids.extent(1);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& elementIndex) {
            double temp[3] = {0.0, 0.0, 0.0};
            const unsigned elemFieldIndex = get_index(bucket[elementIndex]);
            ngp::ConnectedNodesType nodesView = bucket.get_nodes(elementIndex);
            for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
            {
                unsigned idx = get_index(nodesView(nodeIndex));
                for(unsigned k=0; k<dim; ++k) {
                    temp[k] += nodeCoords(idx, k);
                }
            }
            for(unsigned k=0; k<dim; ++k) {
                elementCentroids(elemFieldIndex, k) = temp[k] * 0.125;
            }
        });
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, solo, compact), const int elementIndex) const
    {
        printf("element-solo-compact not implemented!\n");
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, solo, unroll), const int elementIndex) const
    {
        printf("element-solo-unroll not implemented!\n");
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, team, unroll), const TeamHandleType& team) const
    {
        const int elementBucketIndex = team.league_rank();
        const ngp::StaticBucket& bucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, elementBucketIndex);
        const unsigned numElements = bucket.size();
        const unsigned nodesPerElem = bucket.get_num_nodes_per_entity();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& elementIndex) {
            double tempx = 0, tempy = 0, tempz = 0;
            const unsigned elemFieldIndex = get_index(bucket[elementIndex]);
            ngp::ConnectedNodesType elemNodes = bucket.get_nodes(elementIndex);
            for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
            {
                unsigned idx = get_index(elemNodes(nodeIndex));
                tempx += constNodeCoords(idx, 0);
                tempy += constNodeCoords(idx, 1);
                tempz += constNodeCoords(idx, 2);
            }
            elementCentroids(elemFieldIndex, 0) = tempx * 0.125;
            elementCentroids(elemFieldIndex, 1) = tempy * 0.125;
            elementCentroids(elemFieldIndex, 2) = tempz * 0.125;
        });
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, team, unroll), const TeamHandleType& team) const
    {
        printf("element-team-unroll not implemented!\n");
    }

    unsigned bucketCapacity;
    bool combineAllBuckets;

    ngp::StaticMesh ngpMesh;

    DeviceViewMatrixType nodeCoords;
    ConstDeviceViewMatrixType constNodeCoords;
    DeviceViewMatrixType elementCentroids;
    DeviceViewMeshIndicesType meshIndices;
    ConstDeviceViewMeshIndicesType constMeshIndices;

    DeviceViewMatrixType::HostMirror hostNodeCoords;
    DeviceViewMatrixType::HostMirror hostElementCentroids;
    DeviceViewMeshIndicesType::HostMirror hostMeshIndices;
};

void run_bucket_test()
{
    const unsigned numElements = 9;
    const unsigned numNodesPerElement = 4;
    const unsigned bucketId = 0;

    ngp::StaticBucket bucket;
    bucket.initialize(bucketId, stk::topology::ELEM_RANK, numElements, numNodesPerElement);

    unsigned counter = 0;
    for(unsigned elemIndex=0; elemIndex<numElements; ++elemIndex) {
        for(unsigned nodeIndex=0; nodeIndex<numNodesPerElement; ++nodeIndex) {
           bucket.hostConnectivity(elemIndex,nodeIndex) = stk::mesh::Entity(counter);
           ++counter;
        }
    }
    Kokkos::deep_copy(bucket.connectivity, bucket.hostConnectivity);

    double errorCheck = 0;
    Kokkos::parallel_reduce(numElements, KOKKOS_LAMBDA(int elementIndex, double& update) {
        ngp::ConnectedNodesType nodesView = bucket.get_nodes(elementIndex);
        unsigned expectedCounter = elementIndex*numNodesPerElement;
        for(unsigned nodeIndex=0; nodeIndex<numNodesPerElement; ++nodeIndex) {
            if (nodesView(nodeIndex) != expectedCounter) {
                update += 1.0;
            }
            ++expectedCounter;
        }
    }, errorCheck);
    EXPECT_EQ(0.0, errorCheck);
}

TEST_F(MTK_Kokkos, bucket0)
{
    run_bucket_test();
}

TEST_F(MTK_Kokkos, calculate_centroid_field_with_gather_on_device_bucket)
{
    MyApp app;

    GpuGatherBucketScratchData scratch(*app.bulk);
    scratch.initialize(*app.bulk, *app.coords, app.centroid, app.meta.locally_owned_part());

    CentroidCalculator<GpuGatherBucketScratchData> calculator(scratch);

    app.start_timer();
    calculator.calculate_centroids(app.num_repeat, app.choice, app.teamSize);
    app.stop_timer();
    app.report_bandwidth();

    calculator.copy_centroids_to_host();

    for(unsigned bucketIndex=0; bucketIndex<scratch.ngpMesh.num_buckets(stk::topology::ELEM_RANK); ++bucketIndex) {
      const ngp::StaticBucket &bucket = scratch.ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
      unsigned numElements = bucket.size();
      for(unsigned elementIndex=0; elementIndex<numElements; ++elementIndex) {
          unsigned elemFieldIndex = scratch.host_get_index(bucket.host_get_entity(elementIndex));
          calculator.test_centroid_of_element(app.hostCentroid, bucket.host_get_entity(elementIndex), elemFieldIndex);
      }
    }
}


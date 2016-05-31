#include "KokkosBulkDataCentroidCalculation.h"


struct GpuGatherUVMScratchData
{
    GpuGatherUVMScratchData()
    : combineAllBuckets(false)
    {
    }

    unsigned getNumParallelItems() const
    {
        return numParallelItems;
    }

    void setup_elem_entities_and_connectivity_tables(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector, bool combineBuckets)
    {
        const stk::mesh::BucketVector& elementBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);
        unsigned numElementBuckets = elementBuckets.size();

        combineAllBuckets = combineBuckets;

        numParallelItems = numElementBuckets;
        if (combineBuckets) {
            numParallelItems = 0;
            for(const stk::mesh::Bucket* bktptr : elementBuckets) {
                numParallelItems += bktptr->size();
            }
            
            elemEntities = EntityBucketsType("ElemEntities", 1);
            elemEntities(0) = EntityViewType("ElemEntities0", numParallelItems);

            elementNodeConnectivity = ViewOfBucketsType("ElementNodeConnectivity", 1);
            unsigned nodesPerElem = elementBuckets[0]->topology().num_nodes();
            elementNodeConnectivity(0) = BucketConnectivityType("BucketConnectivity0", numParallelItems, nodesPerElem);
        }
        else {
            elemEntities = EntityBucketsType("ElemEntities", numElementBuckets);
            elementNodeConnectivity = ViewOfBucketsType("ElementNodeConnectivity", numElementBuckets);
        }

        unsigned elemOffset = 0;
        for (unsigned elemBucketIndex = 0; elemBucketIndex < numElementBuckets; ++elemBucketIndex)
        {
            const stk::mesh::Bucket& bucket = *elementBuckets[elemBucketIndex];
            unsigned nodesPerElem = bucket.topology().num_nodes();
            if (!combineBuckets) {
                elementNodeConnectivity(bucket.bucket_id()) = BucketConnectivityType("BucketConnectivity"+std::to_string(bucket.bucket_id()), bucket.size(), nodesPerElem);
                elemEntities(bucket.bucket_id()) = EntityViewType("ElemEntities"+std::to_string(bucket.bucket_id()), bucket.size());
                elemOffset = 0;
            }
            unsigned bucket_id = combineBuckets ? 0 : bucket.bucket_id();

            for(unsigned elemIndex = 0; elemIndex < bucket.size(); ++elemIndex)
            {
                stk::mesh::Entity element = bucket[elemIndex];
                elemEntities(bucket_id)(elemOffset) = element;
                const stk::mesh::Entity * elemNodes = bulk.begin_nodes(element);
                for(unsigned iNode = 0; iNode < nodesPerElem; ++iNode)
                {
                    elementNodeConnectivity(bucket_id)(elemOffset, iNode) = elemNodes[iNode];
                }
                ++elemOffset;
            }
        }
    }

    void setup_mesh_indices(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
    {
        const stk::mesh::BucketVector& nodeBuckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);

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
        const stk::mesh::BucketVector& nodeBuckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);
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

    void initialize(const stk::mesh::BulkData& bulk, const CoordFieldType& coords, CoordFieldType& centroid, const stk::mesh::Selector& selector, bool combineBuckets)
    {
        setup_elem_entities_and_connectivity_tables(bulk, selector, combineBuckets);
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
        BucketConnectivityType bucket = elementNodeConnectivity(elementBucketIndex);
        const unsigned numElements = elemEntities(elementBucketIndex).extent(0);
        const unsigned nodesPerElem = bucket.extent(1);
        const unsigned dim = elementCentroids.extent(1);
        double temp[3] = {0.0, 0.0, 0.0};
        for(unsigned elementIndex=0; elementIndex<numElements; ++elementIndex) {
            const unsigned elemFieldIndex = get_index(elemEntities(elementBucketIndex)(elementIndex));
            for(unsigned k=0;k<dim;++k) // loop over x y z coordinates
            {
                temp[k] = 0.0;
            }
            for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
            {
                unsigned idx = get_index(bucket(elementIndex, nodeIndex));
                for(unsigned k=0; k<dim; ++k) {
                    temp[k] += nodeCoords(idx, k);
                }
            }
            for(unsigned k=0; k<dim; ++k) {
                elementCentroids(elemFieldIndex, k) = temp[k] * 0.125;
            }
        }
    }


    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, team, compact), const TeamHandleType& team) const
    {
        const int elementBucketIndex = team.league_rank();
        BucketConnectivityType bucket = elementNodeConnectivity(elementBucketIndex);
        const auto entityView = elemEntities(elementBucketIndex);
        const unsigned numElements = elemEntities(elementBucketIndex).extent(0);
        const unsigned nodesPerElem = bucket.extent(1);
        const unsigned dim = elementCentroids.extent(1);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& elementIndex) {
            double temp[3] = {0.0, 0.0, 0.0};
            const unsigned elemFieldIndex = get_index(entityView(elementIndex));
            for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
            {
                unsigned idx = get_index(bucket(elementIndex, nodeIndex));
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
        BucketConnectivityType bucket = elementNodeConnectivity(0);
        const unsigned nodesPerElem = bucket.extent(1);
        const unsigned dim = elementCentroids.extent(1);
        stk::mesh::Entity elem = elemEntities(0)(elementIndex);
        const unsigned elemFieldIndex = get_index(elem);
        double temp[3] = {0, 0, 0};
        for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
        {
            unsigned idx = get_index(bucket(elementIndex, nodeIndex));
            for(unsigned k=0; k<dim; ++k) {
                temp[k] += nodeCoords(idx, k);
            }
        }
        for(unsigned k=0; k<dim; ++k) {
            elementCentroids(elemFieldIndex, k) = temp[k] * 0.125;
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, solo, unroll), const int elementIndex) const
    {
        BucketConnectivityType bucket = elementNodeConnectivity(0);
        const unsigned nodesPerElem = bucket.extent(1);
        stk::mesh::Entity elem = elemEntities(0)(elementIndex);
        const unsigned elemFieldIndex = get_index(elem);
	double tempx = 0, tempy = 0, tempz = 0;
        for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
        {
            unsigned idx = get_index(bucket(elementIndex, nodeIndex));
	    tempx += constNodeCoords(idx, 0);
	    tempy += constNodeCoords(idx, 1);
	    tempz += constNodeCoords(idx, 2);
        }
	elementCentroids(elemFieldIndex, 0) = tempx * 0.125;
	elementCentroids(elemFieldIndex, 1) = tempy * 0.125;
	elementCentroids(elemFieldIndex, 2) = tempz * 0.125;
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, team, unroll), const TeamHandleType& team) const
    {
        const int elementBucketIndex = team.league_rank();
        BucketConnectivityType bucket = elementNodeConnectivity(elementBucketIndex);
        const auto entityView = elemEntities(elementBucketIndex);
        const unsigned numElements = elemEntities(elementBucketIndex).extent(0);
        const unsigned nodesPerElem = bucket.extent(1);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& elementIndex) {
            double tempx = 0, tempy = 0, tempz = 0;
            const unsigned elemFieldIndex = get_index(entityView(elementIndex));
            for(unsigned nodeIndex=0;nodeIndex<nodesPerElem;++nodeIndex) // loop over every node of this element
            {
                unsigned idx = get_index(bucket(elementIndex, nodeIndex));
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

    }

//    KOKKOS_INLINE_FUNCTION void operator()(const member_type &teamMember) const
//    {
//        unsigned i = teamMember.league_rank();
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, 8), [=] (int j)
//        {
//            const unsigned dim = elementCentroids.extent(1);
//            double temp = elementCentroids(i, j);
//            for(unsigned k=0;k<dim;++k) // loop over x y z coordinates
//                temp += elementNodes(i, j, k); // sum the coordinates
//            elementCentroids(i, j) = temp * 0.125;
//        });
//    }

    unsigned bucketCapacity;
    unsigned numParallelItems;
    bool combineAllBuckets;

    ViewOfBucketsType elementNodeConnectivity;
    EntityBucketsType elemEntities;

    DeviceViewMatrixType nodeCoords;
    ConstDeviceViewMatrixType constNodeCoords;
    DeviceViewMatrixType elementCentroids;
    DeviceViewMeshIndicesType meshIndices;
    ConstDeviceViewMeshIndicesType constMeshIndices;

    DeviceViewMatrixType::HostMirror hostNodeCoords;
    DeviceViewMatrixType::HostMirror hostElementCentroids;
    DeviceViewMeshIndicesType::HostMirror hostMeshIndices;
};

TEST_F(MTK_Kokkos, calculate_centroid_field_with_gather_on_device_uvm)
{
    MyApp app;

    GpuGatherUVMScratchData scratch;
    bool combineBuckets = app.choice==2 ? true : false;
    scratch.initialize(*app.bulk, *app.coords, app.centroid, app.meta.locally_owned_part(), combineBuckets);

    CentroidCalculator<GpuGatherUVMScratchData> calculator(scratch);

    calculator.calculate_centroids(1, app.choice);

    app.start_timer();
    calculator.calculate_centroids(app.num_repeat, app.choice);
    app.stop_timer();
    app.report_bandwidth();

    calculator.copy_centroids_to_host();
    calculator.test_centroid_of_element_1();
}


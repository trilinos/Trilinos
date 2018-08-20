#include "KokkosBulkDataCentroidCalculation.hpp"
#include "Kokkos_Atomic.hpp"

struct GpuGatherFlatScratchData
{
    GpuGatherFlatScratchData(const int choice)
      : appChoice(choice)
    {
    }

    void setup_elem_entities_and_connectivity_tables(const stk::mesh::BulkData& bulk, const stk::mesh::Selector& selector)
    {
        const stk::mesh::BucketVector& elementBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);
        unsigned numElementBuckets = elementBuckets.size();

        connBucketOffsets = DeviceViewIntType("D_connBucketOffsets", numElementBuckets);
        hostConnBucketOffsets = Kokkos::create_mirror_view(connBucketOffsets);

        elemBucketOffsets = DeviceViewIntType("D_elemBucketOffsets", numElementBuckets);
        hostElemBucketOffsets = Kokkos::create_mirror_view(elemBucketOffsets);

        elemsPerBucket = DeviceViewIntType("D_elemsPerBucket", numElementBuckets);
        hostElemsPerBucket = Kokkos::create_mirror_view(elemsPerBucket);

        nodesPerElement = DeviceViewIntType("D_nodesPerElement", numElementBuckets);
        hostNodesPerElement = Kokkos::create_mirror_view(nodesPerElement);

        unsigned numConnectivities = 0;
	unsigned numElements = 0;
	
        for (unsigned elemBucketIndex = 0; elemBucketIndex < numElementBuckets; ++elemBucketIndex)
	{
           const stk::mesh::Bucket& bucket = *elementBuckets[elemBucketIndex];

           unsigned numNodesPerElem = bucket.topology().num_nodes();
	   unsigned numElemsInBucket = bucket.size();

           hostConnBucketOffsets(elemBucketIndex) = numConnectivities;
           hostElemBucketOffsets(elemBucketIndex) = numElements;
	   hostNodesPerElement(elemBucketIndex) = numNodesPerElem;
	   hostElemsPerBucket(elemBucketIndex) = numElemsInBucket;
	   
    	   numConnectivities += numNodesPerElem*numElemsInBucket;
	   numElements += numElemsInBucket;
        }

	Kokkos::deep_copy(elemsPerBucket, hostElemsPerBucket);
	Kokkos::deep_copy(nodesPerElement, hostNodesPerElement);
	Kokkos::deep_copy(connBucketOffsets, hostConnBucketOffsets);
	Kokkos::deep_copy(elemBucketOffsets, hostElemBucketOffsets);

        elementNodeConnectivity = DeviceViewFlatConnectivityType("DElementNodeConnectivity", numConnectivities);
        hostElementNodeConnectivity =  Kokkos::create_mirror_view(elementNodeConnectivity);

        elemEntities = DeviceViewEntitiesType("DElemEntities", numElements);
        hostElemEntities = Kokkos::create_mirror_view(elemEntities);

        for (unsigned elemBucketIndex = 0; elemBucketIndex < numElementBuckets; ++elemBucketIndex)
        {
	    unsigned connBucketOffset = hostConnBucketOffsets(elemBucketIndex);
	    unsigned elemBucketOffset = hostElemBucketOffsets(elemBucketIndex);

            const stk::mesh::Bucket& bucket = *elementBuckets[elemBucketIndex];
            unsigned nodesPerElem = bucket.topology().num_nodes();

            for(unsigned elemIndex = 0; elemIndex < bucket.size(); ++elemIndex)
            {
	        unsigned connElemOffset = elemIndex*nodesPerElem + connBucketOffset;

                stk::mesh::Entity element = bucket[elemIndex];
                hostElemEntities(elemBucketOffset + elemIndex) = element;

                const stk::mesh::Entity * elemNodes = bulk.begin_nodes(element);
                for(unsigned iNode = 0; iNode < nodesPerElem; ++iNode)
                {
		    unsigned nodeOffset = connElemOffset + iNode;
                    hostElementNodeConnectivity(nodeOffset) = elemNodes[iNode];
                }
            }
        }

        Kokkos::deep_copy(elementNodeConnectivity, hostElementNodeConnectivity);
        Kokkos::deep_copy(elemEntities, hostElemEntities);
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
                stk::mesh::FastMeshIndex meshIndex{bucket.bucket_id(), nodeIndex};
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
                stk::mesh::FastMeshIndex meshIndex{bucket.bucket_id(), elemIndex};
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

        elementCentroids = DeviceViewAtomicMatrixType("DElemCentroids", elemAllocSize, bulk.mesh_meta_data().spatial_dimension());
        hostElementCentroids =  Kokkos::create_mirror_view(elementCentroids);
    }

    void initialize(const stk::mesh::BulkData& bulk, const CoordFieldType& coords, CoordFieldType& centroid, const stk::mesh::Selector& selector)
    {
        setup_mesh_indices(bulk, selector);
        setup_node_coords(bulk, coords, selector);
        setup_elem_entities_and_connectivity_tables(bulk, selector);
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

    KOKKOS_INLINE_FUNCTION unsigned get_bucket_id(stk::mesh::Entity entity) const
    {
        const stk::mesh::FastMeshIndex& meshIndex = constMeshIndices(entity.local_offset());
        return meshIndex.bucket_id;
    }

    unsigned getNumParallelItems() const
    {
        if((2 == appChoice) || (4 == appChoice) || (5 == appChoice))
	  return hostElemEntities.extent(0);

        return hostConnBucketOffsets.extent(0);
    }

    KOKKOS_INLINE_FUNCTION void operator()( TYPE_OPERATOR(bucket, solo, compact) , const int elementBucketIndex) const
    {
        const unsigned numElements = elemsPerBucket(elementBucketIndex);
        const unsigned nodesPerElem = nodesPerElement(elementBucketIndex);
        const unsigned dim = elementCentroids.extent(1);
        double temp[3] = {0.0, 0.0, 0.0};
        for(unsigned elementIndex=0; elementIndex<numElements; ++elementIndex) {
   	    int elementOffset = elemBucketOffsets(elementBucketIndex) + elementIndex;	  
	    int connOffset = connBucketOffsets(elementBucketIndex) + elementIndex*nodesPerElem;
	    stk::mesh::Entity element = elemEntities(elementOffset);
	    const unsigned elemFieldIndex = get_index(element);

            for(unsigned k=0;k<dim;++k)  {
                temp[k] = 0.0;
            }
	    for(unsigned nodeIndex=0; nodeIndex<nodesPerElem; ++nodeIndex) // loop over every node of this element
            {
	        const unsigned idx = get_index(elementNodeConnectivity(connOffset + nodeIndex));
		  
		for(unsigned k=0; k<dim; ++k) {
		  temp[k] += nodeCoords(idx, k);
		}
	    }
	    for(unsigned k=0; k<dim; ++k) {
	      elementCentroids(elemFieldIndex, k) = temp[k] * 0.125;
            }
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()( TYPE_OPERATOR(bucket, solo, unroll) , const int elementBucketIndex) const
    {
        const unsigned numElements = elemsPerBucket(elementBucketIndex);
        const unsigned nodesPerElem = nodesPerElement(elementBucketIndex);
        double tempx = 0, tempy = 0, tempz = 0;
        for(unsigned elementIndex=0; elementIndex<numElements; ++elementIndex) {
   	    int elementOffset = elemBucketOffsets(elementBucketIndex) + elementIndex;	  
	    int connOffset = connBucketOffsets(elementBucketIndex) + elementIndex*nodesPerElem;
	    stk::mesh::Entity element = elemEntities(elementOffset);
	    const unsigned elemFieldIndex = get_index(element);

	    tempx = 0;
	    tempy = 0;
	    tempz = 0;

	    for(unsigned nodeIndex=0; nodeIndex<nodesPerElem; ++nodeIndex) // loop over every node of this element
            {
	        const unsigned idx = get_index(elementNodeConnectivity(connOffset + nodeIndex));
		  
		tempx += nodeCoords(idx, 0);
		tempy += nodeCoords(idx, 1);
		tempz += nodeCoords(idx, 2);
	    }

	    elementCentroids(elemFieldIndex, 0) = tempx * 0.125;
	    elementCentroids(elemFieldIndex, 1) = tempy * 0.125;
	    elementCentroids(elemFieldIndex, 2) = tempz * 0.125;
        }
    }
  
    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, team, compact), const TeamHandleType& team) const
    {
        const int elementBucketIndex = team.league_rank();

        const unsigned numElements = elemsPerBucket(elementBucketIndex);
        const unsigned nodesPerElem = nodesPerElement(elementBucketIndex);
        const unsigned dim = elementCentroids.extent(1);

        const int elemBucketOffset = elemBucketOffsets(elementBucketIndex);
	const int connBucketOffset = connBucketOffsets(elementBucketIndex);
	
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& elementIndex) {
            double temp[3] = {0.0, 0.0, 0.0};

   	    int elementOffset = elemBucketOffset + elementIndex;	  
	    int connOffset    = connBucketOffset + elementIndex*nodesPerElem;

	    const unsigned elemFieldIndex = get_index(elemEntities(elementOffset));
	    for(unsigned nodeIndex=0; nodeIndex<nodesPerElem; ++nodeIndex) // loop over every node of this element
            {
	        const unsigned idx = get_index(elementNodeConnectivity(connOffset + nodeIndex));
		for(unsigned k=0; k<dim; ++k) {
		  temp[k] += constNodeCoords(idx, k);
		}		
	    }

            for(unsigned k=0; k<dim; ++k) {
                elementCentroids(elemFieldIndex, k) = temp[k] * 0.125;
            }
       });
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, solo, compact), const int elementIndex) const
    {
        stk::mesh::Entity element = elemEntities(elementIndex);
        const unsigned elementBucketIndex = get_bucket_id(element);
        const unsigned nodesPerElem = nodesPerElement(elementBucketIndex);
        const unsigned dim = elementCentroids.extent(1);
	int elementOffsetInBucket = elementIndex - elemBucketOffsets(elementBucketIndex);	  
	int connOffset = connBucketOffsets(elementBucketIndex) + elementOffsetInBucket*nodesPerElem;
	const unsigned elemFieldIndex = get_index(element);
        double temp[3] = {0, 0, 0};
	for(unsigned nodeIndex=0; nodeIndex<nodesPerElem; ++nodeIndex) // loop over every node of this element
	{
    	    const unsigned idx = get_index(elementNodeConnectivity(connOffset + nodeIndex));
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
        stk::mesh::Entity element = elemEntities(elementIndex);
        const unsigned elementBucketIndex = get_bucket_id(element);
        const unsigned nodesPerElem = nodesPerElement(elementBucketIndex);
	int elementOffsetInBucket = elementIndex - elemBucketOffsets(elementBucketIndex);	  
	int connOffset = connBucketOffsets(elementBucketIndex) + elementOffsetInBucket*nodesPerElem;
	const unsigned elemFieldIndex = get_index(element);
	double tempx = 0, tempy = 0, tempz = 0;
	for(unsigned nodeIndex=0; nodeIndex<nodesPerElem; ++nodeIndex) // loop over every node of this element
	{
    	    const unsigned idx = get_index(elementNodeConnectivity(connOffset + nodeIndex));
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

        const unsigned numElements = elemsPerBucket(elementBucketIndex);
        const unsigned nodesPerElem = nodesPerElement(elementBucketIndex);

        const int elemBucketOffset = elemBucketOffsets(elementBucketIndex);
	const int connBucketOffset = connBucketOffsets(elementBucketIndex);

#if 1
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& elementIndex) {
            double tempx = 0, tempy = 0, tempz = 0;
   	    int elementOffset = elemBucketOffset + elementIndex;	  
	    int connOffset    = connBucketOffset + elementIndex*nodesPerElem;
	    const unsigned elemFieldIndex = get_index(elemEntities(elementOffset));
	    for(unsigned nodeIndex=0; nodeIndex<nodesPerElem; ++nodeIndex) // loop over every node of this element
            {
	        const unsigned idx = get_index(elementNodeConnectivity(connOffset + nodeIndex));
	    
		tempx += constNodeCoords(idx, 0);
		tempy += constNodeCoords(idx, 1);
		tempz += constNodeCoords(idx, 2);
	    }
	    elementCentroids(elemFieldIndex, 0) = tempx * 0.125;
	    elementCentroids(elemFieldIndex, 1) = tempy * 0.125;
	    elementCentroids(elemFieldIndex, 2) = tempz * 0.125;
       });
#else
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, numElements), [&] (const int& elementIndex) {
	    int elementOffset = elemBucketOffset + elementIndex;	  
	    const unsigned elemFieldIndex = get_index(elemEntities(elementOffset));
	
	    elementCentroids(elemFieldIndex, 0) = 0;
	    elementCentroids(elemFieldIndex, 1) = 0;
	    elementCentroids(elemFieldIndex, 2) = 0;
	});

	const unsigned range = numElements*nodesPerElem;
	
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, range), [&] (const int& rangeIndex) {
            const int elementIndex = rangeIndex/nodesPerElem;
	    unsigned nodeIndex = rangeIndex%nodesPerElem;
	    
   	    int elementOffset = elemBucketOffset + elementIndex;	  
	    int connOffset    = connBucketOffset + elementIndex*nodesPerElem;
	    const unsigned elemFieldIndex = get_index(elemEntities(elementOffset));

	    const unsigned idx = get_index(elementNodeConnectivity(connOffset + nodeIndex));

	    const double tempx = constNodeCoords(idx, 0) * 0.125;
	    const double tempy = constNodeCoords(idx, 1) * 0.125;
	    const double tempz = constNodeCoords(idx, 2) * 0.125;

	    Kokkos::atomic_add(&elementCentroids(elemFieldIndex, 0), tempx);
	    Kokkos::atomic_add(&elementCentroids(elemFieldIndex, 1), tempy);
	    Kokkos::atomic_add(&elementCentroids(elemFieldIndex, 2), tempz);
       });

#endif

    }
    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, team, unroll), const TeamHandleType& team) const
    {

    }

    unsigned bucketCapacity;
    int appChoice;
  
    DeviceViewFlatConnectivityType elementNodeConnectivity;
    DeviceViewMatrixType nodeCoords;
    ConstDeviceViewMatrixType constNodeCoords;
    DeviceViewAtomicMatrixType elementCentroids;
    DeviceViewMeshIndicesType meshIndices;
    ConstDeviceViewMeshIndicesType constMeshIndices;
    DeviceViewEntitiesType elemEntities;

    DeviceViewIntType connBucketOffsets;
    DeviceViewIntType elemBucketOffsets;
    DeviceViewIntType elemsPerBucket;
    DeviceViewIntType nodesPerElement;

    DeviceViewFlatConnectivityType::HostMirror hostElementNodeConnectivity;
    DeviceViewMatrixType::HostMirror hostNodeCoords;
    DeviceViewAtomicMatrixType::HostMirror hostElementCentroids;
    DeviceViewMeshIndicesType::HostMirror hostMeshIndices;
    DeviceViewEntitiesType::HostMirror  hostElemEntities;

    DeviceViewIntType::HostMirror hostConnBucketOffsets;
    DeviceViewIntType::HostMirror hostElemBucketOffsets;
    DeviceViewIntType::HostMirror hostElemsPerBucket;
    DeviceViewIntType::HostMirror hostNodesPerElement;
};

TEST_F(NGP_Kokkos, calculate_centroid_field_with_gather_on_device_flat)
{
    MyApp app;

    GpuGatherFlatScratchData scratch(app.choice);
    scratch.initialize(*app.bulk, *app.coords, app.centroid, app.meta.locally_owned_part());

    CentroidCalculator<GpuGatherFlatScratchData> calculator(scratch);
    
    app.start_timer();
    calculator.calculate_centroids(app.num_repeat, app.choice, app.teamSize);
    app.stop_timer();
    app.report_bandwidth();

    calculator.copy_centroids_to_host();
//    calculator.test_centroid_of_element_1();

    for(unsigned elementIndex=0; elementIndex<scratch.hostElemEntities.extent(0); ++elementIndex) {
        calculator.test_centroid_of_element(app.hostCentroid, scratch.hostElemEntities(elementIndex), elementIndex);
    }
}

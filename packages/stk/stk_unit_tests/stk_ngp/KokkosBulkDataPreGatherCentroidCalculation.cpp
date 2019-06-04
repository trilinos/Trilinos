#include "KokkosBulkDataCentroidCalculation.hpp"

struct my_second_operator{};

struct ScratchData
{
    ScratchData()
    {
    }

    void initialize(const stk::mesh::BulkData& bulk, const CoordFieldType& coords, CoordFieldType& centroid, const stk::mesh::Selector& selector)
    {
        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(selector, bulk.buckets(stk::topology::ELEM_RANK), elements);
        unsigned numElements = elements.size();

        elementNodes = DeviceArray3DType("DElementsNumNodes", numElements, 8, bulk.mesh_meta_data().spatial_dimension());
        hostElementNodes =  Kokkos::create_mirror_view(elementNodes);

        for (unsigned elemIndex = 0; elemIndex < numElements; ++elemIndex)
        {
            stk::mesh::Entity element = elements[elemIndex];
            const stk::mesh::Entity * elemNodes = bulk.begin_nodes(element);
            const unsigned numNodesThisElem = bulk.num_nodes(element);
	    EXPECT_EQ(8u, numNodesThisElem);	    
            for(unsigned iNode = 0; iNode < numNodesThisElem; ++iNode)
            {
                stk::mesh::Entity node = elemNodes[iNode];
                double *node_coords = stk::mesh::field_data(coords, node);
                for(unsigned k=0;k<bulk.mesh_meta_data().spatial_dimension();k++)
                {
                    hostElementNodes(elemIndex, iNode, k) = node_coords[k];
                }
            }
        }

        Kokkos::deep_copy(elementNodes, hostElementNodes);

        elementCentroids = DeviceArray2DType("Centroids", numElements, bulk.mesh_meta_data().spatial_dimension());
        hostElementCentroids =  Kokkos::create_mirror_view(elementCentroids);
    }

    unsigned getNumParallelItems() const
    {
        return elementCentroids.extent(0);
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, solo, compact), const int element) const
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, solo, unroll), const int element) const
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, team, compact), const TeamHandleType& team) const
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, solo, unroll) , const int element) const
    {
        const unsigned dim = elementCentroids.extent(1);
	const unsigned numNodes = elementNodes.extent(1);
	
        for(unsigned k=0;k<dim;++k) // loop over x y z coordinates
        {
            double temp = 0.0;
            for(unsigned node=0;node<numNodes;++node) // loop over every node of this element
                temp += elementNodes(element, node, k); // sum the coordinates
            elementCentroids(element, k) = temp * 0.125;
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, team, unroll), const TeamHandleType& team) const
    {
    }

    KOKKOS_INLINE_FUNCTION void operator()(my_second_operator, const int i) const
    {
        const unsigned dim = elementCentroids.extent(1);
        for (unsigned k=0 ; k<dim ; ++k)
          elementCentroids(i, k) *= 0.125; // divide by num nodes to get centroid
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(element, solo, compact), const int element) const
    {
        const unsigned dim = elementCentroids.extent(1);
	const unsigned numNodes = elementNodes.extent(1);
	double temp[3] = {0, 0, 0};
	for(unsigned node=0;node<numNodes;++node) // loop over every node of this element
	{
	  for(unsigned k=0; k<dim; ++k) // loop over x y z coordinates
	    temp[k] += elementNodes(element, node, k); // sum the coordinates

	  for(unsigned k=0; k<dim; ++k)
            elementCentroids(element, k) = temp[k] * 0.125;
        }
    }

    KOKKOS_INLINE_FUNCTION void operator()(TYPE_OPERATOR(bucket, team, unroll), const TeamHandleType& team) const
    {
    }

    DeviceArray3DType elementNodes;
    DeviceArray2DType elementCentroids;

    DeviceArray3DType::HostMirror hostElementNodes;
    DeviceArray2DType::HostMirror hostElementCentroids;
};

TEST_F(NGP_Kokkos, calculate_centroid_field_on_device)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
        MyApp app;
    
        ScratchData scratch;
        scratch.initialize(*app.bulk, *app.coords, app.centroid, app.meta.locally_owned_part());
    
        CentroidCalculator<ScratchData> calculator(scratch);
    
        app.start_timer();
        calculator.calculate_centroids(app.num_repeat, app.choice, app.teamSize );
        app.stop_timer();
        app.report_bandwidth();
    
        calculator.copy_centroids_to_host();
        calculator.test_centroid_of_element_1();
    
        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(app.meta.locally_owned_part(), app.bulk->buckets(stk::topology::ELEM_RANK), elements);
        for(unsigned elementIndex=0; elementIndex<elements.size(); ++elementIndex) {
            calculator.test_centroid_of_element(app.hostCentroid, elements[elementIndex], elementIndex);
        }
    }
}


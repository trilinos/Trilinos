template<typename Scalar , class DeviceType>
struct ForceGather;

template<typename Scalar>
struct ForceGather<Scalar ,KOKKOS_MACRO_DEVICE>{
	
	typedef KOKKOS_MACRO_DEVICE 								device_type;
	typedef device_type::size_type								size_type;

  	typedef Kokkos::MDArrayView<Scalar,device_type> 			scalar_array_d;
  	typedef Kokkos::MDArrayView<int,device_type> 				int_array_d;  	

	int_array_d		node_elemIDs;
	int_array_d		elem_nodeIDs;
	int_array_d		elems_per_node;
	
	scalar_array_d	nodal_force;
	scalar_array_d 	element_force;

  	ForceGather(int_array_d 	& arg_node_elemIDs,
				int_array_d 	& arg_elem_nodeIDs,
				int_array_d 	& arg_elems_per_node,
				scalar_array_d 	& arg_nodal_force,
				scalar_array_d 	& arg_element_force) :
			
				node_elemIDs	(arg_node_elemIDs),
				elem_nodeIDs	(arg_elem_nodeIDs),
				elems_per_node	(arg_elems_per_node),
				nodal_force		(arg_nodal_force),
				element_force	(arg_element_force){}


  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(int inode) const {

		int associated_elements = elems_per_node(inode + 1) - elems_per_node(inode);

		double local_force[] = {0.0, 0.0, 0.0};

	//	for each element that a node belongs to
		for(int i = 0; i < associated_elements; i++){

		//	elems_per_node is a cumulative structure, so 
		//	elems_per_node(inode) should be the index where
		//	a particular row's elem_IDs begin
			int nelem = node_elemIDs(elems_per_node(inode) + i, 0);

		//	find the row in an element's stiffness matrix
		//	that corresponds to inode
			int elem_row_index = node_elemIDs(elems_per_node(inode) + i, 1);

			local_force[0] += element_force(nelem, 0, elem_row_index);
			local_force[1] += element_force(nelem, 1, elem_row_index);
			local_force[2] += element_force(nelem, 2, elem_row_index);

		}

		nodal_force(inode, 0) = local_force[0];
		nodal_force(inode, 1) = local_force[1];
		nodal_force(inode, 2) = local_force[2];

  	}

}; //ForceGather



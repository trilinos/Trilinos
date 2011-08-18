template<typename Scalar , class DeviceType>
struct Gather;

template<typename Scalar>
struct Gather<Scalar ,KOKKOS_MACRO_DEVICE>{
	
	typedef KOKKOS_MACRO_DEVICE 									device_type;
	typedef device_type::size_type									size_type;

	typedef Kokkos::MultiVectorView<Scalar , device_type>			scalar_vector_d;
	typedef Kokkos::MultiVectorView<int , device_type>				int_vector_d;

  	typedef Kokkos::MDArrayView<Scalar,device_type> 				scalar_array_d;
  	typedef Kokkos::MDArrayView<int,device_type> 					int_array_d;  	
	
	scalar_vector_d	A;
  	scalar_vector_d b;
	int_vector_d	ArowIDs;
	int_vector_d 	AcolIDs;

	int_array_d	node_elemIDs;
	int_array_d	elem_nodeIDs;
	int_array_d	elems_per_node;
	
	scalar_array_d	element_stiffness;
	scalar_array_d element_load;

  	Gather(	scalar_vector_d & arg_A,
			scalar_vector_d & arg_b,
			int_vector_d	& arg_ArowIDs,
			int_vector_d	& arg_AcolIDs,
			int_array_d	& arg_node_elemIDs,
			int_array_d	& arg_elem_nodeIDs,
			int_array_d	& arg_elems_per_node,
			scalar_array_d	& arg_element_stiffness,
			scalar_array_d	& arg_element_load) :

			A(arg_A), 
			b(arg_b),
			ArowIDs(arg_ArowIDs),
			AcolIDs(arg_AcolIDs),
			node_elemIDs(arg_node_elemIDs),
			elem_nodeIDs(arg_elem_nodeIDs),
			elems_per_node(arg_elems_per_node),
			element_stiffness(arg_element_stiffness),
			element_load(arg_element_load){}

  	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()(int irow) const {

		int base_index = ArowIDs(irow);
		int last_index = ArowIDs(irow + 1);
		int associated_elements = elems_per_node(irow + 1) - elems_per_node(irow);

	//	for each element that a node belongs to
		for(int i = 0; i < associated_elements; i++){

		//	elems_per_node is a cumulative structure, so 
		//	elems_per_node(irow) should be the index where
		//	a particular row's elem_IDs begin
			int nelem = node_elemIDs(elems_per_node(irow) + i, 0);
			int elem_row_index = node_elemIDs(elems_per_node(irow) + i, 1);

		//	for each node in a particular related element	
			for(int j = 0; j < 8; j++){

			//	gather the contents of the element stiffness
			//	matrix that belong in irow
				int cid = elem_nodeIDs(nelem, j);
				int lower = 0;
				int column_search = 0;
				int upper = last_index - base_index;


				while((lower <= upper) && (AcolIDs(base_index + column_search) != cid)){
				
					column_search = (upper+lower) / 2;				

					if(AcolIDs(base_index + column_search) < cid)
						lower = column_search + 1;
					else
						upper = column_search - 1;

				}

				A(base_index + column_search) += element_stiffness(nelem, elem_row_index, j);
		
			}

			b(irow) += element_load(nelem, elem_row_index);

		}
  	}

}; //Gather



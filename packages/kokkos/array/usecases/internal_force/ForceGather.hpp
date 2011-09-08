template<typename Scalar , class DeviceType>
struct ForceGather;

template<typename Scalar>
struct ForceGather<Scalar ,KOKKOS_MACRO_DEVICE>{

  typedef KOKKOS_MACRO_DEVICE                 device_type;
  typedef device_type::size_type                size_type;

  typedef Kokkos::MDArrayView<Scalar,device_type>       scalar_array_d;
  typedef Kokkos::MDArrayView<int,device_type>         int_array_d;


  const int_array_d    node_elemIDs;   // Node-element connectivity via ids
  const int_array_d    node_elemOffset; // number elements per node

  const scalar_array_d  nodal_force;
  const scalar_array_d   element_force;

  const int        current_state;
  const int        previous_state;

    ForceGather(const int_array_d   & arg_node_elemIDs,
                const int_array_d   & arg_node_elemOffset,
                const scalar_array_d   & arg_nodal_force,
                const scalar_array_d   & arg_element_force,
                const int arg_current_state,
                const int arg_previous_state)
       : node_elemIDs  (arg_node_elemIDs)
       , node_elemOffset  (arg_node_elemOffset)
       , nodal_force    (arg_nodal_force)
       , element_force  (arg_element_force)
       , current_state(arg_current_state)
       , previous_state(arg_previous_state)
      {}


    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()(int inode) const {

      // Getting count as per 'CSR-like' data structure
      const int element_offset = node_elemOffset(inode);
      const int element_count  = node_elemOffset(inode + 1) - element_offset ;

    double local_force[] = {0.0, 0.0, 0.0};

  //  for each element that a node belongs to
    for(int i = 0; i < element_count ; i++){

    //  node_elemOffset is a cumulative structure, so
    //  node_elemOffset(inode) should be the index where
    //  a particular row's elem_IDs begin
      const int nelem = node_elemIDs( element_offset + i, 0);

    //  find the row in an element's stiffness matrix
    //  that corresponds to inode
      const int elem_node_index = node_elemIDs( element_offset + i, 1);

      local_force[0] += element_force(nelem, 0, elem_node_index);
      local_force[1] += element_force(nelem, 1, elem_node_index);
      local_force[2] += element_force(nelem, 2, elem_node_index);
    }

    nodal_force(inode, 0) = local_force[0];
    nodal_force(inode, 1) = local_force[1];
    nodal_force(inode, 2) = local_force[2];

    }

}; //ForceGather



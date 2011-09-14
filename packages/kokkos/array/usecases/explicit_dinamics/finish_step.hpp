template<typename Scalar , class DeviceType>
struct finish_step;

template<typename Scalar>
struct finish_step<Scalar ,KOKKOS_MACRO_DEVICE>{

  typedef KOKKOS_MACRO_DEVICE       device_type;
  typedef device_type::size_type    size_type;

  typedef Kokkos::MDArrayView<Scalar,device_type>   array_type;
  typedef Kokkos::MDArrayView<int,device_type>      int_array;

  typedef Region<Scalar,device_type> MyRegion;


    finish_step(const MyRegion  & region,
                const Scalar    arg_x_bc,
                const Scalar    arg_dt,
                const int       arg_current_state,
                const int       arg_next_state)
       : node_elem_ids(region.node_elem_ids)
       , node_elem_offset(region.node_elem_offset)
       , internal_force(region.internal_force)
       , element_force(region.element_force)
       , nodal_mass(region.nodal_mass)
       , acceleration(region.acceleration)
       , velocity(region.velocity)
       , displacement(region.displacement)
       , model_coords(region.model_coords)
       , x_bc(arg_x_bc)
       , dt(arg_dt)
       , current_state(arg_current_state)
       , next_state(arg_next_state)
      {}


    KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()(int inode) const {

      // Getting count as per 'CSR-like' data structure
      const int element_offset = node_elem_offset(inode);
      const int element_count  = node_elem_offset(inode + 1) - element_offset ;

      double local_force[] = {0.0, 0.0, 0.0};

      //  for each element that a node belongs to
      for(int i = 0; i < element_count ; i++){

        //  node_elem_offset is a cumulative structure, so
        //  node_elem_offset(inode) should be the index where
        //  a particular row's elem_IDs begin
        const int nelem = node_elem_ids( element_offset + i, 0);

        //  find the row in an element's stiffness matrix
        //  that corresponds to inode
        const int elem_node_index = node_elem_ids( element_offset + i, 1);

        local_force[0] += element_force(nelem, 0, elem_node_index);
        local_force[1] += element_force(nelem, 1, elem_node_index);
        local_force[2] += element_force(nelem, 2, elem_node_index);
      }

      internal_force(inode, 0) = local_force[0];
      internal_force(inode, 1) = local_force[1];
      internal_force(inode, 2) = local_force[2];


      Scalar v_new[3];
      Scalar a_current[3];

      if (model_coords(inode,0) != x_bc ) { //not on x boundary
        acceleration(inode,0) = a_current[0] = local_force[0] / nodal_mass(inode);
        acceleration(inode,1) = a_current[1] = local_force[1] / nodal_mass(inode);
        acceleration(inode,2) = a_current[2] = local_force[2] / nodal_mass(inode);
      } else { //enforce fixed BC
        acceleration(inode,0) = a_current[0] = 0;
        acceleration(inode,1) = a_current[1] = 0;
        acceleration(inode,2) = a_current[2] = 0;
      }

      velocity(inode,0,next_state) = v_new[0] = velocity(inode,0,current_state) + dt*a_current[0];
      velocity(inode,1,next_state) = v_new[1] = velocity(inode,1,current_state) + dt*a_current[1];
      velocity(inode,2,next_state) = v_new[2] = velocity(inode,2,current_state) + dt*a_current[2];

      displacement(inode,0,next_state) = displacement(inode,0,current_state) + dt*v_new[0];
      displacement(inode,1,next_state) = displacement(inode,1,current_state) + dt*v_new[1];
      displacement(inode,2,next_state) = displacement(inode,2,current_state) + dt*v_new[2];

    }


    const int_array       node_elem_ids;   // Node-element connectivity via ids
    const int_array       node_elem_offset; // number elements per node

    const array_type    internal_force;
    const array_type    element_force;
    const array_type      nodal_mass;
    const array_type      acceleration;
    const array_type      velocity;
    const array_type      displacement;
    const array_type      model_coords;
    const Scalar          x_bc;
    const Scalar          dt;
    const int             current_state;
    const int             next_state;

}; //finish_step



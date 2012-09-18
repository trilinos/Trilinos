/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

template<typename Scalar , class DeviceType>
struct finish_step;

template<typename Scalar>
struct finish_step<Scalar ,KOKKOSARRAY_MACRO_DEVICE>{

  typedef KOKKOSARRAY_MACRO_DEVICE       device_type;
  typedef device_type::size_type    size_type;

  typedef KokkosArray::MDArray<Scalar,device_type>   array_type;
  typedef KokkosArray::MDArray<int,device_type>      int_array;

  typedef Region<Scalar,device_type> MyRegion;

  typedef KokkosArray::Value<Scalar,device_type>     scalar;


    finish_step(const MyRegion  & region,
                const Scalar    arg_x_bc,
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
       , dt( region.dt)
       , prev_dt( region.prev_dt)
       , x_bc(arg_x_bc)
       , current_state(arg_current_state)
       , next_state(arg_next_state)
      {
        //std::cout << "finish_step dt: " << dt << std::endl;
        //std::cout << "finish_step prev_dt: " << prev_dt << std::endl;
      }


    KOKKOSARRAY_MACRO_DEVICE_FUNCTION
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

      const Scalar tol = 1.0e-7;

      if ( fabs(model_coords(inode,0)-x_bc) > tol ) { //not on x boundary
        acceleration(inode,0) = a_current[0] = -local_force[0] / nodal_mass(inode);
        acceleration(inode,1) = a_current[1] = -local_force[1] / nodal_mass(inode);
        acceleration(inode,2) = a_current[2] = -local_force[2] / nodal_mass(inode);
      } else { //enforce fixed BC
        acceleration(inode,0) = a_current[0] = 0;
        acceleration(inode,1) = a_current[1] = 0;
        acceleration(inode,2) = a_current[2] = 0;
      }

      velocity(inode,0,next_state) = v_new[0] = velocity(inode,0,current_state) + (*prev_dt+*dt)/2.0*a_current[0];
      velocity(inode,1,next_state) = v_new[1] = velocity(inode,1,current_state) + (*prev_dt+*dt)/2.0*a_current[1];
      velocity(inode,2,next_state) = v_new[2] = velocity(inode,2,current_state) + (*prev_dt+*dt)/2.0*a_current[2];

      displacement(inode,0,next_state) = displacement(inode,0,current_state) + *dt*v_new[0];
      displacement(inode,1,next_state) = displacement(inode,1,current_state) + *dt*v_new[1];
      displacement(inode,2,next_state) = displacement(inode,2,current_state) + *dt*v_new[2];

    }


    const int_array       node_elem_ids;   // Node-element connectivity via ids
    const int_array       node_elem_offset; // number elements per node

    const array_type      internal_force;
    const array_type      element_force;
    const array_type      nodal_mass;
    const array_type      acceleration;
    const array_type      velocity;
    const array_type      displacement;
    const array_type      model_coords;
    const scalar          dt;
    const scalar          prev_dt;
    const Scalar          x_bc;
    const int             current_state;
    const int             next_state;

}; //finish_step



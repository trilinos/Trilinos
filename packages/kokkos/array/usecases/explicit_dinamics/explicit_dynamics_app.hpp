#ifndef INTERNALFORCE
#define INTERNALFORCE

#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "BoxMeshFixture.hpp"
#include "region.hpp"
#include "initialize_element.hpp"
#include "initialize_node.hpp"
#include "grad_hgop.hpp"
#include "decomp_rotate.hpp"
#include "divergence.hpp"
#include "finish_step.hpp"

//----------------------------------------------------------------------------

template<typename Scalar, class device_type>
double explicit_dynamics_app( const size_t ex, const size_t ey, const size_t ez )
{
  typedef typename Kokkos::MDArrayView<Scalar,device_type>::HostView  scalar_array_h;
  typedef typename Kokkos::MDArrayView<int,device_type>::HostView     int_array_h;

  typedef typename Kokkos::MDArrayView<Scalar,device_type>            scalar_array_d;
  typedef typename Kokkos::MDArrayView<int,device_type>               int_array_d;

  typedef typename Kokkos::ValueView<Scalar,device_type>::HostView   scalar_h;
  typedef typename Kokkos::ValueView<Scalar,device_type>             scalar_d;

  const int NumStates = 2;

  double compute_time = 0.0;
  double total = 0.0;

  // element block parameters
  const Scalar  lin_bulk_visc = 0.06;
  const Scalar  quad_bulk_visc = 1.2;
  const Scalar  hg_stiffness = 0.0;
  const Scalar  hg_viscosity = 0.0;
  //const Scalar  hg_stiffness = 0.03;
  //const Scalar  hg_viscosity = 0.001;

  // material properties
  const Scalar youngs_modulus=1.0e6;
  const Scalar poissons_ratio=0.0;
  const Scalar  density = 8.0e-4;

  BoxMeshFixture<int_array_h, scalar_array_h> mesh(ex,ey,ez);
  Region<Scalar,device_type>  region( NumStates,
                                      mesh,
                                      lin_bulk_visc,
                                      quad_bulk_visc,
                                      hg_stiffness,
                                      hg_viscosity,
                                      youngs_modulus,
                                      poissons_ratio,
                                      density);

  /*
  scalar_array_h  displacement_h   =  Kokkos::create_mdarray< scalar_array_h >(region.num_nodes, 3, 2); // two state field
  scalar_array_h  acceleration_h   =  Kokkos::create_mdarray< scalar_array_h >(region.num_nodes, 3);
  scalar_array_h  internal_force_h =  Kokkos::create_mdarray< scalar_array_h >(region.num_nodes, 3);
  scalar_array_h  nodal_mass_h     =  Kokkos::create_mdarray< scalar_array_h >(region.num_nodes);
  scalar_array_h  elem_mass_h      =  Kokkos::create_mdarray< scalar_array_h >(region.num_elements);
  scalar_array_h  stress_new_h     =  Kokkos::create_mdarray< scalar_array_h >(region.num_elements,6);
  */

  scalar_array_h  velocity_h     =   Kokkos::create_mdarray< scalar_array_h >(region.num_nodes, 3, 2); // two state field

  //setup the initial condition on velocity
  {
    const unsigned X = 0;
    for (int inode = 0; inode< region.num_nodes; ++inode) {
      if ( region.model_coords(inode,X) == 0) {
        velocity_h(inode,X,0) = 1.0e3;
        velocity_h(inode,X,1) = 1.0e3;
      }
    }
  }

  Kokkos::deep_copy(region.velocity, velocity_h);

  // Parameters required for the internal force computations.
  const Scalar dt = 1.0e-6;
  const Scalar  end_time = 0.0050;


  //--------------------------------------------------------------------------
  // We will call a sequence of functions.  These functions have been
  // grouped into several functors to balance the number of global memory
  // accesses versus requiring too many registers or too much L1 cache.
  // Global memory accees have read/write cost and memory subsystem contention cost.
  //--------------------------------------------------------------------------

  Kokkos::parallel_for( region.num_elements,
      initialize_element<Scalar,device_type>(region)
      );

  Kokkos::parallel_for( region.num_nodes,
      initialize_node<Scalar,device_type>(region)
      );

  int current_state = 0;
  int previous_state = 0;
  int next_state = 0;
  //for (Scalar sim_time = 0.0; sim_time < end_time; sim_time += dt) {
  for (int sim_time = 0; sim_time < 3; ++sim_time) {

    //rotate the states
    previous_state = current_state;
    current_state = next_state;
    ++next_state;
    next_state %= NumStates;

    // First kernel 'grad_hgop' combines three functions:
    // gradient, velocity gradient, and hour glass operator.
    Kokkos::parallel_for( region.num_elements ,
        grad_hgop<Scalar, device_type> ( region,
                                         dt,
                                         current_state,
                                         previous_state
                                       )
        , compute_time );

    total += compute_time;

    // Combine tensor decomposition and rotation functions.
    Kokkos::parallel_for( region.num_elements ,
        decomp_rotate<Scalar, device_type> ( region,
                                             dt,
                                             current_state,
                                             previous_state
                                           )
        , compute_time );

    total += compute_time;

    // Single beastly function in this last functor,
    // did not notice any opportunity for splitting.
    Kokkos::parallel_for( region.num_elements ,
        divergence<Scalar, device_type> ( region,
                                          dt,
                                          current_state,
                                          previous_state
                                        )
        , compute_time );

    total += compute_time;

    // Assembly of elements' contributions to nodal force into
    // a nodal force vector.  Update the accelerations, velocities,
    // displacements.
    //
    // The same pattern can be used for matrix-free residual computations.
    Kokkos::parallel_for( region.num_nodes ,
        finish_step<Scalar, device_type>( region,
                                          ex+1,
                                          dt,
                                          current_state,
                                          next_state
                                        )
        , compute_time );


    total += compute_time;


    std::cout << "Time step = " << sim_time << std::endl << std::endl;

    std::cout << "Element Stress\n";
    for(int i = 0; i<region.num_elements; ++i) {
      std::cout << "(";
      std::cout << region.stress_new(i,0) << ",";
      std::cout << region.stress_new(i,1) << ",";
      std::cout << region.stress_new(i,2) << ",";
      std::cout << region.stress_new(i,3) << ",";
      std::cout << region.stress_new(i,4) << ",";
      std::cout << region.stress_new(i,5) << "), ";
    }
    std::cout << "\n\n";

    for (int inode = 0; inode<region.num_nodes; ++inode) {
      std::cout << "Node = " << inode << std::endl;

      std::cout << "nodal_mass = " << region.nodal_mass(inode) << std::endl;

      std::cout << "internal_force = (";
      std::cout << region.internal_force(inode,0) << ",";
      std::cout << region.internal_force(inode,1) << ",";
      std::cout << region.internal_force(inode,2) << ")" << std::endl;

      std::cout << "acceleration = (";
      std::cout << region.acceleration(inode,0) << ",";
      std::cout << region.acceleration(inode,1) << ",";
      std::cout << region.acceleration(inode,2) << ")" << std::endl;

      std::cout << "velocity (current) = (";
      std::cout << region.velocity(inode,0,current_state) << ",";
      std::cout << region.velocity(inode,1,current_state) << ",";
      std::cout << region.velocity(inode,2,current_state) << ")" << std::endl;

      std::cout << "displacement (current) = (";
      std::cout << region.displacement(inode,0,current_state) << ",";
      std::cout << region.displacement(inode,1,current_state) << ",";
      std::cout << region.displacement(inode,2,current_state) << ")" << std::endl << std::endl;
    }
  }

  return total;
}

#endif

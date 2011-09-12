#ifndef INTERNALFORCE
#define INTERNALFORCE

#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "grad_hgop.hpp"
#include "decomp_rotate.hpp"
#include "divergence.hpp"
#include "BoxMeshFixture.hpp"
#include "ForceGather.hpp"
#include "initialize_element.hpp"
#include "initialize_nodal_mass.hpp"
#include "compute_acceleration_velocity_displacement.hpp"

//----------------------------------------------------------------------------

template<typename Scalar, class device_type>
double internal_force_test( const size_t ex, const size_t ey, const size_t ez )
{
  typedef typename Kokkos::MDArrayView<Scalar,device_type>::HostView  scalar_array_h;
  typedef typename Kokkos::MDArrayView<int,device_type>::HostView     int_array_h;

  typedef typename Kokkos::MDArrayView<Scalar,device_type>            scalar_array_d;
  typedef typename Kokkos::MDArrayView<int,device_type>               int_array_d;

  typedef typename Kokkos::ValueView<Scalar,device_type>::HostView   scalar_h;
  typedef typename Kokkos::ValueView<Scalar,device_type>             scalar_d;

  double compute_time = 0.0;
  double total = 0.0;

  timeval start, stop, result;

  BoxMeshFixture<int_array_h, scalar_array_h> mesh(ex,ey,ez);

  const int NumStates = 2;

  const unsigned nelems = mesh.nelems;;
  const unsigned nnodes = mesh.nnodes;

  //  Initialize Host mesh data structures:
  int_array_h     elem_node_connectivity_h = mesh.elem_node_ids;
  int_array_h     node_elem_offset_h = mesh.node_elem_offset;
  int_array_h     node_elem_ids_h = mesh.node_elem_ids;
  scalar_array_h  model_coords_h = mesh.node_coords;

  scalar_array_h  velocity_h     =  Kokkos::create_mdarray< scalar_array_h >(nnodes, 3, 2); // two state field

  //setup the initial condition on velecity
  {
    const unsigned X = 0;
    for (unsigned inode = 0; inode< nnodes; ++inode) {
      if ( model_coords_h(inode,X) == 0) {
        velocity_h(inode,X,0) = 1.0e3;
        velocity_h(inode,X,1) = 1.0e3;
      }
    }
  }

  //arrays for output
  scalar_array_h  displacement_h   =  Kokkos::create_mdarray< scalar_array_h >(nnodes, 3, 2); // two state field
  scalar_array_h  acceleration_h   =  Kokkos::create_mdarray< scalar_array_h >(nnodes, 3);
  scalar_array_h  internal_force_h =  Kokkos::create_mdarray< scalar_array_h >(nnodes, 3);
  scalar_array_h  nodal_mass_h     =  Kokkos::create_mdarray< scalar_array_h >(nnodes);
  scalar_array_h  elem_mass_h      =  Kokkos::create_mdarray< scalar_array_h >(nelems);
  scalar_array_h  stress_new_h     =  Kokkos::create_mdarray< scalar_array_h >(nelems,6);

  //scalar_h current_time_step_h = Kokkos::create_value< scalar_h >();


  gettimeofday(&start, NULL);


  int_array_d elem_node_connectivity
    = Kokkos::create_mdarray< int_array_d >(  elem_node_connectivity_h.dimension(0),
                                              elem_node_connectivity_h.dimension(1));

  int_array_d node_elem_ids
    = Kokkos::create_mdarray< int_array_d >(  node_elem_ids_h.dimension(0),
                                              node_elem_ids_h.dimension(1));

  int_array_d node_elem_offset
    = Kokkos::create_mdarray< int_array_d >(  node_elem_offset_h.dimension(0),
                                              node_elem_offset_h.dimension(1));


  gettimeofday(&stop, NULL);
  timersub(&stop, &start, &result);
  compute_time = (result.tv_sec + result.tv_usec/1000000.0);

  // Field data required for the internal force computations on the device.

  scalar_array_d  model_coords =  Kokkos::create_mdarray< scalar_array_d >(nnodes, 3);
  scalar_array_d  displacement =  Kokkos::create_mdarray< scalar_array_d >(nnodes, 3, NumStates); // two state field
  scalar_array_d  velocity     =  Kokkos::create_mdarray< scalar_array_d >(nnodes, 3, NumStates); // two state field
  scalar_array_d  acceleration =  Kokkos::create_mdarray< scalar_array_d >(nnodes, 3);
  scalar_array_d  nodal_mass   =  Kokkos::create_mdarray< scalar_array_d >(nnodes);
  scalar_array_d  elem_mass    =  Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  stress_new =    Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  internal_force =   Kokkos::create_mdarray< scalar_array_d >(nnodes, 3);


  scalar_array_d  hg_resist =     Kokkos::create_mdarray< scalar_array_d >(nelems, 12, NumStates); // old and new
  scalar_array_d  rotation =      Kokkos::create_mdarray< scalar_array_d >(nelems, 9, NumStates);  // rotation old and new
  scalar_array_d  gradop12 =      Kokkos::create_mdarray< scalar_array_d >(nelems, 3, 8);
  scalar_array_d  force_new =     Kokkos::create_mdarray< scalar_array_d >(nelems, 3, 8);


  scalar_array_d  hgop =          Kokkos::create_mdarray< scalar_array_d >(nelems, 32, NumStates); // hgop and mid_hgop
  scalar_array_d  vel_grad =      Kokkos::create_mdarray< scalar_array_d >(nelems, 9);
  scalar_array_d  stretch =       Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  vorticity =     Kokkos::create_mdarray< scalar_array_d >(nelems, 3);
  scalar_array_d  rot_stret =     Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  rot_stress =    Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  mid_vol =       Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  shrmod =        Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  dilmod =        Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  elem_t_step =   Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  hg_energy =     Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  intern_energy = Kokkos::create_mdarray< scalar_array_d >(nelems);

  // how do we create in shared memory

  //scalar_h curr_dt_h = Kokkos::create_value< scalar_h >();

  //scalar_d prev_dt = Kokkos::create_value< scalar_d >();
  //scalar_d curr_dt = Kokkos::create_value< scalar_d >();
  //scalar_d next_dt = Kokkos::create_value< scalar_d >();

  // Parameters required for the internal force computations.
  const Scalar dt = 1.0e-6;
  const Scalar  end_time = 0.0050;

  // fudge factors
  const Scalar  lin_bulk_visc = 0.06;
  const Scalar  quad_bulk_visc = 1.2;
  const Scalar  hg_stiffness = 0.03;
  const Scalar  hg_viscosity = 0.001;

  // material properties
  const Scalar youngs_modulus=1.0e6;
  const Scalar poissons_ratio=0.0;

  const Scalar  two_mu = youngs_modulus/(1+poissons_ratio);
  const Scalar  bulk_modulus = youngs_modulus/(3*(1-2*poissons_ratio));

  const Scalar  density = 8.0e-4;

  //  before starting the internal force calculations,
  //  copy the Host generated mesh information to the accelerator device

  Kokkos::deep_copy(elem_node_connectivity, elem_node_connectivity_h);
  Kokkos::deep_copy(node_elem_ids, node_elem_ids_h);
  Kokkos::deep_copy(node_elem_offset, node_elem_offset_h);
  Kokkos::deep_copy(model_coords, model_coords_h);
  Kokkos::deep_copy(velocity, velocity_h);

  //--------------------------------------------------------------------------
  // We will call a sequence of functions.  These functions have been
  // grouped into several functors to balance the number of global memory
  // accesses versus requiring too many registers or too much L1 cache.
  // Global memory accees have read/write cost and memory subsystem contention cost.
  //--------------------------------------------------------------------------

  Kokkos::parallel_for( nelems,
      initialize_element<Scalar,device_type>(
        elem_node_connectivity,
        model_coords,
        elem_mass,
        stretch,
        rotation,
        density
        )
      );

  Kokkos::parallel_for( nnodes,
      initialize_nodal_mass<Scalar,device_type>(
        node_elem_offset,
        node_elem_ids,
        elem_mass,
        nodal_mass
        )
      );

  Kokkos::deep_copy(nodal_mass_h,nodal_mass);
  Kokkos::deep_copy(elem_mass_h,elem_mass);

#if 0
  std::cout << "Element Mass\n";
  for(unsigned i = 0; i<nelems; ++i) {
    std::cout << elem_mass_h(i) << ",";
    if ((i+1)%20 == 0) std::cout << "\n";
  }

  std::cout << "\n\n";

  std::cout << "Nodal Mass\n";
  for(unsigned i = 0; i<nnodes; ++i) {
    std::cout << nodal_mass_h(i) << ",";
    if ((i+1)%20 == 0) std::cout << "\n";
  }

  std::cout << "\n\n";

  std::cout << "Nodal Velocity\n";
  for(unsigned i = 0; i<nnodes; ++i) {
    std::cout << '(';
    std::cout << velocity_h(i,0,0) << ",";
    std::cout << velocity_h(i,1,0) << ",";
    std::cout << velocity_h(i,2,0) << "), ";
    if ((i+1)%10 == 0) std::cout << "\n";
  }

  std::cout << "\n\n";
#endif

  int current_state = 0;
  int previous_state = 1;
  //for (Scalar sim_time = 0.0; sim_time < end_time; sim_time += dt) {
  for (int sim_time = 0; sim_time < 20; ++sim_time) {

    //rotate the states
    previous_state = current_state;
    ++current_state;
    current_state %= NumStates;


    // First kernel 'grad_hgop' combines three functions:
    // gradient, velocity gradient, and hour glass operator.
    Kokkos::parallel_for( nelems , grad_hgop<Scalar, device_type>
        ( elem_node_connectivity,
          model_coords,
          displacement,
          velocity,
          mid_vol,
          vel_grad,
          hgop,
          dt,
          current_state,
          previous_state)     , compute_time );

    total += compute_time;

    // Combine tensor decomposition and rotation functions.

    Kokkos::parallel_for( nelems , decomp_rotate<Scalar, device_type>
        ( rotation,
          vel_grad,
          stretch,
          vorticity,
          rot_stret,
          dt,
          current_state,
          previous_state)     , compute_time );

    total += compute_time;

    // Single beastly function in this last functor,
    // did not notice any opportunity for splitting.

    Kokkos::parallel_for( nelems , divergence<Scalar, device_type>
        ( elem_node_connectivity,
          model_coords,
          displacement,
          velocity,
          force_new,
          vorticity,
          rotation,
          stress_new,
          rot_stress,
          rot_stret,
          gradop12,
          elem_mass,
          dilmod,
          shrmod,
          elem_t_step,
          intern_energy,
          mid_vol,
          hgop,
          hg_resist,
          hg_energy,
          two_mu,
          bulk_modulus,
          hg_stiffness,
          hg_viscosity,
          lin_bulk_visc,
          quad_bulk_visc,
          dt,
          current_state,
          previous_state)   , compute_time );

    total += compute_time;

    // Assembly of elements' contributions to nodal force into
    // a nodal force vector.  The same pattern can be used for
    // matrix-free residual computations.

    Kokkos::parallel_for( nnodes , ForceGather<Scalar, device_type>
        (  node_elem_ids,
           node_elem_offset,
           internal_force,
           force_new,
           current_state,
           previous_state) , compute_time);

    Kokkos::parallel_for( nnodes ,
        compute_acceleration_velocity_displacement<Scalar, device_type>
        (  internal_force,
           nodal_mass,
           acceleration,
           velocity,
           displacement,
           dt,
           current_state,
           previous_state) ,
        compute_time);

    total += compute_time;

  }

  Kokkos::deep_copy(stress_new_h,stress_new);

  std::cout << "Element Stress\n";
  for(unsigned i = 0; i<nelems; ++i) {
    std::cout << "(";
    std::cout << stress_new_h(i,0) << ",";
    std::cout << stress_new_h(i,1) << ",";
    std::cout << stress_new_h(i,2) << ",";
    std::cout << stress_new_h(i,3) << ",";
    std::cout << stress_new_h(i,4) << ",";
    std::cout << stress_new_h(i,5) << "), ";
    if ((i+1)%5 == 0) std::cout << "\n";
  }
  std::cout << "\n\n";

  return total;
}

#endif

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

//----------------------------------------------------------------------------

template<typename Scalar, class device_type>
double internal_force_test( const size_t ex, const size_t ey, const size_t ez )
{
  const int nelems = ex * ey * ez;
  const int nnodes = (ex + 1) * (ey + 1) * (ez + 1);

  double time = 0.0;
  double total = 0.0;

  timeval start, stop, result;

  // Creating the mesh on the Host, so use compatible MDArrays on the Host.
  // The 'HostView' member type is the compatible MDArray on the Host memory space.
  // Use compatible MDArray layouts to avoid the cost of transposing when copying
  // to the device memory.

  typedef typename Kokkos::MDArrayView<Scalar,device_type>           scalar_array_d;
  typedef typename Kokkos::MDArrayView<int,device_type>             int_array_d;

  //  Initialize Host data structures:
  HostView_int     elem_nodeIDs_h ;  // Element-node connectivity array

  typedef typename scalar_array_d::HostView  scalar_array_h;
  typedef typename int_array_d   ::HostView  int_array_h;

  gettimeofday(&start, NULL);

  const BoxMeshFixture< int_array_h , scalar_array_h > box_mesh( ex , ey , ez );

  // Element-nodes coordinates
  scalar_array_h elem_coords_h = Kokkos::create_mdarray< scalar_array_h >( nelems , 3 , 8 );

  // Create and populate the box-mesh connectivity and coordinates.
  // To be considered: to use a node coordinate array and index
  //                   into that array for element-node coordinates.
  mesh_init<HostView_scalar, HostView_int>(  elem_coords_h,
                        elem_nodeIDs_h,
                        node_elemIDs_h,
                        elems_per_node_h,
                        ex, ey, ez);

  //  Copies of the mesh data structures on the device.
  int_array_d   elem_nodeIDs, node_elemIDs, node_elemOffset;

  elem_nodeIDs = Kokkos::create_mdarray< int_array_d >(  elem_nodeIDs_h.dimension(0),
                                                         elem_nodeIDs_h.dimension(1));

  node_elemIDs = Kokkos::create_mdarray< int_array_d >(  node_elemIDs_h.dimension(0),
                                                         node_elemIDs_h.dimension(1));

  elems_per_node = Kokkos::create_mdarray< int_array_d >(  elems_per_node_h.dimension(0),
                                                           elems_per_node_h.dimension(1));


  gettimeofday(&stop, NULL);
  timersub(&stop, &start, &result);
  time = (result.tv_sec + result.tv_usec/1000000.0);

  // Field data required for the internal force computations on the device.

  scalar_array_d  position =      Kokkos::create_mdarray< scalar_array_d >(nelems, 3, 8);
  scalar_array_d  velocity =      Kokkos::create_mdarray< scalar_array_d >(nelems, 3, 8);
  scalar_array_d  hg_resist =     Kokkos::create_mdarray< scalar_array_d >(nelems, 12, 2); // old and new
  scalar_array_d  rotation =      Kokkos::create_mdarray< scalar_array_d >(nelems, 9, 2);  // rotation old and new
  scalar_array_d  gradop12 =      Kokkos::create_mdarray< scalar_array_d >(nelems, 3, 8);
  scalar_array_d  force_new =     Kokkos::create_mdarray< scalar_array_d >(nelems, 3, 8);

  scalar_array_d  nodal_force =   Kokkos::create_mdarray< scalar_array_d >(nnodes, 3);

  scalar_array_d  hgop =          Kokkos::create_mdarray< scalar_array_d >(nelems, 32, 2); // hgop and mid_hgop
  scalar_array_d  vel_grad =      Kokkos::create_mdarray< scalar_array_d >(nelems, 9);
  scalar_array_d  stretch =       Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  s_temp =        Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  vorticity =     Kokkos::create_mdarray< scalar_array_d >(nelems, 3);
  scalar_array_d  rot_stret =     Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  stress_new =    Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  rot_stress =    Kokkos::create_mdarray< scalar_array_d >(nelems, 6);
  scalar_array_d  mid_vol =       Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  shrmod =        Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  dilmod =        Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  elem_mass =     Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  elem_t_step =   Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  hg_energy =     Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  intern_energy = Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  bulk_modulus =  Kokkos::create_mdarray< scalar_array_d >(nelems);
  scalar_array_d  two_mu =        Kokkos::create_mdarray< scalar_array_d >(nelems);

  // Parameters required for the internal force computations.
  const Scalar  dt = 0.25;
  const Scalar  lin_bulk_visc = 0.06;
  const Scalar  quad_bulk_visc = 1.2;
  const Scalar  hg_stiffness = 0.03;
  const Scalar  hg_viscosity = 0.003;
  const bool    scaleHGRotation = false;

  //  before starting the internal force calculations,
  //  copy the Host generated mesh information to the accelerator device

  Kokkos::deep_copy(position,       elem_coords_h);
  Kokkos::deep_copy(elem_nodeIDs,   box_mesh.elem_node_ids);
  Kokkos::deep_copy(node_elemIDs,   box_mesh.node_elem_ids);
  Kokkos::deep_copy(node_elemOffset, box_mesh.node_elem_offset);

  //--------------------------------------------------------------------------
  // We will call a sequence of functions.  These functions have been
  // grouped into several functors to balance the number of global memory
  // accesses versus requiring too many registers or too much L1 cache.
  // Global memory accees have read/write cost and memory subsystem contention cost.
  //--------------------------------------------------------------------------

  // First kernel 'grad_hgop' combines three functions:
  // gradient, velocity gradient, and hour glass operator.

  Kokkos::parallel_for( nelems , grad_hgop<Scalar, device_type>
  (   position,
    velocity,
    mid_vol,
    vel_grad,
    hgop,
    dt )     , time );

  total += time;

  // Combine tensor decomposition and rotation functions.

  Kokkos::parallel_for( nelems , decomp_rotate<Scalar, device_type>
  (    rotation,
    vel_grad,
    s_temp,
    stretch,
    vorticity,
    rot_stret,
    dt)     , time );

  total += time;

  // Single beastly function in this last functor,
  // did not notice any opportunity for splitting.

  Kokkos::parallel_for( nelems , divergence<Scalar, device_type>
  (   position,
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
    scaleHGRotation)   , time );

  total += time;

  // Assembly of elements' contributions to nodal force into
  // a nodal force vector.  The same pattern can be used for
  // matrix-free residual computations.

  Kokkos::parallel_for( nnodes , ForceGather<Scalar, device_type>
  (  node_elemIDs,
    node_elemOffset,
    nodal_force,
    force_new ) , time);

  total += time;


  return total;
}

#endif

/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef INTERNALFORCE
#define INTERNALFORCE

#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <impl/KokkosArray_Timer.hpp>

#include <ExplicitBoxMeshFixture.hpp>
#include <region.hpp>
#include <initialize_element.hpp>
#include <initialize_node.hpp>
#include <grad_hgop.hpp>
#include <decomp_rotate.hpp>
#include <divergence.hpp>
//#include <minimum_stable_time_step.hpp>
#include <finish_step.hpp>

//----------------------------------------------------------------------------

namespace explicit_dynamics {

struct PerformanceData {
  double mesh_time ;
  double init_time ;
  double internal_force_time;
  double central_diff;
  double copy_to_host_time;
  size_t number_of_steps ;

  PerformanceData()
  : mesh_time(0)
  , init_time(0)
  , internal_force_time(0)
  , central_diff(0)
  , copy_to_host_time(0)
  , number_of_steps(0)
  {}

  void best( const PerformanceData & rhs )
  {
    if ( rhs.mesh_time < mesh_time ) mesh_time = rhs.mesh_time ;
    if ( rhs.init_time < init_time ) init_time = rhs.init_time ;
    if ( rhs.internal_force_time < internal_force_time ) internal_force_time = rhs.internal_force_time ;
    if ( rhs.central_diff < central_diff ) central_diff = rhs.central_diff ;
    if ( rhs.copy_to_host_time < copy_to_host_time ) copy_to_host_time = rhs.copy_to_host_time ;
  }
};

template<typename Scalar, class device_type>
void explicit_dynamics_app( const size_t ex,
                            const size_t ey,
                            const size_t ez,
                            const size_t steps ,
                            PerformanceData & perf )
{
  typedef typename KokkosArray::MDArray<Scalar,device_type>::HostMirror  scalar_array_h;
  typedef typename KokkosArray::MDArray<int,device_type>::HostMirror     int_array_h;

  typedef typename KokkosArray::MDArray<Scalar,device_type>            scalar_array_d;
  typedef typename KokkosArray::MDArray<int,device_type>               int_array_d;

  typedef typename KokkosArray::Value<Scalar,device_type>::HostMirror   scalar_h;
  typedef typename KokkosArray::Value<Scalar,device_type>             scalar_d;

  const int NumStates = 2;

  const Scalar user_dt = 1.0e-5;
  //const Scalar  end_time = 0.0050;

  // element block parameters
  const Scalar  lin_bulk_visc = 0.0;
  const Scalar  quad_bulk_visc = 0.0;
  //const Scalar  lin_bulk_visc = 0.06;
  //const Scalar  quad_bulk_visc = 1.2;
  const Scalar  hg_stiffness = 0.0;
  const Scalar  hg_viscosity = 0.0;
  //const Scalar  hg_stiffness = 0.03;
  //const Scalar  hg_viscosity = 0.001;

  // material properties
  const Scalar youngs_modulus=1.0e6;
  const Scalar poissons_ratio=0.0;
  const Scalar  density = 8.0e-4;

  KokkosArray::Impl::Timer wall_clock ;

  BoxMeshFixture<int_array_h, scalar_array_h> mesh(ex,ey,ez);

  scalar_array_h  nodal_mass_h     =  KokkosArray::create_mdarray< scalar_array_h >(mesh.nnodes);
  scalar_array_h  elem_mass_h      =  KokkosArray::create_mdarray< scalar_array_h >(mesh.nelems);


  scalar_array_h  acceleration_h   =  KokkosArray::create_mdarray< scalar_array_h >(mesh.nnodes, 3);
  scalar_array_h  velocity_h     =   KokkosArray::create_mdarray< scalar_array_h >(mesh.nnodes, 3, 2); // two state field
  scalar_array_h  displacement_h   =  KokkosArray::create_mdarray< scalar_array_h >(mesh.nnodes, 3, 2); // two state field
  scalar_array_h  internal_force_h =  KokkosArray::create_mdarray< scalar_array_h >(mesh.nnodes, 3);
  scalar_array_h  stress_new_h     =  KokkosArray::create_mdarray< scalar_array_h >(mesh.nelems,6);


  //setup the initial condition on velocity
  {
    const unsigned X = 0;
    for (int inode = 0; inode< mesh.nnodes; ++inode) {
      if ( mesh.node_coords(inode,X) == 0) {
        velocity_h(inode,X,0) = 1.0e3;
        velocity_h(inode,X,1) = 1.0e3;
      }
    }
  }

  Region<Scalar,device_type>  region( NumStates,
                                      mesh,
                                      lin_bulk_visc,
                                      quad_bulk_visc,
                                      hg_stiffness,
                                      hg_viscosity,
                                      youngs_modulus,
                                      poissons_ratio,
                                      density);

  KokkosArray::deep_copy(region.velocity, velocity_h);


  perf.mesh_time = wall_clock.seconds(); // Mesh and graph allocation and population.
  wall_clock.reset();



  // Parameters required for the internal force computations.


  //--------------------------------------------------------------------------
  // We will call a sequence of functions.  These functions have been
  // grouped into several functors to balance the number of global memory
  // accesses versus requiring too many registers or too much L1 cache.
  // Global memory accees have read/write cost and memory subsystem contention cost.
  //--------------------------------------------------------------------------

  KokkosArray::parallel_for( region.num_elements,
      initialize_element<Scalar,device_type>(region)
      );

  KokkosArray::parallel_for( region.num_nodes,
      initialize_node<Scalar,device_type>(region)
      );

  perf.init_time = wall_clock.seconds(); // Initialization
  wall_clock.reset();

  int current_state = 0;
  int previous_state = 0;
  int next_state = 0;

  const int total_num_steps = steps ;

  perf.number_of_steps = total_num_steps ;

  for (int step = 0; step < total_num_steps; ++step) {

    //rotate the states
    previous_state = current_state;
    current_state = next_state;
    ++next_state;
    next_state %= NumStates;

    wall_clock.reset();

    // First kernel 'grad_hgop' combines three functions:
    // gradient, velocity gradient, and hour glass operator.
    KokkosArray::parallel_for( region.num_elements ,
        grad_hgop<Scalar, device_type> ( region,
                                         current_state,
                                         previous_state
                                       ));

    // Combine tensor decomposition and rotation functions.
    KokkosArray::parallel_for( region.num_elements ,
        decomp_rotate<Scalar, device_type> ( region,
                                             current_state,
                                             previous_state
                                           ));


    // Single beastly function in this last functor,
    // did not notice any opportunity for splitting.
    KokkosArray::parallel_reduce( region.num_elements ,
        divergence<Scalar, device_type> ( region,
                                          user_dt,
                                          current_state,
                                          previous_state
                                        ),
        set_next_time_step<Scalar,device_type>(region));

    device_type::fence();

    perf.internal_force_time += wall_clock.seconds();
    wall_clock.reset();


    // Assembly of elements' contributions to nodal force into
    // a nodal force vector.  Update the accelerations, velocities,
    // displacements.
    // The same pattern can be used for matrix-free residual computations.
    KokkosArray::parallel_for( region.num_nodes ,
        finish_step<Scalar, device_type>( region,
                                          ex,
                                          current_state,
                                          next_state
                                        ));

    device_type::fence();
    perf.central_diff += wall_clock.seconds();
    wall_clock.reset();

#ifdef KOKKOSARRAY_DEVICE_CUDA
    if (step%100 == 0 ) {
      KokkosArray::deep_copy(acceleration_h,region.acceleration);
      KokkosArray::deep_copy(velocity_h,region.velocity);
      KokkosArray::deep_copy(displacement_h,region.displacement);
      KokkosArray::deep_copy(internal_force_h,region.internal_force);
      KokkosArray::deep_copy(stress_new_h,region.stress_new);
    }
#endif

    device_type::fence();
    perf.copy_to_host_time += wall_clock.seconds();
    wall_clock.reset();


  }
}


template <typename Scalar, typename Device>
static void driver( const char * label , int beg , int end , int runs )
{
  int shift = 20;

  std::cout << std::endl ;
  std::cout << "\"MiniExplicitDynamics with KokkosArray " << label << "\"" << std::endl;
  std::cout << std::left << std::setw(shift) << "\"Size\" , ";
  std::cout << std::left << std::setw(shift) << "\"Time Steps\" , ";
  std::cout << std::left << std::setw(shift) << "\"Setup\" , ";
  std::cout << std::left << std::setw(shift) << "\"Initialize\" , ";
  std::cout << std::left << std::setw(shift) << "\"InternalForce\" , ";
  std::cout << std::left << std::setw(shift) << "\"CentralDiff\" , ";
  std::cout << std::left << std::setw(shift) << "\"CopyToHost\" , ";
  std::cout << std::left << std::setw(shift) << "\"TimePerElement\"";

  std::cout << std::endl;

  std::cout << std::left << std::setw(shift) << "\"elements\" , ";
  std::cout << std::left << std::setw(shift) << "\"iterations\" , ";
  std::cout << std::left << std::setw(shift) << "\"microsec\" , ";
  std::cout << std::left << std::setw(shift) << "\"microsec\" , ";
  std::cout << std::left << std::setw(shift) << "\"microsec\" , ";
  std::cout << std::left << std::setw(shift) << "\"microsec\" , ";
  std::cout << std::left << std::setw(shift) << "\"microsec\" , ";
  std::cout << std::left << std::setw(shift) << "\"microsec/element\"";

  std::cout << std::endl;

  const int steps = 1000 ;

  for(int i = beg ; i < end; ++i )
  {
    int two_to_the_i = 1 << i;
    int factor = static_cast<int>(cbrt(static_cast<double>(two_to_the_i)));

    const int ix = 10 * factor;
    const int iy = factor;
    const int iz = factor;
    const int n  = ix * iy * iz ;

    PerformanceData perf , best ;

    for(int j = 0; j < runs; j++){

     explicit_dynamics_app<Scalar,Device>(ix,iy,iz,steps,perf);

     if( j == 0 ) {
       best = perf ;
     }
     else {
       best.best( perf );
     }
   }

   double time_per_element = (best.internal_force_time + best.central_diff)/
           ( n * perf.number_of_steps );

   std::cout << std::setw(shift-3) << n << " , "
             << std::setw(shift-3) << best.number_of_steps << " , "
             << std::setw(shift-3) << best.mesh_time * 1000000 << " , "
             << std::setw(shift-3) << best.init_time * 1000000 << " , "
             << std::setw(shift-3) << ( best.internal_force_time * 1000000 ) / best.number_of_steps << " , "
             << std::setw(shift-3) << ( best.central_diff * 1000000 ) / best.number_of_steps << " , "
             << std::setw(shift-3) << best.copy_to_host_time * 1000000 << " , "
             << std::setw(shift) << time_per_element * 1000000
             << std::endl ;
  }
}


} // namespace explicit_dynamics

#endif

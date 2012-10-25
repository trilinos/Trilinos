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

#ifndef EXPLICIT_DRIVER_HPP
#define EXPLICIT_DRIVER_HPP

#include <sys/time.h>
#include <iostream>
#include <iomanip>

#include <impl/KokkosArray_Timer.hpp>

#include <ParallelDataMap.hpp>
#include <ExplicitFunctors.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Explicit {

struct PerformanceData {
  double mesh_time ;
  double init_time ;
  double internal_force_time ;
  double central_diff ;
  double comm_time ;
  size_t number_of_steps ;

  PerformanceData()
  : mesh_time(0)
  , init_time(0)
  , internal_force_time(0)
  , central_diff(0)
  , comm_time(0)
  , number_of_steps(0)
  {}

  void best( const PerformanceData & rhs )
  {
    if ( rhs.mesh_time < mesh_time ) mesh_time = rhs.mesh_time ;
    if ( rhs.init_time < init_time ) init_time = rhs.init_time ;
    if ( rhs.internal_force_time < internal_force_time ) internal_force_time = rhs.internal_force_time ;
    if ( rhs.central_diff < central_diff ) central_diff = rhs.central_diff ;
    if ( rhs.comm_time < comm_time ) comm_time = rhs.comm_time ;
  }
};

template< typename Scalar , class FixtureType >
PerformanceData run( const typename FixtureType::FEMeshType & mesh ,
                     const int global_max_x ,
                     const int global_max_y ,
                     const int global_max_z ,
                     const unsigned uq_count ,
                     const int steps ,
                     const int print_sample )
{
  typedef Scalar                              scalar_type ;
  typedef FixtureType                         fixture_type ;
  typedef typename fixture_type::device_type  device_type ;

  enum { ElementNodeCount = fixture_type::element_node_count };

  const int total_num_steps = steps ;

  const Scalar user_dt = 5.0e-6;
  //const Scalar  end_time = 0.0050;

  // element block parameters
  const Scalar  lin_bulk_visc = 0.0;
  const Scalar  quad_bulk_visc = 0.0;

  // const Scalar  lin_bulk_visc = 0.06;
  // const Scalar  quad_bulk_visc = 1.2;
  // const Scalar  hg_stiffness = 0.0;
  // const Scalar  hg_viscosity = 0.0;
  // const Scalar  hg_stiffness = 0.03;
  // const Scalar  hg_viscosity = 0.001;

  // material properties
  const Scalar youngs_modulus=1.0e6;
  const Scalar poissons_ratio=0.0;
  const Scalar  density = 8.0e-4;

  const comm::Machine machine = mesh.parallel_data_map.machine ;

  PerformanceData perf_data ;

  KokkosArray::Impl::Timer wall_clock ;

  //------------------------------------
  // Generate fields

  typedef Fields< scalar_type , device_type > fields_type ;

  fields_type mesh_fields( mesh , uq_count ,
                           lin_bulk_visc ,
                           quad_bulk_visc ,
                           youngs_modulus ,
                           poissons_ratio ,
                           density );

  typename fields_type::node_coords_type::HostMirror
    model_coords_h = KokkosArray::create_mirror( mesh_fields.model_coords );

  typename fields_type::geom_array_view::HostMirror
    displacement_h = KokkosArray::create_mirror( mesh_fields.displacement );

  typename fields_type::geom_array_view::HostMirror
    velocity_h = KokkosArray::create_mirror( mesh_fields.velocity );

  KokkosArray::deep_copy( model_coords_h , mesh_fields.model_coords );

  //------------------------------------
  // Initialization

  initialize_element( mesh_fields );
  initialize_node(    mesh_fields );

  const Scalar x_bc = global_max_x ;

  // Initial condition on velocity to initiate a pulse along the X axis
  {
    const unsigned X = 0;
    for ( unsigned inode = 0; inode< mesh_fields.num_nodes; ++inode) {
      for ( unsigned kq = 0 ; kq < uq_count ; ++kq ) {
        if ( model_coords_h(inode,X) == 0 ) {
          velocity_h(inode,kq,X) = 1000 + 100 * kq ;
          velocity_h(inode,kq,X) = 1000 + 100 * kq ;
        }
      }
    }
  }

  KokkosArray::deep_copy( mesh_fields.velocity , velocity_h );
  KokkosArray::deep_copy( mesh_fields.velocity_new , velocity_h );

  //--------------------------------------------------------------------------
  // We will call a sequence of functions.  These functions have been
  // grouped into several functors to balance the number of global memory
  // accesses versus requiring too many registers or too much L1 cache.
  // Global memory accees have read/write cost and memory subsystem contention cost.
  //--------------------------------------------------------------------------

  perf_data.init_time = comm::max( machine , wall_clock.seconds() );

  // Parameters required for the internal force computations.

  perf_data.number_of_steps = total_num_steps ;

  typedef typename
    fields_type::geom_array_view::scalar_type  comm_value_type ;

  const unsigned comm_value_count = 6 ;

  KokkosArray::AsyncExchange< comm_value_type , device_type ,
                              KokkosArray::ParallelDataMap >
    comm_exchange( mesh.parallel_data_map , comm_value_count );

  for ( int step = 0; step < total_num_steps; ++step ) {

    //------------------------------------------------------------------------
    // rotate the state variable views.

    swap( mesh_fields.dt ,           mesh_fields.dt_new );
    swap( mesh_fields.displacement , mesh_fields.displacement_new );
    swap( mesh_fields.velocity ,     mesh_fields.velocity_new );
    swap( mesh_fields.rotation ,     mesh_fields.rotation_new );

    //------------------------------------------------------------------------
    // Communicate "send" nodes' displacement and velocity next_state
    // to the ghosted nodes.
    // buffer packages: { { dx , dy , dz , vx , vy , vz }_node }

    wall_clock.reset();

    pack_state( mesh_fields ,
                comm_exchange.buffer(),
                mesh.parallel_data_map.count_interior ,
                mesh.parallel_data_map.count_send );

    comm_exchange.setup();

    comm_exchange.send_receive();

    unpack_state( mesh_fields ,
                  comm_exchange.buffer() ,
                  mesh.parallel_data_map.count_owned ,
                  mesh.parallel_data_map.count_receive );

    device_type::fence();

    perf_data.comm_time += comm::max( machine , wall_clock.seconds() );

    //------------------------------------------------------------------------

    wall_clock.reset();

    // First kernel 'grad_hgop' combines two functions:
    // gradient, velocity gradient
    gradient( mesh_fields );

    // Combine tensor decomposition and rotation functions.
    decomp_rotate( mesh_fields );

    internal_force( mesh_fields , user_dt );

    device_type::fence();

    perf_data.internal_force_time +=
      comm::max( machine , wall_clock.seconds() );

    //------------------------------------------------------------------------
    // Assembly of elements' contributions to nodal force into
    // a nodal force vector.  Update the accelerations, velocities,
    // displacements.
    // The same pattern can be used for matrix-free residual computations.

    wall_clock.reset();

    nodal_update( mesh_fields , x_bc );

    device_type::fence();

    perf_data.central_diff +=
      comm::max( machine , wall_clock.seconds() );

    if ( print_sample && 0 == step % 100 ) {
      KokkosArray::deep_copy( displacement_h , mesh_fields.displacement_new );
      KokkosArray::deep_copy( velocity_h ,     mesh_fields.velocity_new );

      if ( 1 == print_sample ) {
        for ( unsigned kp = 0 ; kp < uq_count ; ++kp ) {
          std::cout << "step " << step
                    << " : displacement({*,0,0}," << kp << ",0) =" ;
          for ( unsigned i = 0 ; i < mesh_fields.num_nodes_owned ; ++i ) {
            if ( model_coords_h(i,1) == 0 && model_coords_h(i,2) == 0 ) {
                std::cout << " " << displacement_h(i,kp,0);
            }
          }
          std::cout << std::endl ;

          const float tol = 1.0e-6 ;
          const int yb = global_max_y ;
          const int zb = global_max_z ;
          std::cout << "step " << step
                    << " : displacement({*," << yb << "," << zb << "}," << kp << ",0) =" ;
          for ( unsigned i = 0 ; i < mesh_fields.num_nodes_owned ; ++i ) {
            if ( fabs( model_coords_h(i,1) - yb ) < tol &&
                 fabs( model_coords_h(i,2) - zb ) < tol ) {
              std::cout << " " << displacement_h(i,kp,0);
            }
          }
          std::cout << std::endl ;
        }
      }
      else if ( 2 == print_sample ) {
        const unsigned kp = 0 ;

        const float tol = 1.0e-6 ;
        const int xb = global_max_x / 2 ;
        const int yb = global_max_y / 2 ;
        const int zb = global_max_z / 2 ;

        for ( unsigned i = 0 ; i < mesh_fields.num_nodes_owned ; ++i ) {
          if ( fabs( model_coords_h(i,0) - xb ) < tol &&
               fabs( model_coords_h(i,1) - yb ) < tol &&
               fabs( model_coords_h(i,2) - zb ) < tol ) {
            std::cout << "step " << step
                      << " : displacement("
                      << xb << "," << yb << "," << zb << ") = {" 
                      << std::setprecision(6)
                      << " " << displacement_h(i,kp,0)
                      << std::setprecision(2)
                      << " " << displacement_h(i,kp,1)
                      << std::setprecision(2)
                      << " " << displacement_h(i,kp,2)
                      << " }" << std::endl ;
          }
        }
      }
    }
  }

  return perf_data ;
}


template <typename Scalar, typename Device>
static void driver( const char * label , comm::Machine machine ,
                    const size_t elem_count_beg , const size_t elem_count_end ,
                    const size_t uq_count_beg ,   const size_t uq_count_end ,
                    const size_t run_count )
{
  typedef Scalar              scalar_type ;
  typedef Device              device_type ;
  typedef double              coordinate_scalar_type ;
  typedef FixtureElementHex8  fixture_element_type ;

  typedef BoxMeshFixture< coordinate_scalar_type ,
                          device_type ,
                          fixture_element_type > fixture_type ;

  typedef typename fixture_type::FEMeshType mesh_type ;


  const int space = 16 ;

  if ( comm::rank( machine ) == 0 ) {

  std::cout << std::endl ;
  std::cout << "\"MiniExplicitDynamics with KokkosArray " << label
            << "\"" << std::endl;
  std::cout << std::left << std::setw(space) << "\"Element\" , ";
  std::cout << std::left << std::setw(space) << "\"UQ-Count\" , " ;
  std::cout << std::left << std::setw(space) << "\"Time Steps\" , ";
  std::cout << std::left << std::setw(space) << "\"Initialize\" , ";
  std::cout << std::left << std::setw(space) << "\"ElemForce\" , ";
  std::cout << std::left << std::setw(space) << "\"NodeUpdate\" , ";
  std::cout << std::left << std::setw(space) << "\"NodeComm\" , ";
  std::cout << std::left << std::setw(space) << "\"Time/Elem\" , ";
  std::cout << std::left << std::setw(space) << "\"Time/Node\"";

  std::cout << std::endl;

  std::cout << std::left << std::setw(space) << "\"count\" , ";
  std::cout << std::left << std::setw(space) << "\"count\" , ";
  std::cout << std::left << std::setw(space) << "\"iterations\" , ";
  std::cout << std::left << std::setw(space) << "\"microsec\" , ";
  std::cout << std::left << std::setw(space) << "\"microsec\" , ";
  std::cout << std::left << std::setw(space) << "\"microsec\" , ";
  std::cout << std::left << std::setw(space) << "\"microsec\" , ";
  std::cout << std::left << std::setw(space) << "\"microsec\" , ";
  std::cout << std::left << std::setw(space) << "\"microsec\"";

  std::cout << std::endl;

  }

  const int steps = 1000 ;
  const int print_sample = 0 ;

  for( size_t i = elem_count_beg ; i < elem_count_end ; i *= 2 )
  {
    const size_t iz = std::max( (size_t) 1 , (size_t) cbrt( ((double) i) / 2.0 ) );
    const size_t iy = iz + 1 ;
    const size_t ix = 2 * iy ;
    const size_t nelem = ix * iy * iz ;
    const size_t nnode = ( ix + 1 ) * ( iy + 1 ) * ( iz + 1 );

    mesh_type mesh =
      fixture_type::create( comm::size( machine ) , comm::rank( machine ) , ix , iy , iz );

    mesh.parallel_data_map.machine = machine ;

    PerformanceData perf , best ;

    for ( size_t iuq = uq_count_beg ; iuq < uq_count_end ; iuq *= 2 ) {

      for( size_t j = 0; j < run_count ; j++ ) {

       perf = run<Scalar,fixture_type>(mesh,ix,iy,iz,iuq,steps,print_sample);

       if( j == 0 ) {
         best = perf ;
       }
       else {
         best.best( perf );
       }
     }

     if ( comm::rank( machine ) == 0 ) {
       double time_per_element =
         ( best.internal_force_time ) / ( nelem * perf.number_of_steps );
       double time_per_node =
         ( best.comm_time + best.central_diff ) / ( nnode * perf.number_of_steps );

       std::cout << std::setw(space-3) << nelem << " , "
                 << std::setw(space-3) << iuq << " , "
                 << std::setw(space-3) << best.number_of_steps << " , "
                 << std::setw(space-3) << best.init_time * 1000000 << " , "
                 << std::setw(space-3)
                 << ( best.internal_force_time * 1000000 ) / best.number_of_steps << " , "
                 << std::setw(space-3)
                 << ( best.central_diff * 1000000 ) / best.number_of_steps << " , "
                 << std::setw(space-3)
                 << ( best.comm_time * 1000000 ) / best.number_of_steps << " , "
                 << std::setw(space-3) << time_per_element * 1000000 << " , "
                 << std::setw(space-3) << time_per_node * 1000000
                 << std::endl ;
      }
    }
  }
}


} // namespace Explicit

#endif /* #ifndef EXPLICIT_DRIVER_HPP */


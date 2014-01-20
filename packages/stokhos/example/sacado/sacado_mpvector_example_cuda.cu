// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "KokkosCore_config.h"
#ifdef KOKKOS_HAVE_CUDA

#include "stdio.h"
#include "Stokhos_Sacado_Kokkos.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "sacado_mpvector_example.hpp"

#include "Kokkos_Cuda.hpp"

// Partial specialization of vector example runner for CUDA
template <int MaxSize, typename Scalar>
struct MPVectorExample<MaxSize, Scalar, Kokkos::Cuda> {
  typedef Kokkos::Cuda Device;

  static bool run(Storage_Method storage_method,
                  int num_elements, int num_samples,
                  int team_size, int league_size,
                  bool reset, bool print) {
    typedef MPVectorTypes<MaxSize, Scalar, Device> MPT;

    // Initialize Cuda
    const bool init_cuda_mirror = ! Kokkos::Cuda::host_mirror_device_type::is_initialized();

   if ( init_cuda_mirror ) { Kokkos::Cuda::host_mirror_device_type::initialize(1); }
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
    Kokkos::Cuda::print_configuration( std::cout );

    // Setup work request
    if (team_size == -1)
      team_size = num_samples;
    Kokkos::ParallelWorkRequest config(num_elements, team_size);

    int local_vector_size = num_samples / team_size;
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   num_samples % team_size != 0, std::logic_error,
    //   "Number of samples (" << num_samples ") is not a multiple of "
    //   << "thread team size (" << team_size << ")!");

    bool status = false;
    std::string device_name = "Cuda";


    switch ( local_vector_size ) {
    case  1 : status = run_view_kernel<Scalar, 1,Device>( config, num_elements, num_samples, reset, print, device_name); break ;
    case  2 : status = run_view_kernel<Scalar, 2,Device>( config, num_elements, num_samples, reset, print, device_name); break ;
    case  4 : status = run_view_kernel<Scalar, 4,Device>( config, num_elements, num_samples, reset, print, device_name); break ;
    case  8 : status = run_view_kernel<Scalar, 8,Device>( config, num_elements, num_samples, reset, print, device_name); break ;
    case 16 : status = run_view_kernel<Scalar,16,Device>( config, num_elements, num_samples, reset, print, device_name); break ;
    default :
      std::cout << "Local_vector_size (" << local_vector_size << ") NOT IMPLEMENTED" << std::endl ;
      return false ;
    }

    if ( status ) {

    status = run_view_kernel<Scalar, 0,Device>( config, num_elements, num_samples, reset, print, device_name);

    if ( status ) {

    if (storage_method == STATIC)
      status = run_kernels<Scalar, typename MPT::static_vector, typename MPT::static_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    else if (storage_method == STATIC_FIXED) {
      if (local_vector_size == 1)
        status = run_kernels<Scalar, typename MPT::static_fixed_vector_1, typename MPT::static_fixed_vector_1, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 2)
        status = run_kernels<Scalar, typename MPT::static_fixed_vector_2, typename MPT::static_fixed_vector_2, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 4)
        status = run_kernels<Scalar, typename MPT::static_fixed_vector_4, typename MPT::static_fixed_vector_4, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 8)
        status = run_kernels<Scalar, typename MPT::static_fixed_vector_8, typename MPT::static_fixed_vector_8, Device>(config, num_elements, num_samples, reset, print, device_name);
      else {
        std::cout <<  "Storage STATIC Invalid local vector size (" << local_vector_size
                  << ")!" << std::endl;
        status = false;
      }
    }
    else if (storage_method == LOCAL) {
      if (local_vector_size == 1)
        status = run_kernels<Scalar, typename MPT::local_vector_1, typename MPT::local_vector_1, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 2)
        status = run_kernels<Scalar, typename MPT::local_vector_2, typename MPT::local_vector_2, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 4)
        status = run_kernels<Scalar, typename MPT::local_vector_4, typename MPT::local_vector_4, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 8)
        status = run_kernels<Scalar, typename MPT::local_vector_8, typename MPT::local_vector_8, Device>(config, num_elements, num_samples, reset, print, device_name);
      else {
        std::cout <<  "Storage LOCAL Invalid local vector size (" << local_vector_size
                  << ")!" << std::endl;
        status = false;
      }
    }
    else if (storage_method == DYNAMIC)
      status = run_kernels<Scalar, typename MPT::dynamic_vector, typename MPT::dynamic_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    else if (storage_method == DYNAMIC_STRIDED)
      status = run_kernels<Scalar, typename MPT::dynamic_strided_vector, typename MPT::dynamic_strided_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    else if (storage_method == DYNAMIC_THREADED)
      status = run_kernels<Scalar, typename MPT::dynamic_threaded_vector, typename MPT::dynamic_threaded_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    if (storage_method == VIEW_STATIC)
      status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::static_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    else if (storage_method == VIEW_STATIC_FIXED) {
      if (local_vector_size == 1)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::static_fixed_vector_1, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 2)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::static_fixed_vector_2, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 4)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::static_fixed_vector_4, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 8)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::static_fixed_vector_8, Device>(config, num_elements, num_samples, reset, print, device_name);
      else {
        std::cout <<  "Storage VIEW_STATIC_FIXED Invalid local vector size (" << local_vector_size
                  << ")!" << std::endl;
        status = false;
      }
    }
    else if (storage_method == VIEW_LOCAL) {
      if (local_vector_size == 1)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::local_vector_1, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 2)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::local_vector_2, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 4)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::local_vector_4, Device>(config, num_elements, num_samples, reset, print, device_name);
      else if (local_vector_size == 8)
        status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::local_vector_8, Device>(config, num_elements, num_samples, reset, print, device_name);
      else {
        std::cout <<  "Storage VIEW_LOCAL Invalid local vector size (" << local_vector_size
                  << ")!" << std::endl;
        status = false;
      }
    }
    else if (storage_method == VIEW_DYNAMIC)
      status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::dynamic_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    else if (storage_method == VIEW_DYNAMIC_STRIDED)
      status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::dynamic_strided_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    else if (storage_method == VIEW_DYNAMIC_THREADED) {
      status = run_kernels<Scalar, typename MPT::view_vector, typename MPT::dynamic_threaded_vector, Device>(
        config, num_elements, num_samples, reset, print, device_name);
    }

    }}

    // Finalize Cuda
    Kokkos::Cuda::finalize();
    if ( init_cuda_mirror ) { Kokkos::Cuda::host_mirror_device_type::finalize(); }

    return status;
  }

};

template struct MPVectorExample<MaxSize, Scalar, Kokkos::Cuda>;

#endif

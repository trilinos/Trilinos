// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "KokkosCore_config.h"
#ifdef KOKKOS_HAVE_PTHREAD

#include "Stokhos_Sacado_Kokkos.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "sacado_mpvector_example.hpp"
#include <Kokkos_Threads.hpp>

// Partial specialization of vector example runner for Threads
template <int MaxSize, typename Scalar>
struct MPVectorExample<MaxSize, Scalar, Kokkos::Threads> {
  typedef Kokkos::Threads Device;

  static bool run(Storage_Method storage_method,
                  int num_elements, int num_samples,
                  int team_size, int league_size,
                  bool reset, bool print) {
    typedef MPVectorTypes<MaxSize, Scalar, Device> MPT;

    Kokkos::Threads::print_configuration( std::cout );

    // Setup work request
    Kokkos::ParallelWorkRequest config(num_elements, team_size);

    if ( num_samples % team_size ) {
      std::cout << "Number of samples (" << num_samples
                << ") must be a multiple of thread team size (" << team_size
                << ")" << std::endl ;
      return false ;
    }

    const int local_vector_size = num_samples / team_size ;
    bool status = false;
    std::string device_name = "Threads";

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

    if ( ! status ) return status ;

    status = run_view_kernel<Scalar, 0,Device>( config, num_elements, num_samples, reset, print, device_name);

    if ( ! status ) return status ;

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
        std::cout <<  "Storage STATIC_FIXED Invalid local vector size (" << local_vector_size
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
    else if (storage_method == DYNAMIC_THREADED) {
      std::cout << "Threads device doesn't support dynamic-threaded storage!"
                << std::endl;
      status = false;
    }
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
        std::cout <<  "Storage VIEW_STATIC Invalid local vector size (" << local_vector_size
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
      std::cout << "Threads device doesn't support dynamic-threaded storage!"
                << std::endl;
      status = false;
    }

    return status;
  }

};

template struct MPVectorExample<MaxSize, Scalar, Kokkos::Threads>;

#endif

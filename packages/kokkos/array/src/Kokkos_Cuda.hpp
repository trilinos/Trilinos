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

#ifndef KOKKOS_CUDA_HPP
#define KOKKOS_CUDA_HPP

#include <impl/Kokkos_IndexMap.hpp>

#define KOKKOS_CUDA  Kokkos::Cuda

/*--------------------------------------------------------------------------*/

namespace Kokkos {

class Cuda {
public:
  /*--------------------------------*/
  /* required type declarations and functions for a device */
  
  typedef Cuda          memory_space ;
  typedef unsigned int  size_type ;

  template< unsigned Rank = 0 , unsigned N1 = 0 ,
            unsigned N2   = 0 , unsigned N3 = 0 ,
            unsigned N4   = 0 , unsigned N5 = 0 , 
            unsigned N6   = 0 , unsigned N7 = 0 >
  struct IndexMap {
    typedef Impl::IndexMapLeft<memory_space,Rank,N1,N2,N3,N4,N5,N6,N7> type ;
  };
  
  /*--------------------------------*/

  struct SelectDevice {
    const int cuda_device_id ;
    SelectDevice() : cuda_device_id(0) {}
    explicit SelectDevice( int id ) : cuda_device_id( id ) {}
  };

  static void initialize( const SelectDevice = SelectDevice() );

  static void finalize();

  /*--------------------------------*/
  /** \brief  Set the device in a 'sleep' state.
   *
   *  This function sets the device in a "sleep" state
   *  in which it is not ready for work and should
   *  consume less resources.
   *
   *  Return 'true' if the device is in the 'sleep' state.
   *  Return 'false' if the device is actively working and
   *  could not enter the 'sleep' state.
   */
  static bool sleep();

  /** \brief  Wake the device from the 'sleep' state so
   *          it is ready for work.
   *
   *  Return 'true' if the device is in the 'ready' state.
   *  Return 'false' if the device is actively working.
   */
  static bool wake();

  /** \brief  The parallel_for or parallel_reduce dispatch of a
   *          functor may return before the functor completes.
   *          Wait until all dispatched functors complete.
   */
  static void fence();

  /*--------------------------------*/
  /* -- Device specific functions --*/

  static size_type detect_device_count();

  /*--------------------------------*/
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <Cuda/Kokkos_Cuda_MemoryManager.hpp>
#include <Cuda/Kokkos_Cuda_Value.hpp>
#include <Cuda/Kokkos_Cuda_Parallel.hpp>
#include <Cuda/Kokkos_Cuda_ParallelFor.hpp>
#include <Cuda/Kokkos_Cuda_ParallelReduce.hpp>

#endif /* #ifndef KOKKOS_CUDA_HPP */

/* Partial specializations for optional data structures */
#include <Cuda/Kokkos_Cuda_Specialize.hpp>



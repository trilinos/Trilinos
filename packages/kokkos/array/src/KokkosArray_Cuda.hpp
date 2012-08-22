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

#ifndef KOKKOSARRAY_CUDA_HPP
#define KOKKOSARRAY_CUDA_HPP

#include <KokkosArray_Host.hpp>
#include <KokkosArray_Layout.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

class CudaMemorySpace ;

} // namespace Impl
} // namespace KokkosArray

/*--------------------------------------------------------------------------*/

namespace KokkosArray {

/// \class Cuda
/// \brief KokkosArray device that uses CUDA to run on GPUs.
class Cuda {
public:
  //! \name Type declarations that all KokkosArray devices must provide.
  //@{

  typedef Cuda          type ;
  typedef Cuda          layout_type ;
  typedef Cuda          device_type ;
  typedef unsigned int  size_type ;

  typedef Impl::CudaMemorySpace  memory_space ;
  typedef LayoutLeft             array_layout ;

  //--------------------------------------------------------------------------

  struct SelectDevice {
    int cuda_device_id ;
    SelectDevice() : cuda_device_id(0) {}
    explicit SelectDevice( int id ) : cuda_device_id( id ) {}
  };

  //@}
  //! \name Functions that all KokkosArray devices must implement.
  //@{

  /** \brief  Set the device in a "sleep" state.
   *
   * This function sets the device in a "sleep" state in which it is
   * not ready for work.  This may consume less resources than if the
   * device were in an "awake" state, but it may also take time to
   * bring the device from a sleep state to be ready for work.
   *
   * \return True if the device is in the "sleep" state, else false if
   *   the device is actively working and could not enter the "sleep"
   *   state.
   */
  static bool sleep();

  /// \brief Wake the device from the 'sleep' state so it is ready for work.
  ///
  /// \return True if the device is in the "ready" state, else "false"
  ///  if the device is actively working (which also means that it's
  ///  awake).
  static bool wake();

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void fence();

  //! Free any resources being consumed by the device.
  static void finalize();

  //@}
  //! \name Device-specific functions
  //@{

  //! Initialize, telling the CUDA run-time library which device to use.
  static void initialize( const SelectDevice = SelectDevice() );

  static size_type detect_device_count();

  //@}
};

} // namespace KokkosArray

/*--------------------------------------------------------------------------*/

#include <Cuda/KokkosArray_Cuda_MemorySpace.hpp>
#include <Cuda/KokkosArray_Cuda_View.hpp>

#include <Cuda/KokkosArray_Cuda_Parallel.hpp>
#include <Cuda/KokkosArray_Cuda_ParallelFor.hpp>
#include <Cuda/KokkosArray_Cuda_ParallelReduce.hpp>

#endif /* #ifndef KOKKOSARRAY_CUDA_HPP */

//----------------------------------------------------------------------------



/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CUDA_HPP
#define KOKKOS_CUDA_HPP

#include <iosfwd>
#include <vector>

#include <Kokkos_Macros.hpp>
#include <Kokkos_Threads.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class CudaExec ;
}
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/// \class Cuda
/// \brief Kokkos device that uses CUDA to run on GPUs.
class Cuda {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  typedef Cuda                  device_type ;
  typedef CudaSpace             memory_space ;
  typedef CudaSpace::size_type  size_type ;
  typedef LayoutLeft            array_layout ;
  typedef Kokkos::Threads       host_mirror_device_type ;

  //@}
  //--------------------------------------------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

#if defined( __CUDA_ARCH__ ) 
  KOKKOS_INLINE_FUNCTION static bool in_parallel() { return true ; }
#else
  KOKKOS_INLINE_FUNCTION static bool in_parallel() { return false ; }
#endif

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

  /** \brief  Print Cuda configuation */
  static void print_configuration( std::ostream & );

  //@}
  //! \name Device-specific functions
  //@{

  struct SelectDevice {
    int cuda_device_id ;
    SelectDevice() : cuda_device_id(0) {}
    explicit SelectDevice( int id ) : cuda_device_id( id ) {}
  };

  /** \brief  Initialize, telling the CUDA run-time library which device to use. */
  static void initialize( const SelectDevice = SelectDevice() );

  static int is_initialized();

  /** \brief  Cuda device architecture of the selected device.
   *          Matches the __CUDA_ARCH__ specification.
   */
  static size_type device_arch();


  /** \brief  Query device count. */
  static size_type detect_device_count();

  /** \brief  Detect the available devices and their architecture
   *          as defined by the __CUDA_ARCH__ specification.
   */
  static std::vector<unsigned> detect_device_arch();

  //@}
  //--------------------------------------------------------------------------

#if defined( __CUDA_ARCH__ )

  //! \name Functions for the functor device interface
  //@{

  __device__ inline int league_size() const { return gridDim.x ; }
  __device__ inline int league_rank() const { return blockIdx.x ; }

  __device__ inline int team_size() const { return blockDim.x ; }
  __device__ inline int team_rank() const { return threadIdx.x ; }

  __device__ inline void team_barrier() const { __syncthreads(); }
  __device__ inline unsigned int team_barrier_count(bool value) const
             { return __syncthreads_count(value); }

  /* Collectively compute the league-wide unordered exclusive prefix sum.
   * Values are ordered within a team, but not between teams (i.e. the start
   * values of thread 0 in each team are not ordered according to team number).
   * This call does not use a global synchronization. Multiple unordered scans
   * can be in-flight at the same time (using different scratch_views).
   * The scratch-view will hold the complete sum in the end.
   */
  template< class VT >
  __device__ inline typename VT::value_type unordered_scan
             (typename VT::value_type& value, VT& scratch_view);

  /* Collectively compute the team-wide exclusive prefix sum using CUDA Unbound.
   * Values are ordered, the last thread returns the sum of all values
   * in the team less its own value
   */
  template< typename T >
  __device__ inline T team_scan(T& value);

  __device__ inline void * get_shmem( const int size );

  __device__ inline Cuda( Impl::CudaExec & exec ) : m_exec(exec) {}
  __device__ inline Cuda( const Cuda & rhs ) : m_exec(rhs.m_exec) {}

  //@}

private:

  Impl::CudaExec & m_exec ;

  //--------------------------------------------------------------------------
#else

  int league_size() const ;
  int league_rank() const ;

  int team_size() const ;
  int team_rank() const ;

  void team_barrier() const ;
  unsigned int team_barrier_count(bool) const ;

  template< class VT >
    inline typename VT::value_type unordered_scan
             (typename VT::value_type& value, VT& scratch_view);
  template< typename T >
    inline T team_scan(T& value);

  void * get_shmem( const int size );

  Cuda( Impl::CudaExec & );

#endif

};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief Cuda-specific parallel work configuration */

struct CudaWorkConfig {
  Cuda::size_type  grid[3] ;   //< Grid dimensions
  Cuda::size_type  block[3] ;  //< Block dimensions
  Cuda::size_type  shared ;    //< Shared memory size

  CudaWorkConfig()
  {
    enum { WarpSize = 32 };
    grid[0] = grid[1] = grid[2] = 1 ;
    block[1] = block[2] = 1 ;
    block[0] = 8 * WarpSize ;
    shared = 0 ;
  }
};

template< class FunctorType >
inline
void parallel_for( const CudaWorkConfig & work_config ,
                   const FunctorType    & functor )
{
  Impl::ParallelFor< FunctorType , Cuda , CudaWorkConfig >
    ( work_config , functor );
}

template< class FunctorType , class FinalizeType >
inline
void parallel_reduce( const CudaWorkConfig & work_config ,
                      const FunctorType    & functor ,
                      const FinalizeType   & finalize );

template< class FunctorType >
inline
typename FunctorType::value_type
parallel_reduce( const CudaWorkConfig & work_config ,
                 const FunctorType    & functor );

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <Cuda/Kokkos_CudaExec.hpp>

#include <Cuda/Kokkos_Cuda_View.hpp>
#include <Cuda/Kokkos_Cuda_Parallel.hpp>
#include <Cuda/Kokkos_Cuda_ParallelReduce.hpp>

#endif /* #ifndef KOKKOS_CUDA_HPP */

//----------------------------------------------------------------------------



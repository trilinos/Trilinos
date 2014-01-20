/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

/// \file Kokkos_Parallel.hpp
/// \brief Declaration of parallel operators

#ifndef KOKKOS_PARALLEL_HPP
#define KOKKOS_PARALLEL_HPP

#include <cstddef>
#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>
#include <impl/Kokkos_Traits.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/// \class ParallelFor
/// \brief Implementation of the ParallelFor operator that has a
///   partial specialization for the device.
///
/// This is an implementation detail of parallel_for.  Users should
/// skip this and go directly to the nonmember function parallel_for.
template< class FunctorType ,
          class WorkSpec ,
          class DeviceType = typename FunctorType::device_type >
class ParallelFor ;

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

/// \class VectorParallel
/// \brief Request for parallel_for to attempt thread+vector parallelism.
struct VectorParallel
{
  const size_t nwork ;
  VectorParallel( const size_t n ) : nwork(n) {}
  operator size_t () const { return nwork ; }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief Execute \c functor \c work_count times in parallel.
 *
 * A "functor" is a class containing the function to execute in
 * parallel, any data needed for that execution, and a \c device_type
 * typedef.  Here is an example functor for parallel_for:
 *
 * \code
 *  class FunctorType {
 *  public:
 *    typedef  ...  device_type ;
 *    void operator() (IntType iwork) const ;
 *  };
 * \endcode
 *
 * In the above example, \c IntType is any integer type for which a
 * valid conversion from \c size_t to \c IntType exists.  Its
 * <tt>operator()</tt> method defines the operation to parallelize,
 * over the range of integer indices <tt>iwork=[0,work_count-1]</tt>.
 * This compares to a single iteration \c iwork of a \c for loop.
 */
template< class FunctorType >
inline
void parallel_for( const size_t        work_count ,
                   const FunctorType & functor )
{
  Impl::ParallelFor< FunctorType , size_t > tmp( functor , work_count );
}


/** \brief Execute \c functor \c work_count times in parallel, with vectorization.
 *
 * This is like parallel_for, except that it <i>mandates</i>
 * vectorization as well as parallelization of the given functor.  We
 * emphasize "mandates": this means that the user asserts that
 * vectorization is correct, and insists that the compiler vectorize.
 * Mandating vectorization is not always desirable, for example if the
 * body of the functor is complicated.  In some cases, users might
 * want to parallelize over threads, and use vectorization inside the
 * parallel operation.  Furthermore, the compiler might still be able
 * to vectorize through a parallel_for.  Thus, users should take care
 * not to use this execution option arbitrarily.
 */
template< class FunctorType >
inline
void vector_parallel_for( const size_t        work_count ,
                          const FunctorType & functor )
{
  Impl::ParallelFor< FunctorType , VectorParallel > tmp( functor , work_count );
}

template< class DeviceType >
class MultiFunctorParallelFor ;

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/// \class ParallelReduce
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType ,
          class WorkSpec ,
          class DeviceType = typename FunctorType::device_type >
class ParallelReduce ;

/// \class ReduceAdapter
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType ,
          class ValueType = typename FunctorType::value_type >
struct ReduceAdapter ;

} // namespace Impl
} // namespace Kokkos


namespace Kokkos {

/** \brief  Parallel reduction
 *
 * Example of a parallel_reduce functor for a POD (plain old data) value type:
 * \code
 *  class FunctorType { // For POD value type
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type ;
 *    void operator()( <intType> iwork , <podType> & update ) const ;
 *    void init( <podType> & update ) const ;
 *    void join( volatile       <podType> & update ,
 *               volatile const <podType> & input ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> & update ) const ;
 *  };
 * \endcode
 *
 * Example of a parallel_reduce functor for an array of POD (plain old data) values:
 * \code
 *  class FunctorType { // For array of POD value
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type[] ;
 *    void operator()( <intType> , <podType> update[] ) const ;
 *    void init( <podType> update[] ) const ;
 *    void join( volatile       <podType> update[] ,
 *               volatile const <podType> input[] ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> update[] ) const ;
 *  };
 * \endcode
 */
template< class FunctorType >
inline
void parallel_reduce( const size_t        work_count ,
                      const FunctorType & functor )
{
  Impl::ParallelReduce< FunctorType , size_t > reduce( functor , work_count );
}

/** \brief  Parallel reduction and output to host.
 *
 *  If FunctorType::value_type is
 *    - \c PodType,  then \c reference_type is <tt>PodType & </tt>.
 *    - <tt>PodType[]</tt>, then \c reference_type is <tt>PodType * </tt>.
 */
template< class FunctorType >
inline
void parallel_reduce( const size_t work_count ,
                      const FunctorType & functor ,
                      typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type result )
{
  Impl::ParallelReduce< FunctorType, size_t >
    reduce( functor , work_count , Kokkos::Impl::ReduceAdapter< FunctorType >::pointer( result ) );

  reduce.wait();
}

template< class FunctorType >
inline
void parallel_reduce( const VectorParallel & work_count ,
                      const FunctorType & functor ,
                      typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type result )
{
  Impl::ParallelReduce< FunctorType, VectorParallel >
    reduce( functor , work_count , Kokkos::Impl::ReduceAdapter< FunctorType >::pointer( result ) );

  reduce.wait();
}

template< class DeviceType >
class MultiFunctorParallelReduce ;

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/// \class ParallelReduce
/// \brief Implementation detail of parallel_reduce.
///
/// This is an implementation detail of parallel_reduce.  Users should
/// skip this and go directly to the nonmember function parallel_reduce.
template< class FunctorType ,
          class WorkSpec ,
          class DeviceType = typename FunctorType::device_type >
class ParallelScan ;

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

template< class FunctorType >
inline
void parallel_scan( const size_t        work_count ,
                    const FunctorType & functor )
{
  Impl::ParallelScan< FunctorType , size_t > scan( functor , work_count );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Parallel work request for shared memory, league size, and team size.
 *
 *  If the shared size is too large then slow (global) memory will be used.
 *  If the league or team size are too large then they will be reduced.
 */
struct ParallelWorkRequest {
  size_t  league_size ; ///<  Size of league (number of teams in a league)
  size_t  team_size ;   ///<  Size of team (number of threads in a team)

  KOKKOS_INLINE_FUNCTION
  ParallelWorkRequest() : league_size(0), team_size(0) {}

  KOKKOS_INLINE_FUNCTION
  ParallelWorkRequest( size_t s0 , size_t s1 ) : league_size(s0), team_size(s1) {}
};

/** \brief  Execute functor in parallel with work request,
 *          the actual league_size and team_size may be smaller.
 *
 *  class FunctorType {
 *  public:
 *    typedef  ...  device_type ;
 *    void operator()( device_type ) const ;
 *  };
 */
template< class FunctorType >
inline
void parallel_for( const ParallelWorkRequest & request ,
                   const FunctorType         & functor )
{
  Kokkos::Impl::ParallelFor< FunctorType , ParallelWorkRequest >( functor , request );
}

} // namespace Kokkos

namespace Kokkos {

/** \brief  Parallel reduction.
 *
 *  class FunctorType {
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type ; // POD type
 *    void operator()( device_type , <podType> & ) const ;
 *    void init( <podType> & ) const ;
 *    void join( volatile       <podType> & update ,
 *               volatile const <podType> & input ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> & update ) const ;
 *  };
 *
 *  class FunctorType { // For array of POD value
 *  public:
 *    typedef    ...     device_type ;
 *    typedef <podType>  value_type[] ;
 *    void operator()( device_type , <podType> update[] ) const ;
 *    void init( <podType> update[] ) const ;
 *    void join( volatile       <podType> update[] ,
 *               volatile const <podType> input[] ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> update[] ) const ;
 *  };
 */
template< class FunctorType >
inline
void parallel_reduce( const Kokkos::ParallelWorkRequest  & request ,
                      const FunctorType          & functor )
{
  Impl::ParallelReduce< FunctorType , Kokkos::ParallelWorkRequest > reduce( functor , request );
}

template< class FunctorType >
inline
void parallel_reduce( const Kokkos::ParallelWorkRequest  & request ,
                      const FunctorType          & functor ,
                      typename Kokkos::Impl::ReduceAdapter< FunctorType >::reference_type result )
{
  Impl::ParallelReduce< FunctorType , Kokkos::ParallelWorkRequest >
    reduce( functor , request , Kokkos::Impl::ReduceAdapter< FunctorType >::pointer( result ) );

  reduce.wait(); // Wait for reduce to complete and output result
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class Enable = void >
struct FunctorHasJoin : public false_type {};

template< class FunctorType >
struct FunctorHasJoin< FunctorType , typename enable_if< 0 < sizeof( & FunctorType::join ) >::type >
  : public true_type {};

template< class FunctorType , class Enable = void >
struct FunctorHasFinal : public false_type {};

template< class FunctorType >
struct FunctorHasFinal< FunctorType , typename enable_if< 0 < sizeof( & FunctorType::final ) >::type >
  : public true_type {};

template< class FunctorType , class Enable = void >
struct FunctorShmemSize
{
  static inline size_t value( const FunctorType & ) { return 0 ; }
};

template< class FunctorType >
struct FunctorShmemSize< FunctorType , typename enable_if< 0 < sizeof( & FunctorType::shmem_size ) >::type >
{
  static inline size_t value( const FunctorType & f ) { return f.shmem_size() ; }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ScalarType >
struct ReduceAdapter
{
  enum { StaticValueSize = sizeof(ScalarType) };

  typedef ScalarType & reference_type  ;
  typedef ScalarType * pointer_type  ;
  typedef ScalarType   scalar_type  ;

  KOKKOS_INLINE_FUNCTION static
  reference_type reference( void * p ) { return *((ScalarType*) p); }

  KOKKOS_INLINE_FUNCTION static
  reference_type reference( void * p , unsigned i ) { return ((ScalarType*) p)[i]; }

  KOKKOS_INLINE_FUNCTION static
  pointer_type pointer( reference_type p ) { return & p ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_count( const FunctorType & ) { return 1 ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_size( const FunctorType & ) { return sizeof(ScalarType); }

  KOKKOS_INLINE_FUNCTION static
  void copy( const FunctorType & , void * const dst , const void * const src )
    { *((scalar_type*)dst) = *((const scalar_type*)src); }

  KOKKOS_INLINE_FUNCTION static
  void join( const FunctorType & f , volatile void * update , volatile const void * input )
    { f.join( *((volatile ScalarType*)update) , *((volatile const ScalarType*)input) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & f ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    FunctorHasFinal<F>::value )
                                >::type * p )
    { f.final( *((ScalarType *) p ) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    ! FunctorHasFinal<F>::value )
                                >::type * )
    {}
};

template< class FunctorType , class ScalarType >
struct ReduceAdapter< FunctorType , ScalarType[] >
{
  enum { StaticValueSize = 0 };

  typedef ScalarType * reference_type  ;
  typedef ScalarType * pointer_type  ;
  typedef ScalarType   scalar_type  ;

  KOKKOS_INLINE_FUNCTION static
  ScalarType * reference( void * p ) { return (ScalarType*) p ; }

  KOKKOS_INLINE_FUNCTION static
  reference_type reference( void * p , unsigned i ) { return ((ScalarType*) p)+i; }

  KOKKOS_INLINE_FUNCTION static
  pointer_type pointer( reference_type p ) { return p ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_count( const FunctorType & f ) { return f.value_count ; }

  KOKKOS_INLINE_FUNCTION static
  unsigned value_size( const FunctorType & f ) { return f.value_count * sizeof(ScalarType); }

  KOKKOS_INLINE_FUNCTION static
  void copy( const FunctorType & f , void * const dst , const void * const src )
    {
      for ( int i = 0 ; i < int(f.value_count) ; ++i ) {
        ((scalar_type*)dst)[i] = ((const scalar_type*)src)[i];
      }
    }

  KOKKOS_INLINE_FUNCTION static
  void join( const FunctorType & f , volatile void * update , volatile const void * input )
    { f.join( ((volatile ScalarType*)update) , ((volatile const ScalarType*)input) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & f ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    FunctorHasFinal<F>::value )
                                >::type * p )
    { f.final( ((ScalarType *) p ) ); }

  template< class F >
  KOKKOS_INLINE_FUNCTION static
  void final( const F & ,
              typename enable_if< ( is_same<F,FunctorType>::value &&
                                    ! FunctorHasFinal<F>::value )
                                >::type * )
    {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_PARALLEL_HPP */


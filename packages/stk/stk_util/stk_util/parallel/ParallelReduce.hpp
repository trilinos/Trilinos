// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_parallel_ParallelReduce_hpp
#define stk_util_parallel_ParallelReduce_hpp

#include <stk_util/stk_config.h>
#include <stdint.h>                     // for int64_t
#include <cstddef>                      // for size_t
#include <iosfwd>                       // for ostream
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/util/SimpleArrayOps.hpp>  // for BitAnd, BitOr, Copy, etc
#include <string>                       // for string
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/MPI.hpp>
#include "stk_util/environment/ReportHandler.hpp"

//------------------------------------------------------------------------

namespace stk {

#if defined (STK_HAS_MPI)
template<typename T>
void all_reduce_impl( ParallelMachine comm , const T * local , T * global , unsigned count, MPI_Op op )
{
  T * tmp = const_cast<T*>( local );
  BABBLE_STK_PARALLEL_COMM(comm, "                      calling MPI_Allreduce from all_reduce");
  MPI_Allreduce( tmp , global , count , sierra::MPI::Datatype<T>::type() , op , comm );
}

void all_reduce_impl(ParallelMachine comm, const size_t * local, size_t * global, unsigned count, MPI_Op op);

template<typename T>
void all_reduce_max( ParallelMachine comm , const T * local , T * global , unsigned count )
{
  all_reduce_impl(comm, local, global, count, MPI_MAX);
}

template<typename T>
void all_reduce_min( ParallelMachine comm , const T * local , T * global , unsigned count )
{
  all_reduce_impl(comm, local, global, count, MPI_MIN);
}

template<typename T>
void all_reduce_sum( ParallelMachine comm , const T * local , T * global , unsigned count )
{
  all_reduce_impl(comm, local, global, count, MPI_SUM);
}

template<typename T, typename IdType>
void
all_reduce_loc_impl(ParallelMachine comm,
    const T local_extrema[],
    const IdType local_loc[],
    T global_extrema[],
    IdType global_loc[],
    unsigned n,
    MPI_Op mpiOp)
{
    typedef sierra::MPI::Loc<T, IdType> MpiLocType;
    if ( n < 1 ) return;
    MpiLocType * const vin  = new MpiLocType[n] ;
    MpiLocType * const vout = new MpiLocType[n] ;
    for (unsigned i = 0 ; i < n ; ++i ) {
      vin[i].m_value = local_extrema[i] ;
      vin[i].m_loc = local_loc[i] ;
    }
    ThrowRequire(MPI_SUCCESS == MPI_Allreduce( vin, vout, (int) n,
                                       sierra::MPI::Datatype<MpiLocType >::type(),
                                       mpiOp,
                                       comm )
        );

    for (unsigned i = 0 ; i < n ; ++i ) {
      global_extrema[i] = vout[i].m_value ;
      global_loc[i] = vout[i].m_loc ;
    }
    delete[] vin ;
    delete[] vout ;
}



template<typename T, typename IdType>
void
all_reduce_minloc(ParallelMachine comm,
    const T local_extrema[],
    const IdType local_loc[],
    T global_extrema[],
    IdType global_loc[],
    unsigned count)
{
  all_reduce_loc_impl(comm,
                  local_extrema,
                  local_loc,
                  global_extrema,
                  global_loc,
                  count,
                  sierra::MPI::get_mpi_loc_op<T, std::less<T>, IdType>());
}

template<typename T, typename IdType>
void
all_reduce_maxloc(ParallelMachine comm,
    const T local_extrema[],
    const IdType local_loc[],
    T global_extrema[],
    IdType global_loc[],
    unsigned count)
{
  all_reduce_loc_impl(comm,
                  local_extrema,
                  local_loc,
                  global_extrema,
                  global_loc,
                  count,
                  sierra::MPI::get_mpi_loc_op<T, std::greater<T>, IdType>());
}

#else
template<typename T>
void all_reduce_max( ParallelMachine comm , const T * local , T * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

template<typename T>
void all_reduce_min( ParallelMachine comm , const T * local , T * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

template<typename T>
void all_reduce_sum( ParallelMachine comm , const T * local , T * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

template<typename T, typename IdType>
void
all_reduce_minloc(ParallelMachine comm,
    const T local_extrema[],
    const IdType local_loc[],
    T global_extrema[],
    IdType global_loc[],
    unsigned count)
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global_extrema[i] = local_extrema[i] ; }
  for ( unsigned i = 0 ; i < count ; ++i ) { global_loc[i] = local_loc[i] ; }
}

template<typename T, typename IdType>
void
all_reduce_maxloc(ParallelMachine comm,
    const T local_extrema[],
    const IdType local_loc[],
    T global_extrema[],
    IdType global_loc[],
    unsigned count)
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global_extrema[i] = local_extrema[i] ; }
  for ( unsigned i = 0 ; i < count ; ++i ) { global_loc[i] = local_loc[i] ; }
}
#endif

/** \addtogroup parallel_module
 *  \{
 */


/** \brief  Write string from any or all processors
 *          to the ostream on the root processor.
 */
void all_write_string( ParallelMachine ,
                       std::ostream & ,
                       const std::string & );

/** \brief  Parallel bitwise-or to all processors */
void all_reduce_bor( ParallelMachine ,
                     const unsigned * local ,
                     unsigned * global , unsigned count );

/** Aggregated parallel in-place reduce-to-all-processors operations.
 *
 *  example:
 *  <PRE>
 *    ParallelMachine comm = ... ;
 *    double a[5] ;
 *    int    b[3] ;
 *    all_reduce( comm , ReduceSum &lt; 5 &gt;( a ) & ReduceMax &lt; 3 &gt;( b ) );
 *
 *  Reduction options include:
 *    ReduceSum    &lt; N &gt; ( T * )   // Summation
 *    ReduceProd   &lt; N &gt; ( T * )   // Product
 *    ReduceMax    &lt; N &gt; ( T * )   // Maximum
 *    ReduceMin    &lt; N &gt; ( T * )   // Minimum
 *    ReduceBitOr  &lt; N &gt; ( T * )   // Bit-wise OR
 *    ReduceBitAnd &lt; N &gt; ( T * )   // Bit-wise AND
 * </PRE>
 */
template < class ReduceOp >
void all_reduce( ParallelMachine , const ReduceOp & );

/** \} */

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

namespace stk {
namespace {
// Blank namespace so that this class produces local symbols,
// avoiding complaints from a linker of multiple-define symbols.

struct ReduceEnd {
  struct WorkType {};
  void copyin(  WorkType & ) const {}
  void copyout( WorkType & ) const {}
  static void op( WorkType & , WorkType & ) {}
};

// Workhorse class for aggregating reduction operations.

template <class Op, typename T, class Next>
struct Reduce {

  typedef T Type ;
  enum { N = Op::N };

  struct WorkType {
    typename Next::WorkType m_next ;
    Type                    m_value[N];
  };

  Next   m_next ;
  Type * m_value ;

  // Copy values into buffer:
  void copyin( WorkType & w ) const
    { Copy<N>( w.m_value , m_value ); m_next.copyin( w.m_next ); }

  // Copy value out from buffer:
  void copyout( WorkType & w ) const
    { Copy<N>( m_value , w.m_value ); m_next.copyout( w.m_next ); }

  // Reduction function
  static void op( WorkType & out , WorkType & in )
    { Op( out.m_value , in.m_value ); Next::op( out.m_next , in.m_next ); }

  // Aggregate reduction operations, use '&' for left-to-right evaluation
  template<class OpB, typename TB>
  Reduce<OpB, TB, Reduce<Op,T,Next> >
    operator & ( const Reduce<OpB,TB,ReduceEnd> & rhs )
      { return Reduce<OpB, TB, Reduce<Op,T,Next> >( rhs , *this ); }

  // Constructor for aggregation:
  Reduce( const Reduce<Op,T, ReduceEnd> & arg_val , const Next & arg_next )
    : m_next( arg_next ), m_value( arg_val.m_value ) {}

  // Constructor for aggregate member:
  explicit Reduce( Type * arg_value )
   : m_next(), m_value( arg_value ) {}

  static void void_op( void*inv, void*inoutv, int*, ParallelDatatype*);
};

template <class Op, typename T, class Next>
void Reduce<Op,T,Next>::void_op( void*inv, void*inoutv,int*,ParallelDatatype*)
{
  op( * reinterpret_cast<WorkType*>( inoutv ) ,
      * reinterpret_cast<WorkType*>( inv ) );
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {

template<unsigned N, typename T>
inline
Reduce< Sum<N> , T, ReduceEnd> ReduceSum( T * value )
{ return Reduce< Sum<N>, T, ReduceEnd >( value ); }

template<unsigned N, typename T>
inline
Reduce< Prod<N>, T, ReduceEnd > ReduceProd( T * value )
{ return Reduce< Prod<N>, T, ReduceEnd >( value ); }

template<unsigned N, typename T>
inline
Reduce< Max<N>, T, ReduceEnd> ReduceMax( T * value )
{ return Reduce< Max<N>, T, ReduceEnd>( value ); }

template<unsigned N, typename T>
inline
Reduce< Min<N>, T, ReduceEnd> ReduceMin( T * value )
{ return Reduce<Min<N>, T, ReduceEnd>( value ); }

template<unsigned N, typename T>
inline
Reduce< BitOr<N>, T, ReduceEnd> ReduceBitOr( T * value )
{ return Reduce< BitOr<N>, T, ReduceEnd>( value ); }

template<unsigned N, typename T>
inline
Reduce< BitAnd<N>, T, ReduceEnd> ReduceBitAnd( T * value )
{ return Reduce< BitAnd<N>, T, ReduceEnd>( value ); }

//----------------------------------------------------------------------
// all_reduce( comm , ReduceSum<5>( A ) & ReduceMax<3>( B ) );

extern "C" {
typedef void (*ParallelReduceOp)
  ( void * inv , void * outv , int * , ParallelDatatype * );
}

void all_reduce( ParallelMachine  arg_comm ,
                 ParallelReduceOp arg_op ,
                 void           * arg_in ,
                 void           * arg_out ,
                 unsigned         arg_len );

namespace {

template < class ReduceOp >
void all_reduce_driver( ParallelMachine comm , const ReduceOp & op )
{
  typedef typename ReduceOp::WorkType WorkType ;

  WorkType inbuf , outbuf ;

  ParallelReduceOp f =
    reinterpret_cast<ParallelReduceOp>( & ReduceOp::void_op );
  op.copyin( inbuf );
  all_reduce( comm , f , & inbuf, & outbuf, sizeof(WorkType) );
  op.copyout( outbuf );
}

}

template < class ReduceOp >
inline
void all_reduce( ParallelMachine comm , const ReduceOp & op )
{ all_reduce_driver<ReduceOp>( comm , op ); }

}

#endif /* DOXYGEN_COMPILE */

//----------------------------------------------------------------------

#endif


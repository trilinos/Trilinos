#ifndef stk_util_parallel_ParallelReduce_hpp
#define stk_util_parallel_ParallelReduce_hpp

#include <cstddef>
#include <iosfwd>
#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/SimpleArrayOps.hpp>

//------------------------------------------------------------------------

namespace stk {

/** \addtogroup parallel_module
 *  \{
 */

// REFACTOR: Replace ReduceSum with Sum?, etc...  Should be possible

/** \brief  Write string from any or all processors
 *          to the ostream on the root processor.
 */
void all_write_string( ParallelMachine ,
                       std::ostream & ,
                       const std::string & );

/** \brief  Parallel summation to all processors */
void all_reduce_sum( ParallelMachine ,
                     const double * local , double * global , unsigned count );

/** \brief  Parallel summation to all processors */
void all_reduce_sum( ParallelMachine ,
                     const float * local , float * global , unsigned count );

/** \brief  Parallel summation to all processors */
void all_reduce_sum( ParallelMachine ,
                     const int * local , int * global , unsigned count );

/** \brief  Parallel summation to all processors */
void all_reduce_sum( ParallelMachine ,
                     const size_t * local , size_t * global , unsigned count );

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


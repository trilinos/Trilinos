/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/




#include <stk_algsup/AlgorithmRunner.hpp>

#ifdef STK_HAVE_TPI

#include <TPI.h>

#include <stdexcept>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

namespace stk {
namespace {

//----------------------------------------------------------------------

struct RunTPI { 
  const mesh::Selector             & selector ;
  const mesh::PartVector           & union_parts ;
  const std::vector<mesh::Bucket*> & buckets ; 
  const AlgorithmInterface         & alg ;

  RunTPI( const mesh::Selector             & arg_selector ,
          const mesh::PartVector           & arg_union_parts ,
          const std::vector<mesh::Bucket*> & arg_buckets ,
          const AlgorithmInterface         & arg_alg )
    : selector( arg_selector ), 
      union_parts(arg_union_parts), 
      buckets( arg_buckets ), 
      alg( arg_alg ) 
  {}
 
  ~RunTPI() 
  {}
};   

extern "C" {
 
static void RunTPI_join( TPI_Work * work , const void * reduce )
{  
  const RunTPI & myself = * ((const RunTPI *) work->info );

  myself.alg.join( work->reduce , reduce );
}

static void RunTPI_init( TPI_Work * work )
{
  const RunTPI & myself = * ((const RunTPI *) work->info );

  myself.alg.init( work->reduce );
}
   
static void RunTPI_apply( TPI_Work * work )
{
  const RunTPI & myself = * ((const RunTPI *) work->info );

  myself.alg.apply_one( myself.selector ,
                        myself.union_parts ,
                        * myself.buckets[ work->rank ] ,
                        work->reduce 
                        );
}

}
 
class AlgorithmRunnerTPI : public AlgorithmRunnerInterface {
public:

  void run_alg( const mesh::Selector                & selector ,
                const mesh::PartVector              & union_parts ,
                const std::vector< mesh::Bucket * > & buckets ,
                const AlgorithmInterface            & alg ,
                void                                * reduce ) const ;

  AlgorithmRunnerTPI( int nthreads ) : result( 0 <= TPI_Init( nthreads ) ) {}

  const bool result ;
};

void AlgorithmRunnerTPI::run_alg(
  const mesh::Selector                & selector ,
  const mesh::PartVector              & union_parts ,
  const std::vector< mesh::Bucket * > & buckets ,
  const AlgorithmInterface            & alg ,
  void                                * reduce ) const
{
  if ( reduce && ! alg.m_reduce_allocation_size ) {
    std::string msg("AlgorithmRunnerTPI: ERROR reduce value with zero size");
    throw std::invalid_argument(msg);
  }
 
  if ( ! buckets.empty() ) {
 
    RunTPI tmp( selector, union_parts, buckets , alg );
 
    if ( reduce ) {
      TPI_Run_reduce( RunTPI_apply , & tmp ,
                      buckets.size() ,
                      RunTPI_join ,
                      RunTPI_init ,
                      alg.m_reduce_allocation_size ,
                      reduce );
    }
    else {
      TPI_Run( RunTPI_apply , & tmp , buckets.size() , 0 );
    }
  }
}
 
} // namespace

//----------------------------------------------------------------------

AlgorithmRunnerInterface * algorithm_runner_tpi( int nthreads )
{
  static AlgorithmRunnerTPI runner( nthreads );

  return runner.result ? & runner : NULL ;
}

} // namespace stk

#else

namespace stk {

AlgorithmRunnerInterface * algorithm_runner_tpi( int nthreads )
{
  return NULL ;
}

} // namespace stk

#endif


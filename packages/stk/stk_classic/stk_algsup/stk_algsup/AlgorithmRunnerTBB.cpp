/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_algsup/AlgorithmRunner.hpp>

#ifdef STK_HAVE_TBB

#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/scalable_allocator.h>
#include <tbb/partitioner.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>

namespace stk_classic {
namespace {

//----------------------------------------------------------------------

struct RunTBB {
  const mesh::Selector & selector ;
  const mesh::PartVector & union_parts ;
  const std::vector<mesh::Bucket*> & buckets ;
  const AlgorithmInterface         & alg ;

  void operator()(const tbb::blocked_range<int>& r) const;

  RunTBB( const mesh::Selector & arg_selector ,
          const mesh::PartVector & arg_union_parts ,
          const std::vector<mesh::Bucket*> & arg_buckets ,
          const AlgorithmInterface         & arg_alg );

  ~RunTBB();
};

RunTBB::RunTBB(
  const mesh::Selector    & arg_selector ,
  const mesh::PartVector & arg_union_parts ,
  const std::vector<mesh::Bucket*> & arg_buckets ,
  const AlgorithmInterface         & arg_alg )
  : selector( arg_selector ),
    union_parts( arg_union_parts ),
    buckets( arg_buckets ),
    alg( arg_alg )
{}

RunTBB::~RunTBB()
{
}

void RunTBB::operator()( const tbb::blocked_range<int> & r ) const
{
  for ( int i = r.begin() ; i < r.end() ; ++i ) {
    alg.apply_one( selector , union_parts , * buckets[i] , NULL );
  }
}

struct RunTBBreduce {
  const mesh::Selector    & selector ;
  const mesh::PartVector & union_parts ;
  const std::vector<mesh::Bucket*> & buckets ;
  const AlgorithmInterface         & alg ;
  void                             * reduce ;

  void operator()(const tbb::blocked_range<int>& r);

  void join( const RunTBBreduce & rhs ) const ;

  RunTBBreduce( const RunTBBreduce & rhs , tbb::split );

  RunTBBreduce( const mesh::Selector    & arg_selector ,
          const mesh::PartVector & arg_union_parts ,
          const std::vector<mesh::Bucket*> & arg_buckets ,
          const AlgorithmInterface         & arg_alg ,
          void                             * arg_reduce = NULL );

  ~RunTBBreduce();
};

RunTBBreduce::RunTBBreduce( const RunTBBreduce & rhs , tbb::split )
  : selector( rhs.selector ),
    union_parts( rhs.union_parts ),
    buckets(  rhs.buckets ),
    alg(      rhs.alg ),
    reduce( NULL )
{
  if ( rhs.reduce ) {
    reduce = malloc( alg.m_reduce_allocation_size ); //scalable_malloc ?
    alg.init( reduce );
  }
}

RunTBBreduce::~RunTBBreduce()
{
  if ( reduce ) { free( reduce ); /* scalable_free ? */}
}

void RunTBBreduce::join( const RunTBBreduce & rhs ) const
{
  alg.join( reduce , rhs.reduce );
}

void RunTBBreduce::operator()( const tbb::blocked_range<int> & r )
{
  for ( int i = r.begin() ; i < r.end() ; ++i ) {
    alg.apply_one( selector , union_parts, * buckets[i] , reduce );
  }
}

RunTBBreduce::RunTBBreduce(
  const mesh::Selector    & arg_selector ,
  const mesh::PartVector & arg_union_parts ,
  const std::vector<mesh::Bucket*> & arg_buckets ,
  const AlgorithmInterface         & arg_alg ,
  void                              * arg_reduce )
  : selector( arg_selector ),
    union_parts( arg_union_parts ),
    buckets( arg_buckets ),
    alg( arg_alg ),
    reduce( arg_reduce )
{}

//----------------------------------------------------------------------

class AlgorithmRunnerTBB : public AlgorithmRunnerInterface {
public:
  AlgorithmRunnerTBB(int nthreads)
   : tbb_task_init_(NULL)
  {
    tbb_task_init_ = new tbb::task_scheduler_init(nthreads);
  }

  ~AlgorithmRunnerTBB()
  {
    delete tbb_task_init_;
  }

  void run_alg( const mesh::Selector & selector ,
                const mesh::PartVector & union_parts ,
                const std::vector< mesh::Bucket * > & buckets ,
                const AlgorithmInterface      & alg ,
                void                          * reduce ) const ;

private:
  tbb::task_scheduler_init* tbb_task_init_;
};
 
void AlgorithmRunnerTBB::run_alg(
  const mesh::Selector & selector ,
  const mesh::PartVector & union_parts ,
  const std::vector< mesh::Bucket * > & buckets ,
  const AlgorithmInterface      & alg ,
  void * reduce ) const
{
  static tbb::affinity_partitioner ap;

  if ( reduce && ! alg.m_reduce_allocation_size ) {
    std::string msg("AlgorithmRunnerTBB: ERROR reduce value with zero size");
    throw std::invalid_argument(msg);
  }

  if ( ! buckets.empty() ) {

    tbb::blocked_range<int> range( 0 , buckets.size() );

    if ( reduce ) {
      RunTBBreduce tmp( selector , union_parts , buckets , alg, reduce );

      tbb::parallel_reduce( range , tmp , ap );
      tmp.reduce = NULL ; /* Prevent the tbb::scalable_free( tmp.reduce ); */
    }
    else {
      RunTBB tmp( selector , union_parts , buckets , alg );

      tbb::parallel_for( range, tmp , ap);
    }
  }
}

} // namespace

AlgorithmRunnerInterface * algorithm_runner_tbb( int nthreads )
{
  static AlgorithmRunnerTBB runner(nthreads) ;

  return & runner ;
}

} // namespace stk_classic

#else

namespace stk_classic {

AlgorithmRunnerInterface * algorithm_runner_tbb( int nthreads )
{
  return NULL ;
}

} // namespace stk_classic

#endif


/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_algsup/AlgorithmRunner.hpp>

namespace stk {

//----------------------------------------------------------------------

AlgorithmInterface::~AlgorithmInterface(){}

//----------------------------------------------------------------------

void AlgorithmInterface::apply_one(
  const mesh::Selector          & selector ,
  const mesh::PartVector & union_part_vector ,
  const mesh::Bucket            & bucket ,
  void                          * reduce ) const
{
  mesh::PartVector parts ;

  const bool run_it = selector( bucket );
  get_involved_parts( union_part_vector, bucket, parts );

  if ( run_it ) {
    if ( 0 < m_maximum_entity_count ) {
      for ( mesh::Bucket::iterator j = bucket.begin(); j != bucket.end() ; ) {
        mesh::Bucket::iterator e = j ;
        if ( static_cast<ptrdiff_t>( bucket.end() - e ) < static_cast<ptrdiff_t>(m_maximum_entity_count) ) {
          e = bucket.end();
        }
        else {
          e += m_maximum_entity_count ;
        }
        apply( j , e , parts , reduce );
        j = e ;
      }
    }
    else {
      apply( bucket.begin() , bucket.end() , parts , reduce );
    }
  }
}

//void AlgorithmInterface::apply_one(
//  const mesh::Selector          & selector ,
//  const mesh::Bucket            & bucket ,
//  void                          * reduce ) const
//{
//  const mesh::PartVector empty_union_part_vector;
//  apply_one(selector,empty_union_part_vector,bucket,reduce);
//}

//----------------------------------------------------------------------

namespace {

class AlgorithmRunnerNonThread : public AlgorithmRunnerInterface {
public:
 
  void run_alg( const mesh::Selector          & selector,
                const mesh::PartVector & union_parts ,
                const std::vector< mesh::Bucket * > & buckets,
                const AlgorithmInterface & algorithm,
                void * reduce ) const ;
 
  AlgorithmRunnerNonThread() {}
  ~AlgorithmRunnerNonThread() {}
};

void AlgorithmRunnerNonThread::run_alg(
  const mesh::Selector                 & selector ,
  const mesh::PartVector               & union_parts ,
  const std::vector< mesh::Bucket * >  & buckets ,
  const AlgorithmInterface             & algorithm ,
  void                                 * reduce ) const
{
  for ( std::vector< mesh::Bucket * >::const_iterator
        i = buckets.begin() ; i != buckets.end() ; ++i ) {
    algorithm.apply_one( selector , union_parts, **i , reduce );
  }
}

}

//----------------------------------------------------------------------

AlgorithmRunnerInterface * algorithm_runner_non_thread()
{
  static AlgorithmRunnerNonThread runner ;
  return & runner ;
}

} // namespace stk



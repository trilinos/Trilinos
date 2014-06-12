/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_algsup_AlgorithmRunner_hpp
#define stk_algsup_AlgorithmRunner_hpp

#include <utility>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

namespace stk_classic {

//----------------------------------------------------------------------

class AlgorithmInterface ;
class AlgorithmRunnerInterface ;

//----------------------------------------------------------------------

/** A local-serial algorithm runner */
AlgorithmRunnerInterface * algorithm_runner_non_thread();

/** A local-parallel algorithm runner using pthreads */
AlgorithmRunnerInterface * algorithm_runner_tpi( int nthreads );

/** A local-parallel algorithm runner using TBB */
AlgorithmRunnerInterface * algorithm_runner_tbb( int nthreads );

//----------------------------------------------------------------------

class AlgorithmRunnerInterface {
public:
  virtual ~AlgorithmRunnerInterface() {}

  /** \brief  Run an algorithm.
   *
   *  The Algorithm class must satisfy the following interface:
   *
   *  \code
   *  class Algorithm {
   *  public:
   *    <integer_type> maximum_entity_count ;
   *    void apply( mesh::Bucket::iterator i ,
   *                mesh::Bucket::iterator j ,
   *                const mesh::PartVector & selected_parts ) const ;
   *
   *  };
   *  \endcode
   */
  template< class Algorithm >
  void run_parts( const mesh::Selector                & selector ,
                  const mesh::PartVector              & union_parts ,
                  const std::vector< mesh::Bucket * > & buckets ,
                  const Algorithm                     & algorithm ) const ;

  /** \brief  Run an algorithm.
   *
   *  The Algorithm class must satisfy the following interface:
   *
   *  \code
   *  class Algorithm {
   *  public:
   *    <integer_type> maximum_entity_count ;
   *    void apply( mesh::Bucket::iterator i ,
   *                mesh::Bucket::iterator j ) const ;
   *
   *  };
   *  \endcode
   */
  template< class Algorithm >
  void run( const mesh::Selector                & selector ,
            const mesh::PartVector              & union_parts ,
            const std::vector< mesh::Bucket * > & buckets ,
            const Algorithm                     & algorithm ) const ;

  /** \brief  Run an algorithm with a reduction operation.
   *
   *  The Algorithm class must satisfy the following interface:
   *
   *  \code
   *  class Algorithm {
   *  public:
   *    typedef <value_type> reduce_type ;
   *    <integer_type>       reduce_count ;
   *    void init( reduce_type * out ) const ;
   *    void join( reduce_type * inout , const reduce_type * in ) const ;
   *
   *    <integer_type> maximum_entity_count ;
   *    void apply( mesh::Bucket::iterator i ,
   *                mesh::Bucket::iterator j ,
   *                const mesh::PartVector & selected_parts ,
   *                reduce_type * reduce_inout ) const ;
   *
   *  };
   *  \endcode
   *
   *  The apply method is passed a pointer to
   *    'reduce_type value[ reduce_count ]'
   *  into which it must reduce its results and
   *  for which it has thread-safe access.
   *  The join method is passed two pointers to
   *    reduce_type[ reduce_count ]
   *  values such that it must reduce the 'in' value into the 'inout' value;
   *  this method has thread-safe access to both 'in' and 'inout' values.
   */
  template< class Algorithm >
  void run_parts( const mesh::Selector                & selector ,
                  const mesh::PartVector              & union_parts ,
                  const std::vector< mesh::Bucket * > & buckets ,
                  const Algorithm                     & algorithm ,
                  typename Algorithm::reduce_type     * reduce_value ) const ;

  /** \brief  Run an algorithm with a reduction operation.
   *
   *  \code
   *  class Algorithm {
   *  public:
   *    typedef <value_type> reduce_type ;
   *    <integer_type>       reduce_count ;
   *    void init( reduce_type * out ) const ;
   *    void join( reduce_type * inout , const reduce_type * in ) const ;
   *
   *    <integer_type> maximum_entity_count ;
   *    void apply( mesh::Bucket::iterator i ,
   *                mesh::Bucket::iterator j ,
   *                reduce_type * reduce_inout ) const ;
   *
   *  };
   *  \endcode
   */
  template< class Algorithm >
  void run( const mesh::Selector                & selector ,
            const mesh::PartVector              & union_parts ,
            const std::vector< mesh::Bucket * > & buckets ,
            const Algorithm                     & algorithm ,
            typename Algorithm::reduce_type     * reduce_value ) const ;

private:
  AlgorithmRunnerInterface ( const AlgorithmRunnerInterface & );
  AlgorithmRunnerInterface & operator = ( const AlgorithmRunnerInterface & );

protected:

  AlgorithmRunnerInterface() {}

  /** Run many buckets in thread-parallel */
  virtual void run_alg( const mesh::Selector          & selector ,
                        const mesh::PartVector & union_parts ,
                        const std::vector< mesh::Bucket * > & buckets ,
                        const AlgorithmInterface            & algorithm ,
                        void * reduce ) const = 0 ;
};

//----------------------------------------------------------------------

/** Interface for internal wrapper-classes, not part of public API.
*/
class AlgorithmInterface {
public:
  const size_t m_maximum_entity_count ;
  const size_t m_reduce_allocation_size ;

  virtual void init( void * out ) const = 0 ;

  virtual void join( void * inout , const void * in ) const = 0 ;

  virtual void apply( mesh::Bucket::iterator i ,
                      mesh::Bucket::iterator j ,
                      const mesh::PartVector & selected_parts ,
                      void * reduce_inout ) const = 0 ;

  virtual ~AlgorithmInterface();

  //void apply_one( const mesh::Selector          & selector ,
  //                const mesh::Bucket                  & bucket ,
  //                void                                * reduce ) const ;

  void apply_one( const mesh::Selector   & selector ,
                  const mesh::PartVector & union_part_vector ,
                  const mesh::Bucket     & bucket ,
                  void                   * reduce ) const ;

protected:

  explicit AlgorithmInterface( )
    : m_maximum_entity_count( 0 ),
      m_reduce_allocation_size( 0 ) {}

  AlgorithmInterface( size_t count )
    : m_maximum_entity_count( count ),
      m_reduce_allocation_size( 0 ) {}

  AlgorithmInterface( size_t count , size_t size )
    : m_maximum_entity_count( count ),
      m_reduce_allocation_size( size ) {}

private:
  AlgorithmInterface( const AlgorithmInterface & );
  AlgorithmInterface & operator = ( const AlgorithmInterface & );
};

//----------------------------------------------------------------------

namespace {

template< class Algorithm >
class AlgorithmWrapper : public AlgorithmInterface {
private:
  void init( void * ) const {}
  void join( void * , const void * ) const {}
public:
  const Algorithm & m_alg ;

  void apply( mesh::Bucket::iterator i ,
              mesh::Bucket::iterator j ,
              const mesh::PartVector & parts, void * reduce ) const
  { m_alg.apply( i , j ); }

  explicit AlgorithmWrapper( const Algorithm & alg )
    : AlgorithmInterface( alg.maximum_entity_count ), m_alg( alg ) {}
};

template< class Algorithm >
class AlgorithmWrapperParts : public AlgorithmInterface {
private:
  void init( void * ) const {}
  void join( void * , const void * ) const {}
public:
  const Algorithm & m_alg ;

  void apply( mesh::Bucket::iterator i ,
              mesh::Bucket::iterator j ,
              const mesh::PartVector & selected_parts , void * reduce ) const
  { m_alg.apply( i , j , selected_parts ); }

  explicit AlgorithmWrapperParts( const Algorithm & alg )
    : AlgorithmInterface( alg.maximum_entity_count ), m_alg( alg ) {}
};

template< class Algorithm >
class AlgorithmWrapperReduce : public AlgorithmInterface {
public:
  typedef typename Algorithm::reduce_type reduce_type ;

  const Algorithm & m_alg ;

  void init( void * reduce_out ) const
  { m_alg.init( (reduce_type *) reduce_out ); }

  void join( void * reduce_inout , const void * reduce_in ) const
  {
    m_alg.join( (reduce_type *) reduce_inout ,
                (const reduce_type *) reduce_in );
  }

  void apply( mesh::Bucket::iterator i ,
              mesh::Bucket::iterator j ,
              const mesh::PartVector & parts, void * reduce_inout ) const
  {
    m_alg.apply( i , j , (reduce_type*) reduce_inout );
  }

  explicit AlgorithmWrapperReduce( const Algorithm & alg )
    : AlgorithmInterface( alg.maximum_entity_count ,
                          alg.reduce_count * sizeof( reduce_type ) ),
          m_alg( alg ) {}
};

template< class Algorithm >
class AlgorithmWrapperPartsReduce : public AlgorithmInterface {
public:
  typedef typename Algorithm::reduce_type reduce_type ;

  const Algorithm & m_alg ;

  void init( void * reduce_out ) const
  {
    m_alg.init( (reduce_type *) reduce_out );
  }

  void join( void * reduce_inout , const void * reduce_in ) const
  {
    m_alg.join( (reduce_type *) reduce_inout ,
                (const reduce_type *) reduce_in );
  }

  void apply( mesh::Bucket::iterator i ,
              mesh::Bucket::iterator j ,
              const mesh::PartVector & selected_parts ,
              void * reduce_inout ) const
  {
    m_alg.apply( i , j , selected_parts , (reduce_type*) reduce_inout );
  }

  explicit AlgorithmWrapperPartsReduce( const Algorithm & alg )
    : AlgorithmInterface( alg.maximum_entity_count ,
                          alg.reduce_count * sizeof(reduce_type) ),
      m_alg( alg ) {}
};

}

template< class Algorithm >
inline
void AlgorithmRunnerInterface::run(
  const mesh::Selector                & selector ,
  const mesh::PartVector              & union_parts ,
  const std::vector< mesh::Bucket * > & buckets ,
  const Algorithm                     & algorithm ,
  typename Algorithm::reduce_type     * reduce ) const
{
  const AlgorithmWrapperReduce<Algorithm> wrap( algorithm );

  run_alg( selector , union_parts, buckets , wrap , reduce );
}

template< class Algorithm >
inline
void AlgorithmRunnerInterface::run_parts(
  const mesh::Selector                & selector ,
  const mesh::PartVector              & union_parts ,
  const std::vector< mesh::Bucket * > & buckets ,
  const Algorithm                     & algorithm ,
  typename Algorithm::reduce_type     * reduce ) const
{
  const AlgorithmWrapperPartsReduce<Algorithm> wrap( algorithm );

  run_alg( selector , union_parts, buckets , wrap , reduce );
}

template< class Algorithm >
inline
void AlgorithmRunnerInterface::run(
  const mesh::Selector                & selector ,
  const mesh::PartVector              & union_parts ,
  const std::vector< mesh::Bucket * > & buckets ,
  const Algorithm                     & algorithm ) const
{
  const AlgorithmWrapper<Algorithm> wrap( algorithm );

  run_alg( selector , union_parts, buckets , wrap , NULL );
}

template< class Algorithm >
inline
void AlgorithmRunnerInterface::run_parts(
  const mesh::Selector                & selector ,
  const mesh::PartVector              & union_parts ,
  const std::vector< mesh::Bucket * > & buckets ,
  const Algorithm                     & algorithm ) const
{
  const AlgorithmWrapperParts<Algorithm> wrap( algorithm );

  run_alg( selector , union_parts, buckets , wrap , NULL );
}

//----------------------------------------------------------------------

} //namespace stk_classic

#endif


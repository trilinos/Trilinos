/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/DistributedIndex.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

class UnitTestSTKParallelDistributedIndex {
 public:
  static void test_ctor();
  static void test_ctor_bad();
  static void test_update();
  static void test_update_bad();
  static void test_generate();
  static void test_generate_bad();
  static void test_generate_big();
  static void test_update_generate();
};

namespace {

STKUNIT_UNIT_TEST( UnitTestDistributedIndexConstructor , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_ctor(); }

STKUNIT_UNIT_TEST( UnitTestDistributedIndexConstructorBad , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_ctor_bad(); }

STKUNIT_UNIT_TEST( UnitTestDistributedIndexUpdate , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_update(); }

STKUNIT_UNIT_TEST( UnitTestDistributedIndexUpdateBad , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_update_bad(); }
  
STKUNIT_UNIT_TEST( UnitTestDistributedIndexGenerate , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_generate(); }

STKUNIT_UNIT_TEST( UnitTestDistributedIndexGenerateBad , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_generate_bad(); }

STKUNIT_UNIT_TEST( UnitTestDistributedIndexUpdateGenerate , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_update_generate(); }

STKUNIT_UNIT_TEST( UnitTestDistributedIndexGenerateBig , testUnit )
{ UnitTestSTKParallelDistributedIndex::test_generate_big(); }

// Generate spans:
//   [    1..10000]
//   [20001..30000]
//   [40001..50000]
//   [60001..70000]
//   [80001..90000]
//      etc.
void generate_test_spans_10x10000(
  std::vector< stk::parallel::DistributedIndex::KeySpan > & partition_spans )
{
  enum { test_spans_count = 10 };
  enum { test_spans_size  = 10000 };

  partition_spans.resize( test_spans_count );

  for ( unsigned i = 0 ; i < test_spans_count ; ++i ) {
    partition_spans[i].first  = 1 + test_spans_size * i * 2 ;
    partition_spans[i].second = test_spans_size * ( i * 2 + 1 );
  }
}

}

//----------------------------------------------------------------------

void UnitTestSTKParallelDistributedIndex::test_ctor()
{
  typedef stk::parallel::DistributedIndex PDIndex ;
  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  int mpi_rank = stk::parallel_machine_rank(comm);
  int mpi_size = stk::parallel_machine_size(comm);

  std::vector< PDIndex::KeySpan > partition_spans ;

  generate_test_spans_10x10000( partition_spans );

  PDIndex di( comm , partition_spans );

  STKUNIT_EXPECT_EQ(   di.m_comm_rank , mpi_rank );
  STKUNIT_EXPECT_EQ(   di.m_comm_size , mpi_size );
  STKUNIT_EXPECT_TRUE( di.m_key_usage.empty() );
  STKUNIT_ASSERT_EQ(   di.m_key_span.size() , partition_spans.size() );
  for ( size_t i = 0 ; i < di.m_key_span.size() ; ++i ) {
    STKUNIT_EXPECT_EQ( di.m_key_span[i].first , partition_spans[i].first );
    STKUNIT_EXPECT_EQ( di.m_key_span[i].second , partition_spans[i].second );
  }

  // All queries will be empty:

  std::vector<PDIndex::KeyType> keys_to_query ;
  std::vector<PDIndex::KeyProc> sharing_of_local_keys ;

  di.query( sharing_of_local_keys );

  STKUNIT_EXPECT_TRUE( sharing_of_local_keys.empty() );

  di.query( keys_to_query , sharing_of_local_keys );

  STKUNIT_EXPECT_TRUE( sharing_of_local_keys.empty() );

  keys_to_query.push_back( 10 );

  di.query( keys_to_query , sharing_of_local_keys );

  STKUNIT_EXPECT_TRUE( sharing_of_local_keys.empty() );
}

void UnitTestSTKParallelDistributedIndex::test_ctor_bad()
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  {
    // Throw for overlapping span
    std::vector< PDIndex::KeySpan > partition_spans ;

    generate_test_spans_10x10000( partition_spans );

    // Corrupt this span to trigger an error
    partition_spans[5].first = partition_spans[4].second ;

    STKUNIT_ASSERT_THROW( PDIndex di( comm , partition_spans ) , std::runtime_error );
  }

  {
    // Throw for one bad span
    std::vector< PDIndex::KeySpan > partition_spans ;

    generate_test_spans_10x10000( partition_spans );

    // Corrupt this span to trigger an error
    std::swap( partition_spans[5].first , partition_spans[5].second );

    STKUNIT_ASSERT_THROW( PDIndex( comm , partition_spans ) , std::runtime_error );
  }
}

//----------------------------------------------------------------------

void UnitTestSTKParallelDistributedIndex::test_update()
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  int mpi_rank = stk::parallel_machine_rank(comm);
  int mpi_size = stk::parallel_machine_size(comm);

  std::vector< PDIndex::KeySpan > partition_spans ;

  generate_test_spans_10x10000( partition_spans );

  PDIndex di( comm , partition_spans );

  std::vector<PDIndex::KeyType> keys_to_add ;
  std::vector<PDIndex::KeyType> keys_to_remove ;
  std::vector<PDIndex::KeyProc> sharing_of_local_keys ;

  //------------------------------
  // Update nothing:

  di.update_keys( keys_to_add , keys_to_remove );

  STKUNIT_EXPECT_TRUE( di.m_key_usage.empty() );

  //------------------------------
  // Update one key on all processes and
  // one key unique to each process.

  keys_to_add.push_back( partition_spans[0].first + 1 );
  keys_to_add.push_back( partition_spans[1].first + 2 + mpi_rank );

  di.update_keys( keys_to_add , keys_to_remove );

  di.query( sharing_of_local_keys );

  // First key shared by all processes
  // Second key shared just by this process
  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , size_t(mpi_size + 1) );

  di.query( keys_to_add , sharing_of_local_keys );

  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , size_t(mpi_size + 1) );

  //------------------------------
  // Repeat the update, should result in no changes.

  di.update_keys( keys_to_add , keys_to_remove );

  di.query( sharing_of_local_keys );

  // First key shared by all processes
  // Second key shared just by this process
  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , size_t(mpi_size + 1) );

  //------------------------------

  keys_to_remove.clear();
  keys_to_add.clear();

  keys_to_remove.push_back( partition_spans[0].second );

  di.update_keys( keys_to_add , keys_to_remove );

  di.query( sharing_of_local_keys );

  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , size_t(mpi_size + 1) );

  //------------------------------
  // Remove shared key

  keys_to_remove.clear();
  keys_to_add.clear();

  keys_to_remove.push_back( partition_spans[0].first + 1 );

  di.update_keys( keys_to_add , keys_to_remove );

  di.query( sharing_of_local_keys );

  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , size_t(1) );

  //------------------------------
  // Add two shared-by-all

  keys_to_remove.clear();
  keys_to_add.clear();

  keys_to_add.push_back( partition_spans[0].first + 1 );
  keys_to_add.push_back( partition_spans[0].first + 2 );

  di.update_keys( keys_to_add , keys_to_remove );

  di.query( sharing_of_local_keys );

  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , size_t(2*mpi_size + 1) );

  di.query( keys_to_add , sharing_of_local_keys );

  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , size_t(2*mpi_size) );

  //------------------------------

  keys_to_remove.clear();
  keys_to_add.clear();

  // Shared by even rank processes:
  if ( 0 == mpi_rank % 2 ) {
    keys_to_add.push_back( partition_spans[2].first );
  }

  di.update_keys( keys_to_add , keys_to_remove );

  di.query( sharing_of_local_keys );

  {
    size_t expected = 2 * mpi_size + 1 ;
    if ( 0 == mpi_rank % 2 ) { expected += ( mpi_size + 1 ) / 2 ; }
    STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , expected );
  }

}

//----------------------------------------------------------------------

void UnitTestSTKParallelDistributedIndex::test_update_bad()
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  int mpi_rank = stk::parallel_machine_rank(comm);
  int mpi_size = stk::parallel_machine_size(comm);

  std::vector< PDIndex::KeySpan > partition_spans ;

  generate_test_spans_10x10000( partition_spans );

  PDIndex di( comm , partition_spans );

  std::vector<PDIndex::KeyType> keys_to_add ;
  std::vector<PDIndex::KeyType> keys_to_remove ;
  std::vector<PDIndex::KeyProc> sharing_of_local_keys ;

  //------------------------------
  // Invalid key on every process

  keys_to_add.push_back( partition_spans[0].second + 1 + mpi_rank );

  STKUNIT_ASSERT_THROW( di.update_keys( keys_to_add , keys_to_remove ) , std::runtime_error );

  //------------------------------

  keys_to_add.clear();
  if ( mpi_size == mpi_rank + 1 ) {
    keys_to_add.push_back( partition_spans[0].second + 1 );
  }

  STKUNIT_ASSERT_THROW( di.update_keys( keys_to_add , keys_to_remove ) , std::runtime_error );
}

//----------------------------------------------------------------------

void UnitTestSTKParallelDistributedIndex::test_generate()
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  int mpi_rank = stk::parallel_machine_rank(comm);
  // int mpi_size = stk::parallel_machine_size(comm);

  std::vector< PDIndex::KeySpan > partition_spans ;

  generate_test_spans_10x10000( partition_spans );

  PDIndex di( comm , partition_spans );

  std::vector<size_t> requests( partition_spans.size() , size_t(0) );
  std::vector< std::vector<PDIndex::KeyType> > generated_keys ;
  std::vector<PDIndex::KeyProc> sharing_of_local_keys ;

  di.generate_new_keys( requests , generated_keys );

  STKUNIT_EXPECT_TRUE( di.m_key_usage.empty() );

  STKUNIT_ASSERT_EQ( generated_keys.size() , partition_spans.size() );

  for ( size_t i = 0 ; i < generated_keys.size() ; ++i ) {
    STKUNIT_EXPECT_TRUE( generated_keys[i].empty() );
  }

  //----------------------------------------

  size_t total = 0 ;
  for ( size_t i = 0 ; i < requests.size() ; ++i ) {
    requests[i] = i + mpi_rank ;
    total += requests[i] ;
  }

  di.generate_new_keys( requests , generated_keys );

  STKUNIT_ASSERT_EQ( generated_keys.size() , partition_spans.size() );

  for ( size_t i = 0 ; i < generated_keys.size() ; ++i ) {
    STKUNIT_EXPECT_EQ( generated_keys[i].size() , requests[i] );
    for ( size_t j = 0 ; j < generated_keys[i].size() ; ++j ) {
      STKUNIT_EXPECT_TRUE( partition_spans[i].first <= generated_keys[i][j] );
      STKUNIT_EXPECT_TRUE( generated_keys[i][j] <= partition_spans[i].second );
      if ( 0 < j ) {
        STKUNIT_EXPECT_TRUE( generated_keys[i][j-1] < generated_keys[i][j] );
      }
    }
  }

  //----------------------------------------

  di.query( sharing_of_local_keys );

  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , total );

  // Confirm global uniqueness

  for ( size_t i = 0 ; i < generated_keys.size() ; ++i ) {
    di.query( generated_keys[i] , sharing_of_local_keys );
    STKUNIT_EXPECT_EQ( generated_keys[i].size() , sharing_of_local_keys.size() );
    for ( size_t j = 0 ; j < sharing_of_local_keys.size() ; ++j ) {
      STKUNIT_EXPECT_EQ( sharing_of_local_keys[j].second , mpi_rank );
    }
  }

  //----------------------------------------
  // Double the number of keys

  di.generate_new_keys( requests , generated_keys );

  STKUNIT_ASSERT_EQ( generated_keys.size() , partition_spans.size() );

  for ( size_t i = 0 ; i < generated_keys.size() ; ++i ) {
    STKUNIT_EXPECT_EQ( generated_keys[i].size() , requests[i] );
    for ( size_t j = 0 ; j < generated_keys[i].size() ; ++j ) {
      STKUNIT_EXPECT_TRUE( partition_spans[i].first <= generated_keys[i][j] );
      STKUNIT_EXPECT_TRUE( generated_keys[i][j] <= partition_spans[i].second );
      if ( 0 < j ) {
        STKUNIT_EXPECT_TRUE( generated_keys[i][j-1] < generated_keys[i][j] );
      }
    }
  }

  di.query( sharing_of_local_keys );

  STKUNIT_EXPECT_EQ( sharing_of_local_keys.size() , total * 2 );
}

//----------------------------------------------------------------------

void UnitTestSTKParallelDistributedIndex::test_update_generate()
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  int p_rank = stk::parallel_machine_rank(comm);
  int p_size = stk::parallel_machine_size(comm);

  std::vector< PDIndex::KeySpan > partition_spans ;

  generate_test_spans_10x10000( partition_spans );

  PDIndex di( comm , partition_spans );

  std::vector<size_t> requests( partition_spans.size() , size_t(0) );
  std::vector< std::vector<PDIndex::KeyType> > generated_keys ;
  std::vector<PDIndex::KeyProc> sharing_of_local_keys ;

  std::vector<PDIndex::KeyType> keys_to_add ;
  std::vector<PDIndex::KeyType> keys_to_remove ;

  //------------------------------
  // Add ( 5 * j ) odd keys per process
  // starting at the beginning of the partition.

  const size_t old_size_multiplier = 5 ;

  for ( size_t j = 0 ; j < partition_spans.size() ; ++j ) {
    PDIndex::KeyType key_first = partition_spans[j].first ;
    if ( 0 == key_first % 2 ) { ++key_first ; } // Make it odd
    key_first += old_size_multiplier * p_rank ;

    const size_t n = old_size_multiplier * j ;
    for ( size_t i = 0 ; i < n ; ++i ) {
      PDIndex::KeyType key = key_first + 2 * i ;
      keys_to_add.push_back( key );
    }
  }

  di.update_keys( keys_to_add , keys_to_remove );

  //------------------------------
  // Request 20 new keys per process per span
  // The maximum new key will be larger than some spans
  // and within the gaps of other spans.

  const size_t gen_count = 20 ;
  for ( size_t i = 0 ; i < requests.size() ; ++i ) {
    if ( i % 2 ) {
      requests[i] = gen_count ;
    }
    else {
      requests[i] = 0 ;
    }
  }

  di.generate_new_keys( requests , generated_keys );

  for ( size_t i = 0 ; i < requests.size() ; ++i ) {
    STKUNIT_EXPECT_EQ( requests[i] , generated_keys[i].size() );

    const size_t old_count   = p_size * old_size_multiplier * i ;
    const size_t tot_request = p_size * requests[i] ;

    PDIndex::KeyType max_gen_key = partition_spans[i].first ;

    if ( 0 == tot_request ) {
      STKUNIT_EXPECT_TRUE( generated_keys[i].size() == 0 );
    }
    else if ( tot_request < old_count ) {
      // Will only fill in gaps between odd keys
      max_gen_key += 2 * old_count ;

      STKUNIT_EXPECT_TRUE( max_gen_key > generated_keys[i][ requests[i] - 1 ] );
    }
    else {
      // Will fill in gaps contiguously after the old max key
      max_gen_key += old_count + tot_request - 1 ;

      STKUNIT_EXPECT_TRUE( max_gen_key >= generated_keys[i][ requests[i] - 1 ] );
    }

    // Sorted
    for ( size_t j = 0 ; j < generated_keys[i].size() ; ++j ) {
      if ( 0 < j ) {
        STKUNIT_EXPECT_TRUE( generated_keys[i][j-1] < generated_keys[i][j] );
      }
    }
  }
}

//----------------------------------------------------------------------

void UnitTestSTKParallelDistributedIndex::test_generate_bad()
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  int mpi_rank = stk::parallel_machine_rank(comm);
  int mpi_size = stk::parallel_machine_size(comm);

  std::vector< PDIndex::KeySpan > partition_spans ;

  generate_test_spans_10x10000( partition_spans );

  PDIndex di( comm , partition_spans );

  std::vector<size_t> requests( partition_spans.size() , size_t(0) );
  std::vector< std::vector<PDIndex::KeyType> > generated_keys ;
  //------------------------------

  for ( size_t i = 0 ; i < requests.size() ; ++i ) {
    requests[i] = 10 ;
  }

  if( mpi_rank == mpi_size -1 ) {
    requests.clear();
  }

  STKUNIT_ASSERT_THROW( di.generate_new_keys( requests , generated_keys ), std::runtime_error );

  if( mpi_rank == mpi_size -1 ) {
    requests.push_back(2*(partition_spans[0].second - partition_spans[0].first));
  }

  STKUNIT_ASSERT_THROW( di.generate_new_keys( requests , generated_keys ), std::runtime_error );

  for ( size_t i = 0 ; i < requests.size() ; ++i ) {
    requests[i] = partition_spans[i].second - partition_spans[i].first ;
  }

  STKUNIT_ASSERT_THROW( di.generate_new_keys( requests , generated_keys ), std::runtime_error );

}

//----------------------------------------------------------------------

void UnitTestSTKParallelDistributedIndex::test_generate_big()
{
  typedef stk::parallel::DistributedIndex PDIndex ;

  stk::ParallelMachine comm = MPI_COMM_WORLD ;

  const int mpi_rank = stk::parallel_machine_rank(comm);
  const int mpi_size = stk::parallel_machine_size(comm);

  std::vector< PDIndex::KeySpan > partition_spans ;

  generate_test_spans_10x10000( partition_spans );

  PDIndex di( comm , partition_spans );

  std::vector<size_t> requests( partition_spans.size() , size_t(0) );

  std::vector< std::vector<PDIndex::KeyType> > generated_keys ;

  //----------------------------------------

  if ( mpi_rank == mpi_size - 1 ) {
    requests[5] = 5000 ; // Half
  }
  else {
    requests[5] = 0 ;
  }

  di.generate_new_keys( requests , generated_keys );

  STKUNIT_ASSERT_EQ( generated_keys.size() , partition_spans.size() );

  for ( size_t i = 0 ; i < generated_keys.size() ; ++i ) {
    STKUNIT_EXPECT_EQ( generated_keys[i].size() , requests[i] );
  }

  //----------------------------------------

  if ( mpi_rank == mpi_size - 1 ) {
    requests[5] = 4999 ; // Almost half again
  }
  else {
    requests[5] = 0 ;
  }

  di.generate_new_keys( requests , generated_keys );

  STKUNIT_ASSERT_EQ( generated_keys.size() , partition_spans.size() );

  for ( size_t i = 0 ; i < generated_keys.size() ; ++i ) {
    STKUNIT_EXPECT_EQ( generated_keys[i].size() , requests[i] );
  }

  //----------------------------------------
  // This should request generation of one too many keys.

  if ( mpi_rank == 0 ) {
    requests[5] = 2 ; // Exceed the total
  }
  else {
    requests[5] = 0 ;
  }

  STKUNIT_ASSERT_THROW( di.generate_new_keys( requests , generated_keys ) , std::runtime_error );

}


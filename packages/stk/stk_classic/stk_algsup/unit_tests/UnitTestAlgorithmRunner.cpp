/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_mesh/base/GetEntities.hpp>

#include <stk_algsup/AlgorithmRunner.hpp>

#include <unit_tests/UnitTest_helpers.hpp>

static const stk_classic::mesh::EntityRank NODE_RANK = stk_classic::mesh::fem::FEMMetaData::NODE_RANK;

class MyAlg {
public:
  MyAlg() : maximum_entity_count(0), num_applys(0) {}
  virtual ~MyAlg() {}

  void apply(stk_classic::mesh::Bucket::iterator /*ibeg*/,
             stk_classic::mesh::Bucket::iterator /*iend*/) const
  {
    ++num_applys;
  }

  unsigned maximum_entity_count;
  mutable unsigned num_applys;
};//class MyAlg

STKUNIT_UNIT_TEST( UnitTestAlgorithmRunner, UnitTest)
{
#ifdef STK_HAS_MPI
  stk_classic::ParallelMachine comm(MPI_COMM_WORLD);
#else
  stk_classic::ParallelMachine comm(0);
#endif

  const stk_classic::AlgorithmRunnerInterface* alg_runner = NULL;

#ifdef STK_HAVE_TBB
  const int num_threads = 1;
  alg_runner = stk_classic::algorithm_runner_tbb(num_threads);
#elif defined(STK_HAVE_TPI)
  const int num_threads = 1;
  alg_runner = stk_classic::algorithm_runner_tpi(num_threads);
#else
  alg_runner = stk_classic::algorithm_runner_non_thread();
#endif

  STKUNIT_ASSERT( alg_runner != NULL );

  const stk_classic::AlgorithmRunnerInterface* alg_run_nothread = stk_classic::algorithm_runner_non_thread();

  STKUNIT_ASSERT( alg_run_nothread != NULL );

  const unsigned spatial_dim = 3;
  stk_classic::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dim);

  stk_classic::mesh::MetaData & meta_data = stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta);

  stk_classic::mesh::BulkData bulk_data( meta_data, comm );

  fill_utest_mesh_meta_data( fem_meta );
  fill_utest_mesh_bulk_data( bulk_data );

  stk_classic::mesh::Selector selector = meta_data.locally_owned_part() | meta_data.globally_shared_part();
  std::vector<stk_classic::mesh::Bucket*> used_node_buckets;
  stk_classic::mesh::get_buckets(selector, bulk_data.buckets(NODE_RANK), used_node_buckets);

  MyAlg my_alg;

  stk_classic::mesh::PartVector empty_union_vector;
  alg_run_nothread->run(selector, empty_union_vector, bulk_data.buckets(NODE_RANK), my_alg);

  unsigned num_applys = my_alg.num_applys;
  unsigned num_applys_expected = used_node_buckets.size();

  STKUNIT_ASSERT_EQUAL( num_applys, num_applys_expected );

#if defined(STK_HAVE_TBB) || defined(STK_HAVE_TPI)
  my_alg.num_applys = 0;

  alg_runner->run(selector, empty_union_vector, bulk_data.buckets(NODE_RANK), my_alg);

  num_applys = my_alg.num_applys;
  STKUNIT_ASSERT_EQUAL( num_applys, num_applys_expected );
#endif

  my_alg.maximum_entity_count = 1000;
  my_alg.num_applys = 0;

  alg_run_nothread->run(selector, empty_union_vector, bulk_data.buckets(NODE_RANK), my_alg);

  num_applys = my_alg.num_applys;
  num_applys_expected = used_node_buckets.size();

  STKUNIT_ASSERT_EQUAL( num_applys, num_applys_expected );

#if defined(STK_HAVE_TBB) || defined(STK_HAVE_TPI)
  my_alg.num_applys = 0;
  my_alg.num_applys = 0;

  alg_runner->run(selector, empty_union_vector, bulk_data.buckets(NODE_RANK), my_alg);

  num_applys = my_alg.num_applys;
  STKUNIT_ASSERT_EQUAL( num_applys, num_applys_expected );
#endif

}

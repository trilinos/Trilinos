/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/util/perf_util.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>

#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <sstream>

//This very simple test will visit all local elements and traverse the
//element's node relations. It will return the number of nodes visited.
//The purpose of this test is to stress the relation-traversal for a
//performance test.
size_t do_stk_node_rel_test(stk::mesh::BulkData& bulk)
{
  using namespace stk::mesh;

  MetaData& meta = MetaData::get(bulk);

  Selector local = meta.locally_owned_part();

  BucketVector buckets;
  get_buckets(local, bulk.buckets(stk::mesh::MetaData::ELEMENT_RANK), buckets);

  size_t nodes_visited = 0;
  unsigned owner_rank = 0;

  size_t num_elems = 0;
  for(size_t ib=0; ib<buckets.size(); ++ib) {
    const Bucket& b = *buckets[ib];
    num_elems += b.size();

    for(size_t i=0; i<b.size(); ++i) {
      Entity elem = b[i];
      PairIterRelation node_rels = elem.node_relations();

      for(; !node_rels.empty(); ++node_rels) {
        Entity node = node_rels->entity();
        owner_rank += node.owner_rank();
        ++nodes_visited;
      }
    }
  }
  return nodes_visited;
}

STKUNIT_UNIT_TEST(node_rels, node_rels)
{
  //vector of mesh-dimensions holds the number of elements in each dimension.
  //Hard-wired to 3. This test can run with spatial-dimension less than 3,
  //(if generated-mesh can do that) but not greater than 3.
  std::vector<int> mesh_dims(3);
#ifndef NDEBUG
  mesh_dims[0]=50; //num_elems_x
  mesh_dims[1]=50; //num_elems_y
  mesh_dims[2]=50; //num_elems_z
#else
  mesh_dims[0]=100; //num_elems_x
  mesh_dims[1]=100; //num_elems_y
  mesh_dims[2]=100; //num_elems_z
#endif

  std::ostringstream oss;
  oss << mesh_dims[0] << "x" << mesh_dims[1] << "x" << mesh_dims[2] << "|sideset:xXyYzZ";

  double start_time = stk::cpu_time();

  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, oss.str());
  fixture.commit();

  double mesh_create_time = stk::cpu_time() - start_time;

  //compute total number of elements in the mesh:
  const size_t spatial_dim = fixture.getMetaData().spatial_dimension();
  size_t num_elems = mesh_dims[0];
  for(size_t d=1; d<spatial_dim; ++d) {
    num_elems *= mesh_dims[d];
  }

  std::cout << "num_elems: " << num_elems << std::endl;
  std::cout << "sizeof(stk::mesh::Relation): " << sizeof(stk::mesh::Relation) << std::endl;
  std::cout << "sizeof(stk::mesh::Entity): " << sizeof(stk::mesh::Entity) << std::endl;

  start_time = stk::cpu_time();

  const int num_iters = 50;
  for(int t=0; t<num_iters; ++t) {

    size_t nodes_visited = do_stk_node_rel_test(fixture.getBulkData());
    size_t expected_nodes_visited = num_elems*8;

    EXPECT_EQ(nodes_visited, expected_nodes_visited);
  }

  double traverse_time = stk::cpu_time() - start_time;
  double total_time    = mesh_create_time + traverse_time;

  static const int NUM_TIMERS = 3;
  const double timers[NUM_TIMERS] = {mesh_create_time, traverse_time, total_time};
  const char* timer_names[NUM_TIMERS] = {"Create mesh", "Traverse", "Total time"};

  stk::print_timers_and_memory(&timer_names[0], &timers[0], NUM_TIMERS);

}

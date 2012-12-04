#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_topology/topology.hpp>

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

STKUNIT_UNIT_TEST( stk_topology, pretty_print_ranks)
{
  using stk::topology;

  std::vector<std::string> rank_names(topology::NUM_RANKS);
  std::vector<std::string> output_rank_names(topology::NUM_RANKS);

  for(int i=0; i<topology::NUM_RANKS; ++i)
    rank_names[i] = topology::rank_names[i];


  for(topology::rank_t r = topology::BEGIN_RANK; r < topology::END_RANK; ++r)
  {
    std::ostringstream name;
    name << r;

    output_rank_names[r] = name.str();
  }


  STKUNIT_EXPECT_TRUE( std::equal(rank_names.begin(), rank_names.end(), output_rank_names.begin()) );

}

STKUNIT_UNIT_TEST( stk_topology, pretty_print_topologies)
{
  using stk::topology;

  std::vector<std::string> topology_names(topology::NUM_TOPOLOGIES);
  std::vector<std::string> output_topology_names(topology::NUM_TOPOLOGIES);

  for(int i=0; i<topology::NUM_TOPOLOGIES; ++i)
    topology_names[i] = topology::topology_names[i];


  for(topology t = topology::BEGIN_TOPOLOGY; t < topology::END_TOPOLOGY; ++t)
  {
    std::ostringstream name;
    name << t;

    output_topology_names[t] = name.str();
  }

  STKUNIT_EXPECT_TRUE( std::equal(topology_names.begin(), topology_names.end(), output_topology_names.begin()) );

}


STKUNIT_UNIT_TEST( stk_topology, lexicographical_smallest_permutation)
{
  using stk::topology;

  topology t = topology::TRI_3;

  const char nodes[]="bac";

  char permutation_nodes[4] = "";

  int permutation_index = t.lexicographical_smallest_permutation(nodes);
  t.permutation_nodes(nodes,permutation_index,permutation_nodes);

  STKUNIT_EXPECT_EQ( std::string("abc"), std::string(permutation_nodes));

  permutation_index = t.lexicographical_smallest_permutation(nodes,true);
  t.permutation_nodes(nodes,permutation_index,permutation_nodes);

  STKUNIT_EXPECT_EQ( std::string("acb"), std::string(permutation_nodes));
}

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_topology/topology.hpp>
#include <stk_topology/pretty_print.hpp>

#include <iostream>

STKUNIT_MAIN(argc,argv)

STKUNIT_UNIT_TEST( stk_topology, topology)
{
  for (stk::topology t = stk::topology::BEGIN_TOPOLOGY; t < stk::topology::END_TOPOLOGY; ++t)
  {
    stk::verbose_print_topology(std::cout,t);
  }
}

STKUNIT_UNIT_TEST( stk_topology, edge_nodes)
{
  int nodes[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

  int edge[2] = {};

  stk::topology t = stk::topology::HEX_8;

  for (unsigned i=0u, e=t.num_edges(); i<e; ++i)
  {
    t.edge_nodes(nodes,i,edge);

    std::cout << "edge " << i << ": (" << edge[0] << ", " << edge[1] << ")" << std::endl;
  }

}

STKUNIT_UNIT_TEST( stk_topology, face_nodes)
{
  int nodes[] = { 1, 2, 3, 4, 5, 6, 7, 8 };

  int face[4] = {};

  stk::topology t = stk::topology::HEX_8;

  for (unsigned i=0u, e=t.num_faces(); i<e; ++i)
  {
    t.face_nodes(nodes,i,face);

    std::cout << "face " << i << ": (" << face[0] << ", "
                                       << face[1] << ", "
                                       << face[2] << ", "
                                       << face[3] << ")" << std::endl;
  }

}

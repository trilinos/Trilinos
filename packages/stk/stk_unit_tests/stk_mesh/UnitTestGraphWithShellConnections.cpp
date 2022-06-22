#include "gtest/gtest.h"
#include <vector>
#include "mpi.h"
#include <stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemGraphShellConnections.hpp>

namespace
{

struct ProcAndId
{
  ProcAndId(int p, stk::mesh::impl::LocalId i) : proc(p), id(i) { }
  int proc;
  stk::mesh::impl::LocalId id;
};

class GraphWithShells : public ::testing::Test
{
protected:
  GraphWithShells()
    : procRank(stk::parallel_machine_rank(comm)),
      parGraphInfo(procRank),
      graphInfo(graph, parGraphInfo, elementTopologies)
  {
  }

  void run_test_on_num_procs(int numProcs)
  {
    if(stk::parallel_machine_size(comm) == numProcs)
      run_test();
  }

  virtual void run_test() = 0;

  void set_element_topologies_per_proc(const std::vector<std::vector<stk::topology>> &tops)
  {
    elementTopologies = tops[procRank];
  }
  void create_graph_edges_for_hex0_shell1_hex2_on_procs(int hex0Proc, int shell1Proc, int hex2Proc)
  {
    add_edges_on_proc(hex0Proc, graphEdgeHex0Shell1, graphEdgeHex0Hex2);
    add_edges_on_proc(shell1Proc, graphEdgeShell1Hex0, graphEdgeShell1Hex2);
    add_edges_on_proc(hex2Proc, graphEdgeHex2Hex0, graphEdgeHex2Shell1);
  }
  void expect_graph_edges_for_hex0_shell1_hex2_on_procs(ProcAndId hex0, ProcAndId shell1, ProcAndId hex2)
  {
    expect_graph_has_edges_for_element(hex0, {graphEdgeHex0Shell1});
    expect_graph_has_edges_for_element(shell1, {graphEdgeShell1Hex0, graphEdgeShell1Hex2});
    expect_graph_has_edges_for_element(hex2, {graphEdgeHex2Shell1});
  }

private:
  void add_edges_on_proc(int creationProc, const stk::mesh::GraphEdge &graphEdge1, const stk::mesh::GraphEdge &graphEdge2)
  {
    if(procRank == creationProc)
      add_two_edges(graphEdge1, graphEdge2);
  }
  void add_two_edges(const stk::mesh::GraphEdge& graphEdge1, const stk::mesh::GraphEdge& graphEdge2)
  {
    graph.add_new_element();
    std::vector<stk::mesh::GraphEdge> edges = {graphEdge1, graphEdge2};
    graph.add_sorted_edges(edges);
  }
  void expect_graph_has_edges_for_element(ProcAndId procAndId, const std::vector<stk::mesh::GraphEdge> &expectedEdges)
  {
    if(procRank == procAndId.proc)
    {
      ASSERT_EQ(expectedEdges.size(), graph.get_num_edges_for_element(procAndId.id));
      for(size_t i=0; i<expectedEdges.size(); i++)
        EXPECT_EQ(expectedEdges[i], graph.get_edge_for_element(procAndId.id, i));
    }
  }

protected:
  MPI_Comm comm = MPI_COMM_WORLD;
  int procRank;
  stk::mesh::Graph graph;
  stk::mesh::ParallelInfoForGraphEdges parGraphInfo;
  std::vector<stk::topology> elementTopologies;
  stk::mesh::GraphInfo graphInfo;
  stk::mesh::GraphEdge graphEdgeHex0Hex2;
  stk::mesh::GraphEdge graphEdgeHex0Shell1;
  stk::mesh::GraphEdge graphEdgeShell1Hex0;
  stk::mesh::GraphEdge graphEdgeShell1Hex2;
  stk::mesh::GraphEdge graphEdgeHex2Hex0;
  stk::mesh::GraphEdge graphEdgeHex2Shell1;
};

class GraphWithShellsSerial : public GraphWithShells
{
protected:
  GraphWithShellsSerial()
  {
    graphEdgeHex0Hex2 = stk::mesh::GraphEdge(0, 0, 2, 4);
    graphEdgeHex0Shell1 = stk::mesh::GraphEdge(0, 0, 1, 0);
    graphEdgeShell1Hex0 = stk::mesh::GraphEdge(1, 0, 0, 0);
    graphEdgeShell1Hex2 = stk::mesh::GraphEdge(1, 1, 2, 4);
    graphEdgeHex2Hex0 = stk::mesh::GraphEdge(2, 4, 0, 0);
    graphEdgeHex2Shell1 = stk::mesh::GraphEdge(2, 4, 1, 1);
  }

  virtual void run_test()
  {
    set_element_topologies_per_proc({{stk::topology::HEX_8, stk::topology::SHELL_QUAD_4, stk::topology::HEX_8}});
    create_graph_edges_for_hex0_shell1_hex2_on_procs(0, 0, 0);
    remove_graph_edges_blocked_by_shell(graphInfo);
    expect_graph_edges_for_hex0_shell1_hex2_on_procs(ProcAndId(0, 0), ProcAndId(0, 1), ProcAndId(0, 2));
  }
};
TEST_F(GraphWithShellsSerial, HexShellHex)
{
  run_test_on_num_procs(1);
}

class GraphWithShells2Procs : public GraphWithShells
{
protected:
  GraphWithShells2Procs()
  {
    graphEdgeHex0Hex2 = stk::mesh::GraphEdge(0, 0, -3, 4);
    graphEdgeHex0Shell1 = stk::mesh::GraphEdge(0, 0, 1, 0);
    graphEdgeShell1Hex0 = stk::mesh::GraphEdge(1, 0, 0, 0);
    graphEdgeShell1Hex2 = stk::mesh::GraphEdge(1, 1, -3, 4);
    graphEdgeHex2Hex0 = stk::mesh::GraphEdge(0, 4, -1, 0);
    graphEdgeHex2Shell1 = stk::mesh::GraphEdge(0, 4, -2, 1);
  }

  virtual void run_test()
  {
    setup_graph_for_two_procs();
    remove_graph_edges_blocked_by_shell(graphInfo);
    expect_graph_edges_for_hex0_shell1_hex2_on_procs(ProcAndId(0, 0), ProcAndId(0, 1), ProcAndId(1, 0));
  }

  void setup_graph_for_two_procs()
  {
    set_element_topologies_per_proc( { {stk::topology::HEX_8, stk::topology::SHELL_QUAD_4}, {stk::topology::HEX_8}});
    create_graph_edges_for_hex0_shell1_hex2_on_procs(0, 0, 1);
    setup_parallel_graph_info_for_two_procs();
  }

private:
  void setup_parallel_graph_info_for_two_procs()
  {
    setup_parallel_graph_for_proc_0();
    setup_parallel_graph_for_proc_1();
  }
  void setup_parallel_graph_for_proc_0()
  {
    if(procRank == 0)
    {
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeHex0Hex2, stk::mesh::impl::ParallelInfo(1, 0, stk::topology::HEX_8));
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeShell1Hex2, stk::mesh::impl::ParallelInfo(1, 0, stk::topology::HEX_8));
    }
  }
  void setup_parallel_graph_for_proc_1()
  {
    if(procRank == 1)
    {
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeHex2Hex0, stk::mesh::impl::ParallelInfo(0, 0, stk::topology::HEX_8));
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeHex2Shell1, stk::mesh::impl::ParallelInfo(0, 0, stk::topology::SHELL_QUAD_4));
    }
  }
};

TEST_F(GraphWithShells2Procs, HexShellHex)
{
  run_test_on_num_procs(2);
}

class GraphWithShells3Procs : public GraphWithShells
{
protected:
  GraphWithShells3Procs()
  {
    graphEdgeHex0Hex2 = stk::mesh::GraphEdge(0, 0, -3, 4);
    graphEdgeHex0Shell1 = stk::mesh::GraphEdge(0, 0, -2, 0);
    graphEdgeShell1Hex0 = stk::mesh::GraphEdge(0, 0, -1, 0);
    graphEdgeShell1Hex2 = stk::mesh::GraphEdge(0, 1, -3, 4);
    graphEdgeHex2Hex0 = stk::mesh::GraphEdge(0, 4, -1, 0);
    graphEdgeHex2Shell1 = stk::mesh::GraphEdge(0, 4, -2, 1);
  }

  virtual void run_test()
  {
    setup_graph_for_three_procs();
    remove_graph_edges_blocked_by_shell(graphInfo);
    expect_graph_edges_for_hex0_shell1_hex2_on_procs(ProcAndId(0, 0), ProcAndId(1, 0), ProcAndId(2, 0));
  }

private:
  void setup_graph_for_three_procs()
  {
    set_element_topologies_per_proc( { {stk::topology::HEX_8}, {stk::topology::SHELL_QUAD_4}, {stk::topology::HEX_8}});
    create_graph_edges_for_hex0_shell1_hex2_on_procs(0, 1, 2);
    setup_parallel_graph_info_for_three_procs();
  }
  void setup_parallel_graph_info_for_three_procs()
  {
    setup_parallel_graph_for_proc_0();
    setup_parallel_graph_for_proc_1();
    setup_parallel_graph_for_proc_2();
  }
  void setup_parallel_graph_for_proc_0()
  {
    if(procRank == 0)
    {
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeHex0Shell1, stk::mesh::impl::ParallelInfo(1, 0, stk::topology::SHELL_QUAD_4));
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeHex0Hex2, stk::mesh::impl::ParallelInfo(2, 0, stk::topology::HEX_8));
    }
  }
  void setup_parallel_graph_for_proc_1()
  {
    if(procRank == 1)
    {
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeShell1Hex0, stk::mesh::impl::ParallelInfo(0, 0, stk::topology::HEX_8));
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeShell1Hex2, stk::mesh::impl::ParallelInfo(2, 0, stk::topology::HEX_8));
    }
  }
  void setup_parallel_graph_for_proc_2()
  {
    if(procRank == 2)
    {
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeHex2Hex0, stk::mesh::impl::ParallelInfo(0, 0, stk::topology::HEX_8));
      parGraphInfo.insert_parallel_info_for_graph_edge(graphEdgeHex2Shell1, stk::mesh::impl::ParallelInfo(1, 0, stk::topology::SHELL_QUAD_4));
    }
  }
};

TEST_F(GraphWithShells3Procs, HexShellHex)
{
  run_test_on_num_procs(3);
}

}





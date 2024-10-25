// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "MeshFixtureRebalance.hpp"
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_balance/rebalance.hpp>
#include <vector>

namespace {

class Rebalance : public MeshFixtureRebalance
{
public:
  void rebalance_mesh(int numFinalProcs, const std::string & decompMethod = "rcb")
  {
    m_balanceSettings.set_is_rebalancing(true);
    m_balanceSettings.set_output_filename(get_output_file_name());
    m_balanceSettings.set_num_input_processors(stk::parallel_machine_size(get_comm()));
    m_balanceSettings.set_num_output_processors(numFinalProcs);
    m_balanceSettings.setDecompMethod(decompMethod);

    stk::set_outputP0(&stk::outputNull());
    stk::balance::rebalance(m_ioBroker, m_balanceSettings);
    stk::reset_default_output_streams();

    stk::mesh::Part * block1 = m_ioBroker.meta_data().get_part("block_1");
    ASSERT_NE(block1, nullptr);
    ASSERT_TRUE(stk::io::has_original_topology_type(*block1));
  }
};

TEST_F(Rebalance, 1elem_beam2_1procTo1proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  std::vector<OriginalTopology> originalTopologies { {"block_1", "beam2"} };

  setup_initial_mesh_textmesh_override_topology("0,1,BEAM_2,1,2,block_1\n", originalTopologies);
  rebalance_mesh(1);
  test_decomposed_mesh_element_distribution({1});
  test_decomposed_mesh_node_sharing({ {} });
  test_decomposed_mesh_element_topology(originalTopologies);
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elem_assemblies_1procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  std::vector<AssemblyGrouping> assemblyDefinitions { {"assembly1", 100, {"block_1", "block_2"}} };

  setup_initial_mesh_textmesh_add_assemblies("0,1,BEAM_2,1,2,block_1\n"
                                             "0,2,BEAM_2,2,3,block_2\n", assemblyDefinitions);
  rebalance_mesh(2);
  test_decomposed_mesh_element_distribution({1, 1});
  test_decomposed_mesh_node_sharing({ {{2, 1}},
                                      {{2, 0}} });
  test_decomposed_mesh_assemblies(assemblyDefinitions);
  clean_up_temporary_files();
}

TEST_F(Rebalance, 1elem_1procTo1proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("1x1x1");
  rebalance_mesh(1);
  test_decomposed_mesh_element_distribution({1});
  test_decomposed_mesh_node_sharing({ {} });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elems_1procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("1x1x2");
  rebalance_mesh(2);
  test_decomposed_mesh_element_distribution({1, 1});
  test_decomposed_mesh_node_sharing({ {{5,1}, {6,1}, {7,1}, {8,1}},
                                      {{5,0}, {6,0}, {7,0}, {8,0}} });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elems_1disconnectedNode_1procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("1x1x2");
  const int owningProc = 0;
  const stk::mesh::EntityId disconnectedNodeId = add_disconnected_node(owningProc);

  rebalance_mesh(2);
  test_decomposed_mesh_element_distribution({1, 1});
  test_decomposed_mesh_node_sharing({ {{5,1}, {6,1}, {7,1}, {8,1}},
                                      {{5,0}, {6,0}, {7,0}, {8,0}} });
  test_node_existence(disconnectedNodeId);
//  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elems_1disconnectedNode_1procTo2proc_graphBased)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("1x1x2");
  const int owningProc = 0;
  const stk::mesh::EntityId disconnectedNodeId = add_disconnected_node(owningProc);

  rebalance_mesh(2, "parmetis");
  test_node_existence(disconnectedNodeId);
  clean_up_temporary_files();
}

TEST_F(Rebalance, 3elems_1procTo3proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("1x1x3");
  rebalance_mesh(3);
  test_decomposed_mesh_element_distribution({1, 1, 1});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, { 9,2}, {10,2}, {11,2}, {12,2}},
                                      {{ 9,1}, {10,1}, {11,1}, {12,1}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elems_2procTo1proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x2");
  rebalance_mesh(1);
  test_decomposed_mesh_element_distribution({2});
  test_decomposed_mesh_node_sharing({ {} });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elems_3procTo1proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("1x1x2");
  rebalance_mesh(1);
  test_decomposed_mesh_element_distribution({2});
  test_decomposed_mesh_node_sharing({ {} });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elems_2procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x2");
  rebalance_mesh(2);
  test_decomposed_mesh_element_distribution({1, 1});
  test_decomposed_mesh_node_sharing({ {{5,1}, {6,1}, {7,1}, {8,1}},
                                      {{5,0}, {6,0}, {7,0}, {8,0}} });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 4elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x4");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({1, 1, 1, 1});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, { 9,2}, {10,2}, {11,2}, {12,2}},
                                      {{ 9,1}, {10,1}, {11,1}, {12,1}, {13,3}, {14,3}, {15,3}, {16,3}},
                                      {{13,2}, {14,2}, {15,2}, {16,2}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 3elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x3");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({1, 1, 1, 0});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, { 9,2}, {10,2}, {11,2}, {12,2}},
                                      {{ 9,1}, {10,1}, {11,1}, {12,1}                                },
                                      {                                                              } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 5elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x5");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({2, 1, 1, 1});
  test_decomposed_mesh_node_sharing({ {                                { 9,1}, {10,1}, {11,1}, {12,1}},
                                      {{ 9,0}, {10,0}, {11,0}, {12,0}, {13,2}, {14,2}, {15,2}, {16,2}},
                                      {{13,1}, {14,1}, {15,1}, {16,1}, {17,3}, {18,3}, {19,3}, {20,3}},
                                      {{17,2}, {18,2}, {19,2}, {20,2}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 8elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x8");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 2});
  test_decomposed_mesh_node_sharing({ {                                { 9,1}, {10,1}, {11,1}, {12,1}},
                                      {{ 9,0}, {10,0}, {11,0}, {12,0}, {17,2}, {18,2}, {19,2}, {20,2}},
                                      {{17,1}, {18,1}, {19,1}, {20,1}, {25,3}, {26,3}, {27,3}, {28,3}},
                                      {{25,2}, {26,2}, {27,2}, {28,2}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 7elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x7");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 1});
  test_decomposed_mesh_node_sharing({ {                                { 9,1}, {10,1}, {11,1}, {12,1}},
                                      {{ 9,0}, {10,0}, {11,0}, {12,0}, {17,2}, {18,2}, {19,2}, {20,2}},
                                      {{17,1}, {18,1}, {19,1}, {20,1}, {25,3}, {26,3}, {27,3}, {28,3}},
                                      {{25,2}, {26,2}, {27,2}, {28,2}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 9elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x9");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({3, 2, 2, 2});
  test_decomposed_mesh_node_sharing({ {                                {13,1}, {14,1}, {15,1}, {16,1}},
                                      {{13,0}, {14,0}, {15,0}, {16,0}, {21,2}, {22,2}, {23,2}, {24,2}},
                                      {{21,1}, {22,1}, {23,1}, {24,1}, {29,3}, {30,3}, {31,3}, {32,3}},
                                      {{29,2}, {30,2}, {31,2}, {32,2}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 8elems_2procTo8proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x8");
  rebalance_mesh(8);
  test_decomposed_mesh_element_distribution({1, 1, 1, 1, 1, 1, 1, 1});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, { 9,2}, {10,2}, {11,2}, {12,2}},
                                      {{ 9,1}, {10,1}, {11,1}, {12,1}, {13,3}, {14,3}, {15,3}, {16,3}},
                                      {{13,2}, {14,2}, {15,2}, {16,2}, {17,4}, {18,4}, {19,4}, {20,4}},
                                      {{17,3}, {18,3}, {19,3}, {20,3}, {21,5}, {22,5}, {23,5}, {24,5}},
                                      {{21,4}, {22,4}, {23,4}, {24,4}, {25,6}, {26,6}, {27,6}, {28,6}},
                                      {{25,5}, {26,5}, {27,5}, {28,5}, {29,7}, {30,7}, {31,7}, {32,7}},
                                      {{29,6}, {30,6}, {31,6}, {32,6}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 6elems_3procTo6proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("1x1x6");
  rebalance_mesh(6);
  test_decomposed_mesh_element_distribution({1, 1, 1, 1, 1, 1});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, { 9,2}, {10,2}, {11,2}, {12,2}},
                                      {{ 9,1}, {10,1}, {11,1}, {12,1}, {13,3}, {14,3}, {15,3}, {16,3}},
                                      {{13,2}, {14,2}, {15,2}, {16,2}, {17,4}, {18,4}, {19,4}, {20,4}},
                                      {{17,3}, {18,3}, {19,3}, {20,3}, {21,5}, {22,5}, {23,5}, {24,5}},
                                      {{21,4}, {22,4}, {23,4}, {24,4}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 3elems_2procTo3proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x3");
  rebalance_mesh(3);
  test_decomposed_mesh_element_distribution({1, 1, 1});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, { 9,2}, {10,2}, {11,2}, {12,2}},
                                      {{ 9,1}, {10,1}, {11,1}, {12,1}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 4elems_2procTo3proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x4");
  rebalance_mesh(3);
  test_decomposed_mesh_element_distribution({1, 2, 1});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, {13,2}, {14,2}, {15,2}, {16,2}},
                                      {{13,1}, {14,1}, {15,1}, {16,1}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 6elems_2procTo7proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x6");
  rebalance_mesh(7);
  test_decomposed_mesh_element_distribution({1, 1, 1, 1, 1, 1, 0});
  test_decomposed_mesh_node_sharing({ {                                { 5,1}, { 6,1}, { 7,1}, { 8,1}},
                                      {{ 5,0}, { 6,0}, { 7,0}, { 8,0}, { 9,2}, {10,2}, {11,2}, {12,2}},
                                      {{ 9,1}, {10,1}, {11,1}, {12,1}, {13,3}, {14,3}, {15,3}, {16,3}},
                                      {{13,2}, {14,2}, {15,2}, {16,2}, {17,4}, {18,4}, {19,4}, {20,4}},
                                      {{17,3}, {18,3}, {19,3}, {20,3}, {21,5}, {22,5}, {23,5}, {24,5}},
                                      {{21,4}, {22,4}, {23,4}, {24,4}                                },
                                      {                                                              } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2elems_3procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("1x1x2");
  rebalance_mesh(2);
  test_decomposed_mesh_element_distribution({1, 1});
  test_decomposed_mesh_node_sharing({ {                            {5,1}, {6,1}, {7,1}, {8,1}},
                                      {{5,0}, {6,0}, {7,0}, {8,0}                            } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 3elems_3procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("1x1x3");
  rebalance_mesh(2);
  test_decomposed_mesh_element_distribution({2, 1});
  test_decomposed_mesh_node_sharing({ {                                { 9,1}, {10,1}, {11,1}, {12,1}},
                                      {{ 9,0}, {10,0}, {11,0}, {12,0}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 4elems_3procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("1x1x4");
  rebalance_mesh(2);
  test_decomposed_mesh_element_distribution({2, 2});
  test_decomposed_mesh_node_sharing({ {                                { 9,1}, {10,1}, {11,1}, {12,1}},
                                      {{ 9,0}, {10,0}, {11,0}, {12,0}                                } });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2x2elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x2x2");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({1, 1, 1, 1});
  test_decomposed_mesh_node_sharing({ {{ 3,2}, { 4,2}, { 7,1}, { 8,1}, { 9,1}, { 9,2}, { 9,3}, {10,1}, {10,2}, {10,3}},
                                      {{ 7,0}, { 8,0}, { 9,0}, { 9,2}, { 9,3}, {10,0}, {10,2}, {10,3}, {15,3}, {16,3}},
                                      {{ 3,0}, { 4,0}, { 9,0}, { 9,1}, { 9,3}, {10,0}, {10,1}, {10,3}, {11,3}, {12,3}},
                                      {{ 9,0}, { 9,1}, { 9,2}, {10,0}, {10,1}, {10,2}, {11,2}, {12,2}, {15,1}, {16,1}} });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2x2elems_2procTo3proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x2x2");
  rebalance_mesh(3);
  test_decomposed_mesh_element_distribution({1, 1, 2});
  test_decomposed_mesh_node_sharing({ {{ 7,1}, { 8,1}, { 9,1}, { 9,2}, {10,1}, {10,2}, {15,2}, {16,2}},
                                      {{ 3,2}, { 4,2}, { 7,0}, { 8,0}, { 9,0}, { 9,2}, {10,0}, {10,2}},
                                      {{ 3,1}, { 4,1}, { 9,0}, { 9,1}, {10,0}, {10,1}, {15,0}, {16,0}} });
  clean_up_temporary_files();
}

TEST_F(Rebalance, 2x2x2elems_2procTo8proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("2x2x2");
  rebalance_mesh(8);
  test_decomposed_mesh_element_distribution({1, 1, 1, 1, 1, 1, 1, 1});
  test_decomposed_mesh_node_sharing({ {{ 2,4}, { 4,2}, { 5,2}, { 5,4}, { 5,6}, {10,1}, {11,1}, {11,4}, {11,5}, {13,1}, {13,2}, {13,3},
                                               {14,1}, {14,2}, {14,3}, {14,4}, {14,5}, {14,6}, {14,7}},
                                      {{10,0}, {11,0}, {11,4}, {11,5}, {13,0}, {13,2}, {13,3}, {20,5}, {22,3}, {23,3}, {23,5}, {23,7},
                                               {14,0}, {14,2}, {14,3}, {14,4}, {14,5}, {14,6}, {14,7}},
                                      {{ 4,0}, { 5,0}, { 5,4}, { 5,6}, { 8,6}, {13,0}, {13,1}, {13,3}, {16,3}, {17,3}, {17,6}, {17,7},
                                               {14,0}, {14,1}, {14,3}, {14,4}, {14,5}, {14,6}, {14,7}},
                                      {{13,0}, {13,1}, {13,2}, {16,2}, {17,2}, {17,6}, {17,7}, {22,1}, {23,1}, {23,5}, {23,7}, {26,7},
                                               {14,0}, {14,1}, {14,2}, {14,4}, {14,5}, {14,6}, {14,7}},
                                      {{ 2,0}, { 5,0}, { 5,2}, { 5,6}, { 6,6}, {11,0}, {11,1}, {11,5}, {12,5}, {15,5}, {15,6}, {15,7},
                                               {14,0}, {14,1}, {14,2}, {14,3}, {14,5}, {14,6}, {14,7}},
                                      {{11,0}, {11,1}, {11,4}, {12,4}, {15,4}, {15,6}, {15,7}, {20,1}, {23,1}, {23,3}, {23,7}, {24,7},
                                               {14,0}, {14,1}, {14,2}, {14,3}, {14,4}, {14,6}, {14,7}},
                                      {{ 5,0}, { 5,2}, { 5,4}, { 6,4}, { 8,2}, {15,4}, {15,5}, {15,7}, {17,2}, {17,3}, {17,7}, {18,7},
                                               {14,0}, {14,1}, {14,2}, {14,3}, {14,4}, {14,5}, {14,7}},
                                      {{15,4}, {15,5}, {15,6}, {17,2}, {17,3}, {17,6}, {18,6}, {23,1}, {23,3}, {23,5}, {24,5}, {26,3},
                                               {14,0}, {14,1}, {14,2}, {14,3}, {14,4}, {14,5}, {14,6}} });
  clean_up_temporary_files();
}

}

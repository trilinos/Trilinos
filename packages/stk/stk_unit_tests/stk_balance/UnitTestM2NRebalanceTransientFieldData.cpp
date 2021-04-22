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
#include "MeshFixtureM2NRebalance.hpp"
#include <stk_util/environment/EnvData.hpp>
#include <vector>

namespace {

class M2NRebalanceTransientFieldData : public MeshFixtureM2NRebalance
{
public:
  virtual void rebalance_mesh_m2n(int numFinalProcs) override
  {
    m_numFinalProcs = numFinalProcs;
    const bool useNestedDecomp = false;
    stk::balance::M2NParsedOptions parsedOptions{get_output_file_name(), numFinalProcs, useNestedDecomp};
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    stk::balance::internal::rebalanceMtoN(m_ioBroker, *m_targetDecompField, parsedOptions);
    stk::EnvData::instance().m_outputP0 = &std::cout;
  }
};

TEST_F(M2NRebalanceTransientFieldData, 2elems_1procTo1proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh_with_transient_field_data("1x1x2");
  rebalance_mesh_m2n(1);
  test_decomposed_mesh_element_distribution({2});
}

TEST_F(M2NRebalanceTransientFieldData, 4elems_1procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh_with_transient_field_data("1x1x4");
  rebalance_mesh_m2n(2);
  test_decomposed_mesh_element_distribution({2, 2});
}

TEST_F(M2NRebalanceTransientFieldData, 6elems_1procTo3proc)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh_with_transient_field_data("1x1x6");
  rebalance_mesh_m2n(3);
  test_decomposed_mesh_element_distribution({2, 2, 2});
}

TEST_F(M2NRebalanceTransientFieldData, 2elems_2procTo1proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x2");
  rebalance_mesh_m2n(1);
  test_decomposed_mesh_element_distribution({2});
}

TEST_F(M2NRebalanceTransientFieldData, 2elems_3procTo1proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh_with_transient_field_data("1x1x2");
  rebalance_mesh_m2n(1);
  test_decomposed_mesh_element_distribution({2});
}

TEST_F(M2NRebalanceTransientFieldData, 4elems_2procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x4");
  rebalance_mesh_m2n(2);
  test_decomposed_mesh_element_distribution({2, 2});
}

TEST_F(M2NRebalanceTransientFieldData, 8elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x8");
  rebalance_mesh_m2n(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 2});
}

TEST_F(M2NRebalanceTransientFieldData, 7elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x7");
  rebalance_mesh_m2n(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 1});
}

TEST_F(M2NRebalanceTransientFieldData, 9elems_2procTo4proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x9");
  rebalance_mesh_m2n(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 3});
}

TEST_F(M2NRebalanceTransientFieldData, 8elems_2procTo8proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x8");
  rebalance_mesh_m2n(8);
  test_decomposed_mesh_element_distribution({1, 1, 1, 1, 1, 1, 1, 1});
}

TEST_F(M2NRebalanceTransientFieldData, 12elems_3procTo6proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh_with_transient_field_data("1x1x12");
  rebalance_mesh_m2n(6);
  test_decomposed_mesh_element_distribution({2, 2, 2, 2, 2, 2});
}

TEST_F(M2NRebalanceTransientFieldData, 6elems_2procTo3proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x6");
  rebalance_mesh_m2n(3);
  test_decomposed_mesh_element_distribution({2, 2, 2});
}

TEST_F(M2NRebalanceTransientFieldData, 7elems_2procTo3proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x7");
  rebalance_mesh_m2n(3);
  test_decomposed_mesh_element_distribution({2, 2, 3});
}

TEST_F(M2NRebalanceTransientFieldData, 6elems_2procTo7proc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x6");
  rebalance_mesh_m2n(7);
  test_decomposed_mesh_element_distribution({1, 1, 1, 0, 1, 1, 1});
}

TEST_F(M2NRebalanceTransientFieldData, 4elems_3procTo2proc)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh_with_transient_field_data("1x1x4");
  rebalance_mesh_m2n(2);
  test_decomposed_mesh_element_distribution({2, 2});
}

}

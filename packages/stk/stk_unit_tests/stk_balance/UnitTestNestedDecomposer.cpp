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
#include "MeshFixtureDecomposer.hpp"
#include "stk_balance/rebalance.hpp"
#include <vector>

namespace {

class NestedDecomposer : public MeshFixtureDecomposer
{
public:
  virtual void setup_decomposer(int numFinalProcs) override
  {
    m_balanceSettings.set_input_filename("junk.g");
    m_balanceSettings.set_num_output_processors(numFinalProcs);

    m_decomposer = new stk::balance::NestedDecomposer(get_bulk(), m_balanceSettings);
  }
};

TEST_F(NestedDecomposer, SubdomainMapping_FinalProcCountNotEvenlyDivisible)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_decomposer(3), std::logic_error);

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainMapping_FinalProcCountTooSmall)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_decomposer(1), std::logic_error);

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainMapping_EqualInitialAndFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_decomposer(2);
  test_subdomain_mapping({{0}, {1}});

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainMapping_LargerFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_decomposer(6);
  test_subdomain_mapping({{0, 1, 2}, {3, 4, 5}});

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainDecomposition_SingleProcInitialAndFinal)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("generated:1x1x6");
  setup_decomposer(1);
  EXPECT_TRUE(is_nested_decomp());

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainDecomposition_SingleProcInitialToMultipleFinal)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("generated:1x1x6");
  setup_decomposer(3);
  EXPECT_TRUE(is_nested_decomp());

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainDecomposition_FinalProcNotEvenlyDivisible)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_decomposer(3), std::logic_error);

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainDecomposition_FinalProcCountTooSmall)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_decomposer(1), std::logic_error);

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainDecomposition_EqualInitialAndFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_decomposer(2);
  EXPECT_TRUE(is_nested_decomp());

  clean_up_decomposer();
}

TEST_F(NestedDecomposer, SubdomainDecomposition_LargerFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_decomposer(4);
  EXPECT_TRUE(is_nested_decomp());

  clean_up_decomposer();
}


class NestedRebalance : public MeshFixtureRebalance
{
public:
  void rebalance_mesh(int numFinalProcs, const std::string & decompMethod = "rcb")
  {
    const bool useNestedDecomp = true;
    m_balanceSettings.set_is_rebalancing(true);
    m_balanceSettings.set_output_filename(get_output_file_name());
    m_balanceSettings.set_num_input_processors(stk::parallel_machine_size(get_comm()));
    m_balanceSettings.set_num_output_processors(numFinalProcs);
    m_balanceSettings.set_use_nested_decomp(useNestedDecomp);
    m_balanceSettings.setDecompMethod(decompMethod);

    stk::set_outputP0(&stk::outputNull());
    stk::balance::rebalance(m_ioBroker, m_balanceSettings);
    stk::reset_default_output_streams();
  }
};

TEST_F(NestedRebalance, Bar)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x8");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 2});

  clean_up_temporary_files();
}

TEST_F(NestedRebalance, Plate)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x4x4");
  rebalance_mesh(4);
  test_decomposed_mesh_element_distribution({4, 4, 4, 4});

  clean_up_temporary_files();
}

TEST_F(NestedRebalance, Cube)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("4x4x4");
  rebalance_mesh(8);
  test_decomposed_mesh_element_distribution({8, 8, 8, 8, 8, 8, 8, 8});

  clean_up_temporary_files();
}

}

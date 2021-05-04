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
#include "MeshFixtureM2NDecomposer.hpp"
#include <vector>

class M2NDefaultDecomposer : public MeshFixtureM2NDecomposer
{
public:
  virtual void setup_m2n_decomposer(int numFinalProcs) override
  {
    const bool useNestedDecomp = false;
    m_parsedOptions = stk::balance::M2NParsedOptions{"junk.g", numFinalProcs, useNestedDecomp};

    m_decomposer = new stk::balance::internal::M2NDecomposer(get_bulk(), m_balanceSettings, m_parsedOptions);
  }
};

TEST_F(M2NDefaultDecomposer, SubdomainMapping_FinalProcCountNotEvenlyDivisible)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(3);
  test_subdomain_mapping({{0, 2}, {1}});
}

TEST_F(M2NDefaultDecomposer, SubdomainMapping_SmallerFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(1);
  test_subdomain_mapping({{0}, {}});
}

TEST_F(M2NDefaultDecomposer, SubdomainMapping_EqualInitialAndFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(2);
  test_subdomain_mapping({{0}, {1}});
}

TEST_F(M2NDefaultDecomposer, SubdomainMapping_LargerFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(6);
  test_subdomain_mapping({{0, 2, 4}, {1, 3, 5}});
}


TEST_F(M2NDefaultDecomposer, SubdomainDecomposition_LargerFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(4);
  EXPECT_FALSE(is_nested_decomp());
}


class M2NNestedDecomposer : public MeshFixtureM2NDecomposer
{
public:
  virtual void setup_m2n_decomposer(int numFinalProcs) override
  {
    const bool useNestedDecomp = false;
    m_parsedOptions = stk::balance::M2NParsedOptions{"junk.g", numFinalProcs, useNestedDecomp};

    m_decomposer = new stk::balance::internal::M2NDecomposerNested(get_bulk(), m_balanceSettings, m_parsedOptions);
  }
};

TEST_F(M2NNestedDecomposer, SubdomainMapping_FinalProcCountNotEvenlyDivisible)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_m2n_decomposer(3), std::logic_error);
}

TEST_F(M2NNestedDecomposer, SubdomainMapping_FinalProcCountTooSmall)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_m2n_decomposer(1), std::logic_error);
}

TEST_F(M2NNestedDecomposer, SubdomainMapping_EqualInitialAndFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(2);
  test_subdomain_mapping({{0}, {1}});
}

TEST_F(M2NNestedDecomposer, SubdomainMapping_LargerFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(6);
  test_subdomain_mapping({{0, 1, 2}, {3, 4, 5}});
}


TEST_F(M2NNestedDecomposer, SubdomainDecomposition_SingleProcInitialAndFinal)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(1);
  EXPECT_TRUE(is_nested_decomp());
}

TEST_F(M2NNestedDecomposer, SubdomainDecomposition_SingleProcInitialToMultipleFinal)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(3);
  EXPECT_TRUE(is_nested_decomp());
}

TEST_F(M2NNestedDecomposer, SubdomainDecomposition_FinalProcNotEvenlyDivisible)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_m2n_decomposer(3), std::logic_error);
}

TEST_F(M2NNestedDecomposer, SubdomainDecomposition_FinalProcCountTooSmall)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  EXPECT_THROW(setup_m2n_decomposer(1), std::logic_error);
}

TEST_F(M2NNestedDecomposer, SubdomainDecomposition_EqualInitialAndFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(2);
  EXPECT_TRUE(is_nested_decomp());
}

TEST_F(M2NNestedDecomposer, SubdomainDecomposition_LargerFinalProcCount)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  setup_m2n_decomposer(4);
  EXPECT_TRUE(is_nested_decomp());
}


class M2NDefaultRebalance : public MeshFixtureM2NRebalance
{
public:
  virtual void rebalance_mesh_m2n(int numFinalProcs) override
  {
    m_numFinalProcs = numFinalProcs;
    const bool useNestedDecomp = false;
    stk::balance::M2NParsedOptions parsedOptions{get_output_file_name(), numFinalProcs, useNestedDecomp};
    testing::internal::CaptureStdout();
    stk::balance::internal::rebalanceMtoN(m_ioBroker, *m_targetDecompField, parsedOptions);
    testing::internal::GetCapturedStdout();
  }
};

TEST_F(M2NDefaultRebalance, Bar)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x8");
  rebalance_mesh_m2n(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 2});
}

TEST_F(M2NDefaultRebalance, Plate)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x4x4");
  rebalance_mesh_m2n(4);
  test_decomposed_mesh_element_distribution({4, 4, 4, 4});
}

TEST_F(M2NDefaultRebalance, Cube)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("4x4x4");
  rebalance_mesh_m2n(8);
  test_decomposed_mesh_element_distribution({8, 8, 8, 8, 8, 8, 8, 8});
}


class M2NNestedRebalance : public MeshFixtureM2NRebalance
{
public:
  virtual void rebalance_mesh_m2n(int numFinalProcs) override
  {
    m_numFinalProcs = numFinalProcs;
    const bool useNestedDecomp = true;
    stk::balance::M2NParsedOptions parsedOptions{get_output_file_name(), numFinalProcs, useNestedDecomp};
    testing::internal::CaptureStdout();
    stk::balance::internal::rebalanceMtoN(m_ioBroker, *m_targetDecompField, parsedOptions);
    testing::internal::GetCapturedStdout();
  }
};

TEST_F(M2NNestedRebalance, Bar)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x1x8");
  rebalance_mesh_m2n(4);
  test_decomposed_mesh_element_distribution({2, 2, 2, 2});
}

TEST_F(M2NNestedRebalance, Plate)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("1x4x4");
  rebalance_mesh_m2n(4);
  test_decomposed_mesh_element_distribution({4, 4, 4, 4});
}

TEST_F(M2NNestedRebalance, Cube)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("4x4x4");
  rebalance_mesh_m2n(8);
  test_decomposed_mesh_element_distribution({8, 8, 8, 8, 8, 8, 8, 8});
}

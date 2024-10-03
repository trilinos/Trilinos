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
#ifndef MESHFIXTUREDECOMPOSER_HPP
#define MESHFIXTUREDECOMPOSER_HPP

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Decomposer.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <vector>
#include <unistd.h>

class MeshFixtureDecomposer : public stk::unit_test_util::MeshFixture
{
protected:
  MeshFixtureDecomposer()
    : MeshFixture(),
      m_decomposer(nullptr)
  { }

  ~MeshFixtureDecomposer() override {
    delete m_decomposer;
  }

  void setup_initial_mesh(const std::string & inputMeshFile)
  {
    setup_mesh(inputMeshFile, stk::mesh::BulkData::AUTO_AURA);
  }

  virtual void setup_decomposer(int numFinalProcs) = 0;

  void test_subdomain_mapping(const std::vector<std::vector<unsigned>> & expectedInitialToFinalSubdomainMapping)
  {
    std::vector<unsigned> finalToInitialMapping = m_decomposer->map_new_subdomains_to_original_processors();

    unsigned numInitialSubdomains = static_cast<unsigned>(get_bulk().parallel_size());
    ASSERT_EQ(expectedInitialToFinalSubdomainMapping.size(), numInitialSubdomains) << "Wrong number of initial subdomains";

    for (unsigned initialSubdomain = 0; initialSubdomain < numInitialSubdomains; ++initialSubdomain) {
      std::vector<unsigned> mappedSubdomains;
      for (unsigned finalSubdomain = 0; finalSubdomain < finalToInitialMapping.size(); ++finalSubdomain) {
        if (finalToInitialMapping[finalSubdomain] == initialSubdomain) {
          mappedSubdomains.push_back(finalSubdomain);
        }
      }

      EXPECT_EQ(expectedInitialToFinalSubdomainMapping[initialSubdomain], mappedSubdomains);
    }
  }

  bool is_nested_decomp()
  {
    std::vector<unsigned> finalToInitialMapping = m_decomposer->map_new_subdomains_to_original_processors();

    std::set<unsigned> myNestedSubdomains;
    const unsigned myRank = static_cast<unsigned>(get_bulk().parallel_rank());
    for (unsigned newSubdomain = 0; newSubdomain < finalToInitialMapping.size(); ++newSubdomain) {
      if (finalToInitialMapping[newSubdomain] == myRank) {
        myNestedSubdomains.insert(newSubdomain);
      }
    }


    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    stk::balance::DecompositionChangeList decomp = m_decomposer->get_partition();
    stk::EnvData::instance().m_outputP0 = &std::cout;

    bool isNested = true;
    for (const auto & entityProc : decomp) {
      const unsigned targetProc = static_cast<unsigned>(entityProc.second);
      if (myNestedSubdomains.find(targetProc) == myNestedSubdomains.end()) {
        isNested = false;
        break;
      }
    }

    return stk::is_true_on_all_procs(get_comm(), isNested);
  }

  stk::balance::Decomposer * m_decomposer;
  stk::balance::BalanceSettings m_balanceSettings;
};


#endif

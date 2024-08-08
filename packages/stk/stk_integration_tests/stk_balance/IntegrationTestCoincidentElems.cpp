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

#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/balance.hpp>
#include <gtest/gtest.h>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>
#include <test_utils/StkBalanceRunner.hpp>

namespace {

class CoincidentElems : public stk::unit_test_util::MeshFixture
{
protected:
  CoincidentElems()
    : balanceRunner(get_comm()),
      outputDir("output")
  {
    balanceRunner.set_output_dir(outputDir);
  }

  void expect_coincidents_on_same_proc(stk::mesh::EntityId elem1Id, stk::mesh::EntityId elem2Id)
  {
    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, elem1Id);
    if(get_bulk().is_valid(elem1))
    {
      stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, elem2Id);
      EXPECT_TRUE(get_bulk().is_valid(elem2));
    }
  }

  void test_decomp_and_balance(const std::string& fileName,
                               stk::mesh::EntityId elem1Id,
                               stk::mesh::EntityId elem2Id,
                               const size_t numGlobalSides,
                               const std::vector<SideTestUtil::Side>& expectedElemSides)
  {
    balanceRunner.set_filename(fileName);

    EXPECT_NO_THROW(balanceRunner.run_end_to_end());

    setup_mesh(outputDir + "/" + fileName, stk::mesh::BulkData::NO_AUTO_AURA);

    expect_coincidents_on_same_proc(elem1Id, elem2Id);
    SideTestUtil::expect_global_num_sides_in_part(get_bulk(), numGlobalSides, get_meta().universal_part());
    SideTestUtil::expect_all_sides_exist_for_elem_side(get_bulk(), fileName, expectedElemSides);
  }

  stk::integration_test_utils::StkBalanceRunner balanceRunner;
  stk::balance::GraphCreationSettings graphOptions;

private:
  const std::string outputDir;
};

TEST_F(CoincidentElems, balance_coincidentsNotSplit)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  std::string meshDesc =
      "0,1,HEX_8,1,2,3,4,5,6,7,8\n"
      "1,2,HEX_8,5,6,7,8,9,10,11,12\n"
      "1,3,HEX_8,5,6,7,8,9,10,11,12";
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

  expect_coincidents_on_same_proc(2, 3);
  stk::balance::balanceStkMesh(graphOptions, get_bulk());
  expect_coincidents_on_same_proc(2, 3);
}

TEST_F(CoincidentElems, linearDecompJL_noThrow)
{
  test_decomp_and_balance("JL.e", 1, 2, 1, std::vector<SideTestUtil::Side>{{1,5}, {2,5}});
}

TEST_F(CoincidentElems, linearDecompALJ_noThrow)
{
  test_decomp_and_balance("ALJ.e", 2, 3, 1, std::vector<SideTestUtil::Side>{{1,5}, {2,4}, {3,4}});
}

TEST_F(CoincidentElems, linearDecompARefLA_noThrow)
{
  test_decomp_and_balance("ARefLA.e", 3, 4, 2, std::vector<SideTestUtil::Side>{{1,5}, {2,4}, {3,0}, {3,1}, {4,0}, {4,1}});
}

}

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

#include "stk_unit_test_utils/stk_mesh_fixtures/SelectorFixture.hpp"  // for SelectorFixture
#include <algorithm>                    // for sort
#include <gtest/gtest.h>                // for AssertHelper, ASSERT_FALSE, etc
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <vector>                       // for vector
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }

namespace {

using stk::mesh::Bucket;
using stk::mesh::Part;
using stk::mesh::Selector;
using stk::mesh::BulkData;
using stk::mesh::Entity;

// Just make copies so we don't have to worry about constness
template <class T>
void sort_and_compare_eq(std::vector<T*> results,
                         std::vector<T*> expected_results)
{
  std::sort(expected_results.begin(), expected_results.end());
  std::sort(results.begin(), results.end());
  ASSERT_EQ(results.size(), expected_results.size());
  for (unsigned i = 0; i < results.size(); ++i) {
    ASSERT_TRUE(results[i] == expected_results[i]);
  }
}

TEST( UnitTestGetBuckets, ExampleFixture )
{
  // Using the SelectorFixture, test for correct bucket membership and
  // correct results from get_buckets.

  // Generate mesh

  stk::mesh::fixtures::SelectorFixture fix ;
  fix.m_meta_data.commit();

  fix.m_bulk_data.modification_begin();
  fix.generate_mesh();
  ASSERT_TRUE(fix.m_bulk_data.modification_end());
  const stk::mesh::BulkData& mesh = fix.m_bulk_data;

  // Check bucket membership correctness

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity1);
    ASSERT_TRUE(  bucket.member( fix.m_partA ) );
    ASSERT_FALSE( bucket.member( fix.m_partB ) );
    ASSERT_FALSE( bucket.member( fix.m_partC ) );
    ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity2);
    ASSERT_TRUE(  bucket.member( fix.m_partA ) );
    ASSERT_TRUE(  bucket.member( fix.m_partB ) );
    ASSERT_FALSE( bucket.member( fix.m_partC ) );
    ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity3);
    ASSERT_FALSE( bucket.member( fix.m_partA ) );
    ASSERT_TRUE(  bucket.member( fix.m_partB ) );
    ASSERT_TRUE(  bucket.member( fix.m_partC ) );
    ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity4);
    ASSERT_FALSE( bucket.member( fix.m_partA ) );
    ASSERT_FALSE( bucket.member( fix.m_partB ) );
    ASSERT_TRUE(  bucket.member( fix.m_partC ) );
    ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

  {
    const Bucket & bucket = mesh.bucket(fix.m_entity5);
    ASSERT_FALSE( bucket.member( fix.m_partA ) );
    ASSERT_FALSE( bucket.member( fix.m_partB ) );
    ASSERT_FALSE( bucket.member( fix.m_partC ) );
    ASSERT_FALSE( bucket.member( fix.m_partD ) );
  }

}

} // empty namespace

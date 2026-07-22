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

#include <stddef.h>                     // for NULL
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/baseImpl/PartRepository.hpp>  // for PartRepository
#include <gtest/gtest.h>
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Part;
using stk::mesh::impl::PartRepository;

class UnitTestPartRepository
{
public:
  UnitTestPartRepository();
  ~UnitTestPartRepository() {}

  const int spatial_dimension;
  MetaData meta;
  stk::mesh::impl::PartRepository partRepo;
  stk::mesh::impl::PartRepository partRepo_1;
  stk::mesh::impl::PartRepository partRepo_2;

  stk::mesh::Part * universal_part;
  stk::mesh::Part * part_A;
  stk::mesh::Part * part_B;
  stk::mesh::Part * part_C;
  stk::mesh::Part * part_D;
  stk::mesh::Part * part_1_A;
  stk::mesh::Part * part_1_B;
  stk::mesh::Part * part_2_A;
  const stk::topology * singleton;
};

UnitTestPartRepository::UnitTestPartRepository()
  : spatial_dimension(3),
    meta(spatial_dimension),
    partRepo(&meta),
    partRepo_1(&meta),
    partRepo_2(&meta),
    universal_part(partRepo.universal_part()),
    part_A   (partRepo.declare_part("A",stk::topology::NODE_RANK) ),
    part_B   (partRepo.declare_part("B",stk::topology::NODE_RANK) ),
    part_C   (partRepo.declare_part("C",stk::topology::NODE_RANK) ),
    part_D   (partRepo.declare_part("D",stk::topology::EDGE_RANK) ),
    part_1_A (partRepo_1.declare_part("A",stk::topology::NODE_RANK) ),
    part_1_B (partRepo_1.declare_part("B",stk::topology::NODE_RANK) ),
    part_2_A (partRepo_2.declare_part("A",stk::topology::NODE_RANK) ),
    singleton(nullptr)
{
  meta.commit();
}

namespace {

TEST( UnitTestPartRepository, construct )
{
  UnitTestPartRepository upr;
  EXPECT_TRUE( upr.universal_part != NULL );
}

TEST( UnitTestPartRepository, universal_in_subset )
{
  UnitTestPartRepository upr;
  ASSERT_THROW(
        upr.partRepo.declare_subset(*upr.part_A,*upr.universal_part),
        std::runtime_error
        );
}

TEST( UnitTestPartRepository, subset_equal_superset )
{
  UnitTestPartRepository upr;
  ASSERT_THROW(
        upr.partRepo.declare_subset(*upr.part_A,*upr.part_A),
        std::runtime_error
        );
}


TEST( UnitTestPartRepository, circular_subset )
{
  UnitTestPartRepository upr;
  upr.partRepo.declare_subset(*upr.part_A,*upr.part_B);
  upr.partRepo.declare_subset(*upr.part_B,*upr.part_C);
  ASSERT_THROW(
        upr.partRepo.declare_subset(*upr.part_C,*upr.part_A),
        std::runtime_error
        );
}

TEST( UnitTestPartRepository, inconsistent_rank_subset )
{
  UnitTestPartRepository upr;
  // lower rank cannot contain higher rank:
  ASSERT_THROW(
        upr.partRepo.declare_subset(*upr.part_A,*upr.part_D),
        std::runtime_error
        );
  // higher rank can contain lower rank:
  ASSERT_NO_THROW(
        upr.partRepo.declare_subset(*upr.part_D,*upr.part_A)
        );
}

TEST( UnitTestPartRepository, two_part_repositories )
{
  UnitTestPartRepository upr;
  // subset/superset parts must come from same part repository
  ASSERT_THROW(
        upr.partRepo_1.declare_subset(*upr.part_1_A,*upr.part_2_A),
        std::runtime_error
        );
}

//Test covers declare_attribute_no_delete in PartRepository.hpp and PartImpl.hpp
TEST( UnitTestPartRepository, declare_attribute_no_delete )
{
  UnitTestPartRepository upr;
  upr.partRepo.declare_attribute_no_delete(*upr.part_A, upr.singleton);

}

TEST( UnitTestPartRepository, get_all_parts )
{
  UnitTestPartRepository upr;
  const stk::mesh::PartVector & pv = upr.partRepo.get_all_parts();
  ASSERT_EQ( pv.size(), 5u );
}

TEST( UnitTestPartRepository, get_mesh_parts )
{
  UnitTestPartRepository upr;
  stk::mesh::PartVector mpv = upr.partRepo.get_mesh_parts();
  ASSERT_EQ( mpv.size(), 4u );
}


}

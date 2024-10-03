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

#include <string>                       // for to_string
#include <Ioss_Utils.h>                 // for Utils
#include <stddef.h>                     // for size_t
#include <iosfwd>                       // for ostream
#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <string>                       // for allocator, operator+, etc
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper
#include "stk_mesh/base/Types.hpp"      // for PartVector

enum { SpaceDim = 3 };

TEST(UnitTestGmeshFixture, testUnit)
{
  const size_t num_x = 1;
  const size_t num_y = 2;
  const size_t num_z = 3;
  const size_t num_surf = 6;
  std::string config_mesh = std::to_string(num_x) + "x" +
                            std::to_string(num_y) + "x" +
                            std::to_string(num_z) + "|sideset:xXyYzZ";
  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, config_mesh);

  fixture.commit();

  const std::vector<std::string> & sideset_names = fixture.getSidesetNames();
  ASSERT_EQ( num_surf, sideset_names.size() );

  for( size_t i = 0; i < num_surf; ++i ) {
    std::string surf_name =  (std::string)"surface_" + std::to_string(i+1);
    ASSERT_TRUE(surf_name == sideset_names[i]);
  }

  // Needed to test field data
  stk::mesh::Field<double> * coord_field = fixture.getMetaData().get_field<double>(stk::topology::NODE_RANK, "coordinates");
  ASSERT_TRUE( coord_field );

  const stk::mesh::PartVector & side_parts = fixture.getSideParts();
  ASSERT_EQ( sideset_names.size(), side_parts.size() );
}

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_MortonIndex.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <array>

namespace krino
{

TEST(MortonIndex, specified_index)
{
  const std::vector<int> valid_indices {0, 1, 2, 4, 3434, 23, 1414415, 76285, MAX_VALID_MORTON3D_INDEX};
  for (int valid_index : valid_indices)
  {
    EXPECT_EQ(valid_index, morton3d_decode_index(morton3d_encode_index(valid_index)));
  }
#ifndef NDEBUG
  const std::vector<int> invalid_indices {MAX_VALID_MORTON3D_INDEX+1, 3000000, 2147483647};
  for (int invalid_index : invalid_indices)
  {
    EXPECT_ANY_THROW(morton3d_encode_index(invalid_index));
  }
#endif
}

TEST(MortonIndex, specified_indices)
{
  const std::vector<int> valid_indices {0, 1, 2, 4, 3434, 23, 1414415, 76285, MAX_VALID_MORTON3D_INDEX};
  for (int ix : valid_indices)
  {
    for (int iy : valid_indices)
    {
      for (int iz : valid_indices)
      {
        std::array<int,3> out_indices = morton3d_decode_indices(morton3d_encode_indices({{ix,iy,iz}}));
        EXPECT_EQ(ix, out_indices[0]);
        EXPECT_EQ(iy, out_indices[1]);
        EXPECT_EQ(iz, out_indices[2]);
      }
    }
  }
}

} // namespace krino

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MortonIndex_h
#define Akri_MortonIndex_h

#include <stk_util/util/ReportHandler.hpp>
#include <array>

namespace krino {

constexpr int MAX_VALID_MORTON3D_INDEX = 0x1fffff;

uint64_t morton3d_encode_index(const int i)
{
  STK_ThrowAssertMsg(i >= 0, "Index must be >= 0 to encode: " << i);
  STK_ThrowAssertMsg(i <= MAX_VALID_MORTON3D_INDEX, "Index must be <= " << MAX_VALID_MORTON3D_INDEX << " to encode: " << i);
  uint64_t m = i & 0x1fffff;
  m = (m | m << 32) & 0x001f00000000ffff;
  m = (m | m << 16) & 0x001f0000ff0000ff;
  m = (m | m << 8) & 0x100f00f00f00f00f;
  m = (m | m << 4) & 0x10c30c30c30c30c3;
  m = (m | m << 2) & 0x1249249249249249;
  return m;
}

int morton3d_decode_index(uint64_t m)
{
  m &= 0x1249249249249249;
  m = (m ^ (m >> 2)) & 0x10c30c30c30c30c3;
  m = (m ^ (m >> 4)) & 0x100f00f00f00f00f;
  m = (m ^ (m >> 8)) & 0x001f0000ff0000ff;
  m = (m ^ (m >> 16)) & 0x001f00000000ffff;
  m = (m ^ (m >> 32)) & 0x1fffff;
  return static_cast<int>(m);
}

uint64_t morton3d_encode_indices(const std::array<int,3> & indices)
{
  return morton3d_encode_index(indices[2]) * 4 + morton3d_encode_index(indices[1]) * 2 + morton3d_encode_index(indices[0]);
}

std::array<int,3> morton3d_decode_indices(const uint64_t m)
{
  const int ix = morton3d_decode_index(m);
  const int iy = morton3d_decode_index(m >> 1);
  const int iz = morton3d_decode_index(m >> 2);
  return {{ix,iy,iz}};
}

} // namespace krino

#endif // Akri_MortonIndex_h

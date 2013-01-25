/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#include <KokkosArray_Host.hpp>

#include <KokkosArray_View.hpp>

#include <iostream>
#include <iomanip>


namespace {

template< typename TiledArray>
struct TestTileArray
{
  static void run(int dimension)
  {
    TiledArray a("",dimension,dimension);

    for (int tj = 0, tile_cols = a.tiles_in_dimension_1(); tj < tile_cols; ++tj) {
    for (int ti = 0, tile_rows = a.tiles_in_dimension_0(); ti < tile_rows; ++ti) {
      typename TiledArray::tile_type tile = a.tile(ti,tj);
      for (int j=0, cols = tile.dimension_1(); j < cols; ++j) {
      for (int i=0, rows = tile.dimension_0(); i < rows; ++i) {
        tile(i,j) = rows*cols*(ti + a.tiles_in_dimension_0()*tj) + i + j*rows;
      }}
    }}

    for (int j=0; j<dimension; ++j) {
    for (int i=0; i<dimension; ++i) {
      ptrdiff_t offset = &a(i,j) - &a(0,0);
      EXPECT_EQ( offset, a(i,j) );
    }}
  }

};


} // namespace



TEST( host_tile, tile_1x1)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<1,1>, Host> >::run(9);
}

TEST( host_tile, tile_2x2)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<2,2>, Host> >::run(9);
}

TEST( host_tile, tile_3x3)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<3,3>, Host> >::run(9);
}

TEST( host_tile, tile_4x4)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<4,4>, Host> >::run(9);
}

TEST( host_tile, tile_5x5)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<5,5>, Host> >::run(9);
}

TEST( host_tile, tile_6x6)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<6,6>, Host> >::run(9);
}

TEST( host_tile, tile_7x7)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<7,7>, Host> >::run(9);
}

TEST( host_tile, tile_8x8)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<8,8>, Host> >::run(9);
}

TEST( host_tile, tile_9x9)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<9,9>, Host> >::run(9);
}

TEST( host_tile, tile_10x10)
{
  using namespace KokkosArray;
  TestTileArray< View< ptrdiff_t**, LayoutTileLeft<10,10>, Host> >::run(9);
}


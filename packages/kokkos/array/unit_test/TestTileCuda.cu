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
#include "TestTile.hpp"

#include <KokkosArray_Cuda.hpp>

TEST( host_tile, tile_1x1)
{
  static const size_t dim = 9;
  typedef KokkosArray::LayoutTileLeft<1,1> tile_layout;
  typedef ReduceTileErrors< KokkosArray::Cuda, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = KokkosArray::parallel_reduce(dim, functor_type(array) );
  EXPECT_EQ( errors, 0u);
}

TEST( host_tile, tile_2x2)
{
  static const size_t dim = 9;
  typedef KokkosArray::LayoutTileLeft<2,2> tile_layout;
  typedef ReduceTileErrors< KokkosArray::Cuda, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = KokkosArray::parallel_reduce(dim, functor_type(array) );
  EXPECT_EQ( errors, 0u);
}

TEST( host_tile, tile_4x4)
{
  static const size_t dim = 9;
  typedef KokkosArray::LayoutTileLeft<4,4> tile_layout;
  typedef ReduceTileErrors< KokkosArray::Cuda, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = KokkosArray::parallel_reduce(dim, functor_type(array) );
  EXPECT_EQ( errors, 0u);
}

TEST( host_tile, tile_8x8)
{
  static const size_t dim = 9;
  typedef KokkosArray::LayoutTileLeft<8,8> tile_layout;
  typedef ReduceTileErrors< KokkosArray::Cuda, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = KokkosArray::parallel_reduce(dim, functor_type(array) );
  EXPECT_EQ( errors, 0u);
}

TEST( host_tile, tile_16x16)
{
  static const size_t dim = 9;
  typedef KokkosArray::LayoutTileLeft<16,16> tile_layout;
  typedef ReduceTileErrors< KokkosArray::Cuda, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = KokkosArray::parallel_reduce(dim, functor_type(array) );
  EXPECT_EQ( errors, 0u);
}

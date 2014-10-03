/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#include <iostream>

#include <Kokkos_Core.hpp>

#include <Kokkos_Bitset.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <Kokkos_Vector.hpp>


//----------------------------------------------------------------------------
#include <TestBitset.hpp>
#include <TestUnorderedMap.hpp>
#include <TestVector.hpp>
#include <TestDualView.hpp>
#include <TestSegmentedView.hpp>

namespace Test {


void cuda_test_insert_close(  uint32_t num_nodes
                            , uint32_t num_inserts
                            , uint32_t num_duplicates
                           )
{
  test_insert< Kokkos::Cuda >( num_nodes, num_inserts, num_duplicates, true);
}

void cuda_test_insert_far(  uint32_t num_nodes
                          , uint32_t num_inserts
                          , uint32_t num_duplicates
                         )
{
  test_insert< Kokkos::Cuda >( num_nodes, num_inserts, num_duplicates, false);
}

void cuda_test_failed_insert(  uint32_t num_nodes )
{
  test_failed_insert< Kokkos::Cuda >( num_nodes );
}

void cuda_test_deep_copy(  uint32_t num_nodes )
{
  test_deep_copy< Kokkos::Cuda >( num_nodes );
}

void cuda_test_vector_combinations(unsigned int size)
{
  test_vector_combinations<int,Kokkos::Cuda>(size);
}

void cuda_test_dualview_combinations(unsigned int size)
{
  test_dualview_combinations<int,Kokkos::Cuda>(size);
}

void cuda_test_segmented_view(unsigned int size)
{
  test_segmented_view<double,Kokkos::Cuda>(size);
}

void cuda_test_bitset()
{
  test_bitset<Kokkos::Cuda>();
}

}

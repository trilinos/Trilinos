/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

// include for ThrustGPUNode method implementations
#include "Kokkos_ThrustGPUNode.cuh"

// includes for all operators
#include "TestOps.hpp"
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

template void Kokkos::ThrustGPUNode::parallel_for<InitOp<int> >(int, int, InitOp<int>);
template SumOp<float>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<SumOp<float> >(int, int, SumOp<float>);
template void Kokkos::ThrustGPUNode::parallel_for<InitOp<float> >(int, int, InitOp<float>);
template NullOp::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<NullOp>(int, int, NullOp);
template SumOp<int>::ReductionType Kokkos::ThrustGPUNode::parallel_reduce<SumOp<int> >(int, int, SumOp<int>);

struct thrust_test_constant_float
{
  __host__ __device__
  float operator()() {return 1.0f;}
};
struct thrust_test_constant_int
{
  __host__ __device__
  int operator()() {return 1;}
};

void thrust_float_alloc(int N, thrust::device_vector<float> &buff) {
  buff.resize(N);
}
void thrust_int_alloc(int N, thrust::device_vector<int> &buff) {
  buff.resize(N);
}

void thrust_float_init(thrust::device_vector<float> &buff) {
  thrust::generate( buff.begin(), buff.end(), thrust_test_constant_float() );
}
void thrust_int_init(thrust::device_vector<int> &buff) {
  thrust::generate( buff.begin(), buff.end(), thrust_test_constant_int() );
}

float thrust_float_sum(const thrust::device_vector<float> &buff) {
  return thrust::reduce( buff.begin(), buff.end(), 0.0f, thrust::plus<float>() );
}
int thrust_int_sum(const thrust::device_vector<int> &buff) {
  return thrust::reduce( buff.begin(), buff.end(), 0,    thrust::plus<int>() );
}

/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_NODAPI_SER_NATIVETESTS_HPP_
#define KOKKOS_NODAPI_SER_NATIVETESTS_HPP_

#include <Teuchos_ScalarTraits.hpp>

namespace {

  using Kokkos::SerialNode;

  template <>
  std::pair<double,double> nativeTimings<float,SerialNode>(int N, int numIters, float &result, const RCP<SerialNode> &node) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,SerialNode init"), sTime("float,SerialNode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      float *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        for (int i=0; i < N; ++i) {
          ptr[i] = 1;
        }
      }
    }
    float sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      const float *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        sum = ptr[0];
        for (int i=1; i < N; ++i) {
          sum += ptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,SerialNode>(int N, int numIters, int &result, const RCP<SerialNode> &node) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,SerialNode init"), sTime("int,SerialNode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      int *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        for (int i=0; i < N; ++i) {
          ptr[i] = 1;
        }
      }
    }
    int sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      const int *ptr = buff.getRawPtr();
      for (int t=0; t < numIters; ++t) {
        sum = ptr[0];
        for (int i=1; i < N; ++i) {
          sum += ptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

}

#endif // KOKKOS_NODAPI_SER_NATIVETESTS_HPP_

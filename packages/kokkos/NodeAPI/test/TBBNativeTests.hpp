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

#ifndef KOKKOS_NODAPI_TBB_NATIVETESTS_HPP_
#define KOKKOS_NODAPI_TBB_NATIVETESTS_HPP_

#include <Teuchos_ScalarTraits.hpp>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

namespace {

  using Kokkos::TBBNode;

  template <class T>
  class TBBNodeTestInit {
  private:
    T * const a_;
  public:
    TBBNodeTestInit(T* a) : a_(a) {}
    void operator()(const tbb::blocked_range<int> &r) const {
      T * const my_a = a_;
      for (int i=r.begin(); i != r.end(); ++i) {
        my_a[i] = Teuchos::ScalarTraits<T>::one();
      }
    }
  };

  template <class T>
  class TBBNodeTestSum {
  private:
    const T * const a_;
  public:
    T sum;
    TBBNodeTestSum(const T* a) : a_(a), sum(Teuchos::ScalarTraits<T>::zero()) {}
    TBBNodeTestSum(TBBNodeTestSum<T> &other, tbb::split) : a_(other.a_), sum(Teuchos::ScalarTraits<T>::zero()) {}
    void join(const TBBNodeTestSum<T> &other) { sum += other.sum; }
    void operator()(const tbb::blocked_range<int> &r) {
      const T* const my_a = a_;
      for (int i=r.begin(); i != r.end(); ++i) {
        sum += my_a[i];
      }
    }
  };

  template <>
  std::pair<double,double> nativeTimings<float,TBBNode>(int N, int numIters, float &result, const RCP<TBBNode> &node) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,TBBNode init"), sTime("float,TBBNode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        tbb::parallel_for( tbb::blocked_range<int>(0,N), TBBNodeTestInit<float>(buff.getRawPtr()), tbb::auto_partitioner() );
      }
    }
    float sum = 0.0f;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TBBNodeTestSum<float> op(buff.getRawPtr());
        tbb::parallel_reduce( tbb::blocked_range<int>(0,N), op, tbb::auto_partitioner() );
        sum = op.sum;
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,TBBNode>(int N, int numIters, int &result, const RCP<TBBNode> &node) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,TBBNode init"), sTime("int,TBBNode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        tbb::parallel_for( tbb::blocked_range<int>(0,N), TBBNodeTestInit<int>(buff.getRawPtr()), tbb::auto_partitioner() );
      }
    }
    int sum = 0;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TBBNodeTestSum<int> op(buff.getRawPtr());
        tbb::parallel_reduce( tbb::blocked_range<int>(0,N), op, tbb::auto_partitioner() );
        sum = op.sum;
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

}

#endif // KOKKOS_NODAPI_TBB_NATIVETESTS_HPP_

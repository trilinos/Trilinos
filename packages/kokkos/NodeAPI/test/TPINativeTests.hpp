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

#ifndef KOKKOS_NODAPI_TPI_NATIVETESTS_HPP_
#define KOKKOS_NODAPI_TPI_NATIVETESTS_HPP_

#include <Teuchos_ScalarTraits.hpp>

namespace {

  using Kokkos::TPINode;

  template <class T>
  struct TPIInit {
    T * x;
    unsigned int N;
    static void work( TPI_Work * work ) {
      struct TPIInit<T> * const w = (TPIInit<T> *) work->info ;
      unsigned int begin;
      const unsigned int max = ( w->N + work->count - 1 ) / work->count ;
      const unsigned int end = w->N - max * ( work->count - ( work->rank + 1 ) );
      if ( work->rank ) {
        begin  = end - max ;
      }
      else {
        begin  = 0 ;
      }
      T * const my_x = w->x;
      for (unsigned int i=begin; i != end; ++i) {my_x[i] = Teuchos::ScalarTraits<T>::one();}
    }
  };

  template <class T>
  struct TPISum {
    const T * x;
    unsigned int N;
    // reduction
    static void work( TPI_Work * work ) {
      struct TPIInit<T> * const w = (TPIInit<T> *) work->info ;
      T * const dst = (T *) work->reduce;
      unsigned int begin;
      const unsigned int max = ( w->N + work->count - 1 ) / work->count ;
      const unsigned int end = w->N - max * ( work->count - ( work->rank + 1 ) );
      if ( work->rank ) {
        begin  = end - max ;
      }
      else {
        begin  = 0 ;
      }
      T * const my_x = w->x;
      for (unsigned int i=begin; i != end; ++i) {*dst += my_x[i];}
    }
    // initialization
    static void init( TPI_Work * work ) {
      T * const dst = (T *) work->reduce ;
      (*dst) = Teuchos::ScalarTraits<T>::zero();
    }
    // combination
    static void join( TPI_Work * work , const void * arg_src ) {
      T * const dst = (T *) work->reduce ;
      const T * const src = (const T *) arg_src ;
      (*dst) += (*src);
    }
  };

  template <>
  std::pair<double,double> nativeTimings<float,TPINode>(int N, int numIters, float &result, const RCP<TPINode> &node) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,TPINode init"), sTime("float,TPINode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        TPIInit<float> data;
        data.x = buff.getRawPtr();
        data.N = N;
        TPI_Run_threads( &TPIInit<float>::work, &data, 0 );
      }
    }
    float sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TPISum<float> data;
        data.x = buff.getRawPtr();
        data.N = N;
        sum = 0.0f;
        TPI_Run_threads_reduce( &TPISum<float>::work, &data, &TPISum<float>::join, &TPISum<float>::init, sizeof(sum), &sum );
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,TPINode>(int N, int numIters, int &result, const RCP<TPINode> &node) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,TPINode init"), sTime("int,TPINode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        TPIInit<int> data;
        data.x = buff.getRawPtr();
        data.N = N;
        TPI_Run_threads( &TPIInit<int>::work, &data, 0 );
      }
    }
    int sum;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        TPISum<int> data;
        data.x = buff.getRawPtr();
        data.N = N;
        sum = 0;
        TPI_Run_threads_reduce( &TPISum<int>::work, &data, &TPISum<int>::join, &TPISum<int>::init, sizeof(sum), &sum );
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

}

#endif // KOKKOS_NODAPI_TPI_NATIVETESTS_HPP_

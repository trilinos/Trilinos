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

#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <vector>
#include <iostream>

#include <stk_util/environment/CPUTime.hpp>

namespace {

struct Timer {
  Timer() : startTime(stk::cpu_time()) {}

  double elapsed() const { return stk::cpu_time() - startTime; }

private:
  double startTime;
};

}

void force_calculation( long long sum )
{
  static int count = 0;
  //volatile to prevent the compiler from optimizing the loops away
  static volatile long long value;

  EXPECT_TRUE( count==0 || value == sum );

  value = sum;
  ++count;
}

struct sum_func
{
  sum_func() : sum(0) {}
  long long sum;
  void operator()( unsigned i )
  {
    sum += i;
  }
};

TEST( loop_iteration, loop_iteration)
{
  Timer  total_time;

  const size_t vector_size = 100000000;
  const size_t num_repetitions = 20;

  std::vector<unsigned> data;
  data.reserve(vector_size);

  {
    Timer timer;
    std::srand(10);
    std::generate_n(std::back_inserter(data), vector_size, &std::rand);
    double elapsedtime = timer.elapsed();
    std::cout << "Construction time: " << elapsedtime << " seconds." << std::endl;
  }
  {
    Timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::accumulate( data.begin(), data.end(), 0ll );
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Iterator accumulate took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    unsigned* begin = &data.front();
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::accumulate( begin, begin+data.size(), 0ll );
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Pointer accumulate took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      for( std::vector<unsigned>::iterator it = data.begin(), iend = data.end(); it != iend; ++it ) {
        sum += *it;
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Iterator for loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      unsigned* it = &data.front();
      for( unsigned* iend = it+data.size(); it != iend; ++it ) {
        sum += *it;
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Pointer for loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      for( std::size_t i = 0, size = data.size(); i != size; ++i ) {
        sum += data[i];
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Index for (size_t) loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      for( int i = 0, size = data.size(); i != size; ++i ) {
        sum += data[i];
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Index for (int) loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::for_each( data.begin(), data.end(), sum_func() ).sum;
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Index for_each took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    unsigned* begin = &data.front();
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::for_each( begin, begin+data.size(), sum_func() ).sum;
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Pointer for_each took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    Timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      for( unsigned i : data ) {
        sum += i;
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "C++11 range-for took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }

  double elapsedtime = total_time.elapsed();

  std::cout << "Total time: " << elapsedtime << " seconds." << std::endl;
}

// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stk_util/stk_config.h>

// restrict this file to only build if KokkosCore is enabled
#ifdef HAVE_STK_KokkosCore

#include <Kokkos_Serial.hpp>
#include <Kokkos_Threads.hpp>
#include <Kokkos_OpenMP.hpp>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_View.hpp>

#include <iostream>

namespace {

#if defined(KOKKOS_HAVE_PTHREAD)
#define KOKKOS_THREAD_DEVICE KOKKOS_THREAD_DEVICE
#elif defined(KOKKOS_HAVE_OPENMP)
#define KOKKOS_THREAD_DEVICE Kokkos::OpenMP
#else
#define KOKKOS_THREAD_DEVICE Kokkos::Serial
#endif

// setup and tear down for the KokkosThreads unit tests
class KokkosThreads : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    unsigned num_threads = 8;

    // if hwloc is present we will get better thread placement and
    // numa aware allocation.
    // Currently sierra does not have the hwloc TPL enabled
    if (Kokkos::hwloc::available()) {
      std::cout<<"Getting num_threads from Kokkos::hwloc..."<<std::endl;
      num_threads = Kokkos::hwloc::get_available_numa_count()
                    * Kokkos::hwloc::get_available_cores_per_numa()
                 // * Kokkos::hwloc::get_available_threads_per_core()
                    ;

    }
    else {
      std::cout<<"Kokkos::hwloc not available."<<std::endl;
    }

    KOKKOS_THREAD_DEVICE::initialize( num_threads );

    std::cout << "Kokkos thread device 'print_configuration' output:\n";
    KOKKOS_THREAD_DEVICE::print_configuration(std::cout, true);

    std::cout << "\nNumber of Threads: " << num_threads << std::endl;
  }

  static void TearDownTestCase()
  {
    KOKKOS_THREAD_DEVICE::finalize();
  }
};

const size_t RUN_TIME_DIMENSION = 4000000;
const size_t COMPILE_TIME_DIMENSION = 3;

} // unnamed namespace


TEST_F( KokkosThreads, SerialInitialize)
{
  // allocate a rank 2 array witn that is RUN_TIME_DIMENSION x COMPILE_TIME_DIMENSION

  // View will default initialize all the values unless it is explicitly disabled, ie,
  // Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> a("node views", RUN_TIME_DIMENSION);
  // zero fills the array, but
  // Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> a( Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);
  // will allocate without initializing the array

  Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> a( Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);

  for (size_t i=0; i < a.dimension_0(); ++i) {
    for (size_t x=0; x < a.dimension_1(); ++x) {
      a(i,x) = i;
    }
  }

  // get a const view to the same array
  // this view shares the same memory as a, but cannot modify the values
  Kokkos::View<const unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> b = a;

  for (size_t i=0; i < b.dimension_0(); ++i) {
    for (size_t x=0; x < b.dimension_1(); ++x) {
      EXPECT_EQ(i, b(i,x));
    }
  }
}

// Not available until c++11 support is enable in Sierra and Kokkos
#if defined (KOKKOS_HAVE_C_PLUS_PLUS_11_LAMBDA)
TEST_F( KokkosThreads, LambdaInitialize)
{
  Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> a( Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);

  Kokkos::parallel_for<KOKKOS_THREAD_DEVICE>(
    a.dimension_0() ,
    [=](size_t i) {
      for (size_t x=0; x < a.dimension_1(); ++x) {
        a(i,x) = i;
      }
    }
  );

  Kokkos::View<const unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> b = a;

  int num_error = 0;
  // Cannot portably call a GTEST macro in parallel
  // count the errors and test that they are equal to zero
  Kokkos::parallel_reduce<KOKKOS_THREAD_DEVICE, int /*reduction value type */>(
    b.dimension_0() ,
    [](int & local_errors)                                    // init lambda
    { local_errors = 0; } ,
    [=](size_t i, int & local_errors) {                       // operator() lambda
      for (size_t x=0; x < b.dimension_1(); ++x)
        local_errors += i == b(i,x) ? 0 : 1;
    } ,
    [](volatile int & dst_err, volatile int const& src_err)   // join lambda
    { dst_err += src_err; } ,
    num_errors                                                // where to store the result
  );
  EXPECT_EQ( 0, num_errors);

}
#endif


namespace {

// Functors need to initialize and check the view with out c++11 support

template <typename View>
struct InitializeView
{
  // need a device_type typedef for all parallel functors
  typedef typename View::device_type device_type;

  View a;

  // get a view to the a
  template <typename RhsView>
  InitializeView( RhsView const& arg_a )
    : a(arg_a)
  {}

  void apply()
  {
    // call parallel_for on this functor
    Kokkos::parallel_for( a.dimension_0(), *this);
  }

  // initialize the a
  KOKKOS_INLINE_FUNCTION
  void operator()(size_t i) const
  {
    for (size_t x=0; x < a.dimension_1(); ++x) {
      a(i,x) = i;
    }
  }
};

template <typename View>
struct CheckView
{
  // need a device_type typedef for all parallel functors
  typedef typename View::device_type device_type;

  // need a value_type typedef for the reduction
  typedef int value_type;

  View a;

  // get a view to the a
  template <typename RhsView>
  CheckView( RhsView const& arg_a )
    : a(arg_a)
  {}

  // return the number of errors found
  value_type apply()
  {
    int num_errors = 0;
    // call a parallel_reduce to count the errors
    Kokkos::parallel_reduce( a.dimension_0(), *this, num_errors);
    return num_errors;
  }

  // initialize the reduction type
  KOKKOS_INLINE_FUNCTION
  void init(value_type & v) const
  { v = 0; }

  // this threads contribution to the reduction type
  // check that the value is equal to the expected value
  // otherwise increment the error count
  KOKKOS_INLINE_FUNCTION
  void operator()(size_t i, value_type & error) const
  {
    for (size_t x=0; x < a.dimension_1(); ++x) {
      error += i == a(i,x) ? 0 : 1;
    }
  }

  // join two threads together
  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, volatile value_type const& src) const
  { dst += src; }
};

struct ForceFunctor
{
    typedef typename KOKKOS_THREAD_DEVICE device_type;

    ForceFunctor(const double* disp, const double* vel, const double* acc, double *force) :
            mDisp(disp), mVel(vel), mAcc(acc), mForce(force), alpha(-1.4), beta(0.3333333), gamma(3.14159)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t rowIndex) const
    {
        mForce[rowIndex] = alpha * mDisp[rowIndex] + beta * mVel[rowIndex] + gamma * mAcc[rowIndex];
    }

private:
    const double *mDisp;
    const double *mVel;
    const double *mAcc;
    double *mForce;
    const double alpha;
    const double beta;
    const double gamma;
};

} // unnameed namespace

TEST_F( KokkosThreads, ParallelInitialize)
{
  typedef Kokkos::View<unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> view_type;
  typedef Kokkos::View<const unsigned*[COMPILE_TIME_DIMENSION], KOKKOS_THREAD_DEVICE> const_view_type;

  view_type a(Kokkos::allocate_without_initializing, "node views", RUN_TIME_DIMENSION);

  // call the InitializeView functor
  {
    InitializeView<view_type> f(a);
    f.apply();
  }

  // call the CheckView functor
  // and expect no errors
  {
    const_view_type b = a;
    CheckView<view_type> f(a);
    const int num_errors = f.apply();
    EXPECT_EQ( 0, num_errors);
  }
}

#endif


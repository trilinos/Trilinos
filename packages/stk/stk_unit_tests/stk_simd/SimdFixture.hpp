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

#ifndef UNITTEST_SIMDFIXTURE_HPP
#define UNITTEST_SIMDFIXTURE_HPP

#include <gtest/gtest.h>
#include <stk_simd/Simd.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <algorithm>
#include <functional>

namespace stk
{
namespace unit_test_simd
{

typedef std::vector<double, non_std::AlignedAllocator<double,64> > DoubleVector;
typedef std::vector<stk::simd::Double> SimdDoubleVector;

typedef std::vector<bool> BoolVector;
typedef std::vector<stk::simd::Bool> SimdBoolVector;

typedef stk::simd::Double* DoublesPtr;

struct Range {
  double left;
  double right;
};

double rand_in(Range range) {
  return range.left + (range.right-range.left) * (double(rand()) / (double)RAND_MAX);
}

// #define TIME_SIMD_OPERATIONS

int size() { 
#ifdef TIME_SIMD_OPERATIONS
  return 8192;
#else
  return 32;
#endif
}
int simd_size() { return size() / stk::simd::ndoubles; }


class TestSimdMathFunction : public ::testing::Test
{
 protected:
  TestSimdMathFunction() { }

  template <typename F1, typename F2>
  void test_operator(F1 scalarOperation, F2 simdOperation)
  {
    DoubleVector scalarSolution = apply_generic_operator_on_scalar_array( scalarOperation );
    SimdDoubleVector simdSolution = apply_generic_operator_on_simd_array( simdOperation );
#ifndef TIME_SIMD_OPERATIONS
    verify_arrays_equal(scalarSolution, simdSolution);
#endif
  }

  template <typename ArgFunc>
  DoubleVector apply_generic_operator_on_scalar_array(const ArgFunc& argFunc) const {
    DoubleVector solution(size());
    const int sz = size();
#ifdef TIME_SIMD_OPERATIONS
    Kokkos::Impl::Timer scalarTimer;
#endif
#if defined(__INTEL_COMPILER)
#pragma novector
#endif
    for (int i=0; i < sz; ++i)
      solution[i] = argFunc(i);
#ifdef TIME_SIMD_OPERATIONS
    printf("Sclr loop took: %g seconds.\n", scalarTimer.seconds());
#endif
    return solution;
  }

  template <typename ArgFunc>
  SimdDoubleVector apply_generic_operator_on_simd_array(const ArgFunc& argFunc) const {
    SimdDoubleVector solution(simd_size());
    const int sz = simd_size();
#ifdef TIME_SIMD_OPERATIONS
    Kokkos::Impl::Timer simdTimer;
#endif
    for (int i=0; i < sz; ++i)
      solution[i] = argFunc(i);
#ifdef TIME_SIMD_OPERATIONS
    printf("Simd loop took: %g seconds.\n", simdTimer.seconds());
#endif
    return solution;
  }

  virtual std::string get_input_vals(int i) const = 0;

  void verify_arrays_equal(const DoubleVector& scalarVec, const SimdDoubleVector& simdVec) const {
    const double* simdSol = reinterpret_cast<const double*>(simdVec.data());
    for (int i=0; i < size(); ++i) {
      EXPECT_DOUBLE_EQ(scalarVec[i], simdSol[i]) << "input values = " + get_input_vals(i);
    }
  }
};


class SimdScalarDoubleHybrid
{
 public:

  SimdScalarDoubleHybrid()
  {
    initialize_array();
  }

  inline double scalar(int i) const { return array[i]; }
  inline stk::simd::Double simd(int i) const { return arraySimd[i]; }

  void set_values_in(Range range) {
    for (int i=0; i<size(); ++i)
      array[i] = rand_in(range);
  }

 private:

  void initialize_array() {
    array.resize(size());
    set_values_in(Range{0.0, 0.5});
    arraySimd = stk::simd::simd_ptr_cast(array.data());
  }

  DoubleVector array;
  stk::simd::Double* arraySimd;
};

class SimdScalarBoolHybrid
{
 public:

  SimdScalarBoolHybrid()
  {
    initialize_array();
  }

  inline bool scalar(int i) const { return array[i]; }
  inline stk::simd::Bool simd(int i) const { return arraySimd[i]; }

  void set_alternating_bools() {
    for (int i=0; i < size(); ++i)
      array[i] = i%2==0 ? true : false;

    DoubleVector onesAndZeros(size());
    for (int i=0; i < size(); ++i)
      onesAndZeros[i] = array[i] ? 1.0 : 0.0;

    stk::simd::Double* onesAndZerosSimd = stk::simd::simd_ptr_cast(onesAndZeros.data());

    for (int i=0; i < simd_size(); ++i) {
      arraySimd[i] = onesAndZerosSimd[i] > 0.5;
    }
  }

 private:

  void initialize_array() {
    array.resize(size());
    arraySimd.resize(simd_size());
    set_alternating_bools();
  }

  BoolVector array;
  std::vector<stk::simd::Bool> arraySimd;
};


class MathFunctionWithOneDoubleArg : public TestSimdMathFunction
{
 protected:
  SimdScalarDoubleHybrid x;

  template <typename ScalarFunc, typename SimdFunc>
  void test_simd_operator(ScalarFunc scalarFunc, SimdFunc simdFunc,
                          Range range) {
    x.set_values_in(range);
    test_operator([&](int i) { return scalarFunc(x.scalar(i)); },
                  [&](int i) { return simdFunc(x.simd(i)); });
  }

  virtual std::string get_input_vals(int i) const override {
    return std::to_string(x.scalar(i));
  }
};

class MathFunctionWithTwoDoubleArg : public TestSimdMathFunction
{
 protected:
  SimdScalarDoubleHybrid x;
  SimdScalarDoubleHybrid y;

  template <typename ScalarFunc, typename SimdFunc>
  void test_simd_operator(ScalarFunc scalarFunc, SimdFunc simdFunc, 
                          Range rangeX, Range rangeY) {
    x.set_values_in(rangeX);
    y.set_values_in(rangeY);
    test_operator([&](int i) { return scalarFunc(x.scalar(i), y.scalar(i)); },
                  [&](int i) { return simdFunc(x.simd(i), y.simd(i)); });
  }

  virtual std::string get_input_vals(int i) const override {
    return std::to_string(x.scalar(i))+" "+std::to_string(y.scalar(i));
  }
};

class MathFunctionWithThreeDoubleArg : public TestSimdMathFunction
{
 protected:
  SimdScalarDoubleHybrid x;
  SimdScalarDoubleHybrid y;
  SimdScalarDoubleHybrid z;

  template <typename ScalarFunc, typename SimdFunc>
  void test_simd_operator(ScalarFunc scalarFunc, SimdFunc simdFunc, 
                          Range rangeX, Range rangeY, Range rangeZ) {
    x.set_values_in(rangeX);
    y.set_values_in(rangeY);
    z.set_values_in(rangeZ);
    test_operator([&](int i) { return scalarFunc(x.scalar(i), y.scalar(i), z.scalar(i)); },
                  [&](int i) { return simdFunc(x.simd(i), y.simd(i), z.simd(i)); });
  }

  virtual std::string get_input_vals(int i) const override {
    return std::to_string(x.scalar(i))+" "+std::to_string(y.scalar(i))+" "+std::to_string(z.scalar(i));
  }
};

class MathFunctionWithBoolAndDoubleArg : public TestSimdMathFunction
{
 protected:
  SimdScalarBoolHybrid condition;
  SimdScalarDoubleHybrid x;

  template <typename ScalarFunc, typename SimdFunc>
  void test_simd_operator(ScalarFunc scalarFunc, SimdFunc simdFunc, 
                          Range rangeX) {
    x.set_values_in(rangeX);
    test_operator([&](int i) { return scalarFunc(condition.scalar(i), x.scalar(i)); },
                  [&](int i) { return simdFunc(condition.simd(i), x.simd(i)); });
  }

  virtual std::string get_input_vals(int i) const override {
    return std::to_string(x.scalar(i));
  }
};

class MathFunctionWithBoolAndTwoDoubleArg : public TestSimdMathFunction
{
 protected:
  SimdScalarBoolHybrid condition;
  SimdScalarDoubleHybrid x;
  SimdScalarDoubleHybrid y;

  template <typename ScalarFunc, typename SimdFunc>
  void test_simd_operator(ScalarFunc scalarFunc, SimdFunc simdFunc, 
                          Range rangeX, Range rangeY) {
    x.set_values_in(rangeX);
    y.set_values_in(rangeY);
    test_operator([&](int i) { return scalarFunc(condition.scalar(i), x.scalar(i), y.scalar(i)); },
                  [&](int i) { return simdFunc(condition.simd(i), x.simd(i), y.simd(i)); });
  }

  virtual std::string get_input_vals(int i) const override {
    return std::to_string(x.scalar(i))+" "+std::to_string(y.scalar(i));
  }
};

}
}

#endif

#ifndef UNITTEST_SIMDFIXTURE_HPP
#define UNITTEST_SIMDFIXTURE_HPP

#include <gtest/gtest.h>
#include <stk_simd/Simd.hpp>
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

int size() { return 64; }
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
    verify_arrays_equal(scalarSolution, simdSolution);
  }

  template <typename ArgFunc>
  DoubleVector apply_generic_operator_on_scalar_array(ArgFunc argFunc) const {
    DoubleVector solution(size());
    for (int i=0; i < size(); ++i)
      solution[i] = argFunc(i);
    return solution;
  }

  template <typename ArgFunc>
  SimdDoubleVector apply_generic_operator_on_simd_array(ArgFunc argFunc) const {
    SimdDoubleVector solution(simd_size());
    for (int i=0; i < simd_size(); ++i)
      solution[i] = argFunc(i);
    return solution;
  }

  virtual std::string get_input_vals(int i) const = 0;

  void verify_arrays_equal(const DoubleVector& scalarVec, const SimdDoubleVector& simdVec) const {
    const double* simdSol = reinterpret_cast<const double*>(simdVec.data());
    for (int i=0; i < size(); ++i) {
      EXPECT_EQ(scalarVec[i], simdSol[i]) << "input values = " + get_input_vals(i);
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

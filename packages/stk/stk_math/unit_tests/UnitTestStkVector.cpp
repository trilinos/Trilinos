#include <stk_math/StkVector.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <cmath>
#include <string>
#include "UnitTestStkVectorUtils.hpp"

namespace 
{
const size_t ZERO=0;
const size_t ONE=1;
const size_t TWO=2;
const size_t THREE=3;
const size_t FOUR=4;
const size_t TEN=10;

double tol = 1.e-10;

using namespace stk::math::unitTestStkVectorUtils;

TEST(stk_math_stkVector,Construction_from_data_double_into_double)
{
  double values[THREE] = { 2.0, 3.0, 4.0};
  stk::math::Vec<double, THREE> vec(values);
  expect_equal(values, vec, tol);
}

TEST(stk_math_stkVector,Construction_from_data_double_into_float)
{
  double values[THREE] = { 2.0, 3.0, 4.0};
  stk::math::Vec<float, THREE> vec(values);
  expect_equal(values, vec, tol);
}

TEST(stk_math_stkVector, Construction_from_data_with_buffer1)
{
  // Create a vector longer than the data we pass in.
  // The remainder should be padded with 0.0.
  double values[FOUR] = { 2.0, 3.0, 4.0, 5.0};
  stk::math::Vec<double, TEN> vec(values, FOUR);

  size_t i=0;
  for ( ; i< FOUR; i++)
    EXPECT_EQ(values[i], vec[i]);

  for ( ; i< TEN; i++)
    EXPECT_EQ(0.0, vec[i]);
}

#ifndef NDEBUG
TEST(stk_math_stkVector, Construction_from_data_with_buffer2)
{
  // Create a vector shorter than the data we pass in.
  // The remainder should be padded with 0.0.
  double values[TEN] = { 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
  ASSERT_DEATH((stk::math::Vec<double, FOUR>{values, TEN}),  ".*Assertion.*");
}

TEST(stk_math_stkVector, Construction_from_data_with_null_buffer1)
{
  // Create a vector shorter than the data we pass in.
  // The remainder should be padded with 0.0.
  ASSERT_DEATH((stk::math::Vec<double, FOUR>{(double*)nullptr, FOUR}),  ".*Assertion.*");
}
#endif

TEST(stk_math_stkVector, Construction_from_data_with_null_buffer2)
{
  stk::math::Vec<double, TEN> vec((double*)nullptr, ZERO);
  expect_equal(std::vector<double>(10, 0.0), vec, tol);
}

TEST(stk_math_stkVector,Copy_Constructor)
{
  double values[THREE] = { 2.0, 3.0, 4.0};
  stk::math::Vec<double, THREE> vec(values);
  stk::math::Vec<double, THREE> copy(vec);
  expect_equal(values, vec, tol);
}

TEST(stk_math_stkVector,Const_Copy_Constructor)
{
  double values[THREE] = { 2.0, 3.0, 4.0};
  const stk::math::Vec<double, THREE> const_vec(values);
  const stk::math::Vec<double, THREE> const_copy(const_vec);
  expect_equal(values, const_copy, tol);
}

TEST(stk_math_stkVector, Construction_doublevec_from_doublexyz)
{
  double xyz[THREE] = {2.0, 3.0, 4.0};
  stk::math::Vec<double, THREE> vec(xyz[0], xyz[1], xyz[2]);
  expect_equal(xyz, vec, tol);
}

TEST(stk_math_stkVector, Construction_floatvec_from_doublexyz)
{
  double xyz[THREE] = {2.0, 3.0, 4.0};
  stk::math::Vec<float, THREE> vec(xyz[0], xyz[1], xyz[2]);
  expect_equal(xyz, vec, tol);
}

TEST(stk_math_stkVector, Construction_doublevec_from_doublexy)
{
  double xyz[TWO] = {2.0, 3.0};
  stk::math::Vec<double, TWO> vec(xyz[0], xyz[1]);
  expect_equal(xyz, vec, tol);
}

TEST(stk_math_stkVector, Construction_floatvec_from_doublexy)
{
  double xyz[TWO] = {2.0, 3.0};
  stk::math::Vec<float, TWO> vec(xyz[0], xyz[1]);
  expect_equal(xyz, vec, tol);
}

TEST(stk_math_stkVector, Construction_doublevec_from_doublex)
{
  double xyz[ONE] = {2.0};
  stk::math::Vec<double, ONE> vec(xyz[0]);
  expect_equal(xyz, vec, tol);
}

TEST(stk_math_stkVector, Construction_floatvec_from_doublex)
{
  double xyz[ONE] = {2.0};
  stk::math::Vec<float, ONE> vec(xyz[0]);
  expect_equal(xyz, vec, tol);
}

TEST(stk_math_stkVector, Construction_default_constructor)
{
  stk::math::Vec<float, TEN> vec;
  expect_equal(std::vector<double>(10, 0.0), vec, tol);
}

TEST(stk_math_stkVector, Construction_MemberInit_constructor)
{
  EXPECT_NO_THROW((stk::math::Vec<float, TEN>{stk::math::MemberInit::NONE}));
  //There is no way to test the contents, because the constructor above
  //does NOT initialize the contents of the vector.  But we can check its length
  stk::math::Vec<float, TEN> vec(stk::math::MemberInit::NONE);
  EXPECT_EQ( TEN, vec.dimension());
}

TEST(stk_math_stkVector, operator_equals)
{
  double values[THREE] = { 2.0, 3.0, 4.0};
  stk::math::Vec<double, THREE> vec(values);
  stk::math::Vec<double, THREE> copy;

  copy = vec;
  expect_equal(values, vec, tol);
}

TEST(stk_math_stkVector,const_operator_equals)
{
  double values[THREE] = { 2.0, 3.0, 4.0};
  const stk::math::Vec<double, THREE> vec(values);
  stk::math::Vec<double, THREE> copy;

  copy = vec;
  expect_equal(values, copy, tol);
}

TEST(stk_math_stkVector, const_operator_minus_equals)
{
  const std::array<double, THREE> gold_vec1 = {{1., 0., -2}};
  const std::array<double, THREE> gold_vec2 = {{1., 3., 6.}};
  stk::math::Vec<double, THREE> vec1(2., 3., 4.);
  const stk::math::Vec<double, THREE> vec2(gold_vec2[0], gold_vec2[1], gold_vec2[2]);

  vec1 -= vec2; // This is what we want to test.

  expect_equal(gold_vec1, vec1, tol);
  expect_equal(gold_vec2, vec2, tol);
}

TEST(stk_math_stkVector, const_operator_plus_equals)
{
  const std::array<double, THREE> gold_vec1 = {{1.1, 5.2, -3.}};
  const std::array<double, THREE> gold_vec2 = {{-1., 2.2, -7.}};
  stk::math::Vec<double, THREE> vec1(2.1, 3., 4.);
  const stk::math::Vec<double, THREE> vec2(gold_vec2[0], gold_vec2[1], gold_vec2[2]);

  vec1 += vec2; // This is what we want to test.

  expect_equal(gold_vec1, vec1, tol);
  expect_equal(gold_vec2, vec2, tol);
}

TEST(stk_math_stkVector, const_operator_multiply_equals)
{
  stk::math::Vec<double, THREE> orig(2.1, 3., -4.);

  const std::array<double, THREE> gold1 = {{0., 0., 0.}};
  stk::math::Vec<double, THREE> vec1(orig);
  vec1*= 0.;

  const std::array<double, THREE> gold2 = {{4.2, 6., -8.}};
  stk::math::Vec<double, THREE> vec2(orig);
  vec2*= 2.;

  const std::array<double, THREE> gold3 = {{-.21, -.3, .4}};
  stk::math::Vec<double, THREE> vec3(orig);
  vec3*= -.1;

  expect_equal(gold1, vec1, tol);
  expect_equal(gold2, vec2, tol);
  expect_equal(gold3, vec3, tol);
}

TEST(stk_math_stkVector, const_operator_divide_equals)
{
  stk::math::Vec<double, THREE> orig(2.1, 3., -4.);

  const std::array<double, THREE> gold1 = {{6.3, 9., -12.}};
  stk::math::Vec<double, THREE> vec1(orig);
  vec1 /= .333333333333333333333;

  const std::array<double, THREE> gold2 = {{-1.05, -1.5, 2.}};
  stk::math::Vec<double, THREE> vec2(orig);
  vec2 /= -2.;

  expect_equal(gold1, vec1, tol);
  expect_equal(gold2, vec2, tol);
}

TEST(stk_math_stkVector, const_operator_unary_negative)
{
  const std::array<double, THREE> gold = {{-2.1, -3., 4.}};
  stk::math::Vec<double, THREE> vec(2.1, 3., -4.);
  stk::math::Vec<double, THREE> answer;

  answer = -vec;

  expect_equal(gold, answer, tol);
}

TEST(stk_math_stkVector, const_operator_subtract)
{
  stk::math::Vec<double, THREE> vec1(2.1, -3.14, -4.);
  stk::math::Vec<double, THREE> vec2(1.1, -3., 10.);
  const std::array<double, THREE> gold1 = {{1., -.14, -14.}};
  const std::array<double, THREE> gold2 = {{-1., .14, 14.}};
  stk::math::Vec<double, THREE> answer1, answer2;

  answer1 = vec1 - vec2;
  answer2 = vec2 - vec1;

  expect_equal(gold1, answer1, tol);
  expect_equal(gold2, answer2, tol);
}

TEST(stk_math_stkVector, const_operator_add)
{
  stk::math::Vec<double, THREE> vec1(2.1, -3.14, -4.);
  stk::math::Vec<double, THREE> vec2(1.1, -3., 10.);
  const std::array<double, THREE> gold = {{3.2, -6.14, 6.}};
  stk::math::Vec<double, THREE> answer1, answer2;

  answer1 = vec1 + vec2;
  answer2 = vec2 + vec1;

  expect_equal(gold, answer1, tol);
  expect_equal(gold, answer2, tol);
}

TEST(stk_math_stkVector, zero_length)
{
    stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    EXPECT_TRUE(zeroVec.zero_length());
    stk::math::Vec<double, THREE> unitVec(1., 0., 0.);
    EXPECT_TRUE(!unitVec.zero_length());
}

TEST(stk_math_stkVector, equality)
{
    stk::math::Vec<double, THREE> xVec1(1., 0., 0.);
    stk::math::Vec<double, THREE> xVec2(1., 0., 0.);
    stk::math::Vec<double, THREE> yVec(0., 1., 0.);
    stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    EXPECT_TRUE(xVec1 == xVec1);
    EXPECT_TRUE(xVec1 == xVec2);
    EXPECT_FALSE(xVec1 == yVec);
    EXPECT_FALSE(zeroVec == yVec);
}

TEST(stk_math_stkVector, nonEquality)
{
    stk::math::Vec<double, THREE> xVec1(1., 0., 0.);
    stk::math::Vec<double, THREE> xVec2(1., 0., 0.);
    stk::math::Vec<double, THREE> yVec(0., 1., 0.);
    stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    EXPECT_FALSE(xVec1 != xVec1);
    EXPECT_FALSE(xVec1 != xVec2);
    EXPECT_TRUE(xVec1 != yVec);
    EXPECT_TRUE(zeroVec != yVec);
}

TEST(stk_math_stkVector, less_componentwise_comparison)
{
    stk::math::Vec<double, THREE> smaller(0., 1., 0.);
    stk::math::Vec<double, THREE> larger(1., 0., 0.);
    EXPECT_TRUE(smaller < larger);
    EXPECT_FALSE(smaller < smaller);
    EXPECT_FALSE(larger < smaller);
}

void expect_vec_is_nan(const stk::math::Vec<double, THREE> &vec)
{
    EXPECT_TRUE(std::isnan(vec[0]));
    EXPECT_TRUE(std::isnan(vec[1]));
    EXPECT_TRUE(std::isnan(vec[2]));
}

TEST(stk_math_stkVector, unitize)
{
    stk::math::Vec<double, THREE> vec(0., 2., 0.);
    vec.unitize();
    expect_equal(std::vector<double>{0., 1., 0.}, vec, tol);

    stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    zeroVec.unitize();
    expect_vec_is_nan(zeroVec);
}

TEST(stk_math_stkVector, unitVector)
{
    stk::math::Vec<double, THREE> vec(0., 2., 0.);
    expect_equal(std::vector<double>{0., 1., 0.}, vec.unit_vector(), tol);
    expect_equal(std::vector<double>{0., 2., 0.}, vec, tol);

    stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    expect_vec_is_nan(zeroVec.unit_vector());
    expect_equal(std::vector<double>{0., 0., 0.}, zeroVec, tol);
}

TEST(stk_math_stkVector, length)
{
    stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    EXPECT_NEAR(0., zeroVec.length(), tol);
    stk::math::Vec<double, THREE> vec(3., 4., 0.);
    EXPECT_NEAR(5., vec.length(), tol);
    stk::math::Vec<double, THREE> vec2(0., 3., 4.);
    EXPECT_NEAR(5., vec2.length(), tol);
}

TEST(stk_math_stkVector, lengthSquared)
{
    stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    EXPECT_NEAR(0., zeroVec.length_squared(), tol);
    stk::math::Vec<double, THREE> vec(1., 2., 3.);
    EXPECT_NEAR(14., vec.length_squared(), tol);
}

TEST(stk_math_stkVector, bracket_operator)
{
    stk::math::Vec<double, THREE> vec(0., 1., 2.);
    for(size_t i=0; i<THREE; i++)
        EXPECT_EQ(double(i), vec[i]);
}

TEST(stk_math_stkVector, const_bracket_operator)
{
    const stk::math::Vec<double, THREE> vec(0., 1., 2.);
    for(size_t i=0; i<THREE; i++)
        EXPECT_EQ(double(i), vec[i]);
}

TEST(stk_math_stkVector, pointer_to_data)
{
    stk::math::Vec<double, THREE> vec(0., 1., 2.);
    double *vecPtr = vec.data();
    double gold[THREE] = {0., 1., 2.};
    for(size_t i=0; i<THREE; i++)
        EXPECT_EQ(gold[i], vecPtr[i]);
}

TEST(stk_math_stkVector, const_pointer_to_data)
{
    const stk::math::Vec<double, THREE> vec(0., 1., 2.);
    const double *vecPtr = vec.data();
    double gold[THREE] = {0., 1., 2.};
    for(size_t i=0; i<THREE; i++)
        EXPECT_EQ(gold[i], vecPtr[i]);
}

TEST(stk_math_stkVector, extraction_operator)
{
    const stk::math::Vec<double, THREE> vec(0., 1., 2.);
    std::ostringstream oss;
    oss << vec;
    std::string output = oss.str();
    EXPECT_TRUE(output.find("0"));
    EXPECT_TRUE(output.find("1"));
    EXPECT_TRUE(output.find("2"));
}

TEST(stk_math_stkVector, dimension)
{
    const stk::math::Vec<double, 1> vec;
    EXPECT_EQ(1u, vec.dimension());
    const stk::math::Vec<double, 2> vec2;
    EXPECT_EQ(2u, vec2.dimension());
}

template <typename Vec>
void test_begin_end(double gold, Vec &vec)
{
    EXPECT_EQ(3u, vec.end() - vec.begin());
    for(auto iter=vec.begin(); iter!=vec.end(); ++iter)
        EXPECT_NEAR(gold, *iter, tol);
}

TEST(stk_math_stkVector, begin_end)
{
    stk::math::Vec<double, THREE> vec(1., 1., 1.);
    test_begin_end(1., vec);
    for(auto iter=vec.begin(); iter!=vec.end(); ++iter)
        *iter = 2.;
    test_begin_end(2., vec);
}

TEST(stk_math_stkVector, const_begin_end)
{
    const stk::math::Vec<double, THREE> vec(1., 1., 1.);
    test_begin_end(1., vec);
}

TEST(stk_math_stkVector, scalar_times_vector)
{
    const stk::math::Vec<double, THREE> vec(0., 1., 2.);
    const stk::math::Vec<double, THREE> result = 2 * vec;
    expect_equal(std::vector<double>{0., 2., 4.}, result, tol);
}

TEST(stk_math_stkVector, vector_times_scalar)
{
    const stk::math::Vec<double, THREE> vec(0., 1., 2.);
    const stk::math::Vec<double, THREE> result = vec * -2;
    expect_equal(std::vector<double>{0., -2., -4.}, result, tol);
}

TEST(stk_math_stkVector, vector_divided_by_scalar)
{
    const stk::math::Vec<double, THREE> vec(0., 1., 2.);
    const stk::math::Vec<double, THREE> result = vec / 2;
    expect_equal(std::vector<double>{0., 0.5, 1.}, result, tol);
}

TEST(stk_math_stkVector, dot_product)
{
    const stk::math::Vec<double, THREE> a(0., 1., 2.);
    const stk::math::Vec<double, THREE> b(2., 2., 2.);
    EXPECT_NEAR(6., Dot(a,b), tol);
}

TEST(stk_math_stkVector, cross_product)
{
    const stk::math::Vec<double, THREE> a(1., 0., 0.);
    const stk::math::Vec<double, THREE> b(0., 1., 0.);
    expect_equal(std::vector<double>{0., 0., 1.}, Cross(a,b), tol);
    const stk::math::Vec<double, THREE> c(1., 1., 1.);
    const stk::math::Vec<double, THREE> d(2., 2., 2.);
    expect_equal(std::vector<double>{0., 0., 0.}, Cross(c,d), tol);
    const stk::math::Vec<double, THREE> zeroVec(0., 0., 0.);
    expect_equal(std::vector<double>{0., 0., 0.}, Cross(zeroVec,d), tol);
}

TEST(stk_math_stkVector, cross_product_with_unit_x)
{
    expect_equal(std::vector<double>{0., 0., -1.}, crossX(stk::math::Vec<double, THREE>(0., 1., 0.)), tol);
    expect_equal(std::vector<double>{0., -1., 0.}, crossX(stk::math::Vec<double, THREE>(0., 0., -1.)), tol);
    expect_equal(std::vector<double>{0., 0., 0.}, crossX(stk::math::Vec<double, THREE>(-1., 0., 0.)), tol);
}

TEST(stk_math_stkVector, cross_product_with_unit_y)
{
    expect_equal(std::vector<double>{0., 0., 1.}, crossY(stk::math::Vec<double, THREE>(1., 0., 0.)), tol);
    expect_equal(std::vector<double>{1., 0., 0.}, crossY(stk::math::Vec<double, THREE>(0., 0., -1.)), tol);
    expect_equal(std::vector<double>{0., 0., 0.}, crossY(stk::math::Vec<double, THREE>(0., 6., 0.)), tol);
}

TEST(stk_math_stkVector, cross_product_with_unit_z)
{
    expect_equal(std::vector<double>{0., -1., 0.}, crossZ(stk::math::Vec<double, THREE>(1., 0., 0.)), tol);
    expect_equal(std::vector<double>{1., 0., 0.}, crossZ(stk::math::Vec<double, THREE>(0., 1., 0.)), tol);
    expect_equal(std::vector<double>{0., 0., 0.}, crossZ(stk::math::Vec<double, THREE>(0., 0., -13.)), tol);
}

TEST(stk_math_stkVector, to_string)
{
    stk::math::Vec<int,THREE> someVec{3,5,7};
    EXPECT_EQ("3 5 7", someVec.to_string());
}

TEST(stk_math_stkVector, whenSettingVectorToInvalid_makeSureItIsInvalid)
{
    stk::math::Vector3d someVec;
    EXPECT_TRUE(someVec.is_valid());
    someVec.set_invalid();
    EXPECT_FALSE(someVec.is_valid());
    someVec = {1,2,3};
    EXPECT_TRUE(someVec.is_valid());
}
}


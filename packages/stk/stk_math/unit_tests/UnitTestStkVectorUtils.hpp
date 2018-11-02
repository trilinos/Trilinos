/*
 * UnitTestStkVectorUtils.hpp
 *
 *  Created on: Oct 18, 2018
 *      Author: mlstate
 */

#ifndef STK_STK_MATH_UNIT_TESTS_UNITTESTSTKVECTORUTILS_HPP_
#define STK_STK_MATH_UNIT_TESTS_UNITTESTSTKVECTORUTILS_HPP_

namespace stk {
namespace math {
namespace unitTestStkVectorUtils{

template <typename GoldT, typename VecT, unsigned Dim>
void expect_equal(GoldT *goldValues, const stk::math::Vec<VecT, Dim> &vec, const double tol)
{
    size_t i = 0;
    for(VecT val : vec)
    {
        EXPECT_LT(i, Dim);
        EXPECT_NEAR(goldValues[i], val, tol);
        i++;
    }
    EXPECT_EQ(Dim, i);
}

template <typename GoldArray, typename VecT, unsigned Dim>
void expect_equal(const GoldArray &goldValues, const stk::math::Vec<VecT, Dim> &vec, const double tol)
{
    expect_equal(goldValues.data(), vec, tol);
}


}
}
}


#endif /* STK_STK_MATH_UNIT_TESTS_UNITTESTSTKVECTORUTILS_HPP_ */

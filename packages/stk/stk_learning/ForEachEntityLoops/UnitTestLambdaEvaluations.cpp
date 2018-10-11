#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <functional>

namespace
{


template <typename FUNC>
void transformOverVector(const std::vector<int> &vecToSum, const FUNC &myFunc)
{
    for(size_t i=0; i<vecToSum.size(); i++)
    {
        myFunc(vecToSum[i]);
    }
}
template <typename FUNC>
void transformOverVectorWithSum(const std::vector<int> &vecToSum, const FUNC &myFunc, int &sum)
{
    for(size_t i=0; i<vecToSum.size(); i++)
    {
        myFunc(vecToSum[i], sum);
    }
}



TEST(LambdaEvaluation, capture_sum)
{
    int sum = 0;
    auto myLambda =
            [&sum](int x)
            {
                sum += x;
            };

    std::vector<int> vec(10, 2);
    transformOverVector(vec, myLambda);

    EXPECT_EQ(20, sum);
}


TEST(LambdaEvaluation, pass_in_sum)
{
    int sum = 0;
    auto myLambda =
            [](int x, int &sum)
            {
                sum += x;
            };

    std::vector<int> vec(10, 2);
    transformOverVectorWithSum(vec, myLambda, sum);

    EXPECT_EQ(20, sum);
}


class MySumFunctor
{
public:
    MySumFunctor(int &sum):
        mSum(sum)
    {
    }
    void operator()(int x) const
    {
        mSum += x;
    }
private:
    int &mSum;
};
TEST(LambdaEvaluation, functor_class)
{
    int sum = 0;
    MySumFunctor mySum(sum);

    std::vector<int> vec(10, 2);
    transformOverVector(vec, mySum);

    EXPECT_EQ(20, sum);
}


void doSum(int x, int &sum)
{
    sum += x;
}
TEST(LambdaEvaluation, function_pass_in_sum)
{
    int sum = 0;
    std::vector<int> vec(10, 2);
    transformOverVectorWithSum(vec, doSum, sum);

    EXPECT_EQ(20, sum);
}


std::function<void(int)> get_my_lambda(int &sum)
{
    return [&sum](int x)
    {
        sum += x;
    };
}
TEST(LambdaEvaluation, try_returning_lambda)
{
    int sum = 0;
    std::vector<int> vec(10, 2);
    transformOverVector(vec, get_my_lambda(sum));

    EXPECT_EQ(20, sum);
}


#define STK_LAMBDA []
#define STK_LAMBDA_USING(...) [__VA_ARGS__]
#define STK_LAMBDA_USING_ALL_LOCALS_BY_VALUE [=]
#define STK_LAMBDA_USING_ALL_LOCALS_BY_REFERENCE [&]
TEST(LambdaEvaluation, pass_in_sum_with_macro)
{
    int sum = 0;
    auto myLambda =
            STK_LAMBDA(int x, int &sum)
            {
                sum += x;
            };

    std::vector<int> vec(10, 2);
    transformOverVectorWithSum(vec, myLambda, sum);

    EXPECT_EQ(20, sum);
}
TEST(LambdaEvaluation, capture_sum_with_macro)
{
    int sum = 0;
    const int extra = 1;

    std::vector<int> vec(10, 2);
    transformOverVector(vec,
        STK_LAMBDA_USING(&sum)(int x)
        {
            sum += x + extra;
        }
    );

    EXPECT_EQ(30, sum);
}
TEST(LambdaEvaluation, capture_sum_with_capture_all_macro)
{
    int sum = 0;
    auto myLambda =
            STK_LAMBDA_USING_ALL_LOCALS_BY_REFERENCE(int x)
            {
                sum += x;
            };

    std::vector<int> vec(10, 2);
    transformOverVector(vec, myLambda);

    EXPECT_EQ(20, sum);
}


}

#include <gtest/gtest.h>
#include <stk_util/util/StkVector.hpp>

namespace {

TEST(StkVectorTest, create_sizeIsZero)
{
    stk::Vector<int> vec;
    EXPECT_EQ(0u, vec.size());
}
TEST(StkVectorTest, createSizeOne_sizeIsOne)
{
    stk::Vector<int> vec(1);
    EXPECT_EQ(1u, vec.size());
}
TEST(StkVectorTest, setSingleInt_getInt)
{
    stk::Vector<int> vec(1);
    vec[0] = 2;
    EXPECT_EQ(2, vec[0]);
}
TEST(StkVectorTest, setSingleDouble_getDouble)
{
    stk::Vector<double> vec(1);
    vec[0] = 0.5;
    EXPECT_EQ(0.5, vec[0]);
}
TEST(StkVectorTest, setTwoInts_getInts)
{
    stk::Vector<int> vec(2);
    vec[0] = 1;
    vec[1] = 2;
    EXPECT_EQ(1, vec[0]);
    EXPECT_EQ(2, vec[1]);
}
TEST(StkVectorTest, pushBackOnEmtpyVector_sizeIsOne)
{
    stk::Vector<double> vec;
    vec.push_back(2);
    EXPECT_EQ(1u, vec.size());
}
TEST(StkVectorTest, pushBackOnEmtpyVector_getValue)
{
    stk::Vector<double> vec;
    vec.push_back(2);
    EXPECT_EQ(2, vec[0]);
}
TEST(StkVectorTest, pushTwoValues_sizeIsTwo)
{
    stk::Vector<double> vec;
    vec.push_back(2);
    vec.push_back(3);
    EXPECT_EQ(2u, vec.size());
}
TEST(StkVectorTest, pushTwoValues_getTwoValues)
{
    stk::Vector<double> vec;
    vec.push_back(2);
    vec.push_back(3);
    EXPECT_EQ(2, vec[0]);
    EXPECT_EQ(3, vec[1]);
}
TEST(StkVectorTest, createWithName_getName)
{
    stk::Vector<double> vec("vec");
    EXPECT_EQ("vec", vec.name());
}
TEST(StkVectorTest, createSizedWithName_getName)
{
    stk::Vector<double> vec("vec", 1);
    EXPECT_EQ("vec", vec.name());
    EXPECT_EQ(1u, vec.size());
    EXPECT_EQ(0.0, vec[0]);
}
TEST(StkVectorTest, createSizedAndInitWithName_getName)
{
    stk::Vector<double> vec("vec", 1, 2.0);
    EXPECT_EQ("vec", vec.name());
    EXPECT_EQ(1u, vec.size());
    EXPECT_EQ(2.0, vec[0]);
}
TEST(StkVectorTest, createNamedPushBack_getName)
{
    stk::Vector<double> vec("vec");
    vec.push_back(3.0);
    EXPECT_EQ("vec", vec.name());
    EXPECT_EQ(1u, vec.size());
    EXPECT_EQ(3.0, vec[0]);
}

class StkVectorSpy : public stk::Vector<double>
{
public:
    StkVectorSpy() : stk::Vector<double>("spy"), numAllocations(0)
    {
    }
    size_t num_allocations() const { return numAllocations; }
protected:
    virtual Kokkos::View<double*>::HostMirror get_new_vals_of_size(size_t s)
    {
        numAllocations++;
        return stk::Vector<double>::get_new_vals_of_size(s);
    }
private:
    size_t numAllocations;
};
void push_back_num_things(stk::Vector<double> &vec, size_t num)
{
    for(size_t i=0; i<num; i++)
        vec.push_back(2);
}
TEST(StkVectorTest, eightPushBacks_fourAllocations)
{
    StkVectorSpy vec;
    push_back_num_things(vec, 8);
    EXPECT_EQ(4u, vec.num_allocations());
}
TEST(StkVectorTest, ninePushBacks_fiveAllocations)
{
    StkVectorSpy vec;
    push_back_num_things(vec, 9);
    EXPECT_EQ(5u, vec.num_allocations());
}
TEST(StkVectorTest, eightPushBacks_sizeIsEightCapacityIsEight)
{
    stk::Vector<double> vec;
    push_back_num_things(vec, 8);
    EXPECT_EQ(8u, vec.size());
    EXPECT_EQ(8u, vec.capacity());
}
TEST(StkVectorTest, ninePushBacks_sizeIsNineCapacityIsSixteen)
{
    stk::Vector<double> vec;
    push_back_num_things(vec, 9);
    EXPECT_EQ(9u, vec.size());
    EXPECT_EQ(16u, vec.capacity());
}

TEST(StkVectorTest, specifySizeAndValue_haveNumValues)
{
    stk::Vector<double> vec(4, 3.0);
    ASSERT_EQ(4u, vec.size());
    for(size_t i=0; i<vec.size(); i++)
        EXPECT_EQ(3.0, vec[i]);
}
TEST(StkVectorTest, resizeEmptyVector_hasNumEntries)
{
    stk::Vector<double> vec;
    vec.resize(4);
    ASSERT_EQ(4u, vec.size());
}
TEST(StkVectorTest, resizeEmptyVector_hasZeroes)
{
    stk::Vector<double> vec;
    vec.resize(4);
    for(size_t i=0; i<vec.size(); i++)
        EXPECT_EQ(0, vec[i]);
}
TEST(StkVectorTest, resizeEmptyVectorWithInit_hasValue)
{
    stk::Vector<double> vec;
    vec.resize(4, 2);
    for(size_t i=0; i<vec.size(); i++)
        EXPECT_EQ(2, vec[i]);
}
TEST(StkVectorTest, resizeFullVector_hasNewSize)
{
    stk::Vector<double> vec(3);
    vec.resize(5);
    ASSERT_EQ(5u, vec.size());
}
TEST(StkVectorTest, resizeFullVectorWithInit_newEntriesHaveValue)
{
    stk::Vector<double> vec(3, 1);
    vec.resize(5, 2);
    for(size_t i=0; i<3; i++)
        EXPECT_EQ(1, vec[i]);
    for(size_t i=3; i<5; i++)
        EXPECT_EQ(2, vec[i]);
}
TEST(StkVectorTest, resizeToSameSize_doesNothing)
{
    stk::Vector<double> vec(3, 2);
    vec.resize(3);
    ASSERT_EQ(3u, vec.size());
    for(size_t i=0; i<3; i++)
        EXPECT_EQ(2, vec[i]);
}
TEST(StkVectorTest, resizeToSmaller_decreasesSize)
{
    stk::Vector<double> vec(6);
    vec.resize(4);
    ASSERT_EQ(4u, vec.size());
}
TEST(StkVectorTest, resizeToSmaller_noAllocation)
{
    StkVectorSpy vec;
    vec.resize(6);
    EXPECT_EQ(1u, vec.num_allocations());
    vec.resize(4);
    EXPECT_EQ(1u, vec.num_allocations());
}
TEST(StkVectorTest, resizeToSmallerThenLarger_noAllocation)
{
    StkVectorSpy vec;
    vec.resize(6);
    EXPECT_EQ(1u, vec.num_allocations());
    vec.resize(4);
    EXPECT_EQ(1u, vec.num_allocations());
    vec.resize(6);
    EXPECT_EQ(6u, vec.size());
    EXPECT_EQ(1u, vec.num_allocations());
}
TEST(StkVectorTest, sizeZero_isEmpty)
{
    stk::Vector<double> vec;
    EXPECT_TRUE(vec.empty());
}
TEST(StkVectorTest, sizeOne_isNotEmpty)
{
    stk::Vector<double> vec(1);
    EXPECT_TRUE(!vec.empty());
}
TEST(StkVectorTest, clear_empty)
{
    stk::Vector<double> vec(1);
    vec.clear();
    EXPECT_TRUE(vec.empty());
}

}

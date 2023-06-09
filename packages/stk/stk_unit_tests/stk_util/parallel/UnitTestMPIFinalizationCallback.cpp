#include "gtest/gtest.h"
#include "stk_util/parallel/MPIFinalizationCallback.hpp"

namespace {

class DummyClass
{
  public:
    DummyClass(int& destructorCount) :
      m_destructorCallback([&](){destructorCount++;})
    {}

  private:
    stk::MPIFinalizationCallback m_destructorCallback;
};

class MPIFinalizationCallbackTester : public ::testing::Test
{
  protected:
    int destructorCount = 0;
};

}

TEST_F(MPIFinalizationCallbackTester, CalledByDestructor)
{
  EXPECT_EQ(destructorCount, 0);
  {
    DummyClass tmp(destructorCount);
    EXPECT_EQ(destructorCount, 0);
  }
  EXPECT_EQ(destructorCount, 1);
}

TEST_F(MPIFinalizationCallbackTester, DefaultConstructor)
{
  EXPECT_EQ(destructorCount, 0);
  {
    stk::MPIFinalizationCallback finalizer;
    EXPECT_EQ(destructorCount, 0);
  }
  EXPECT_EQ(destructorCount, 0);
}

TEST_F(MPIFinalizationCallbackTester, EmptyFunction)
{
  EXPECT_EQ(destructorCount, 0);
  {
    std::function<void()> emptyFunction;
    stk::MPIFinalizationCallback finalizer(emptyFunction);
    EXPECT_EQ(destructorCount, 0);
  }
  EXPECT_EQ(destructorCount, 0);
}
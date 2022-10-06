#include "gtest/gtest.h"
#include "stk_util/parallel/MPIDatatypeGenerator.hpp"
#include "stk_util/parallel/MPITagManager.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace {

struct OneChar
{
  explicit OneChar(int offset=0) : one(1 + offset) {}
  unsigned char one;
};

bool operator==(const OneChar& lhs, const OneChar& rhs)
{
  return lhs.one == rhs.one;
}

struct TwoChar
{
  explicit TwoChar(int offset=0) : one(1 + offset), two(2 + offset) {}
  unsigned char one;
  unsigned char two;
};

bool operator==(const TwoChar& lhs, const TwoChar& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two;
}

struct ThreeChar
{
  explicit ThreeChar(int offset=0) : one(1 + offset), two(2 + offset), three(3 + offset) {}
  unsigned char one;
  unsigned char two;
  unsigned char three;
};

bool operator==(const ThreeChar& lhs, const ThreeChar& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two && lhs.three == rhs.three;
}

struct FourChar
{
  explicit FourChar(int offset=0) : one(1 + offset), two(2 + offset), three(3 + offset), four(4 + offset) {}
  unsigned char one;
  unsigned char two;
  unsigned char three;
  unsigned char four;
};

bool operator==(const FourChar& lhs, const FourChar& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two && lhs.three == rhs.three &&
         lhs.four == rhs.four;
}

struct FiveChar
{
  explicit FiveChar(int offset=0) : one(1 + offset), two(2 + offset), three(3 + offset), four(4 + offset), five(5 + offset) {}
  unsigned char one;
  unsigned char two;
  unsigned char three;
  unsigned char four;
  unsigned char five;
};

struct ThirtySevenChar
{
  explicit ThirtySevenChar(int offset=0)
  {
    for (int i=0; i < 37; ++i)
    {
      values[i] = i + 1 + offset;
    }
  }
  unsigned char values[37];

};

bool operator==(const ThirtySevenChar& lhs, const ThirtySevenChar& rhs)
{
  for (int i=0; i < 37; ++i)
  {
    if (lhs.values[i] != rhs.values[i])
      return false;
  }

  return true;
}


bool operator==(const FiveChar& lhs, const FiveChar& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two && lhs.three == rhs.three &&
         lhs.four == rhs.four && lhs.five == rhs.five;
}

struct CharInt
{
  explicit CharInt(int offset=0) : one(1 + offset), two(2 + offset) {}
  unsigned char one;
  int two;
};

bool operator==(const CharInt& lhs, const CharInt& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two;
}

struct IntChar
{
  explicit IntChar(int offset=0) : one(1 + offset), two(2 + offset) {}
  int one;
  unsigned char two;
};

bool operator==(const IntChar& lhs, const IntChar& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two;
}

struct CharIntChar
{
  explicit CharIntChar(int offset=0) : one(1 + offset), two(2 + offset), three(3 + offset) {}

  unsigned char one;
  int two;
  unsigned char three;
};

bool operator==(const CharIntChar& lhs, const CharIntChar& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two && lhs.three == rhs.three;
}

struct LongLongChar
{
  explicit LongLongChar(int offset=0) : one(1 + offset), two(2 + offset) {}
  long long one;
  unsigned char two;
};

bool operator==(const LongLongChar& lhs, const LongLongChar& rhs)
{
  return lhs.one == rhs.one && lhs.two == rhs.two;
}

template <typename T>
void test_datatype_size()
{
  MPI_Datatype dtype = stk::generate_mpi_datatype<T>();

  int size;
  MPI_Type_size(dtype, &size);
  EXPECT_EQ(sizeof(T), size_t(size));

  MPI_Aint lowerBound, extent;
  MPI_Type_get_extent(dtype, &lowerBound, &extent);
  EXPECT_EQ(0, lowerBound);
  EXPECT_EQ(sizeof(T), size_t(extent));
}

template <typename T>
void send_datatype(size_t numElements)
{
  const int src = 0, dest=1;
  MPI_Datatype dtype = stk::generate_mpi_datatype<T>();

  int commrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
  auto tag = stk::get_mpi_tag_manager().get_tag(MPI_COMM_WORLD);

  std::vector<T> sendbuf(numElements), recvbuf(numElements);
  for (size_t i=0; i < numElements; ++i)
    sendbuf[i] = T(i);

  if (commrank == src)
  {
    MPI_Send(sendbuf.data(), numElements, dtype, dest, 
             tag, MPI_COMM_WORLD);

  } else
  {
    MPI_Recv(recvbuf.data(), numElements, dtype, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (commrank == dest)
  {
    for (size_t i=0; i < numElements; ++i)
    {
      EXPECT_EQ(recvbuf[i], T(i));
    }
  }
}

}


TEST(MPIDatatypeGenerator, Sizes)
{
  test_datatype_size<OneChar>();
  test_datatype_size<TwoChar>();
  test_datatype_size<ThreeChar>();
  test_datatype_size<FourChar>();
  test_datatype_size<FiveChar>();
  test_datatype_size<ThirtySevenChar>();


  test_datatype_size<CharInt>();
  test_datatype_size<IntChar>();
  test_datatype_size<CharIntChar>();
  test_datatype_size<LongLongChar>();
}

TEST(MPIDatatypeGenerator, SendSmallMsg)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  send_datatype<OneChar>(1024);
  send_datatype<TwoChar>(1024);
  send_datatype<ThreeChar>(1024);
  send_datatype<FourChar>(1024);
  send_datatype<FiveChar>(1024);
  send_datatype<ThirtySevenChar>(1024);


  send_datatype<CharInt>(1024);
  send_datatype<IntChar>(1024);
  send_datatype<CharIntChar>(1024);
  send_datatype<LongLongChar>(1024);
}

TEST(MPIDatatypeGenerator, SendLargeMsg)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2)
  {
    GTEST_SKIP();
  }

  send_datatype<OneChar>(1024*256);
  send_datatype<TwoChar>(1024*256);
  send_datatype<ThreeChar>(1024*256);
  send_datatype<FourChar>(1024*256);
  send_datatype<FiveChar>(1024*256);
  send_datatype<ThirtySevenChar>(1024*256);


  send_datatype<CharInt>(1024*256);
  send_datatype<IntChar>(1024*256);
  send_datatype<CharIntChar>(1024*256);
  send_datatype<LongLongChar>(1024*256);
}


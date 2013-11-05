#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <algorithm>
#include <vector>
#include <iterator>
#include <cstdlib>

namespace std {
template <typename Key, typename Proc>
std::ostream & operator<<(std::ostream & out, std::pair<stk::search::ident::IdentProc<Key,Proc>,stk::search::ident::IdentProc<Key,Proc> > const& ident)
{
  return out << "[" << ident.first << ":" << ident.second << "]";
}
} // namespace std


namespace {

// For testing purposes


STKUNIT_UNIT_TEST(stk_search, coarse_search_3D)
{
  typedef stk::search::ident::IdentProc<uint64_t, unsigned> Ident;
  typedef stk::search::box::AxisAlignedBoundingBox<Ident, double, 3> Box;
  typedef std::vector<Box> BoxVector;
  typedef std::vector<std::pair<Ident,Ident> > Output;

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  int size = stk::parallel_machine_size(comm);
  int rank = stk::parallel_machine_rank(comm);

  double data[6];

  BoxVector local_domain, local_range;
  data[0] = rank + 0.1; data[1] = 0.0; data[2] = 0.0;
  data[3] = rank + 0.9; data[4] = 1.0; data[5] = 1.0;
  local_domain.push_back(Box(data, Ident(rank*4, rank)));

  data[0] = rank + 0.1; data[1] = 2.0; data[2] = 0.0;
  data[3] = rank + 0.9; data[4] = 3.0; data[5] = 1.0;
  local_domain.push_back(Box(data, Ident(rank*4+1, rank)));

  data[0] = rank + 0.6; data[1] = 0.5; data[2] = 0.0;
  data[3] = rank + 1.4; data[4] = 1.5; data[5] = 1.0;
  local_range.push_back(Box(data, Ident(rank*4+2, rank)));

  data[0] = rank + 0.6; data[1] = 2.5; data[2] = 0.0;
  data[3] = rank + 1.4; data[4] = 3.5; data[5] = 1.0;
  local_range.push_back(Box(data, Ident(rank*4+3, rank)));

  Output output;
  stk::search::FactoryOrder order;
  order.m_communicator = comm;
  stk::search::coarse_search(output, local_range, local_domain, order);

  if (size == 1) {
    STKUNIT_EXPECT_EQ(output.size(), 2u);

    STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(0, 0), Ident(2, 0)));
    STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(1, 0), Ident(3, 0)));
  }
  else {
    if (rank == 0) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(0, 0), Ident(2, 0)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(1, 0), Ident(3, 0)));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(Ident(4, 1), Ident(2, 0)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(Ident(5, 1), Ident(3, 0)));
    }
    else if (rank == size - 1) {
      STKUNIT_EXPECT_EQ(output.size(), 4u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(rank*4,     rank), Ident(rank*4 - 2, rank - 1)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(rank*4,     rank), Ident(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(Ident(rank*4 + 1, rank), Ident(rank*4 - 1, rank - 1)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(Ident(rank*4 + 1, rank), Ident(rank*4 + 3, rank   )));
    }
    else {
      STKUNIT_EXPECT_EQ(output.size(), 6u);

      STKUNIT_EXPECT_EQ(output[0], std::make_pair(Ident(rank*4,         rank    ), Ident(rank*4 - 2, rank - 1)));
      STKUNIT_EXPECT_EQ(output[1], std::make_pair(Ident(rank*4,         rank    ), Ident(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[2], std::make_pair(Ident(rank*4 + 1,     rank    ), Ident(rank*4 - 1, rank - 1)));
      STKUNIT_EXPECT_EQ(output[3], std::make_pair(Ident(rank*4 + 1,     rank    ), Ident(rank*4 + 3, rank   )));
      STKUNIT_EXPECT_EQ(output[4], std::make_pair(Ident((rank+1)*4,     rank + 1), Ident(rank*4 + 2, rank   )));
      STKUNIT_EXPECT_EQ(output[5], std::make_pair(Ident((rank+1)*4 + 1, rank + 1), Ident(rank*4 + 3, rank   )));
    }
  }
}

} //namespace

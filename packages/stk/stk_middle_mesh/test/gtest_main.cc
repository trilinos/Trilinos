#ifndef STK_BUILT_IN_SIERRA

#include "mpi.h"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  auto result = RUN_ALL_TESTS();

  MPI_Finalize();
  return result;
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
#endif

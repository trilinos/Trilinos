#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>

int main(int argc, char **argv)
{
  boost::mpi::environment env(argc,argv);
  ::testing::InitGoogleTest(&argc,argv);

  int error = RUN_ALL_TESTS();

  // Make stk_samba unit-tests compatible with sierra test framework
  // *without* adding dependency to stk_util.
  if (error == 0) {
    std::cout << "STKUNIT_ALL_PASS" << std::endl;
  }

  return error;
}


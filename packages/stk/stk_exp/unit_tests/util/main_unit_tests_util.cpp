#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>

int main( int argc, char **argv )
{
  boost::mpi::environment env(argc,argv);
  ::testing::InitGoogleTest(&argc,argv);

  return RUN_ALL_TESTS();
}

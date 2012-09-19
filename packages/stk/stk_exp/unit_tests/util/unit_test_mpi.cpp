#include <gtest/gtest.h>
#include <iostream>

#include <boost/mpi/communicator.hpp>

TEST( MPI, communicator )
{
  boost::mpi::communicator world;
}

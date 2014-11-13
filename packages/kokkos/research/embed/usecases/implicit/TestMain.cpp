
#include <Kokkos_Macros.hpp>
#include <ParallelComm.hpp>
#include <sstream>
#include <iostream>

namespace Test {

int test_host( comm::Machine , std::istream & );
int test_cuda( comm::Machine , std::istream & );

}

int main( int argc , char ** argv )
{
  comm::Machine machine = comm::Machine::init( & argc , & argv );

  // Broadcast command line to all processes as a string
  const std::string cmd = comm::command_line( machine , argc , argv );

  std::istringstream input( cmd );

  std::string which ; input >> which ;

  if ( which == std::string("host") ) {
    Test::test_host( machine , input );
  }
#if defined(KOKKOS_HAVE_CUDA)
  else if ( which == std::string("cuda") ) {
    Test::test_cuda( machine , input );
  }
#endif
  else {
    std::cerr << "Expected host OR cuda" << std::endl ;
  }

  comm::Machine::finalize();

  return 0 ;
}


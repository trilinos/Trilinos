
#include <string>
#include <sstream>
#include <iostream>
#include <ParallelComm.hpp>

void test_host( comm::Machine machine , std::istream & );
void test_cuda( comm::Machine machine , std::istream & );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  comm::Machine machine = comm::Machine::init( & argc , & argv );

  std::string argline ;
  if ( 0 == comm::rank( machine ) ) {
    for ( int i = 1 ; i < argc ; ++i ) {
      argline.append(" ").append( argv[i] );
    }
  }

#ifdef HAVE_MPI
  {
    int length = argline.length();
    MPI_Bcast( & length , 1 , MPI_INT , 0 , machine.mpi_comm );
    argline.resize( length , ' ' );
    MPI_Bcast( (void*) argline.data() , length , MPI_CHAR , 0 , machine.mpi_comm );
  }
#endif /* HAVE_MPI */

  std::istringstream input( argline );

  std::string which ; input >> which ;

  if ( which == std::string("host") ) {
    test_host( machine , input );
  }
  else if ( which == std::string("cuda") ) {
#if HAVE_CUDA
    test_cuda( machine , input );
#endif
  }
  else {
    std::cout << "command line = (host|cuda)" << std::endl ;
  }

  comm::Machine::finalize();

  return 0 ;
}


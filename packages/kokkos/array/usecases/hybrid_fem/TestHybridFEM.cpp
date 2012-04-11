
#include <ParallelComm.hpp>

void test_host( comm::Machine machine );
void test_cuda( comm::Machine machine );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  comm::Machine machine = comm::Machine::init( & argc , & argv );

  test_host( machine );

#if HAVE_CUDA
  test_cuda( machine );
#endif

  comm::Machine::finalize();

  return 0 ;
}


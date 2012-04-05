
#include <ParallelDistributedComm.hpp>

void test_host( comm::Machine machine );
void test_cuda( comm::Machine machine );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  comm::Machine machine = comm::Machine::init( & argc , & argv );

#if HAVE_CUDA
  test_cuda( machine );
#endif

  test_host( machine );

  comm::Machine::finalize();

  return 0 ;
}


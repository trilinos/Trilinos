
#include <string>
#include <sstream>
#include <iostream>
#include <ParallelComm.hpp>

//----------------------------------------------------------------------------

void test_host_explicit( comm::Machine machine , 
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t uq_count_begin ,
                         size_t uq_count_end ,
                         size_t count_run );

//----------------------------------------------------------------------------

void test_cuda_explicit( comm::Machine machine ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t uq_count_begin ,
                         size_t uq_count_end ,
                         size_t count_run );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace {

void run( const std::string & argline , comm::Machine machine )
{
  std::istringstream input( argline );

  bool cmd_error = false ;

  std::string which ; input >> which ;

  if ( which == std::string("host") ) {
    size_t host_node_count = 0 ;
    size_t host_node_thread_count = 1 ;

    input >> host_node_count ;
    input >> host_node_thread_count ;
    input >> which ;

    if ( which == std::string("explicit") ) {
 
      size_t mesh_node_begin = 100 ;
      size_t mesh_node_end   = 300 ;
      size_t uq_count_begin  = 1 ;
      size_t uq_count_end    = 1 ;
      size_t run             =   1 ;
      input >> mesh_node_begin >> mesh_node_end
            >> uq_count_begin >> uq_count_end
            >> run ;
      test_host_explicit( machine ,
                          host_node_count ,
                          host_node_thread_count ,
                          mesh_node_begin ,
                          mesh_node_end ,
                          uq_count_begin ,
                          uq_count_end ,
                          run );

    }
    else {
      cmd_error = true ;
    }
  }
#if HAVE_CUDA
  else if ( which == std::string("cuda") ) {
  
    input >> which ;

    if ( which == std::string("explicit") ) {
  
      size_t mesh_node_begin = 100 ;
      size_t mesh_node_end   = 300 ;
      size_t uq_count_begin  = 1 ;
      size_t uq_count_end    = 1 ;
      size_t run             =   1 ;
      input >> mesh_node_begin >> mesh_node_end
            >> uq_count_begin >> uq_count_end
            >> run ;
      test_cuda_explicit( machine ,
                          mesh_node_begin ,
                          mesh_node_end ,
                          uq_count_begin ,
                          uq_count_end ,
                          run );
   
    }
  }
#endif
  else {
    cmd_error = true ;
  }

  if ( cmd_error && 0 == comm::rank( machine ) ) {
    std::cout << "Expecting command line with" << std::endl
              << "    host NumNumaNode NumThreadPerNode <test>" << std::endl
              << "    cuda <test>" << std::endl
              << "where <test> is" << std::endl
              << "    explicit  NumElem_begin NumElem_end NumUQ_begin NumUQ_end NumRun"
              << std::endl ;
  }
}

} // namespace

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  comm::Machine machine = comm::Machine::init( & argc , & argv );

  const unsigned comm_rank = comm::rank( machine );

  // Turn command line into a string,
  // broadcast it to all processes,
  // turn string into an input stream for parsing.

  std::string argline ;

  if ( 0 == comm_rank ) {
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

  try {
    run( argline , machine );
  }
  catch( const std::exception & x ) {
    std::cerr << "P" << comm_rank << " throw: " << x.what() << std::endl ;
  }
  catch( ... ) {
    std::cerr << "P" << comm_rank << " throw: unknown exception" << std::endl ;
  }

  comm::Machine::finalize();

  return 0 ;
}


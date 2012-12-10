
#include <string>
#include <sstream>
#include <iostream>
#include <ParallelComm.hpp>

//----------------------------------------------------------------------------

void test_box_partition( bool print );

//----------------------------------------------------------------------------

void test_host_query( comm::Machine );

void test_host_fixture( comm::Machine machine ,
                        size_t gang_count ,
                        size_t gang_worker_count ,
                        size_t nx , size_t ny , size_t nz );

void test_host_implicit( comm::Machine machine , 
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_host_explicit( comm::Machine machine , 
                         size_t gang_count ,
                         size_t gang_worker_count ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_host_nonlinear( comm::Machine machine , 
                          size_t gang_count ,
                          size_t gang_worker_count ,
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run );

void test_host_nonlinear_quadratic( comm::Machine machine , 
                                    size_t gang_count ,
                                    size_t gang_worker_count ,
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run );


//----------------------------------------------------------------------------

void test_cuda_query( comm::Machine );

void test_cuda_fixture( comm::Machine machine ,
                        size_t nx , size_t ny , size_t nz );

void test_cuda_implicit( comm::Machine machine ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_cuda_explicit( comm::Machine machine ,
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t count_run );

void test_cuda_nonlinear( comm:: Machine machine ,
                          size_t elem_count_begin ,
                          size_t elem_count_end ,
                          size_t count_run );

void test_cuda_nonlinear_quadratic( comm::Machine machine , 
                                    size_t elem_count_begin ,
                                    size_t elem_count_end ,
                                    size_t count_run );


//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace {

void run( const std::string & argline , comm::Machine machine )
{
  std::istringstream input( argline );

  bool cmd_error = false ;

  std::string which ; input >> which ;

  if ( which == std::string("query") ) {
    test_host_query( machine );
#if HAVE_CUDA
    test_cuda_query( machine );
#endif
  }
  else if ( which == std::string("partition") ) {
    if ( 0 == comm::rank( machine ) ) {
      test_box_partition( false /* print flag */ );
    }
  }
  else {
    if ( which == std::string("host") ) {
      size_t host_gang_count = 0 ;
      size_t host_gang_worker_count = 1 ;

      input >> host_gang_count ;
      input >> host_gang_worker_count ;
      input >> which ;
      if ( which == std::string("fixture") ) {
 
        size_t nx = 0 , ny = 0 , nz = 0 ;
        input >> nx >> ny >> nz ;
        test_host_fixture( machine , host_gang_count , host_gang_worker_count , nx , ny , nz );
 
      }
      else if ( which == std::string("explicit") ) {
 
        size_t mesh_node_begin = 100 ;
        size_t mesh_node_end   = 300 ;
        size_t run             =   1 ;
        input >> mesh_node_begin >> mesh_node_end >> run ;
        test_host_explicit( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );
 
      }
      else if ( which == std::string("implicit") ) {
 
        size_t mesh_node_begin = 100 ;
        size_t mesh_node_end   = 300 ;
        size_t run             =   1 ;
        input >> mesh_node_begin >> mesh_node_end >> run ;
        test_host_implicit( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );
 
      }
      else if ( which == std::string("nonlinear") ) {
 
        size_t mesh_node_begin = 100 ;
        size_t mesh_node_end   = 300 ;
        size_t run             =   1 ;
        input >> mesh_node_begin >> mesh_node_end >> run ;
        test_host_nonlinear( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );
 
      }   
      else if ( which == std::string("nonlinear_quadratic") ) {
 
        size_t mesh_node_begin = 100 ;
        size_t mesh_node_end   = 300 ;
        size_t run             =   1 ;
        input >> mesh_node_begin >> mesh_node_end >> run ;
        test_host_nonlinear_quadratic( machine , host_gang_count , host_gang_worker_count , mesh_node_begin , mesh_node_end , run );
 
      }   
      else {
        cmd_error = true ;
      }
    }
#if HAVE_CUDA
    else if ( which == std::string("cuda") ) {
  
      input >> which ;
  
      if ( which == std::string("fixture") ) {
  
        size_t nx = 0 , ny = 0 , nz = 0 ;
        input >> nx >> ny >> nz ;
        test_cuda_fixture( machine , nx , ny , nz );
  
      }
      else if ( which == std::string("implicit") ) {
  
        size_t mesh_node_begin = 100 ;
        size_t mesh_node_end   = 300 ;
        size_t run             =   1 ;
        input >> mesh_node_begin >> mesh_node_end >> run ;
        test_cuda_implicit( machine , mesh_node_begin , mesh_node_end , run );
     
      }
      else if ( which == std::string("explicit") ) {
  
        size_t mesh_node_begin = 100 ;
        size_t mesh_node_end   = 300 ;
        size_t run             =   1 ;
        input >> mesh_node_begin >> mesh_node_end >> run ;
        test_cuda_explicit( machine , mesh_node_begin , mesh_node_end , run );
     
      }
      else if ( which == std::string("nonlinear") ) {

	size_t mesh_node_begin = 100;
	size_t mesh_node_end   = 300;
	size_t run             =   1;
        input >> mesh_node_begin >> mesh_node_end >> run ;
	test_cuda_nonlinear( machine , mesh_node_begin, mesh_node_end, run ); 

      }
      else if ( which == std::string("nonlinear_quadratic") ) {

	size_t mesh_node_begin = 100;
	size_t mesh_node_end   = 300;
	size_t run             =   1;
        input >> mesh_node_begin >> mesh_node_end >> run ;
	test_cuda_nonlinear_quadratic( machine , mesh_node_begin, mesh_node_end, run ); 

      }
    }
#endif
    else {
      cmd_error = true ;
    }
  }

  if ( cmd_error && 0 == comm::rank( machine ) ) {
    std::cout << "Expecting command line with" << std::endl
              << "    query" << std::endl
              << "    partition" << std::endl
              << "    host NumNumaNode NumThreadPerNode <test>" << std::endl
              << "    cuda <test>" << std::endl
              << "where <test> is" << std::endl
              << "    fixture   NumElemX NumElemY NumElemZ" << std::endl
              << "    implicit  NumElemBegin NumElemEnd NumRun" << std::endl 
              << "    explicit  NumElemBegin NumElemEnd NumRun" << std::endl
              << "    nonlinear NumElemBegin NumElemEnd NumRun" << std::endl 
              << "    nonlinear_quadratic NumElemBegin NumElemEnd NumRun" << std::endl ;

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


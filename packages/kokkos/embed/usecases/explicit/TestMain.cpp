
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iostream>

//----------------------------------------------------------------------------

void test_host_explicit( size_t numa_node_count ,
                         size_t numa_node_thread_count ,
                         size_t elem_count ,
                         size_t iter_count );

//----------------------------------------------------------------------------

void test_cuda_explicit( size_t elem_count ,
                         size_t iter_count );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  size_t elem_count = 1 ; 
  size_t iter_count = 1 ; 
  size_t numa_node_thread_count = 1 ;
  size_t numa_node_count = 1 ;

  if ( 1 < argc ) { elem_count = atoi( argv[1] ); }
  if ( 2 < argc ) { iter_count = atoi( argv[2] ); }
  if ( 3 < argc ) { numa_node_thread_count = atoi( argv[3] ); }
  if ( 4 < argc ) { numa_node_count = atoi( argv[4] ); }

  std::cout << argv[0]
            << " elem( " << elem_count
            << " ) iter( " << iter_count
            << " ) thread( " << numa_node_thread_count
            << " ) numa( " << numa_node_count
            << " )" << std::endl ;


  test_host_explicit( numa_node_count ,
                      numa_node_thread_count ,
                      elem_count ,
                      iter_count );

  test_cuda_explicit( elem_count ,
                      iter_count );

  return 0 ;
}


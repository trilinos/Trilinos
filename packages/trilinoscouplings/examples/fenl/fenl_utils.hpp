#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include <utility>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>

//----------------------------------------------------------------------------

#include <Teuchos_Comm.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Command line processing:

enum { CMD_USE_THREADS = 0
       , CMD_USE_NUMA
       , CMD_USE_CORE_PER_NUMA
       , CMD_USE_CUDA
       , CMD_USE_OPENMP
       , CMD_USE_CUDA_DEV
       , CMD_USE_FIXTURE_X
       , CMD_USE_FIXTURE_Y
       , CMD_USE_FIXTURE_Z
       , CMD_USE_FIXTURE_BEGIN
       , CMD_USE_FIXTURE_END
       , CMD_USE_FIXTURE_QUADRATIC
       , CMD_USE_ATOMIC
       , CMD_USE_TRIALS
       , CMD_USE_UQ_ENSEMBLE
       , CMD_USE_UQ_DIM
       , CMD_USE_UQ_ORDER
       , CMD_PRINT
       , CMD_ECHO
       , CMD_ERROR
       , CMD_COUNT };

void print_cmdline( std::ostream & s , const int cmd[] )
{
  if ( cmd[ CMD_USE_THREADS ] ) {
    s << " Threads(" << cmd[ CMD_USE_THREADS ]
      << ") NUMA(" << cmd[ CMD_USE_NUMA ]
      << ") CORE_PER_NUMA(" << cmd[ CMD_USE_CORE_PER_NUMA ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_OPENMP ] ) {
    s << " OpenMP(" << cmd[ CMD_USE_OPENMP ]
      << ") NUMA(" << cmd[ CMD_USE_NUMA ]
      << ") CORE_PER_NUMA(" << cmd[ CMD_USE_CORE_PER_NUMA ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_X ] ) {
    s << " Fixture(" << cmd[ CMD_USE_FIXTURE_X ]
      << "x" << cmd[ CMD_USE_FIXTURE_Y ]
      << "x" << cmd[ CMD_USE_FIXTURE_Z ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_BEGIN ] ) {
    s << " Fixture( " << cmd[ CMD_USE_FIXTURE_BEGIN ]
      << " .. " << cmd[ CMD_USE_FIXTURE_END ]
      << " )" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) {
    s << " Quadratic-Element" ;
  }
  if ( cmd[ CMD_USE_UQ_ENSEMBLE ] ) {
    s << " UQ ensemble" ;
  }
  if ( cmd[ CMD_USE_UQ_DIM ] ) {
    s << " UQ dimension(" << cmd[ CMD_USE_UQ_DIM ] << ")" << ")" ;
  }
  if ( cmd[ CMD_USE_UQ_DIM ] ) {
    s << " UQ order(" << cmd[ CMD_USE_UQ_ORDER ] << ")" << ")" ;
  }
  if ( cmd[ CMD_USE_CUDA ] ) {
    s << " CUDA(" << cmd[ CMD_USE_CUDA_DEV ] << ")" ;
  }
  if ( cmd[ CMD_USE_ATOMIC ] ) {
    s << " ATOMIC" ;
  }
  if ( cmd[ CMD_USE_TRIALS ] ) {
    s << " TRIALS(" << cmd[ CMD_USE_TRIALS ] << ")" ;
  }
  if ( cmd[ CMD_PRINT ] ) {
    s << " PRINT" ;
  }
  s << std::endl ;
}

std::vector< size_t >
print_headers( std::ostream & s , const int cmd[] , const int comm_rank )
{
  if ( 0 == comm_rank ) {
    if ( cmd[ CMD_USE_THREADS ] ) { s << "THREADS , " << cmd[ CMD_USE_THREADS ] ; }
    else if ( cmd[ CMD_USE_OPENMP ] ) { s << "OPENMP , " << cmd[ CMD_USE_OPENMP ] ; }
    else if ( cmd[ CMD_USE_CUDA ] ) { s << "CUDA" ; }

    if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) { s << " , QUADRATIC-ELEMENT" ; }
    else { s << " , LINEAR-ELEMENT" ; }

    if ( cmd[ CMD_USE_ATOMIC ] ) { s << " , USING ATOMICS" ; }
    if ( cmd[ CMD_USE_UQ_ENSEMBLE ] ) { s << " , USING UQ ENSEMBLE" ; }
    if ( cmd[ CMD_USE_UQ_DIM ] ) { s << " , UQ DIM , " << cmd[ CMD_USE_UQ_DIM ]; }
    if ( cmd[ CMD_USE_UQ_ORDER ] ) { s << " , UQ ORDER , " << cmd[ CMD_USE_UQ_ORDER ]; }
  }

  std::vector< std::pair<std::string,std::string> > headers;

  headers.push_back(std::make_pair("ELEMS","count"));
  headers.push_back(std::make_pair("NODES","count"));
  headers.push_back(std::make_pair("NEWTON","iter"));
  headers.push_back(std::make_pair("CG","iter"));
  headers.push_back(std::make_pair("MAP_RATIO","ratio"));
  headers.push_back(std::make_pair("SET_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("SCAN/NODE","millisec"));
  headers.push_back(std::make_pair("GRAPH_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("SORT/NODE","millisec"));
  headers.push_back(std::make_pair("ELEM_GRAPH_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("MATRIX_CREATE/NODE","millisec"));
  headers.push_back(std::make_pair("MATRIX_FILL/NODE","millisec"));
  headers.push_back(std::make_pair("BOUNDARY/NODE","millisec"));
  headers.push_back(std::make_pair("CG/ITER/ROW","millisec"));
  headers.push_back(std::make_pair("ERROR","ratio"));
  headers.push_back(std::make_pair("RESPONSE","mean"));
  headers.push_back(std::make_pair("RESPONSE","std.dev."));

  // find print widths
  size_t min_width = 10;
  std::vector< size_t > widths(headers.size());
  for (size_t i=0, ie=headers.size(); i<ie; ++i)
    widths[i] = std::max(min_width, headers[i].first.size()+1);

  // print column headers
  if ( 0 == comm_rank ) {
    s << std::endl ;
    for (size_t i=0; i<headers.size(); ++i)
      s << std::setw(widths[i]) << headers[i].first << " ,";
    s << "\b\b  " << std::endl;
    for (size_t i=0; i<headers.size(); ++i)
      s << std::setw(widths[i]) << headers[i].second << " ,";
    s << "\b\b  " << std::endl;

    s << std::scientific;
    s.precision(3);
  }

  return widths;
}

void print_perf_value( std::ostream & s ,
                       const std::vector<size_t> & widths ,
                       const Kokkos::Example::FENL::Perf & perf )
{
  int i=0;
  s << std::setw(widths[i++]) << perf.global_elem_count << " ,";
  s << std::setw(widths[i++]) << perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << double(perf.newton_iter_count) / perf.uq_count << " ,";
  s << std::setw(widths[i++]) << double(perf.cg_iter_count) / perf.uq_count << " ,";
  s << std::setw(widths[i++]) << perf.map_ratio << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_node_set * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( perf.scan_node_count * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_graph_entries * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( perf.sort_graph_entries * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_element_graph * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( perf.create_sparse_matrix * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( perf.bc_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << ( ( perf.cg_time * 1000.0 ) / perf.cg_iter_count ) / (perf.global_node_count*perf.uq_count) << " ,";
  s << std::setw(widths[i++]) << perf.error_max / perf.uq_count << " ,";
  s << std::setw(widths[i++]) << perf.response_mean << " ,";
  s << std::setw(widths[i])   << perf.response_std_dev;
  s << std::endl ;
}

void parse_cmdline( int argc , char ** argv, int cmdline[],
                    const Teuchos::Comm<int>& comm )
{
  const int comm_rank = comm.getRank();

  for ( int i = 0 ; i < CMD_COUNT ; ++i ) cmdline[i] = 0 ;

  if ( 0 == comm_rank ) {
    for ( int i = 1 ; i < argc ; ++i ) {
      if ( 0 == strcasecmp( argv[i] , "threads" ) ) {
        cmdline[ CMD_USE_THREADS ] = atoi( argv[++i] );
      }
      else if ( 0 == strcasecmp( argv[i] , "openmp" ) ) {
        cmdline[ CMD_USE_OPENMP ] = atoi( argv[++i] );
      }
      else if ( 0 == strcasecmp( argv[i] , "cores" ) ) {
        sscanf( argv[++i] , "%dx%d" ,
                cmdline + CMD_USE_NUMA ,
                cmdline + CMD_USE_CORE_PER_NUMA );
      }
      else if ( 0 == strcasecmp( argv[i] , "cuda" ) ) {
        cmdline[ CMD_USE_CUDA ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "cuda-dev" ) ) {
        cmdline[ CMD_USE_CUDA ] = 1 ;
        cmdline[ CMD_USE_CUDA_DEV ] = atoi( argv[++i] ) ;
      }
      else if ( 0 == strcasecmp( argv[i] , "fixture" ) ) {
        sscanf( argv[++i] , "%dx%dx%d" ,
                cmdline + CMD_USE_FIXTURE_X ,
                cmdline + CMD_USE_FIXTURE_Y ,
                cmdline + CMD_USE_FIXTURE_Z );
      }
      else if ( 0 == strcasecmp( argv[i] , "fixture-range" ) ) {
        sscanf( argv[++i] , "%d..%d" ,
                cmdline + CMD_USE_FIXTURE_BEGIN ,
                cmdline + CMD_USE_FIXTURE_END );
      }
      else if ( 0 == strcasecmp( argv[i] , "fixture-quadratic" ) ) {
        cmdline[ CMD_USE_FIXTURE_QUADRATIC ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "ensemble" ) ) {
        cmdline[ CMD_USE_UQ_ENSEMBLE ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "uq-dim" ) ) {
        cmdline[ CMD_USE_UQ_DIM ] = atoi( argv[++i] ) ;
      }
      else if ( 0 == strcasecmp( argv[i] , "uq-order" ) ) {
        cmdline[ CMD_USE_UQ_ORDER ] = atoi( argv[++i] ) ;
      }
      else if ( 0 == strcasecmp( argv[i] , "atomic" ) ) {
        cmdline[ CMD_USE_ATOMIC ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "trials" ) ) {
        cmdline[ CMD_USE_TRIALS ] = atoi( argv[++i] ) ;
      }
      else if ( 0 == strcasecmp( argv[i] , "print" ) ) {
        cmdline[ CMD_PRINT ] = 1 ;
      }
      else if ( 0 == strcasecmp( argv[i] , "echo" ) ) {
        cmdline[ CMD_ECHO ] = 1 ;
      }
      else {
        cmdline[ CMD_ERROR ] = 1 ;

        std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl ;
      }
    }

    if ( ! cmdline[ CMD_USE_TRIALS ] ) { cmdline[ CMD_USE_TRIALS ] = 1 ; }

    if ( ! cmdline[ CMD_USE_FIXTURE_X ] &&
         ! cmdline[ CMD_USE_FIXTURE_BEGIN ] ) {
      cmdline[ CMD_USE_FIXTURE_X ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Y ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Z ] = 2 ;
    }

    if ( cmdline[ CMD_ECHO ] ) { print_cmdline( std::cout , cmdline ); }
  }

  //--------------------------------------------------------------------------

  comm.broadcast( int(0) , int(CMD_COUNT * sizeof(int)) , (char *) cmdline );

}


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

#include <KokkosCore_config.h>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Threads.hpp>

#if defined( KOKKOS_HAVE_CUDA )
#include <Kokkos_Cuda.hpp>
#endif

#if defined( KOKKOS_HAVE_OPENMP )
#include <Kokkos_OpenMP.hpp>
#endif

#include <fenl.hpp>

//----------------------------------------------------------------------------

#include <Tpetra_Version.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// #include <Tpetra_Vector.hpp>

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

void print_perf_value( std::ostream & s , const std::vector<size_t> & widths,  const Kokkos::Example::FENL::Perf & perf )
{
  int i=0;
  s << std::setw(widths[i++]) << perf.global_elem_count << " ,";
  s << std::setw(widths[i++]) << perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << perf.newton_iter_count << " ,";
  s << std::setw(widths[i++]) << perf.cg_iter_count << " ,";
  s << std::setw(widths[i++]) << perf.map_ratio << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_node_set * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.scan_node_count * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_graph_entries * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.sort_graph_entries * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_element_graph * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.create_sparse_matrix * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.fill_time * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( perf.bc_time * 1000.0 ) / perf.global_node_count << " ,";
  s << std::setw(widths[i++]) << ( ( perf.cg_time * 1000.0 ) / perf.cg_iter_count ) / perf.global_node_count << " ,";
  s << std::setw(widths[i])   << perf.error_max;
  s << std::endl ;
}

//----------------------------------------------------------------------------

template< class Device , Kokkos::Example::BoxElemPart::ElemOrder ElemOrder >
bool run( const Teuchos::RCP<const Teuchos::Comm<int> > & comm , const int cmd[] )
{
  typedef typename ::Kokkos::Compat::KokkosDeviceWrapperNode<Device> NodeType;
  bool success = true;
  try {

  const int comm_rank = comm->getRank();

  if ( 0 == comm_rank ) {
    if ( cmd[ CMD_USE_THREADS ] ) { std::cout << "THREADS , " << cmd[ CMD_USE_THREADS ] ; }
    else if ( cmd[ CMD_USE_OPENMP ] ) { std::cout << "OPENMP , " << cmd[ CMD_USE_OPENMP ] ; }
    else if ( cmd[ CMD_USE_CUDA ] ) { std::cout << "CUDA" ; }

    if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) { std::cout << " , QUADRATIC-ELEMENT" ; }
    else { std::cout << " , LINEAR-ELEMENT" ; }

    if ( cmd[ CMD_USE_ATOMIC ] ) { std::cout << " , USING ATOMICS" ; }
  }

  std::vector< std::pair<std::string,std::string> > headers;

  Teuchos::ParameterList params;
  params.set("Verbose",     0);
  if ( cmd[ CMD_USE_THREADS ] )
    params.set("Num Threads", cmd[CMD_USE_THREADS]);
  if ( cmd[ CMD_USE_NUMA ] && cmd[ CMD_USE_CORE_PER_NUMA ] ) {
    params.set("Num NUMA", cmd[ CMD_USE_NUMA ]);
    params.set("Num CoresPerNUMA", cmd[ CMD_USE_CORE_PER_NUMA ]);
  }
  if ( cmd[ CMD_USE_CUDA_DEV ] )
    params.set("Device", cmd[ CMD_USE_CUDA_DEV ] );

  // Create Tpetra Node
  ::Teuchos::RCP<NodeType> node = ::Teuchos::rcp (new NodeType(params));

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

  // find print widths
  size_t min_width = 10;
  std::vector< size_t > widths(headers.size());
  for (size_t i=0, ie=headers.size(); i<ie; ++i)
    widths[i] = std::max(min_width, headers[i].first.size()+1);

  // print column headers
  if ( 0 == comm_rank ) {
    std::cout << std::endl ;
    for (size_t i=0; i<headers.size(); ++i)
      std::cout << std::setw(widths[i]) << headers[i].first << " ,";
    std::cout << "\b\b  " << std::endl;
    for (size_t i=0; i<headers.size(); ++i)
      std::cout << std::setw(widths[i]) << headers[i].second << " ,";
    std::cout << "\b\b  " << std::endl;

    std::cout << std::scientific;
    std::cout.precision(3);
  }

  if ( cmd[ CMD_USE_FIXTURE_BEGIN ] ) {
    for ( int i = cmd[CMD_USE_FIXTURE_BEGIN] ; i < cmd[CMD_USE_FIXTURE_END] * 2 ; i *= 2 ) {
      int nelem[3] ;
      nelem[0] = std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
      nelem[1] = 1 + nelem[0] ;
      nelem[2] = 2 * nelem[0] ;

      const Kokkos::Example::FENL::Perf perf =
        cmd[ CMD_USE_FIXTURE_QUADRATIC ]
        ? Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemQuadratic >
            ( comm , node , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
        : Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemLinear >
            ( comm , node , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
        ;

      if ( 0 == comm_rank ) { print_perf_value( std::cout , widths, perf ); }
    }
  }
  else {
    int nelem[3] = { cmd[ CMD_USE_FIXTURE_X ] ,
                     cmd[ CMD_USE_FIXTURE_Y ] ,
                     cmd[ CMD_USE_FIXTURE_Z ] };

    const Kokkos::Example::FENL::Perf perf =
      cmd[ CMD_USE_FIXTURE_QUADRATIC ]
      ? Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemQuadratic >
          ( comm , node , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
      : Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemLinear >
          ( comm , node , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
      ;

    if ( 0 == comm_rank ) { print_perf_value( std::cout , widths, perf ); }
  }


  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return success;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);

  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  const int comm_rank = comm->getRank();
  Teuchos::ParameterList params;

  //--------------------------------------------------------------------------

  int cmdline[ CMD_COUNT ] ;

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

    if ( cmdline[ CMD_ECHO ] ) { print_cmdline( std::cout , cmdline ); }
  }

  //--------------------------------------------------------------------------

  comm->broadcast( int(0) , int(CMD_COUNT * sizeof(int)) , (char *) cmdline );

  if ( ! cmdline[ CMD_ERROR ] && ! cmdline[ CMD_ECHO ] ) {

    if ( ! cmdline[ CMD_USE_TRIALS ] ) { cmdline[ CMD_USE_TRIALS ] = 1 ; }

    if ( ! cmdline[ CMD_USE_FIXTURE_X ] && ! cmdline[ CMD_USE_FIXTURE_BEGIN ] ) {
      cmdline[ CMD_USE_FIXTURE_X ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Y ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Z ] = 2 ;
    }

#if defined( KOKKOS_HAVE_PTHREAD )

    if ( cmdline[ CMD_USE_THREADS ] ) {
      run< Kokkos::Threads , Kokkos::Example::BoxElemPart::ElemLinear >( comm , cmdline );
    }

#endif

#if defined( KOKKOS_HAVE_OPENMP )

    if ( cmdline[ CMD_USE_OPENMP ] ) {
      run< Kokkos::OpenMP , Kokkos::Example::BoxElemPart::ElemLinear >( comm , cmdline );
    }

#endif

#if defined( KOKKOS_HAVE_CUDA )
    if ( cmdline[ CMD_USE_CUDA ] ) {
      run< Kokkos::Cuda , Kokkos::Example::BoxElemPart::ElemLinear >( comm , cmdline );
    }

#endif

  }

  //--------------------------------------------------------------------------

  return cmdline[ CMD_ERROR ] ? -1 : 0 ;
}


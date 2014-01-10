
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include <utility>
#include <string>
#include <sstream>
#include <iostream>

#include <KokkosCore_config.h>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Threads.hpp>

#if defined( KOKKOS_HAVE_CUDA )
#include <Kokkos_Cuda.hpp>
#endif

#if defined( KOKKOS_HAVE_OPENMP )
#include <Kokkos_OpenMP.hpp>
#endif

#include <WrapMPI.hpp>
#include <fenl.hpp>

//----------------------------------------------------------------------------

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

void print_perf_value( std::ostream & s , const Kokkos::Example::FENL::Perf & perf )
{
  s << perf.global_elem_count << " , "
    << perf.global_node_count << " , "
    << perf.newton_iter_count << " , "
    << perf.cg_iter_count << " , "
    << ( perf.graph_time * 1000.0 ) / perf.global_node_count << " , "
    << ( perf.fill_time * 1000.0 ) / perf.global_node_count << " , "
    << ( perf.bc_time * 1000.0 ) / perf.global_node_count << " , "
    << ( ( perf.cg_time * 1000.0 ) / perf.cg_iter_count ) / perf.global_node_count << " , "
    << perf.error_max
    << std::endl ;
}

template< class Device , Kokkos::Example::BoxElemPart::ElemOrder ElemOrder >
void run( MPI_Comm comm , const int cmd[] )
{
  if ( cmd[ CMD_USE_THREADS ] ) { std::cout << "THREADS , " << cmd[ CMD_USE_THREADS ] ; }
  else if ( cmd[ CMD_USE_OPENMP ] ) { std::cout << "OPENMP , " << cmd[ CMD_USE_OPENMP ] ; }
  else if ( cmd[ CMD_USE_CUDA ] ) { std::cout << "CUDA" ; }

  if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) { std::cout << " , QUADRATIC-ELEMENT" ; }
  else { std::cout << " , LINEAR-ELEMENT" ; }

  if ( cmd[ CMD_USE_ATOMIC ] ) { std::cout << " , USING ATOMICS" ; }

  std::cout << std::endl ;
  std::cout << "ELEMS , NODES , NEWTON , CG   , GRAPH/NODE , FILL/NODE , BOUNDARY/NODE , CG/ITER/ROW , ERROR" << std::endl ;
  std::cout << "count , count , iter   , iter , millisec ,   millisec  , millisec      , millisec    , ratio" << std::endl ;

  if ( cmd[ CMD_USE_FIXTURE_BEGIN ] ) {
    for ( int i = cmd[CMD_USE_FIXTURE_BEGIN] ; i < cmd[CMD_USE_FIXTURE_END] * 2 ; i *= 2 ) {
      int nelem[3] ;
      nelem[0] = std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
      nelem[1] = 1 + nelem[0] ;
      nelem[2] = 2 * nelem[0] ;

      const Kokkos::Example::FENL::Perf perf =
        cmd[ CMD_USE_FIXTURE_QUADRATIC ]
        ? Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemQuadratic >
            ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
        : Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemLinear >
            ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
        ;

      print_perf_value( std::cout , perf );
    }
  }
  else {
    int nelem[3] = { cmd[ CMD_USE_FIXTURE_X ] ,
                     cmd[ CMD_USE_FIXTURE_Y ] ,
                     cmd[ CMD_USE_FIXTURE_Z ] };

    const Kokkos::Example::FENL::Perf perf =
      cmd[ CMD_USE_FIXTURE_QUADRATIC ]
      ? Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemQuadratic >
          ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
      : Kokkos::Example::FENL::fenl< Device , Kokkos::Example::BoxElemPart::ElemLinear >
          ( comm , cmd[CMD_PRINT], cmd[CMD_USE_TRIALS], cmd[CMD_USE_ATOMIC], nelem )
      ;

    print_perf_value( std::cout , perf );
  }
}

//----------------------------------------------------------------------------

int main( int argc , char ** argv )
{
  int comm_rank = 0 ;
  int comm_size = 1 ;

#if defined( KOKKOS_HAVE_MPI )
  MPI_Init( & argc , & argv );
  MPI_Comm comm = MPI_COMM_WORLD ;
  MPI_Comm_rank( comm , & comm_rank );
  MPI_Comm_size( comm , & comm_size );
#else
  MPI_Comm comm = 0 ;
#endif

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

#if defined( KOKKOS_HAVE_MPI )
  MPI_Bcast( cmdline , CMD_COUNT , MPI_INT , 0 , comm );
#endif

  if ( ! cmdline[ CMD_ERROR ] && ! cmdline[ CMD_ECHO ] ) {

    if ( ! cmdline[ CMD_USE_TRIALS ] ) { cmdline[ CMD_USE_TRIALS ] = 1 ; }

    if ( ! cmdline[ CMD_USE_FIXTURE_X ] && ! cmdline[ CMD_USE_FIXTURE_BEGIN ] ) {
      cmdline[ CMD_USE_FIXTURE_X ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Y ] = 2 ;
      cmdline[ CMD_USE_FIXTURE_Z ] = 2 ;
    }

#if defined( KOKKOS_HAVE_PTHREAD )

    if ( cmdline[ CMD_USE_THREADS ] ) {

      if ( cmdline[ CMD_USE_NUMA ] && cmdline[ CMD_USE_CORE_PER_NUMA ] ) {
        Kokkos::Threads::initialize( cmdline[ CMD_USE_THREADS ] ,
                                     cmdline[ CMD_USE_NUMA ] ,
                                     cmdline[ CMD_USE_CORE_PER_NUMA ] );
      }
      else {
        Kokkos::Threads::initialize( cmdline[ CMD_USE_THREADS ] );
      }

      run< Kokkos::Threads , Kokkos::Example::BoxElemPart::ElemLinear >( comm , cmdline );

      Kokkos::Threads::finalize();
    }

#endif

#if defined( KOKKOS_HAVE_OPENMP )

    if ( cmdline[ CMD_USE_OPENMP ] ) {

      if ( cmdline[ CMD_USE_NUMA ] && cmdline[ CMD_USE_CORE_PER_NUMA ] ) {
        Kokkos::OpenMP::initialize( cmdline[ CMD_USE_OPENMP ] ,
                                     cmdline[ CMD_USE_NUMA ] ,
                                     cmdline[ CMD_USE_CORE_PER_NUMA ] );
      }
      else {
        Kokkos::OpenMP::initialize( cmdline[ CMD_USE_THREADS ] );
      }

      run< Kokkos::OpenMP , Kokkos::Example::BoxElemPart::ElemLinear >( comm , cmdline );

      Kokkos::OpenMP::finalize();
    }

#endif

#if defined( KOKKOS_HAVE_CUDA )
    if ( cmdline[ CMD_USE_CUDA ] ) {
      // Use the last device:

      Kokkos::Cuda::host_mirror_device_type::initialize();
      Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice( cmdline[ CMD_USE_CUDA_DEV ] ) );

      run< Kokkos::Cuda , Kokkos::Example::BoxElemPart::ElemLinear >( comm , cmdline );

      Kokkos::Cuda::finalize();
      Kokkos::Cuda::host_mirror_device_type::finalize();
    }

#endif

  }

#if defined( KOKKOS_HAVE_MPI )
  MPI_Finalize();
#endif

  return cmdline[ CMD_ERROR ] ? -1 : 0 ;
}


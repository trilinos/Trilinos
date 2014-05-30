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
#include <fenl_utils.hpp>

//----------------------------------------------------------------------------

#include <Tpetra_Version.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

//----------------------------------------------------------------------------

template< class Device , Kokkos::Example::BoxElemPart::ElemOrder ElemOrder >
bool run( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ,
          const int cmd[] )
{
  typedef typename Kokkos::Compat::KokkosDeviceWrapperNode<Device> NodeType;
  bool success = true;
  try {

  // Print output headers
  const int comm_rank = comm->getRank();
  const std::vector< size_t > widths =
    print_headers( std::cout , cmd , comm_rank );

  using Kokkos::Example::FENL::ManufacturedSolution;
  using Kokkos::Example::FENL::ElementComputationConstantCoefficient;
  using Kokkos::Example::BoxElemPart;
  using Kokkos::Example::FENL::fenl;
  using Kokkos::Example::FENL::Perf;

  const double bc_lower_value = 1 ;
  const double bc_upper_value = 2 ;
  const ManufacturedSolution
    manufactured_solution( 0 , 1 , bc_lower_value , bc_upper_value  );

  ElementComputationConstantCoefficient
    diffusion_coefficient( manufactured_solution.K );

  // Create Tpetra Node
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
  Teuchos::RCP<NodeType> node = Teuchos::rcp (new NodeType(params));

  double response = 0;
  if ( cmd[ CMD_USE_FIXTURE_BEGIN ] ) {
    for ( int i = cmd[CMD_USE_FIXTURE_BEGIN] ; i < cmd[CMD_USE_FIXTURE_END] * 2 ; i *= 2 ) {
      int nelem[3] ;
      nelem[0] = std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
      nelem[1] = 1 + nelem[0] ;
      nelem[2] = 2 * nelem[0] ;

      Perf perf;
      if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] )
        perf = fenl< double , Device , BoxElemPart::ElemQuadratic >
          ( comm , node , cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] ,
            cmd[CMD_USE_ATOMIC] , cmd[CMD_USE_BELOS] , cmd[CMD_USE_MUELU] ,
            nelem , diffusion_coefficient , manufactured_solution ,
            manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
            true , response);
      else
        perf = fenl< double , Device , BoxElemPart::ElemLinear >
          ( comm , node , cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] ,
            cmd[CMD_USE_ATOMIC] , cmd[CMD_USE_BELOS] , cmd[CMD_USE_MUELU] ,
            nelem , diffusion_coefficient , manufactured_solution ,
            manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
            true , response);

      perf.response_mean = response;
      perf.response_std_dev = 0.0;

      if ( 0 == comm_rank ) { print_perf_value( std::cout , widths, perf ); }
    }
  }
  else {
    int nelem[3] = { cmd[ CMD_USE_FIXTURE_X ] ,
                     cmd[ CMD_USE_FIXTURE_Y ] ,
                     cmd[ CMD_USE_FIXTURE_Z ] };

    Perf perf;
    if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] )
      perf = fenl< double , Device , BoxElemPart::ElemQuadratic >
        ( comm , node , cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] ,
          cmd[CMD_USE_ATOMIC] , cmd[CMD_USE_BELOS] , cmd[CMD_USE_MUELU] ,
          nelem , diffusion_coefficient , manufactured_solution ,
          manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
          true , response);
    else
      perf = fenl< double , Device , BoxElemPart::ElemLinear >
        ( comm , node , cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] ,
          cmd[CMD_USE_ATOMIC] , cmd[CMD_USE_BELOS] , cmd[CMD_USE_MUELU] ,
          nelem , diffusion_coefficient , manufactured_solution ,
          manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
          true , response);

    perf.response_mean = response;
    perf.response_std_dev = 0.0;

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

  //--------------------------------------------------------------------------

  int cmdline[ CMD_COUNT ] ;
  clp_return_type rv = parse_cmdline( argc, argv, cmdline, *comm );
  if (rv==CLP_HELP)
    return(EXIT_SUCCESS);
  else if (rv==CLP_ERROR)
    return(EXIT_FAILURE);

  if ( ! cmdline[ CMD_ERROR ] && ! cmdline[ CMD_ECHO ] ) {

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

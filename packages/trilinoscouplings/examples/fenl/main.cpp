#include <Kokkos_Core.hpp>

#include <fenl.hpp>
#include <fenl_utils.hpp>

//----------------------------------------------------------------------------

#include <Tpetra_Version.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

//----------------------------------------------------------------------------

template< class Device >
bool run( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ,
          const CMD & cmd )
{
  typedef typename Kokkos::Compat::KokkosDeviceWrapperNode<Device> NodeType;
  bool success = true;
  try {

  const int comm_rank = comm->getRank();

  // Create Tpetra Node
  Teuchos::RCP<NodeType> node = createKokkosNode<NodeType>( cmd , *comm );

  // Print output headers
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

  double response = 0;
  if ( cmd.CMD_USE_FIXTURE_BEGIN  ) {
    for ( int i = cmd.CMD_USE_FIXTURE_BEGIN ; i < cmd.CMD_USE_FIXTURE_END * 2 ; i *= 2 ) {
      int nelem[3] ;
      nelem[0] = std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
      nelem[1] = 1 + nelem[0] ;
      nelem[2] = 2 * nelem[0] ;

      Perf perf;
      if ( cmd.CMD_USE_FIXTURE_QUADRATIC  )
        perf = fenl< double , Device , BoxElemPart::ElemQuadratic >
          ( comm , node , cmd.CMD_PRINT , cmd.CMD_USE_TRIALS ,
            cmd.CMD_USE_ATOMIC , cmd.CMD_USE_BELOS , cmd.CMD_USE_MUELU , cmd.CMD_USE_MEANBASED , 
            nelem , diffusion_coefficient , manufactured_solution ,
            manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
            true , response);
      else
        perf = fenl< double , Device , BoxElemPart::ElemLinear >
          ( comm , node , cmd.CMD_PRINT , cmd.CMD_USE_TRIALS ,
            cmd.CMD_USE_ATOMIC , cmd.CMD_USE_BELOS , cmd.CMD_USE_MUELU , cmd.CMD_USE_MEANBASED ,
            nelem , diffusion_coefficient , manufactured_solution ,
            manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
            true , response);

      perf.response_mean = response;
      perf.response_std_dev = 0.0;

      if ( 0 == comm_rank ) {
        print_perf_value( std::cout , cmd , widths , perf );
      }
    }
  }
  else {
    int nelem[3] = { cmd.CMD_USE_FIXTURE_X  ,
                     cmd.CMD_USE_FIXTURE_Y  ,
                     cmd.CMD_USE_FIXTURE_Z  };

    Perf perf;
    if ( cmd.CMD_USE_FIXTURE_QUADRATIC  )
      perf = fenl< double , Device , BoxElemPart::ElemQuadratic >
        ( comm , node , cmd.CMD_PRINT , cmd.CMD_USE_TRIALS ,
          cmd.CMD_USE_ATOMIC , cmd.CMD_USE_BELOS , cmd.CMD_USE_MUELU ,
          cmd.CMD_USE_MEANBASED ,
          nelem , diffusion_coefficient , manufactured_solution ,
          manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
          true , response);
    else
      perf = fenl< double , Device , BoxElemPart::ElemLinear >
        ( comm , node , cmd.CMD_PRINT , cmd.CMD_USE_TRIALS ,
          cmd.CMD_USE_ATOMIC , cmd.CMD_USE_BELOS , cmd.CMD_USE_MUELU ,
          cmd.CMD_USE_MEANBASED ,
          nelem , diffusion_coefficient , manufactured_solution ,
          manufactured_solution.T_zmin , manufactured_solution.T_zmax ,
          true , response);

    perf.response_mean = response;
    perf.response_std_dev = 0.0;

    if ( 0 == comm_rank ) {
      print_perf_value( std::cout , cmd , widths , perf );
    }
  }

  if ( cmd.CMD_SUMMARIZE  ) {
    Teuchos::TimeMonitor::report (comm.ptr (), std::cout);
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
    Teuchos::DefaultComm<int>::getComm();

  //--------------------------------------------------------------------------
  CMD cmdline;
  clp_return_type rv = parse_cmdline( argc, argv, cmdline, *comm, false );
  if (rv==CLP_HELP)
    return(EXIT_SUCCESS);
  else if (rv==CLP_ERROR)
    return(EXIT_FAILURE);

  if ( cmdline.CMD_VTUNE  ) {
    connect_vtune(comm->getRank());
  }

  if ( ! cmdline.CMD_ERROR  && ! cmdline.CMD_ECHO  ) {

#if defined( HAVE_TPETRA_INST_PTHREAD )
    if ( cmdline.CMD_USE_THREADS ) {
      run< Kokkos::Threads >( comm , cmdline );
    }
#endif

#if defined( HAVE_TPETRA_INST_OPENMP )
    if ( cmdline.CMD_USE_OPENMP ) {
      run< Kokkos::OpenMP >( comm , cmdline );
    }
#endif

#if defined( HAVE_TPETRA_INST_CUDA )
    if ( cmdline.CMD_USE_CUDA ) {
      run< Kokkos::Cuda >( comm , cmdline );
    }
#endif

  }

  //--------------------------------------------------------------------------

  return cmdline.CMD_ERROR  ? -1 : 0 ;
}

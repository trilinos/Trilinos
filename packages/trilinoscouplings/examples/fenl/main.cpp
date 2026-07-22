// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Kokkos_Core.hpp>

#include <fenl.hpp>
#include <fenl_impl.hpp>
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
  bool success = true;
  try {

  const int comm_rank = comm->getRank();

  // Print output headers
  const std::vector< size_t > widths =
    print_headers( std::cout , cmd , comm_rank );

  using Kokkos::Example::FENL::ElementComputationLinearCoefficient;
  using Kokkos::Example::BoxElemPart;
  using Kokkos::Example::FENL::fenl;
  using Kokkos::Example::FENL::Perf;

  const double bc_lower_value = 1 ;
  const double bc_upper_value = 2 ;

  ElementComputationLinearCoefficient
    linear_diffusion_coefficient( cmd.USE_DIFF_COEFF_LINEAR ,
                                  cmd.USE_DIFF_COEFF_CONSTANT );

  int nelem[3] = { cmd.USE_FIXTURE_X  ,
                   cmd.USE_FIXTURE_Y  ,
                   cmd.USE_FIXTURE_Z  };

  Perf perf;
  double response = 0;
  Teuchos::Array<double> response_gradient;
  if ( cmd.USE_FIXTURE_QUADRATIC  ) {
    perf = fenl< double , Device , BoxElemPart::ElemQuadratic >
      ( comm , cmd.USE_FENL_XML_FILE ,
        cmd.PRINT , cmd.USE_TRIALS ,
        cmd.USE_ATOMIC , cmd.USE_BELOS , cmd.USE_MUELU ,
        cmd.USE_MEANBASED ,
        nelem , linear_diffusion_coefficient, cmd.USE_ISOTROPIC , cmd.USE_COEFF_SRC ,
        cmd.USE_COEFF_ADV , bc_lower_value , bc_upper_value ,
        response , response_gradient );
  }
  else {
    perf = fenl< double , Device , BoxElemPart::ElemLinear >
      ( comm , cmd.USE_FENL_XML_FILE ,
        cmd.PRINT , cmd.USE_TRIALS ,
        cmd.USE_ATOMIC , cmd.USE_BELOS , cmd.USE_MUELU ,
        cmd.USE_MEANBASED ,
        nelem , linear_diffusion_coefficient, cmd.USE_ISOTROPIC , cmd.USE_COEFF_SRC ,
        cmd.USE_COEFF_ADV , bc_lower_value , bc_upper_value ,
        response , response_gradient );
  }

  perf.response_mean = response;
  perf.response_std_dev = 0.0;

  if ( 0 == comm_rank ) {
    print_perf_value( std::cout , cmd , widths , perf );
  }

  if ( cmd.SUMMARIZE  ) {
    Teuchos::TimeMonitor::report (comm.ptr (), std::cout);
    print_memory_usage(std::cout, *comm);
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

  {
  Kokkos::initialize(argc, argv);

  if ( cmdline.VTUNE  ) {
    connect_vtune(comm->getRank());
  }

  if ( ! cmdline.ERROR  && ! cmdline.ECHO  ) {

#if defined( HAVE_TPETRA_SERIAL )
    if ( cmdline.USE_SERIAL ) {
      run< Kokkos::Serial >( comm , cmdline );
    }
#endif

#if defined( HAVE_TPETRA_PTHREAD )
    if ( cmdline.USE_THREADS ) {
      run< Kokkos::Threads >( comm , cmdline );
    }
#endif

#if defined( HAVE_TPETRA_OPENMP )
    if ( cmdline.USE_OPENMP ) {
      run< Kokkos::OpenMP >( comm , cmdline );
    }
#endif

#if defined( HAVE_TPETRA_CUDA )
    if ( cmdline.USE_CUDA ) {
      run< Kokkos::Cuda >( comm , cmdline );
    }
#endif

  }

  }
  Kokkos::finalize();

  //--------------------------------------------------------------------------

  return cmdline.ERROR  ? -1 : 0 ;
}

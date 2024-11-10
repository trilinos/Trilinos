// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_View_MP_Vector_Fwd.hpp"
#include "Kokkos_View_UQ_PCE_Fwd.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Stokhos.hpp"

//----------------------------------------------------------------------------

#include <Kokkos_Core.hpp>

#include <fenl_pce.hpp>
#include <fenl_utils.hpp>

//----------------------------------------------------------------------------

#include <Tpetra_Version.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_MultiVector.hpp>

//----------------------------------------------------------------------------

// For unit-testing
#include "Teuchos_TestingHelpers.hpp"

template< class Device >
bool run( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ,
          const CMD & cmd)
{
  bool success = true;
  try {

  const int comm_rank = comm->getRank();

  // Set up stochastic discretization
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using one_d_basis    = Stokhos::OneDOrthogPolyBasis<int,double>;
  using legendre_basis = Stokhos::LegendreBasis<int,double>;
  using order_type     = Stokhos::LexographicLess< Stokhos::MultiIndex<int> >;
  using product_basis  = Stokhos::TotalOrderBasis<int,double,order_type>;
  using Cijk           = Stokhos::Sparse3Tensor<int,double>;
  using  quadrature    = Stokhos::Quadrature<int,double>;
  const int dim = cmd.USE_UQ_DIM;
  const int order = cmd.USE_UQ_ORDER ;
  Array< RCP<const one_d_basis> > bases(dim);
  for (int i=0; i<dim; i++)
    bases[i] = rcp(new legendre_basis(order, true));
  RCP<const product_basis> basis = rcp(new product_basis(bases));
  RCP<Cijk> cijk = basis->computeTripleProductTensor();

  using Storage = Stokhos::DynamicStorage<int,double,Device>;
  using Scalar  = Sacado::UQ::PCE<Storage>;
  typename Scalar::cijk_type kokkos_cijk =
    Stokhos::create_product_tensor<Device>(*basis, *cijk);
  Kokkos::setGlobalCijkTensor(kokkos_cijk);

  // Create quadrature data used by assembly
  RCP<const quadrature> quad;
  if ( cmd.USE_SPARSE  ) {
    Stokhos::TotalOrderIndexSet<int> index_set(dim, order);
    quad =
      rcp(new Stokhos::SmolyakSparseGridQuadrature<int,double>(basis,
                                                                 index_set));
  }
  else
    quad = rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
  const int num_pce                         = basis->size();
  const int num_quad_points                 = quad->size();
  const Array<double>& quad_weights         = quad->getQuadWeights();
  const Array< Array<double> >& quad_points = quad->getQuadPoints();
  const Array< Array<double> >& quad_values = quad->getBasisAtQuadPoints();

  // Align number of quadrature points to ensemble size used in assembly
  const int align = 32;
  const int mask = align-1;
  const int num_quad_points_aligned = (num_quad_points + mask) & ~mask;

  // Copy quadrature data to view's for assembly kernels
  using QD                = Kokkos::Example::FENL::QuadratureData<Device>;
  using quad_weights_type = typename QD::quad_weights_type;
  using quad_values_type  = typename QD::quad_values_type;
  QD qd;
  qd.weights_view =
    quad_weights_type( "quad weights", num_quad_points_aligned );
  qd.points_view =
    quad_values_type( "quad points", num_quad_points_aligned, dim );
  qd.values_view =
    quad_values_type( "quad values", num_quad_points_aligned, num_pce );
  typename quad_weights_type::HostMirror host_weights_view =
    Kokkos::create_mirror_view( qd.weights_view );
  typename quad_values_type::HostMirror host_points_view =
    Kokkos::create_mirror_view( qd.points_view );
  typename quad_values_type::HostMirror host_values_view =
    Kokkos::create_mirror_view( qd.values_view );
  for (int qp=0; qp<num_quad_points; ++qp) {
    host_weights_view(qp) = quad_weights[qp];
    for (int i=0; i<dim; ++i)
      host_points_view(qp,i) = quad_points[qp][i];
    for (int i=0; i<num_pce; ++i)
      host_values_view(qp,i) = quad_values[qp][i];
  }
  for (int qp=num_quad_points; qp<num_quad_points_aligned; ++qp) {
    host_weights_view(qp) = 0.0;
    for (int i=0; i<dim; ++i)
      host_points_view(qp,i) = quad_points[num_quad_points-1][i];
    for (int i=0; i<num_pce; ++i)
      host_values_view(qp,i) = quad_values[num_quad_points-1][i];
  }
  Kokkos::deep_copy( qd.weights_view, host_weights_view );
  Kokkos::deep_copy( qd.points_view, host_points_view );
  Kokkos::deep_copy( qd.values_view, host_values_view );

  // Print output headers
  const std::vector< size_t > widths =
    print_headers( std::cout , cmd , comm_rank );

  using Kokkos::Example::FENL::ElementComputationKLCoefficient;
  using Kokkos::Example::FENL::ExponentialKLCoefficient;
  using Kokkos::Example::BoxElemPart;
  using Kokkos::Example::FENL::fenl;
  using Kokkos::Example::FENL::Perf;

  const double bc_lower_value = 1 ;
  const double bc_upper_value = 2 ;

  int nelem[3] = { cmd.USE_FIXTURE_X  ,
                   cmd.USE_FIXTURE_Y  ,
                   cmd.USE_FIXTURE_Z  };

  // Create KL diffusion coefficient
  const double kl_mean = cmd.USE_MEAN;
  const double kl_variance = cmd.USE_VAR;
  const double kl_correlation = cmd.USE_COR;
  const bool kl_exp = cmd.USE_EXPONENTIAL;
  const double kl_exp_shift = cmd.USE_EXP_SHIFT;
  const double kl_exp_scale = cmd.USE_EXP_SCALE;
  const bool kl_disc_exp_scale = cmd.USE_DISC_EXP_SCALE;
  //typedef ElementComputationKLCoefficient< Scalar, double, Device > KL;
  typedef ExponentialKLCoefficient< Scalar, double, Device > KL;
  KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, dim,
                            kl_exp, kl_exp_shift, kl_exp_scale,
                            kl_disc_exp_scale );
  typedef typename KL::RandomVariableView RV;
  typedef typename RV::HostMirror HRV;
  RV rv("KL Random Variables", dim);
  HRV hrv = Kokkos::create_mirror_view(rv);

  // Set random variables
  // ith random variable \xi_i = \psi_I(\xi) / \psi_I(1.0)
  // where I is determined by the basis ordering (since the component basis
  // functions have unit two-norm, \psi_I(1.0) might not be 1.0).  We compute
  // this by finding the index of the multivariate term that is first order in
  // the ith slot, all other orders 0
  Teuchos::Array<double> point(dim, 1.0);
  Teuchos::Array<double> basis_vals(num_pce);
  basis->evaluateBases(point, basis_vals);
  for (int i=0; i<dim; ++i) {
    Stokhos::MultiIndex<int> term(dim, 0);
    term[i] = 1;
    int index = basis->index(term);
    hrv(i).fastAccessCoeff(index) = 1.0 / basis_vals[index];
  }
  Kokkos::deep_copy( rv, hrv );
  diffusion_coefficient.setRandomVariables(rv);

  // Compute stochastic response using stochastic Galerkin method
  Scalar response = 0;
  Teuchos::Array<Scalar> response_gradient;
  Perf perf;
  if ( cmd.USE_FIXTURE_QUADRATIC  )
    perf = fenl< Scalar , Device , BoxElemPart::ElemQuadratic >
      ( comm , cmd.USE_FENL_XML_FILE ,
        cmd.PRINT , cmd.USE_TRIALS ,
        cmd.USE_ATOMIC , cmd.USE_BELOS , cmd.USE_MUELU ,
        cmd.USE_MEANBASED ,
        nelem , diffusion_coefficient , cmd.USE_ISOTROPIC , cmd.USE_COEFF_SRC ,
        cmd.USE_COEFF_ADV , bc_lower_value , bc_upper_value ,
        response, response_gradient, qd );
  else
    perf = fenl< Scalar , Device , BoxElemPart::ElemLinear >
      ( comm , cmd.USE_FENL_XML_FILE ,
        cmd.PRINT , cmd.USE_TRIALS ,
        cmd.USE_ATOMIC , cmd.USE_BELOS , cmd.USE_MUELU ,
        cmd.USE_MEANBASED ,
        nelem , diffusion_coefficient , cmd.USE_ISOTROPIC , cmd.USE_COEFF_SRC ,
        cmd.USE_COEFF_ADV , bc_lower_value , bc_upper_value ,
        response , response_gradient, qd );

  // std::cout << "newton count = " << perf.newton_iter_count
  //           << " cg count = " << perf.cg_iter_count << std::endl;
  perf.uq_count = num_quad_points;
  perf.newton_iter_count *= num_quad_points;
  perf.cg_iter_count *= num_pce;
  perf.map_ratio *= num_pce;
  perf.fill_node_set *= num_pce;
  perf.scan_node_count *= num_pce;
  perf.fill_graph_entries *= num_pce;
  perf.sort_graph_entries *= num_pce;
  perf.fill_element_graph *= num_pce;

  // Compute response mean, variance
  perf.response_mean = response.mean();
  perf.response_std_dev = response.standard_deviation();

  //std::cout << std::endl << response << std::endl;

  Kokkos::setGlobalCijkTensor(typename Scalar::cijk_type());

  if ( 0 == comm_rank ) {
    print_perf_value( std::cout , cmd , widths , perf );
  }

  if ( cmd.SUMMARIZE  ) {
    Teuchos::TimeMonitor::report (comm.ptr (), std::cout);
    print_memory_usage(std::cout, *comm);
  }

  // If we are running as a unit-test, check mean and variance
  if (cmd.UNIT_TEST) {
    TEUCHOS_TEST_FLOATING_EQUALITY( perf.response_mean, cmd.TEST_MEAN, cmd.TEST_TOL, std::cout, success );
    TEUCHOS_TEST_FLOATING_EQUALITY( perf.response_std_dev, cmd.TEST_STD_DEV, cmd.TEST_TOL, std::cout, success );
    if (success)
      std::cout << "Test Passed!" << std::endl;
    else
      std::cout << "Test Failed!" << std::endl;
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
  clp_return_type rv = parse_cmdline( argc, argv, cmdline, *comm, true );

  {
  Kokkos::initialize(argc, argv);

  // Print a warning if we are using the non-mean-based preconditioner
  if (cmdline.USE_MUELU && !cmdline.USE_MEANBASED &&
      comm->getRank() == 0) {
    std::cout << "Warning:  The non-mean-based preconditioner for PCE is "
              << "not correctly applying the inverse diagonal in the smoother."
              << std::endl;
  }
#ifndef HAVE_TPETRA_EXPLICIT_INSTANTIATION
  if (cmdline.USE_MUELU && cmdline.USE_MEANBASED &&
      comm->getRank() == 0) {
    std::cout << "Warning:  The mean-based preconditioner for PCE requires "
              << "specializations introduced through ETI, however it is not "
              << "enabled.  Memory errors are likely!"
              << std::endl;
  }
#endif
  if (rv==CLP_HELP)
    return(EXIT_SUCCESS);
  else if (rv==CLP_ERROR)
    return(EXIT_FAILURE);

  if ( cmdline.VTUNE  ) {
    connect_vtune(comm->getRank());
  }

  if ( ! cmdline.ERROR  && ! cmdline.ECHO  ) {

    // If execution space not specified, use the default
    if (!cmdline.USE_SERIAL && !cmdline.USE_THREADS && !cmdline.USE_OPENMP &&
        !cmdline.USE_CUDA)
      run< Kokkos::DefaultExecutionSpace >( comm , cmdline );

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

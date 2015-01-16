#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos.hpp"

//----------------------------------------------------------------------------

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
#include <Tpetra_MultiVector.hpp>

//----------------------------------------------------------------------------

template< class Device >
bool run( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ,
          const CMD & cmd)
{
  typedef typename Kokkos::Compat::KokkosDeviceWrapperNode<Device> NodeType;
  bool success = true;
  try {

  const int comm_rank = comm->getRank();

  // Create Tpetra Node -- do this first as it initializes host/device
  Teuchos::RCP<NodeType> node = createKokkosNode<NodeType>( cmd , *comm );

  // Set up stochastic discretization
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Stokhos::OneDOrthogPolyBasis<int,double> one_d_basis;
  typedef Stokhos::LegendreBasis<int,double> legendre_basis;
  typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,double,order_type> product_basis;
  typedef Stokhos::Sparse3Tensor<int,double> Cijk;
  typedef Stokhos::Quadrature<int,double> quadrature;
  const int dim = cmd.CMD_USE_UQ_DIM;
  const int order = cmd.CMD_USE_UQ_ORDER ;
  Array< RCP<const one_d_basis> > bases(dim);
  for (int i=0; i<dim; i++)
    bases[i] = rcp(new legendre_basis(order, true));
  RCP<const product_basis> basis = rcp(new product_basis(bases));
  RCP<Cijk> cijk = basis->computeTripleProductTensor();

  typedef Stokhos::DynamicStorage<int,double,Device> Storage;
  typedef Sacado::UQ::PCE<Storage> Scalar;
  typename Scalar::cijk_type kokkos_cijk =
    Stokhos::create_product_tensor<Device>(*basis, *cijk);
  Kokkos::setGlobalCijkTensor(kokkos_cijk);

  // Create quadrature data used by assembly
  RCP<const quadrature> quad;
  if ( cmd.CMD_USE_SPARSE  ) {
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
  typedef Kokkos::Example::FENL::QuadratureData<Device> QD;
  typedef typename QD::quad_weights_type quad_weights_type;
  typedef typename QD::quad_values_type quad_values_type;
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

  using Kokkos::Example::FENL::TrivialManufacturedSolution;
  using Kokkos::Example::FENL::ElementComputationKLCoefficient;
  using Kokkos::Example::BoxElemPart;
  using Kokkos::Example::FENL::fenl;
  using Kokkos::Example::FENL::Perf;

  const double bc_lower_value = 1 ;
  const double bc_upper_value = 2 ;
  const TrivialManufacturedSolution manufactured_solution;

  int nelem[3] = { cmd.CMD_USE_FIXTURE_X  ,
                   cmd.CMD_USE_FIXTURE_Y  ,
                   cmd.CMD_USE_FIXTURE_Z  };

  // Create KL diffusion coefficient
  const double kl_mean = cmd.CMD_USE_MEAN;
  const double kl_variance = cmd.CMD_USE_VAR;
  const double kl_correlation = cmd.CMD_USE_COR;
  typedef ElementComputationKLCoefficient< Scalar, double, Device > KL;
  KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, dim );
  typedef typename KL::RandomVariableView RV;
  typedef typename RV::HostMirror HRV;
  RV rv = diffusion_coefficient.getRandomVariables();
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

  // Compute stochastic response using stochastic Galerkin method
  Scalar response = 0;
  Perf perf;
  if ( cmd.CMD_USE_FIXTURE_QUADRATIC  )
    perf = fenl< Scalar , Device , BoxElemPart::ElemQuadratic >
      ( comm , node , cmd.CMD_PRINT , cmd.CMD_USE_TRIALS ,
        cmd.CMD_USE_ATOMIC , cmd.CMD_USE_BELOS , cmd.CMD_USE_MUELU ,
        cmd.CMD_USE_MEANBASED ,
        nelem , diffusion_coefficient , manufactured_solution ,
        bc_lower_value , bc_upper_value ,
        false , response, qd );
  else
    perf = fenl< Scalar , Device , BoxElemPart::ElemLinear >
      ( comm , node , cmd.CMD_PRINT , cmd.CMD_USE_TRIALS ,
        cmd.CMD_USE_ATOMIC , cmd.CMD_USE_BELOS , cmd.CMD_USE_MUELU ,
        cmd.CMD_USE_MEANBASED ,
        nelem , diffusion_coefficient , manufactured_solution ,
        bc_lower_value , bc_upper_value ,
        false , response , qd );

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

  if ( 0 == comm_rank ) {
    print_perf_value( std::cout , cmd , widths , perf );
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
  clp_return_type rv = parse_cmdline( argc, argv, cmdline, *comm, true );

  // Print a warning if we are using the non-mean-based preconditioner
  if (cmdline.CMD_USE_MUELU && !cmdline.CMD_USE_MEANBASED &&
      comm->getRank() == 0) {
    std::cout << "Warning:  The non-mean-based preconditioner for PCE is "
              << "not correctly applying the inverse diagonal in the smoother."
              << std::endl;
  }

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

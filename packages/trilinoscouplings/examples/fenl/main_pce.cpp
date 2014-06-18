#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos.hpp"

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
#include <fenl_utils.hpp>

//----------------------------------------------------------------------------

#include <Tpetra_Version.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_Comm.hpp>
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
  Teuchos::RCP<NodeType> node = createKokkosNode<NodeType>( cmd , comm_rank );

  // Set up stochastic discretization
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Stokhos::OneDOrthogPolyBasis<int,double> one_d_basis;
  typedef Stokhos::LegendreBasis<int,double> legendre_basis;
  typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,double,order_type> product_basis;
  typedef Stokhos::Sparse3Tensor<int,double> Cijk;
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

  // typedef Stokhos::TensorProductQuadrature<int,double> quadrature;
  // RCP<const quadrature> quad     = rcp(new quadrature(basis));
  // const int num_quad_points                 = quad->size();
  // const Array<double>& quad_weights         = quad->getQuadWeights();
  // const Array< Array<double> >& quad_points = quad->getQuadPoints();
  // const Array< Array<double> >& quad_values = quad->getBasisAtQuadPoints();

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
  const double kl_correlation = cmd.CMD_USE_VAR;
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
  Teuchos::Array<double> basis_vals(basis->size());
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
        false , response);
  else
    perf = fenl< Scalar , Device , BoxElemPart::ElemLinear >
      ( comm , node , cmd.CMD_PRINT , cmd.CMD_USE_TRIALS ,
        cmd.CMD_USE_ATOMIC , cmd.CMD_USE_BELOS , cmd.CMD_USE_MUELU ,
        cmd.CMD_USE_MEANBASED ,
        nelem , diffusion_coefficient , manufactured_solution ,
        bc_lower_value , bc_upper_value ,
        false , response);

  // std::cout << "newton count = " << perf.newton_iter_count
  //           << " cg count = " << perf.cg_iter_count << std::endl;
  int pce_size = basis->size();
  perf.uq_count = pce_size;
  perf.newton_iter_count *= pce_size;
  perf.cg_iter_count *= pce_size;
  perf.map_ratio *= pce_size;
  perf.fill_node_set *= pce_size;
  perf.scan_node_count *= pce_size;
  perf.fill_graph_entries *= pce_size;
  perf.sort_graph_entries *= pce_size;
  perf.fill_element_graph *= pce_size;

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
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  //--------------------------------------------------------------------------
  CMD cmdline; 
  parse_cmdline( argc, argv, cmdline, *comm );
  if ( ! cmdline.CMD_USE_UQ_DIM  ) cmdline.CMD_USE_UQ_DIM  = 3 ;
  if ( ! cmdline.CMD_USE_UQ_ORDER  ) cmdline.CMD_USE_UQ_ORDER  = 2 ;


  if ( cmdline.CMD_VTUNE  ) {
    connect_vtune(comm->getRank());
  }

  if ( ! cmdline.CMD_ERROR  && ! cmdline.CMD_ECHO  ) {

#if defined( KOKKOS_HAVE_PTHREAD )
    if ( cmdline.CMD_USE_THREADS ) {
      run< Kokkos::Threads >( comm , cmdline );
    }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
    if ( cmdline.CMD_USE_OPENMP ) {
      run< Kokkos::OpenMP >( comm , cmdline );
    }
#endif

#if defined( KOKKOS_HAVE_CUDA )
    if ( cmdline.CMD_USE_CUDA ) {
      run< Kokkos::Cuda >( comm , cmdline );
    }
#endif

  }

  //--------------------------------------------------------------------------

  return cmdline.CMD_ERROR  ? -1 : 0 ;
}

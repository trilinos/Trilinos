#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos.hpp"

//----------------------------------------------------------------------------

#include <Kokkos_Core.hpp>

#include <fenl_ensemble.hpp>
#include <fenl_utils.hpp>

//----------------------------------------------------------------------------

#include <Tpetra_Version.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

//----------------------------------------------------------------------------

template< class Device , int VectorSize >
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
  typedef Stokhos::Quadrature<int,double> quadrature;
  const int dim = cmd.USE_UQ_DIM ;
  const int order = cmd.USE_UQ_ORDER ;
  Array< RCP<const one_d_basis> > bases(dim);
  for (int i=0; i<dim; i++)
    bases[i] = rcp(new legendre_basis(order, true));
  RCP<const product_basis> basis = rcp(new product_basis(bases));
  RCP<const quadrature> quad;
  int num_quad_points;
  Array<double> quad_weights;
  Array< Array<double> > quad_points;
  Array< Array<double> > quad_values;
  if ( cmd.USE_UQ_FAKE > 0 ) {
    // Create fake UQ problem of size cmd.USE_UQ_FAKE, initializing
    // points, weights, values to 0
    num_quad_points = cmd.USE_UQ_FAKE;
    quad_weights.resize(num_quad_points);
    quad_points.resize(num_quad_points);
    quad_values.resize(num_quad_points);
    for (int i=0; i<num_quad_points; ++i) {
      quad_points[i].resize(dim);
      quad_values[i].resize(basis->size());
    }
  }
  else {
    if ( cmd.USE_SPARSE  ) {
      Stokhos::TotalOrderIndexSet<int> index_set(dim, order);
      quad =
        rcp(new Stokhos::SmolyakSparseGridQuadrature<int,double>(basis,
                                                                 index_set));
    }
    else
      quad = rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    num_quad_points = quad->size();
    quad_weights    = quad->getQuadWeights();
    quad_points     = quad->getQuadPoints();
    quad_values     = quad->getBasisAtQuadPoints();
  }

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

  const double kl_mean = cmd.USE_MEAN;
  const double kl_variance = cmd.USE_VAR;
  const double kl_correlation = cmd.USE_COR;

  int nelem[3] = { cmd.USE_FIXTURE_X,
                   cmd.USE_FIXTURE_Y,
                   cmd.USE_FIXTURE_Z};
  Perf perf_total;
  perf_total.uq_count = num_quad_points;

  typedef Stokhos::DynamicStorage<int,double,Device> PCEStorage;
  typedef Sacado::UQ::PCE<PCEStorage> PCEScalar;
  PCEScalar response_pce( typename PCEScalar::cijk_type(), basis->size() );

  // Compute PCE of response propagating blocks of quadrature
  // points at a time
  if ( cmd.USE_UQ_ENSEMBLE > 0 ) {

    typedef Stokhos::StaticFixedStorage<int,double,VectorSize,Device> Storage;
    typedef Sacado::MP::Vector<Storage> Scalar;

    // Set global vector size -- this is mandatory
    Kokkos::global_sacado_mp_vector_size = VectorSize;

    //typedef ElementComputationKLCoefficient< Scalar, double, Device > KL;
    typedef ExponentialKLCoefficient< Scalar, double, Device > KL;
    KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, dim );
    typedef typename KL::RandomVariableView RV;
    typedef typename RV::HostMirror HRV;
    RV rv = diffusion_coefficient.getRandomVariables();
    HRV hrv = Kokkos::create_mirror_view(rv);

    const int num_qp_blocks = ( num_quad_points + VectorSize-1 ) / VectorSize;

    // Loop over quadrature points
    for (int qp_block=0; qp_block<num_qp_blocks; ++qp_block) {
      const int qp_begin = qp_block * VectorSize;
      const int qp_end = qp_begin + VectorSize <= num_quad_points ?
        qp_begin+VectorSize : num_quad_points;

      // Set random variables
      for (int qp=qp_begin, j=0; qp<qp_end; ++qp, ++j)
        for (int i=0; i<dim; ++i)
          hrv(i).fastAccessCoeff(j) = quad_points[qp][i];
      if (qp_end - qp_begin < VectorSize)
        for (int j=qp_end-qp_begin; j<VectorSize; ++j)
          for (int i=0; i<dim; ++i)
            hrv(i).fastAccessCoeff(j) = quad_points[qp_end-1][i];
      Kokkos::deep_copy( rv, hrv );

      // Evaluate response on qp block
      Scalar response = 0;
      Perf perf =
        fenl< Scalar , Device , BoxElemPart::ElemLinear >
        ( comm , node , cmd.USE_FENL_XML_FILE ,
          cmd.PRINT , cmd.USE_TRIALS ,
          cmd.USE_ATOMIC , cmd.USE_BELOS , cmd.USE_MUELU ,
          cmd.USE_MEANBASED ,
          nelem , diffusion_coefficient , cmd.USE_COEFF_SRC ,
          cmd.USE_COEFF_ADV , bc_lower_value , bc_upper_value ,
          response);

      perf.newton_iter_count *= VectorSize;
      perf.cg_iter_count *= VectorSize;
      perf.map_ratio *= VectorSize;
      perf.fill_node_set *= VectorSize;
      perf.scan_node_count *= VectorSize;
      perf.fill_graph_entries *= VectorSize;
      perf.sort_graph_entries *= VectorSize;
      perf.fill_element_graph *= VectorSize;
      perf_total.increment(perf, !cmd.USE_BELOS);

      // Sum response into integral computing response PCE coefficients
      for (int qp=qp_begin, j=0; qp<qp_end; ++qp, ++j) {
        double r = response.coeff(j);
        double w = quad_weights[qp];
        for (int i=0; i<basis->size(); ++i)
          response_pce.fastAccessCoeff(i) += r*w*quad_values[qp][i];
      }

    }

  }

  // Compute PCE of response propagating one quadrature point at a time
  else {

    typedef double Scalar;
    //typedef ElementComputationKLCoefficient< Scalar, double, Device > KL;
    typedef ExponentialKLCoefficient< Scalar, double, Device > KL;
    KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, dim );
    typedef typename KL::RandomVariableView RV;
    typedef typename RV::HostMirror HRV;
    RV rv = diffusion_coefficient.getRandomVariables();
    HRV hrv = Kokkos::create_mirror_view(rv);

    // Loop over quadrature points
    for (int qp=0; qp<num_quad_points; ++qp) {

      // Set random variables
      for (int i=0; i<dim; ++i)
        hrv(i) = quad_points[qp][i];
      Kokkos::deep_copy( rv, hrv );

      // Evaluate response on qp block
      Scalar response = 0;
      Perf perf =
        fenl< Scalar , Device , BoxElemPart::ElemLinear >
        ( comm , node , cmd.USE_FENL_XML_FILE ,
          cmd.PRINT , cmd.USE_TRIALS ,
          cmd.USE_ATOMIC , cmd.USE_BELOS , cmd.USE_MUELU ,
          cmd.USE_MEANBASED ,
          nelem , diffusion_coefficient , cmd.USE_COEFF_SRC ,
          cmd.USE_COEFF_ADV , bc_lower_value , bc_upper_value ,
          response);

      perf_total.increment(perf, !cmd.USE_BELOS);

      // Sum response into integral computing response PCE coefficients
      double r = response;
      double w = quad_weights[qp];
      for (int i=0; i<basis->size(); ++i)
        response_pce.fastAccessCoeff(i) += r*w*quad_values[qp][i];

    }

  }

  //std::cout << std::endl << response_pce << std::endl;

  perf_total.response_mean = response_pce.mean();
  perf_total.response_std_dev = response_pce.standard_deviation();

  if ( 0 == comm_rank ) {
    print_perf_value( std::cout , cmd , widths , perf_total );
  }

  if ( cmd.SUMMARIZE  ) {
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
  if (rv==CLP_HELP)
    return(EXIT_SUCCESS);
  else if (rv==CLP_ERROR)
    return(EXIT_FAILURE);

  if ( cmdline.VTUNE  ) {
    connect_vtune(comm->getRank());
  }

  if ( ! cmdline.ERROR  && ! cmdline.ECHO  ) {

#if defined( HAVE_TPETRA_PTHREAD )
    if ( cmdline.USE_THREADS ) {
#if defined(__MIC__)
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::Threads >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::Threads , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
#else
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 4 )
        run< Kokkos::Threads ,  4 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::Threads , 16 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::Threads , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
#endif
    }
#endif

#if defined( HAVE_TPETRA_OPENMP )
    if ( cmdline.USE_OPENMP ) {
#if defined(__MIC__)
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::OpenMP , 16 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::OpenMP , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
#else
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 4 )
        run< Kokkos::OpenMP ,  4 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::OpenMP , 16 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::OpenMP , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
#endif
    }
#endif

#if defined( HAVE_TPETRA_CUDA )
    if ( cmdline.USE_CUDA ) {
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::Cuda , 16 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::Cuda , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
    }
#endif

  }

  //--------------------------------------------------------------------------

  return cmdline.ERROR ? -1 : 0 ;
}

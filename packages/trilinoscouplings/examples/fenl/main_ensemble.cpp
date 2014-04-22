#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
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

//----------------------------------------------------------------------------

template< class Device , Kokkos::Example::BoxElemPart::ElemOrder ElemOrder >
bool run( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ,
          const int cmd[] )
{
  typedef typename Kokkos::Compat::KokkosDeviceWrapperNode<Device> NodeType;
  bool success = true;
  try {

  // Set up stochastic discretization
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Stokhos::OneDOrthogPolyBasis<int,double> one_d_basis;
  typedef Stokhos::LegendreBasis<int,double> legendre_basis;
  typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,double,order_type> product_basis;
  typedef Stokhos::TensorProductQuadrature<int,double> quadrature;
  const int dim = cmd[ CMD_USE_UQ_DIM ];
  const int order = cmd[ CMD_USE_UQ_ORDER ];
  Array< RCP<const one_d_basis> > bases(dim);
  for (int i=0; i<dim; i++)
    bases[i] = rcp(new legendre_basis(order, true));
  RCP<const product_basis> basis = rcp(new product_basis(bases));
  RCP<const quadrature> quad     = rcp(new quadrature(basis));
  const int num_quad_points                 = quad->size();
  const Array<double>& quad_weights         = quad->getQuadWeights();
  const Array< Array<double> >& quad_points = quad->getQuadPoints();
  //const Array< Array<double> >& quad_values = quad->getBasisAtQuadPoints();

  // Print output headers
  const int comm_rank = comm->getRank();
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

  const double kl_mean = 1.0;
  const double kl_variance = 0.1;
  const double kl_correlation = 0.25;

  int nelem[3] = { cmd[ CMD_USE_FIXTURE_X ] ,
                   cmd[ CMD_USE_FIXTURE_Y ] ,
                   cmd[ CMD_USE_FIXTURE_Z ] };
  Perf perf_total;
  perf_total.uq_count = num_quad_points;

  double mean = 0.0;
  double std_dev = 0.0;

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

  // Compute mean/std. dev. of response propagating blocks of quadrature
  // points at a time
  if ( cmd[ CMD_USE_UQ_ENSEMBLE ] ) {

    static const bool is_cuda =
      Kokkos::Impl::is_same<Device,Kokkos::Cuda>::value;
    static const int VectorSize = is_cuda ? 16 : 4;
    typedef Stokhos::StaticFixedStorage<int,double,VectorSize,Device> Storage;
    typedef Sacado::MP::Vector<Storage> Scalar;

    // Set global vector size -- this is mandatory
    Kokkos::global_sacado_mp_vector_size = VectorSize;

    typedef ElementComputationKLCoefficient< Scalar, double, Device > KL;
    KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, dim );
    typedef typename KL::RandomVariableView RV;
    RV rv = diffusion_coefficient.getRandomVariables();

    const int num_qp_blocks = ( num_quad_points + VectorSize-1 ) / VectorSize;

    // Loop over quadrature points
    for (int qp_block=0; qp_block<num_qp_blocks; ++qp_block) {
      const int qp_begin = qp_block * VectorSize;
      const int qp_end = qp_begin + VectorSize <= num_quad_points ?
        qp_begin+VectorSize : num_quad_points;

      // Set random variables -- using UVM here
      for (int i=0; i<dim; ++i)
        for (int qp=qp_begin, j=0; qp<qp_end; ++qp, ++j)
          rv(i).fastAccessCoeff(j) = quad_points[qp][i];

      // Evaluate response on qp block
      Scalar response = 0;
      Perf perf;
      if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] )
        perf = fenl< Scalar , Device , BoxElemPart::ElemQuadratic >
          ( comm , node ,
            cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] , cmd[CMD_USE_ATOMIC] ,
            nelem , diffusion_coefficient , manufactured_solution ,
            bc_lower_value , bc_upper_value ,
            false , response);
      else
        perf = fenl< Scalar , Device , BoxElemPart::ElemLinear >
          ( comm , node ,
            cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] , cmd[CMD_USE_ATOMIC] ,
            nelem , diffusion_coefficient , manufactured_solution ,
            bc_lower_value , bc_upper_value ,
            false , response);

      // std::cout << "newton count = " << perf.newton_iter_count
      //           << " cg count = " << perf.cg_iter_count << std::endl;
      perf.newton_iter_count *= VectorSize;
      perf.cg_iter_count *= VectorSize;
      perf.map_ratio *= VectorSize;
      perf.fill_node_set *= VectorSize;
      perf.scan_node_count *= VectorSize;
      perf.fill_graph_entries *= VectorSize;
      perf.sort_graph_entries *= VectorSize;
      perf.fill_element_graph *= VectorSize;
      perf_total.increment(perf);

      // Sum response into mean, variance
      for (int qp=qp_begin, j=0; qp<qp_end; ++qp, ++j) {
        double r = response.coeff(j);
        double w = quad_weights[qp];
        mean    += r     * w;
        std_dev += r * r * w;
      }

    }

  }

  // Compute mean/std. dev. of response propagating one quadrature point
  // at a time
  else {

    typedef double Scalar;
    typedef ElementComputationKLCoefficient< Scalar, double, Device > KL;
    KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, dim );
    typedef typename KL::RandomVariableView RV;
    RV rv = diffusion_coefficient.getRandomVariables();

    // Loop over quadrature points
    for (int qp=0; qp<num_quad_points; ++qp) {

      // Set random variables -- using UVM here
      for (int i=0; i<dim; ++i)
        rv(i) = quad_points[qp][i];

      // Evaluate response on qp block
      Scalar response = 0;
      Perf perf;
      if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] )
        perf = fenl< Scalar , Device , BoxElemPart::ElemQuadratic >
          ( comm , node ,
            cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] , cmd[CMD_USE_ATOMIC] ,
            nelem , diffusion_coefficient , manufactured_solution ,
            bc_lower_value , bc_upper_value ,
            false , response);
      else
        perf = fenl< Scalar , Device , BoxElemPart::ElemLinear >
          ( comm , node ,
            cmd[CMD_PRINT] , cmd[CMD_USE_TRIALS] , cmd[CMD_USE_ATOMIC] ,
            nelem , diffusion_coefficient , manufactured_solution ,
            bc_lower_value , bc_upper_value ,
            false , response);

      // std::cout << "newton count = " << perf.newton_iter_count
      //           << " cg count = " << perf.cg_iter_count << std::endl;
      perf_total.increment(perf);

      // Sum response into mean, variance
      double r = response;
      double w = quad_weights[qp];
      mean    += r     * w;
      std_dev += r * r * w;

    }

  }

  perf_total.response_mean = mean;
  perf_total.response_std_dev = std::sqrt(std_dev - mean*mean);

  if ( 0 == comm_rank ) { print_perf_value( std::cout , widths, perf_total ); }

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
  parse_cmdline( argc, argv, cmdline, *comm );
  if ( ! cmdline[ CMD_USE_UQ_DIM ] ) cmdline[ CMD_USE_UQ_DIM ] = 3 ;
  if ( ! cmdline[ CMD_USE_UQ_ORDER ] ) cmdline[ CMD_USE_UQ_ORDER ] = 2 ;

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

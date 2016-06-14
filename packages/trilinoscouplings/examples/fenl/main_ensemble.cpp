#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos.hpp"

//----------------------------------------------------------------------------

#include <Kokkos_Core.hpp>

#include <fenl_ensemble.hpp>
#include <fenl_utils.hpp>
#include <SampleGrouping.hpp>

//----------------------------------------------------------------------------

#include <Tpetra_Version.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#ifdef HAVE_TRILINOSCOUPLINGS_TASMANIAN
#include <TasmanianSparseGrid.hpp>
#endif

#include <fstream>

//----------------------------------------------------------------------------

template < class ProblemType, class CoeffFunctionType >
void run_samples(
  const Teuchos::Comm<int>& comm ,
  ProblemType& problem,
  const CoeffFunctionType& coeff_function,
  const Teuchos::RCP<Kokkos::Example::FENL::SampleGrouping<double> >& grouper,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const CMD & cmd ,
  const double bc_lower_value,
  const double bc_upper_value,
  const Teuchos::Array< Teuchos::Array<double> >& points,
  Teuchos::Array<double>& responses,
  Teuchos::Array<int>& iterations,
  Kokkos::Example::FENL::Perf& perf_total)
{
  typedef typename CoeffFunctionType::RandomVariableView RV;
  typedef typename RV::HostMirror HRV;
  RV rv = coeff_function.getRandomVariables();
  HRV hrv = Kokkos::create_mirror_view(rv);

  const int num_samples = points.size();
  const int dim = rv.dimension_0();;
  for (int sample=0; sample<num_samples; ++sample) {

    // Set random variable values to this sample
    for (int i=0; i<dim; ++i)
      hrv(i) = points[sample][i];
    Kokkos::deep_copy( rv, hrv );

    // Evaluate response at quadrature point
    double response = 0;
    Kokkos::Example::FENL::Perf perf =
      fenl( problem , fenlParams ,
            cmd.PRINT , cmd.USE_TRIALS , cmd.USE_ATOMIC ,
            cmd.USE_BELOS , cmd.USE_MUELU , cmd.USE_MEANBASED ,
            coeff_function , cmd.USE_ISOTROPIC ,
            cmd.USE_COEFF_SRC , cmd.USE_COEFF_ADV ,
            bc_lower_value , bc_upper_value ,
            response);

    responses[sample] = response;
    iterations[sample] = perf.cg_iter_count;

    if (cmd.PRINT_ITS && 0 == comm.getRank()) {
      std::cout << sample << " : " << perf.cg_iter_count << " ( ";
      for (int i=0; i<dim; ++i)
        std::cout << hrv(i) << " ";
      std::cout << ")" << std::endl;
    }

     // Increment timing statistics
    perf_total.increment(perf, !cmd.USE_BELOS);

  }
}

template < class Storage,
           class Device ,
           Kokkos::Example::BoxElemPart::ElemOrder ElemOrder,
           class CoeffFunctionType >
void run_samples(
  const Teuchos::Comm<int>& comm ,
  Kokkos::Example::FENL::Problem< Sacado::MP::Vector<Storage>, Device, ElemOrder>& problem ,
  const CoeffFunctionType & coeff_function,
  const Teuchos::RCP<Kokkos::Example::FENL::SampleGrouping<double> >& grouper,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const CMD & cmd ,
  const double bc_lower_value,
  const double bc_upper_value,
  const Teuchos::Array< Teuchos::Array<double> >& points,
  Teuchos::Array<double>& responses,
  Teuchos::Array<int>& iterations,
  Kokkos::Example::FENL::Perf& perf_total)
{
  using Teuchos::Array;
  using Teuchos::Ordinal;

  typedef typename Sacado::MP::Vector<Storage> Scalar;
  typedef typename CoeffFunctionType::RandomVariableView RV;
  typedef typename RV::HostMirror HRV;
  static const int VectorSize = Storage::static_size;

  // Group points into ensembles
  Array< Array<Ordinal> > groups;
  Ordinal num_duplicate = 0;
  grouper->group(VectorSize, points, groups, num_duplicate);

  const int num_groups = groups.size();
  RV rv = coeff_function.getRandomVariables();
  HRV hrv = Kokkos::create_mirror_view(rv);
  const int dim = rv.dimension_0();

  // Loop over quadrature point groups
  for (int group=0; group<num_groups; ++group) {

    // Set random variables
    for (int qp=0; qp<VectorSize; ++qp)
      for (int i=0; i<dim; ++i)
        hrv(i).fastAccessCoeff(qp) = points[groups[group][qp]][i];
    Kokkos::deep_copy( rv, hrv );

    // Evaluate response at quadrature point
    Scalar response = 0;
    Kokkos::Example::FENL::Perf perf =
      fenl( problem , fenlParams ,
            cmd.PRINT , cmd.USE_TRIALS , cmd.USE_ATOMIC ,
            cmd.USE_BELOS , cmd.USE_MUELU , cmd.USE_MEANBASED ,
            coeff_function , cmd.USE_ISOTROPIC ,
            cmd.USE_COEFF_SRC , cmd.USE_COEFF_ADV ,
            bc_lower_value , bc_upper_value ,
            response);

    // Save response -- note currently all samples within an ensemble
    // get the same number of iterations
    for (int qp=0; qp<VectorSize; ++qp) {
      responses[groups[group][qp]] = response.coeff(qp);
      iterations[groups[group][qp]] = perf.cg_iter_count;
    }

    if (cmd.PRINT_ITS && 0 == comm.getRank()) {
      std::cout << group << " : " << perf.cg_iter_count << " ( ";
      for (int qp=0; qp<VectorSize; ++qp)
        std::cout << groups[group][qp] << " ";
      std::cout << ")";
      std::cout << " ( ";
      for (int i=0; i<dim; ++i)
        std::cout << hrv(i) << " ";
      std::cout << ")" << std::endl;
    }

    // Adjust timing statistics for ensemble size
    perf.newton_iter_count *= VectorSize;
    perf.cg_iter_count *= VectorSize;
    perf.map_ratio *= VectorSize;
    perf.fill_node_set *= VectorSize;
    perf.scan_node_count *= VectorSize;
    perf.fill_graph_entries *= VectorSize;
    perf.sort_graph_entries *= VectorSize;
    perf.fill_element_graph *= VectorSize;

    // Increment timing statistics
    perf_total.increment(perf, !cmd.USE_BELOS);

  }
}

template< class ProblemType, class CoeffFunctionType >
void run_stokhos(
  const Teuchos::Comm<int>& comm ,
  ProblemType& problem ,
  const CoeffFunctionType & coeff_function,
  const Teuchos::RCP<Kokkos::Example::FENL::SampleGrouping<double> >& grouper,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const CMD & cmd ,
  const double bc_lower_value,
  const double bc_upper_value,
  Kokkos::Example::FENL::Perf& perf_total)
{
  // Set up stochastic discretization
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Ordinal;
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

  perf_total.uq_count = num_quad_points;

  typedef typename ProblemType::DeviceType Device;
  typedef Stokhos::DynamicStorage<int,double,Device> PCEStorage;
  typedef Sacado::UQ::PCE<PCEStorage> PCEScalar;
  PCEScalar response_pce( typename PCEScalar::cijk_type(), basis->size() );

  // Evaluate response at each quadrature point
  Array<double> responses(num_quad_points);
  Array<int> iterations(num_quad_points);
  run_samples(comm, problem, coeff_function, grouper,
              fenlParams, cmd,
              bc_lower_value, bc_upper_value,
              quad_points, responses, iterations, perf_total);

  // Integrate responses into PCE
  for (int qp=0; qp<num_quad_points; ++qp) {
    double r = responses[qp];
    double w = quad_weights[qp];
    for (int i=0; i<basis->size(); ++i)
      response_pce.fastAccessCoeff(i) += r*w*quad_values[qp][i];
  }

  //std::cout << std::endl << response_pce << std::endl;

  perf_total.response_mean = response_pce.mean();
  perf_total.response_std_dev = response_pce.standard_deviation();
}

template< class ProblemType, class CoeffFunctionType >
void run_tasmanian(
  const Teuchos::Comm<int>& comm ,
  ProblemType& problem ,
  const CoeffFunctionType & coeff_function,
  const Teuchos::RCP<Kokkos::Example::FENL::SampleGrouping<double> >& grouper,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const CMD & cmd ,
  const double bc_lower_value,
  const double bc_upper_value,
  Kokkos::Example::FENL::Perf& perf_total)
{
#ifdef HAVE_TRILINOSCOUPLINGS_TASMANIAN

  using Teuchos::Array;

  // Start up Tasmanian
  TasGrid::TasmanianSparseGrid sparseGrid;

  // Algorithmic parameters
  const int dim = cmd.USE_UQ_DIM;
  const int qoi = 2;
  const int initial_level = cmd.USE_UQ_INIT_LEVEL;
  const int max_level = cmd.USE_UQ_MAX_LEVEL;
  const int max_order = 1;
  const double tol = cmd.USE_UQ_TOL;
  const TasGrid::TypeOneDRule rule = TasGrid::rule_localp;
  const TasGrid::TypeRefinement refinement = TasGrid::refine_classic;
  const int qoi_to_refine = 0;

  // Create the initial grid
  sparseGrid.makeLocalPolynomialGrid(dim, qoi, initial_level, max_order, rule);
  int num_new_points = sparseGrid.getNumNeeded();

  perf_total.uq_count = num_new_points;
  int level = initial_level;
  while (num_new_points > 0 && level <= max_level) {

    if (cmd.PRINT_ITS && 0 == comm.getRank()) {
      std::cout << "Tasmanian grid level " << level
                << ", " << num_new_points << " points"
                << std::endl;
    }

    // Get the sample points
    const double *points = sparseGrid.getNeededPoints();

    // Copy points into Teuchos arrays
    Array< Array<double> > quad_points(num_new_points);
    for (int i=0; i<num_new_points; ++i) {
      quad_points[i].resize(dim);
      for (int j=0; j<dim; ++j)
        quad_points[i][j] = points[dim*i+j];
    }

    // Evaluate response on those points
    Array<double> responses(num_new_points);
    Array<int> iterations(num_new_points);
    run_samples(comm, problem, coeff_function, grouper,
                fenlParams, cmd,
                bc_lower_value, bc_upper_value,
                quad_points, responses, iterations, perf_total);

    // Load responses back into Tasmanian
    Array<double> tas_responses(qoi*num_new_points);
    for (int i=0; i<num_new_points; ++i) {
      tas_responses[i*qoi]   = responses[i];              // for mean
      tas_responses[i*qoi+1] = responses[i]*responses[i]; // for variance
    }
    sparseGrid.loadNeededPoints(&tas_responses[0]);

    // Refine the grid
    sparseGrid.setSurplusRefinement(tol, refinement, qoi_to_refine);

    // Get the number of new points
    num_new_points = sparseGrid.getNumNeeded();
    perf_total.uq_count += num_new_points;
    ++level;
  }

  if (level > max_level && comm.getRank() == 0)
    std::cout << "Warning:  Tasmanian did not achieve refinement tolerance "
              << tol << std::endl;

  // Compute mean and standard deviation of response
  double s[qoi];
  sparseGrid.integrate(s);
  const double weight = std::pow(0.5, dim); // uniform measure in dim dimensions
  s[0] *= weight; s[1] *= weight;
  perf_total.response_mean = s[0];
  perf_total.response_std_dev = std::sqrt(s[1]-s[0]*s[0]);

#else

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TASMANIAN not available.  Please re-configure with TASMANIAN TPL enabled.");

#endif
}


template< class ProblemType, class CoeffFunctionType >
void run_file(
  const Teuchos::Comm<int>& comm ,
  ProblemType& problem ,
  const CoeffFunctionType & coeff_function,
  const Teuchos::RCP<Kokkos::Example::FENL::SampleGrouping<double> >& grouper,
  const Teuchos::RCP<Teuchos::ParameterList>& fenlParams,
  const CMD & cmd ,
  const double bc_lower_value,
  const double bc_upper_value,
  Kokkos::Example::FENL::Perf& perf_total)
{
  using Teuchos::Array;

  const int dim = cmd.USE_UQ_DIM;
  int num_quad_points;
  Array< Array<double > > quad_points;

  // Open and read sample points
  std::ifstream fin("samples.txt");
  fin >> num_quad_points;
  quad_points.resize(num_quad_points);
  for (int i=0; i<num_quad_points; ++i) {
    quad_points[i].resize(dim);
    for (int j=0; j<dim; ++j)
      fin >> quad_points[i][j];
  }
  fin.close();

  // Evaluate response at each quadrature point
  Array<double> responses(num_quad_points);
  Array<int> iterations(num_quad_points);
  run_samples(comm, problem, coeff_function, grouper,
              fenlParams, cmd,
              bc_lower_value, bc_upper_value,
              quad_points, responses, iterations, perf_total);

  // Write responses to file, including solver iterations
  if (comm.getRank() == 0) {
    std::ofstream fout("responses.txt");
    fout << num_quad_points << std::endl;
    for (int i=0; i<num_quad_points; ++i) {
      fout << responses[i] << " " << iterations[i] << std::endl;
    }
    fout.close();
  }

  perf_total.response_mean = 0.0;
  perf_total.response_std_dev = 0.0;
}

template< class Device , int VectorSize >
bool run( const Teuchos::RCP<const Teuchos::Comm<int> > & comm ,
          const CMD & cmd)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef typename Kokkos::Compat::KokkosDeviceWrapperNode<Device> NodeType;
  bool success = true;
  try {

  const int comm_rank = comm->getRank();

  // Create Tpetra Node -- do this first as it initializes host/device
  Teuchos::RCP<NodeType> node = createKokkosNode<NodeType>( cmd , *comm );

  // Print output headers
  const std::vector< size_t > widths =
    print_headers( std::cout , cmd , comm_rank );

  using Kokkos::Example::FENL::ElementComputationKLCoefficient;
  using Kokkos::Example::FENL::ExponentialKLCoefficient;
  using Kokkos::Example::FENL::Problem;
  using Kokkos::Example::BoxElemPart;
  using Kokkos::Example::FENL::fenl;
  using Kokkos::Example::FENL::Perf;

  const double bc_lower_value = 1 ;
  const double bc_upper_value = 2 ;
  const double geom_bubble[3] = { 1.0 , 1.0 , 1.0 };

  const double kl_mean = cmd.USE_MEAN;
  const double kl_variance = cmd.USE_VAR;
  const double kl_correlation = cmd.USE_COR;
  const int kl_dim = cmd.USE_UQ_DIM ;
  const bool kl_exp = cmd.USE_EXPONENTIAL;
  const double kl_exp_shift = cmd.USE_EXP_SHIFT;
  const double kl_exp_scale = cmd.USE_EXP_SCALE;

  int nelem[3] = { cmd.USE_FIXTURE_X,
                   cmd.USE_FIXTURE_Y,
                   cmd.USE_FIXTURE_Z};
  Perf perf_total;

  // Read in any params from xml file
  Teuchos::RCP<Teuchos::ParameterList> fenlParams = Teuchos::parameterList();
  Teuchos::updateParametersFromXmlFileAndBroadcast(
    cmd.USE_FENL_XML_FILE, fenlParams.ptr(), *comm);

  // Compute PCE of response propagating blocks of quadrature
  // points at a time
  if ( cmd.USE_UQ_ENSEMBLE > 0 ) {

    typedef Stokhos::StaticFixedStorage<int,double,VectorSize,Device> Storage;
    typedef Sacado::MP::Vector<Storage> Scalar;

    // Set global vector size -- this is mandatory
    Kokkos::global_sacado_mp_vector_size = VectorSize;

    typedef ExponentialKLCoefficient< Scalar, double, Device > KL;
    KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, kl_dim,
                              kl_exp, kl_exp_shift, kl_exp_scale );

    // Problem setup
    typedef Problem< Scalar, Device , BoxElemPart::ElemLinear > ProblemType;
    ProblemType problem( comm , node , nelem , geom_bubble , cmd.PRINT );

    // Grouping method
    RCP< Kokkos::Example::FENL::SampleGrouping<double> > grouper;
    if (cmd.USE_GROUPING == GROUPING_NATURAL)
      grouper = rcp(new Kokkos::Example::FENL::NaturalGrouping<double>);
    else if (cmd.USE_GROUPING == GROUPING_MAX_ANISOTROPY) {
      typedef ExponentialKLCoefficient< double, double, Device > DKL;
      DKL diff_coeff( kl_mean, kl_variance, kl_correlation, kl_dim,
                      kl_exp, kl_exp_shift, kl_exp_scale);
      typedef typename ProblemType::FixtureType Mesh;
      grouper =
        rcp(new Kokkos::Example::FENL::MaxAnisotropyGrouping<double,Mesh,DKL>(
              comm, problem.fixture, diff_coeff));
    }

    if (cmd.USE_UQ_SAMPLING == SAMPLING_STOKHOS)
      run_stokhos(*comm, problem, diffusion_coefficient, grouper,
                  fenlParams, cmd, bc_lower_value, bc_upper_value, perf_total);
    else if (cmd.USE_UQ_SAMPLING == SAMPLING_TASMANIAN)
      run_tasmanian(*comm, problem, diffusion_coefficient, grouper,
                    fenlParams, cmd, bc_lower_value, bc_upper_value, perf_total);
    else if (cmd.USE_UQ_SAMPLING == SAMPLING_FILE)
      run_file(*comm, problem, diffusion_coefficient, grouper,
               fenlParams, cmd, bc_lower_value, bc_upper_value, perf_total);

  }

  // Compute PCE of response propagating one quadrature point at a time
  else {

    typedef double Scalar;
    typedef ExponentialKLCoefficient< Scalar, double, Device > KL;
    KL diffusion_coefficient( kl_mean, kl_variance, kl_correlation, kl_dim,
                              kl_exp, kl_exp_shift, kl_exp_scale );

    // Problem setup
    Problem< Scalar, Device , BoxElemPart::ElemLinear > problem(
      comm , node , nelem , geom_bubble , cmd.PRINT );

    if (cmd.USE_UQ_SAMPLING == SAMPLING_STOKHOS)
      run_stokhos(*comm, problem, diffusion_coefficient, Teuchos::null,
                  fenlParams, cmd, bc_lower_value, bc_upper_value, perf_total);
    else if (cmd.USE_UQ_SAMPLING == SAMPLING_TASMANIAN)
      run_tasmanian(*comm, problem, diffusion_coefficient, Teuchos::null,
                    fenlParams, cmd, bc_lower_value, bc_upper_value, perf_total);
    else if (cmd.USE_UQ_SAMPLING == SAMPLING_FILE)
      run_file(*comm, problem, diffusion_coefficient, Teuchos::null,
               fenlParams, cmd, bc_lower_value, bc_upper_value, perf_total);

  }

  if ( 0 == comm_rank ) {
    print_perf_value( std::cout , cmd , widths , perf_total );
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
  clp_return_type rv = parse_cmdline( argc, argv, cmdline, *comm, true );
  if (rv==CLP_HELP)
    return(EXIT_SUCCESS);
  else if (rv==CLP_ERROR)
    return(EXIT_FAILURE);

  if ( cmdline.VTUNE  ) {
    connect_vtune(comm->getRank());
  }

  if ( ! cmdline.ERROR  && ! cmdline.ECHO  ) {

#if defined( HAVE_TPETRA_SERIAL )
    if ( cmdline.USE_SERIAL ) {
#if defined(__MIC__)
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::Serial , 16 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::Serial , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
#else
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 4 )
        run< Kokkos::Serial ,  4 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 8 )
        run< Kokkos::Serial ,  8 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::Serial , 16 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::Serial , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
#endif
    }
#endif

#if defined( HAVE_TPETRA_PTHREAD )
    if ( cmdline.USE_THREADS ) {
#if defined(__MIC__)
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 16 )
        run< Kokkos::Threads , 16 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 32 )
        run< Kokkos::Threads , 32 >( comm , cmdline );
      else
        std::cout << "Invalid ensemble size!" << std::endl;
#else
      if ( cmdline.USE_UQ_ENSEMBLE == 0 ||
           cmdline.USE_UQ_ENSEMBLE == 4 )
        run< Kokkos::Threads ,  4 >( comm , cmdline );
      else if ( cmdline.USE_UQ_ENSEMBLE == 8 )
        run< Kokkos::Threads ,  8 >( comm , cmdline );
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
      else if ( cmdline.USE_UQ_ENSEMBLE == 8 )
        run< Kokkos::OpenMP ,  8 >( comm , cmdline );
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

// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <fenl_utils.hpp>

#include <utility>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#if defined( KOKKOS_ENABLE_CUDA )
#include <cuda_runtime_api.h>
#endif

// For vtune
#include <sys/types.h>
#include <unistd.h>

// For memory proviling
#include <sys/time.h>
#include <sys/resource.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_CommHelpers.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Command line processing:

// Parse command line
clp_return_type parse_cmdline( int argc , char ** argv, CMD & cmdline,
                               const Teuchos::Comm<int>& comm,
                               const bool uq )
{
  Teuchos::ParameterList params;
  Teuchos::CommandLineProcessor clp(false);

  const int num_grouping_types = 4;
  const GroupingType grouping_values[] = {
    GROUPING_NATURAL, GROUPING_MAX_ANISOTROPY, GROUPING_MORTAN_Z,
    GROUPING_TASMANIAN_SURROGATE };
  const char *grouping_names[] = { "natural", "max-anisotropy", "mortan-z", "tasmanian-surrogate" };

  const int num_sampling_types = 4;
  const SamplingType sampling_values[] = {
    SAMPLING_STOKHOS, SAMPLING_TASMANIAN, SAMPLING_FILE, SAMPLING_VPS };
  const char *sampling_names[] = { "stokhos", "tasmanian", "file", "vps" };

  clp.setOption("serial", "no-serial",     &cmdline.USE_SERIAL, "use the serial device");
  clp.setOption("threads",                 &cmdline.USE_THREADS, "number of pthreads threads");
  clp.setOption("openmp",                  &cmdline.USE_OPENMP,  "number of openmp threads");
  clp.setOption("numa",                    &cmdline.USE_NUMA,  "number of numa nodes");
  clp.setOption("cores",                   &cmdline.USE_CORE_PER_NUMA, "cores per numa node");
  clp.setOption("cuda", "no-cuda",         &cmdline.USE_CUDA,  "use the CUDA device");
  clp.setOption("device",                  &cmdline.USE_CUDA_DEV,  "CUDA device ID.  Set to default of -1 to use the default device as determined by the local node MPI rank and --ngpus");
  clp.setOption("ngpus",                   &cmdline.USE_NGPUS, "Number of GPUs per node for multi-GPU runs via MPI");
  std::string fixtureSpec="2x2x2";
  clp.setOption("fixture",                 &fixtureSpec,  "fixture string: \"XxYxZ\"");
  clp.setOption("fixture-x",               &cmdline.USE_FIXTURE_X,  "fixture");
  clp.setOption("fixture-y",               &cmdline.USE_FIXTURE_Y,  "fixture");
  clp.setOption("fixture-z",               &cmdline.USE_FIXTURE_Z,  "fixture");
  clp.setOption("fixture-quadratic", "no-fixture-quadratic", &cmdline.USE_FIXTURE_QUADRATIC,  "quadratic");

  clp.setOption("atomic", "no-atomic",      &cmdline.USE_ATOMIC ,  "atomic");
  clp.setOption("trials",                   &cmdline.USE_TRIALS,  "trials");
  clp.setOption("xml-file",                 &cmdline.USE_FENL_XML_FILE, "XML file containing solver parameters");
  clp.setOption("belos", "no-belos",        &cmdline.USE_BELOS ,  "use Belos solver");
  clp.setOption("muelu", "no-muelu",        &cmdline.USE_MUELU,  "use MueLu preconditioner");
  clp.setOption("mean-based", "no-mean-based", &cmdline.USE_MEANBASED,  "use mean-based preconditioner");
  if(cmdline.USE_MUELU || cmdline.USE_MEANBASED)
    cmdline.USE_BELOS = true;

  clp.setOption("sampling", &cmdline.USE_UQ_SAMPLING, num_sampling_types, sampling_values, sampling_names, "UQ sampling method");
  clp.setOption("uq-fake",                  &cmdline.USE_UQ_FAKE,  "setup a fake UQ problem of this size");
  clp.setOption("uq-dim",                   &cmdline.USE_UQ_DIM,  "UQ dimension");
  clp.setOption("uq-order",                 &cmdline.USE_UQ_ORDER,  "UQ order");
  clp.setOption("uq-init-level",            &cmdline.USE_UQ_INIT_LEVEL,  "Initial adaptive sparse grid level");
  clp.setOption("uq-max-level",             &cmdline.USE_UQ_MAX_LEVEL,  "Max adaptive sparse grid level");
  clp.setOption("uq-max-samples",           &cmdline.USE_UQ_MAX_SAMPLES,  "Max number of samples to run");
  clp.setOption("uq-tol",                   &cmdline.USE_UQ_TOL,  "Adaptive sparse grid tolerance");
  clp.setOption("diff-coeff-linear",        &cmdline.USE_DIFF_COEFF_LINEAR,  "Linear term in diffusion coefficient");
  clp.setOption("diff-coeff-constant",      &cmdline.USE_DIFF_COEFF_CONSTANT,  "Constant term in diffusion coefficient");
  clp.setOption("mean",                     &cmdline.USE_MEAN,  "KL diffusion mean");
  clp.setOption("var",                      &cmdline.USE_VAR,  "KL diffusion variance");
  clp.setOption("cor",                      &cmdline.USE_COR,  "KL diffusion correlation");
  clp.setOption("exponential", "no-exponential", &cmdline.USE_EXPONENTIAL,  "take exponential of KL diffusion coefficient");
  clp.setOption("exp-shift",                &cmdline.USE_EXP_SHIFT,  "Linear shift of exponential of KL diffusion coefficient");
  clp.setOption("exp-scale",                &cmdline.USE_EXP_SCALE,  "Multiplicative scale of exponential of KL diffusion coefficient");
  clp.setOption("discontinuous-exp-scale", "continuous-exp-scale", &cmdline.USE_DISC_EXP_SCALE,  "use discontinuous scale factor on exponential");
  clp.setOption("isotropic", "anisotropic", &cmdline.USE_ISOTROPIC,  "use isotropic or anisotropic diffusion coefficient");
  clp.setOption("coeff-src",                &cmdline.USE_COEFF_SRC,  "Coefficient for source term");
  clp.setOption("coeff-adv",                &cmdline.USE_COEFF_ADV,  "Coefficient for advection term");
  clp.setOption("sparse", "tensor",         &cmdline.USE_SPARSE ,  "use sparse or tensor grid");
  clp.setOption("ensemble",                 &cmdline.USE_UQ_ENSEMBLE,  "UQ ensemble size.  This needs to be a valid choice based on available instantiations.");
  clp.setOption("grouping", &cmdline.USE_GROUPING, num_grouping_types, grouping_values, grouping_names, "Sample grouping method for ensemble propagation");
  clp.setOption("surrogate-grouping-level", &cmdline.TAS_GROUPING_INITIAL_LEVEL,  "Starting level for surrogate-based grouping");

  clp.setOption("vtune", "no-vtune",       &cmdline.VTUNE ,  "connect to vtune");
  clp.setOption("verbose", "no-verbose",   &cmdline.VERBOSE, "print verbose intialization info");
  clp.setOption("print", "no-print",        &cmdline.PRINT,  "print detailed test output");
  clp.setOption("print-its", "no-print-its",&cmdline.PRINT_ITS,  "print solver iterations after each sample");
  clp.setOption("summarize", "no-summarize",&cmdline.SUMMARIZE,  "summarize Teuchos timers at end of run");
  clp.setOption("unit-test", "no-unit-test",&cmdline.UNIT_TEST,  "whether code is running as a unit-test and should check values");
  clp.setOption("test-mean",                &cmdline.TEST_MEAN,  "test value for mean when running as unit-test");
  clp.setOption("test-std-dev",             &cmdline.TEST_STD_DEV,  "test value for standard deviation when running as unit-test");
  clp.setOption("test-tol",                 &cmdline.TEST_TOL,  "tolerance for unit testing");

  bool doDryRun = false;
  clp.setOption("echo", "no-echo",          &doDryRun,  "dry-run only");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return CLP_HELP;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return CLP_ERROR;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

#if defined( KOKKOS_ENABLE_CUDA )
  // Set CUDA device based on local node rank
  if (cmdline.USE_CUDA && cmdline.USE_CUDA_DEV == -1) {
    int local_rank = 0;
    char *str;
    if ((str = std::getenv("MV2_COMM_WORLD_LOCAL_RANK")))
      local_rank = std::atoi(str);
    else if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK")))
      local_rank = std::atoi(str);
    else if ((str = std::getenv("SLURM_LOCALID")))
      local_rank = std::atoi(str);
    cmdline.USE_CUDA_DEV = local_rank % cmdline.USE_NGPUS;

    // Check device is valid
    int num_device; cudaGetDeviceCount(&num_device);
    TEUCHOS_TEST_FOR_EXCEPTION(
      cmdline.USE_CUDA_DEV >= cmdline.USE_NGPUS, std::logic_error,
      "Invalid device ID " << cmdline.USE_CUDA_DEV << ".  You probably are trying" <<
      " to run with too many GPUs per node");
  }
#endif

  sscanf( fixtureSpec.c_str() , "%dx%dx%d" ,
          &cmdline.USE_FIXTURE_X ,
          &cmdline.USE_FIXTURE_Y ,
          &cmdline.USE_FIXTURE_Z );

  cmdline.USE_UQ = uq;

  if (doDryRun) {
    print_cmdline( std::cout , cmdline );
    cmdline.ECHO  = 1;
  } else {
    cmdline.ECHO  = 0;
  }
  cmdline.ERROR  = 0 ;

  return CLP_OK;

}

// Print command line
void print_cmdline( std::ostream & s , const CMD & cmd )
{
  if ( cmd.USE_SERIAL  ) {
    s << " Serial" ;
  }
  if ( cmd.USE_THREADS  ) {
    s << " Threads(" << cmd.USE_THREADS
      << ") NUMA(" << cmd.USE_NUMA
      << ") CORE_PER_NUMA(" << cmd.USE_CORE_PER_NUMA
      << ")" ;
  }
  if ( cmd.USE_OPENMP  ) {
    s << " OpenMP(" << cmd.USE_OPENMP
      << ") NUMA(" << cmd.USE_NUMA
      << ") CORE_PER_NUMA(" << cmd.USE_CORE_PER_NUMA
      << ")" ;
  }
  if ( cmd.USE_FIXTURE_X  ) {
    s << " Fixture(" << cmd.USE_FIXTURE_X
      << "x" << cmd.USE_FIXTURE_Y
      << "x" << cmd.USE_FIXTURE_Z
      << ")" ;
  }
  if ( cmd.USE_FIXTURE_QUADRATIC  ) {
    s << " Quadratic-Element" ;
  }
  if ( cmd.USE_UQ_ENSEMBLE  ) {
    s << " UQ ensemble(" << cmd.USE_UQ_ENSEMBLE << ")" ;
  }
  if ( cmd.USE_UQ_FAKE ) {
    s << " UQ fake(" << cmd.USE_UQ_FAKE  << ")" ;
  }
  if ( cmd.USE_UQ_DIM  ) {
    s << " UQ dimension(" << cmd.USE_UQ_DIM  << ")" ;
  }
  if ( cmd.USE_UQ_ORDER  ) {
    s << " UQ order(" << cmd.USE_UQ_ORDER  << ")" ;
  }
  if ( cmd.USE_DIFF_COEFF_LINEAR ) {
    s << " Diffusion Coefficient A(" << cmd.USE_DIFF_COEFF_LINEAR << ")" ;
  }
  if ( cmd.USE_DIFF_COEFF_CONSTANT ) {
    s << " Diffusion Coefficient B(" << cmd.USE_DIFF_COEFF_CONSTANT << ")" ;
  }
  if ( cmd.USE_VAR  ) {
    s << " KL variance(" << cmd.USE_VAR << ")" ;
  }
  if ( cmd.USE_MEAN  ) {
    s << " KL mean(" << cmd.USE_MEAN << ")" ;
  }
  if ( cmd.USE_COR  ) {
    s << " KL correlation(" << cmd.USE_COR << ")" ;
  }
  if ( cmd.USE_EXPONENTIAL ) {
    s << " KL exponential(" << cmd.USE_EXPONENTIAL << ")" ;
  }
  if ( cmd.USE_EXP_SHIFT ) {
    s << " KL exponential shift(" << cmd.USE_EXP_SHIFT << ")" ;
  }
  if ( cmd.USE_EXP_SCALE ) {
    s << " KL exponential scale(" << cmd.USE_EXP_SCALE << ")" ;
  }
  if ( cmd.USE_ISOTROPIC ) {
    s << " isotropic" ;
  }
  if ( !cmd.USE_ISOTROPIC ) {
    s << " anisotropic" ;
  }
  if ( cmd.USE_COEFF_SRC  ) {
    s << " Source coefficient(" << cmd.USE_COEFF_SRC  << ")" << ")" ;
  }
  if ( cmd.USE_COEFF_ADV  ) {
    s << " Advection coefficient" << cmd.USE_COEFF_ADV  << ")" << ")" ;
  }
  if ( cmd.USE_SPARSE  ) {
    s << " Sparse grid" ;
  }
  if ( !cmd.USE_SPARSE  ) {
    s << " Tensor grid" ;
  }
  if ( cmd.USE_CUDA  ) {
    s << " CUDA(" << cmd.USE_CUDA_DEV  << ")" ;
  }
  if ( cmd.USE_ATOMIC  ) {
    s << " ATOMIC" ;
  }
  if ( cmd.USE_TRIALS  ) {
    s << " TRIALS(" << cmd.USE_TRIALS  << ")" ;
  }
  if ( cmd.USE_BELOS  ) {
    s << " BELOS" ;
  }
  if ( cmd.USE_MUELU  ) {
    s << " MUELU" ;
  }
  if ( cmd.USE_MEANBASED  ) {
    s << " MEAN-BASED" ;
  }
  if ( cmd.VTUNE  ) {
    s << " VTUNE" ;
  }
  if ( cmd.VERBOSE  ) {
    s << " VERBOSE" ;
  }
  if ( cmd.PRINT  ) {
    s << " PRINT" ;
  }
  if ( cmd.PRINT_ITS  ) {
    s << " PRINT_ITS" ;
  }
  if ( cmd.SUMMARIZE  ) {
    s << " SUMMARIZE" ;
  }
  s << std::endl ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Display performance:

// Print timing headers
std::vector< size_t >
print_headers( std::ostream & s , const CMD & cmd , const int comm_rank )
{
  if ( 0 == comm_rank ) {
   if ( cmd.USE_SERIAL  ) { s << "SERIAL"  ; }
   if ( cmd.USE_THREADS  ) { s << "THREADS , " << cmd.USE_THREADS  ; }
   else if ( cmd.USE_OPENMP  ) { s << "OPENMP , " << cmd.USE_OPENMP  ; }
   else if ( cmd.USE_CUDA  ) { s << "CUDA" ; }

   if ( cmd.USE_FIXTURE_QUADRATIC  ) { s << " , QUADRATIC-ELEMENT" ; }
   else { s << " , LINEAR-ELEMENT" ; }

   s << " , FIXTURE , "
     << cmd.USE_FIXTURE_X << "x"
     << cmd.USE_FIXTURE_Y << "x"
     << cmd.USE_FIXTURE_Z ;

   s << " , SOURCE COEFFICIENT " << cmd.USE_COEFF_SRC
     << " , ADVECTION COEFFICIENT " << cmd.USE_COEFF_ADV ;

   if ( !cmd.USE_UQ ) {
     s << " , DIFFUSION COEFFICIENT( "
       << cmd.USE_DIFF_COEFF_LINEAR << " , "
       << cmd.USE_DIFF_COEFF_CONSTANT << " )" ;
     if ( cmd.USE_ISOTROPIC )
       s << " ISOTROPIC" ;
     else
       s << " ANISOTROPIC" ;
   }

   if ( cmd.USE_ATOMIC ) { s << " , USING ATOMICS" ; }
   if ( cmd.USE_BELOS  ) { s << " , USING BELOS" ; }
   if ( cmd.USE_MUELU  ) { s << " , USING MUELU" ; }

   if ( cmd.USE_UQ ) {
     s << " , KL MEAN , " << cmd.USE_MEAN ;
     s << " , KL VAR , " << cmd.USE_VAR ;
     s << " , KL COR , " << cmd.USE_COR ;
     s << " , KL EXP , " << cmd.USE_EXPONENTIAL ;
     s << " , KL EXP SHIFT, " << cmd.USE_EXP_SHIFT ;
     s << " , KL EXP SCALE, " << cmd.USE_EXP_SCALE ;
     if ( cmd.USE_ISOTROPIC )
       s << " ISOTROPIC" ;
     else
       s << " ANISOTROPIC" ;
     if ( cmd.USE_UQ_FAKE ) {
       s << " , UQ FAKE , " << cmd.USE_UQ_FAKE ;
     }
     else {
       s << " , UQ DIM , " << cmd.USE_UQ_DIM ;
       s << " , UQ ORDER , " << cmd.USE_UQ_ORDER;
       if ( cmd.USE_SPARSE  ) { s << " , USING SPARSE GRID" ; }
       else { s << " , USING TENSOR GRID" ; }
     }
     if ( cmd.USE_UQ_ENSEMBLE  ) {
       s << " , USING UQ ENSEMBLE , " << cmd.USE_UQ_ENSEMBLE ;
     }
     if ( cmd.USE_MEANBASED  ) { s << " , MEAN-BASED PREC" ; }
   }

  }

  std::vector< std::pair<std::string,std::string> > headers;

  headers.push_back(std::make_pair("ELEMS","count"));
  headers.push_back(std::make_pair("NODES","count"));
  if ( cmd.USE_UQ  ) {
    headers.push_back(std::make_pair("SAMPLES","count"));
    headers.push_back(std::make_pair("NEWTON","iter"));
    headers.push_back(std::make_pair("SOLVER","iter"));
    headers.push_back(std::make_pair("IMPORT/NODE","millisec"));
    headers.push_back(std::make_pair("MATRIX_FILL/NODE","millisec"));
    headers.push_back(std::make_pair("BOUNDARY/NODE","millisec"));
    headers.push_back(std::make_pair("MAT_VEC/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("SOLVER/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("PREC/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("PREC SETUP","millisec"));
    headers.push_back(std::make_pair("SOLVER TOTAL","millisec"));
    headers.push_back(std::make_pair("RESPONSE","mean"));
    headers.push_back(std::make_pair("RESPONSE","std.dev."));
  }
  else {
    headers.push_back(std::make_pair("NEWTON","iter"));
    headers.push_back(std::make_pair("SOLVER","iter"));
    headers.push_back(std::make_pair("MAP_RATIO","ratio"));
    headers.push_back(std::make_pair("SET_FILL/NODE","millisec"));
    headers.push_back(std::make_pair("SCAN/NODE","millisec"));
    headers.push_back(std::make_pair("GRAPH_FILL/NODE","millisec"));
    headers.push_back(std::make_pair("SORT/NODE","millisec"));
    headers.push_back(std::make_pair("ELEM_GRAPH_FILL/NODE","millisec"));
    headers.push_back(std::make_pair("MATRIX_CREATE/NODE","millisec"));
    headers.push_back(std::make_pair("MATRIX_FILL/NODE","millisec"));
    headers.push_back(std::make_pair("BOUNDARY/NODE","millisec"));
    headers.push_back(std::make_pair("MAT_VEC/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("SOLVER/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("RESIDUAL ERROR","ratio"));
  }

  // find print widths
  size_t min_width = 10;
  std::vector< size_t > widths(headers.size());
  for (size_t i=0, ie=headers.size(); i<ie; ++i)
    widths[i] = std::max(min_width, headers[i].first.size()+1);

  // print column headers
  if ( 0 == comm_rank ) {
    s << std::endl ;
    for (size_t i=0; i<headers.size(); ++i)
      s << std::setw(widths[i]) << headers[i].first << " ,";
    s << "\b\b  " << std::endl;
    for (size_t i=0; i<headers.size(); ++i)
      s << std::setw(widths[i]) << headers[i].second << " ,";
    s << "\b\b  " << std::endl;
  }

  return widths;
}

// Print times
void print_perf_value( std::ostream & s ,
                       const CMD & cmd ,
                       const std::vector<size_t> & widths ,
                       const Kokkos::Example::FENL::Perf & perf )
{
  int i=0;
  s << std::scientific;
  s.precision(3);
  s << std::setw(widths[i++]) << perf.global_elem_count << " ,";
  s << std::setw(widths[i++]) << perf.global_node_count << " ,";

  if ( cmd.USE_UQ  ) {
    // Note:  cg_iter_count is already a sum across all samples,
    // so don't scale cg times by uq_count
    s << std::setw(widths[i++]) << perf.uq_count << " ,";
    s << std::setw(widths[i++]) << double(perf.newton_iter_count) / perf.uq_count << " ,";
    s << std::setw(widths[i++]) << double(perf.cg_iter_count) / perf.uq_count << " ,";
    s << std::setw(widths[i++]) << ( perf.import_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( perf.fill_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( perf.bc_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( ( perf.mat_vec_time * 1000.0 ) / perf.cg_iter_count ) / (perf.global_node_count) << " ,";
    s << std::setw(widths[i++]) << ( ( perf.cg_iter_time * 1000.0 ) / perf.cg_iter_count ) / (perf.global_node_count) << " ,";
    s << std::setw(widths[i++]) << ( ( perf.prec_apply_time * 1000.0 ) / perf.cg_iter_count ) / (perf.global_node_count) << " ,";
    s << std::setw(widths[i++]) << ( perf.prec_setup_time * 1000.0 ) / (perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( perf.cg_total_time * 1000.0 ) / (perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << perf.response_mean << " ,";
    s << std::setw(widths[i])   << perf.response_std_dev;
  }
  else {
    s << std::setw(widths[i++]) << double(perf.newton_iter_count) / perf.uq_count << " ,";
  s << std::setw(widths[i++]) << double(perf.cg_iter_count) / perf.uq_count << " ,";
    s << std::setw(widths[i++]) << perf.map_ratio << " ,";
    s << std::setw(widths[i++]) << ( perf.fill_node_set * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( perf.scan_node_count * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( perf.fill_graph_entries * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( perf.sort_graph_entries * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( perf.fill_element_graph * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( perf.create_sparse_matrix * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( perf.fill_time * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( perf.bc_time * 1000.0 ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( ( perf.mat_vec_time * 1000.0 ) / perf.cg_iter_count ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << ( ( perf.cg_iter_time * 1000.0 ) / perf.cg_iter_count ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << perf.error_max << " ,";
  }
    s << std::endl ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Profiling:

// Connect executable to vtune for profiling
void connect_vtune(const int p_rank) {
  std::stringstream cmd;
  pid_t my_os_pid=getpid();
  const std::string vtune_loc =
    "amplxe-cl";
  const std::string output_dir = "./vtune/vtune.";
  cmd << vtune_loc
      << " -collect hotspots -result-dir " << output_dir << p_rank
      << " -target-pid " << my_os_pid << " &";
  if (p_rank == 0)
    std::cout << cmd.str() << std::endl;
  system(cmd.str().c_str());
  system("sleep 10");
}

// Get memory usage in MB
MemUsage get_memory_usage(const Teuchos::Comm<int>& comm) {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  size_t mem = usage.ru_maxrss;
#if defined(__APPLE__)
  mem /= 1024; // Apple returns bytes instead of kilobytes
#endif

  size_t max_mem = 0;
  size_t min_mem = 0;
  size_t tot_mem = 0;
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &mem, &max_mem);
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, 1, &mem, &min_mem);
  Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &mem, &tot_mem);

  MemUsage mem_usage;
  mem_usage.max_mem = static_cast<double>(max_mem) / 1024.0;
  mem_usage.min_mem = static_cast<double>(min_mem) / 1024.0;
  mem_usage.tot_mem = static_cast<double>(tot_mem) / 1024.0;

  return mem_usage;
}

// Print memory usage to stream
void print_memory_usage(std::ostream& s, const Teuchos::Comm<int>& comm) {
  MemUsage mem =  get_memory_usage(comm);
  if ( 0 == comm.getRank() ) {
    s << std::fixed;
    s.precision(3);
    s << "Memory usage across all processors (MB):" << std::endl
      << "\t Max:  " << mem.max_mem << std::endl
      << "\t Min:  " << mem.min_mem << std::endl
      << "\t Tot:  " << mem.tot_mem << std::endl;
  }
}

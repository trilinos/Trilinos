#include <fenl_utils.hpp>

#include <utility>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

#if defined( KOKKOS_HAVE_CUDA )
#include <cuda_runtime_api.h>
#endif

// For vtune
#include <sys/types.h>
#include <unistd.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TestForException.hpp>

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

  clp.setOption("threads",                 &cmdline.CMD_USE_THREADS, "number of pthreads threads");
  clp.setOption("openmp",                  &cmdline.CMD_USE_OPENMP,  "number of openmp threads");
  clp.setOption("numa",                    &cmdline.CMD_USE_NUMA,  "number of numa nodes");
  clp.setOption("cores",                   &cmdline.CMD_USE_CORE_PER_NUMA, "cores per numa node");
  clp.setOption("cuda", "no-cuda",         &cmdline.CMD_USE_CUDA,  "use CUDA");
  clp.setOption("device",                  &cmdline.CMD_USE_CUDA_DEV,  "CUDA device ID.  Set to default of -1 to use the default device as determined by the local node MPI rank and --ngpus");
  clp.setOption("ngpus",                   &cmdline.CMD_USE_NGPUS, "Number of GPUs per node for multi-GPU runs via MPI");
  std::string fixtureSpec="2x2x2";
  clp.setOption("fixture",                 &fixtureSpec,  "fixture string: \"XxYxZ\"");
  clp.setOption("fixture-x",               &cmdline.CMD_USE_FIXTURE_X,  "fixture");
  clp.setOption("fixture-y",               &cmdline.CMD_USE_FIXTURE_Y,  "fixture");
  clp.setOption("fixture-z",               &cmdline.CMD_USE_FIXTURE_Z,  "fixture");

  std::string fixtureRange;
  clp.setOption("fixture-range",           &fixtureRange,  "fixture range: \"x..y\"");
  clp.setOption("fixture-begin",           &cmdline.CMD_USE_FIXTURE_BEGIN,  "fixture begin");
  clp.setOption("fixture-end",             &cmdline.CMD_USE_FIXTURE_END,  "fixture end");
  clp.setOption("fixture-quadratic", "no-fixture-quadratic", &cmdline.CMD_USE_FIXTURE_QUADRATIC,  "quadratic");

  clp.setOption("atomic", "no-atomic",      &cmdline.CMD_USE_ATOMIC ,  "atomic");
  clp.setOption("trials",                   &cmdline.CMD_USE_TRIALS,  "trials");
  clp.setOption("belos", "no-belos",        &cmdline.CMD_USE_BELOS ,  "use Belos solver");
  clp.setOption("muelu", "no-muelu",        &cmdline.CMD_USE_MUELU,  "use MueLu preconditioner");
  clp.setOption("mean-based", "no-mean-based", &cmdline.CMD_USE_MEANBASED,  "use mean-based preconditioner");
  if(cmdline.CMD_USE_MUELU || cmdline.CMD_USE_MEANBASED)
    cmdline.CMD_USE_BELOS = true;

  clp.setOption("uq-fake",                     &cmdline.CMD_USE_UQ_FAKE,  "setup a fake UQ problem of this size");
  clp.setOption("uq-dim",                   &cmdline.CMD_USE_UQ_DIM,  "UQ dimension");
  clp.setOption("uq-order",                 &cmdline.CMD_USE_UQ_ORDER,  "UQ order");
  clp.setOption("mean",                     &cmdline.CMD_USE_MEAN,  "KL diffusion mean");
  clp.setOption("var",                      &cmdline.CMD_USE_VAR,  "KL diffusion variance");
  clp.setOption("cor",                      &cmdline.CMD_USE_COR,  "KL diffusion correlation");
  clp.setOption("sparse", "tensor",         &cmdline.CMD_USE_SPARSE ,  "use sparse or tensor grid");
  clp.setOption("ensemble",                 &cmdline.CMD_USE_UQ_ENSEMBLE,  "UQ ensemble size.  This needs to be a valid choice based on available instantiations.");

  clp.setOption("vtune", "no-vtune",       &cmdline.CMD_VTUNE ,  "connect to vtune");
  clp.setOption("verbose", "no-verbose",   &cmdline.CMD_VERBOSE, "print verbose intialization info");
  clp.setOption("print", "no-print",        &cmdline.CMD_PRINT,  "print detailed test output");
  clp.setOption("summarize", "no-summarize",&cmdline.CMD_SUMMARIZE,  "summarize Teuchos timers at end of run");

  bool doDryRun = false;
  clp.setOption("echo", "no-echo",          &doDryRun,  "dry-run only");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return CLP_HELP;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return CLP_ERROR;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

#if defined( KOKKOS_HAVE_CUDA )
  // Set CUDA device based on local node rank
  if (cmdline.CMD_USE_CUDA && cmdline.CMD_USE_CUDA_DEV == -1) {
    int local_rank = 0;
    char *str;
    if ((str = std::getenv("MV2_COMM_WORLD_LOCAL_RANK")))
      local_rank = std::atoi(str);
    else if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK")))
      local_rank = std::atoi(str);
    else if ((str = std::getenv("SLURM_LOCALID")))
      local_rank = std::atoi(str);
    cmdline.CMD_USE_CUDA_DEV = local_rank % cmdline.CMD_USE_NGPUS;

    // Check device is valid
    int num_device; cudaGetDeviceCount(&num_device);
    TEUCHOS_TEST_FOR_EXCEPTION(
      cmdline.CMD_USE_CUDA_DEV >= cmdline.CMD_USE_NGPUS, std::logic_error,
      "Invalid device ID " << cmdline.CMD_USE_CUDA_DEV << ".  You probably are trying" <<
      " to run with too many GPUs per node");
  }
#endif

  sscanf( fixtureSpec.c_str() , "%dx%dx%d" ,
          &cmdline.CMD_USE_FIXTURE_X ,
          &cmdline.CMD_USE_FIXTURE_Y ,
          &cmdline.CMD_USE_FIXTURE_Z );
  sscanf( fixtureRange.c_str(), "%d..%d" ,
          &cmdline.CMD_USE_FIXTURE_BEGIN ,
          &cmdline.CMD_USE_FIXTURE_END );

  cmdline.CMD_USE_UQ = uq;

  if (doDryRun) {
    print_cmdline( std::cout , cmdline );
    cmdline.CMD_ECHO  = 1;
  } else {
    cmdline.CMD_ECHO  = 0;
  }
  cmdline.CMD_ERROR  = 0 ;

  return CLP_OK;

}

// Print command line
void print_cmdline( std::ostream & s , const CMD & cmd )
{
  if ( cmd.CMD_USE_THREADS  ) {
    s << " Threads(" << cmd.CMD_USE_THREADS
      << ") NUMA(" << cmd.CMD_USE_NUMA
      << ") CORE_PER_NUMA(" << cmd.CMD_USE_CORE_PER_NUMA
      << ")" ;
  }
  if ( cmd.CMD_USE_OPENMP  ) {
    s << " OpenMP(" << cmd.CMD_USE_OPENMP
      << ") NUMA(" << cmd.CMD_USE_NUMA
      << ") CORE_PER_NUMA(" << cmd.CMD_USE_CORE_PER_NUMA
      << ")" ;
  }
  if ( cmd.CMD_USE_FIXTURE_X  ) {
    s << " Fixture(" << cmd.CMD_USE_FIXTURE_X
      << "x" << cmd.CMD_USE_FIXTURE_Y
      << "x" << cmd.CMD_USE_FIXTURE_Z
      << ")" ;
  }
  if ( cmd.CMD_USE_FIXTURE_BEGIN  ) {
    s << " Fixture( " << cmd.CMD_USE_FIXTURE_BEGIN
      << " .. " << cmd.CMD_USE_FIXTURE_END
      << " )" ;
  }
  if ( cmd.CMD_USE_FIXTURE_QUADRATIC  ) {
    s << " Quadratic-Element" ;
  }
  if ( cmd.CMD_USE_UQ_ENSEMBLE  ) {
    s << " UQ ensemble(" << cmd.CMD_USE_UQ_ENSEMBLE << ")" ;
  }
  if ( cmd.CMD_USE_UQ_FAKE ) {
    s << " UQ fake(" << cmd.CMD_USE_UQ_FAKE  << ")" ;
  }
  if ( cmd.CMD_USE_UQ_DIM  ) {
    s << " UQ dimension(" << cmd.CMD_USE_UQ_DIM  << ")" ;
  }
  if ( cmd.CMD_USE_UQ_ORDER  ) {
    s << " UQ order(" << cmd.CMD_USE_UQ_ORDER  << ")" ;
  }
  if ( cmd.CMD_USE_VAR  ) {
    s << " KL variance(" << cmd.CMD_USE_VAR  << ")" ;
  }
  if ( cmd.CMD_USE_MEAN  ) {
    s << " KL mean(" << cmd.CMD_USE_MEAN  << ")" ;
  }
  if ( cmd.CMD_USE_COR  ) {
    s << " KL correlation(" << cmd.CMD_USE_COR  << ")" << ")" ;
  }
  if ( cmd.CMD_USE_SPARSE  ) {
    s << " Sparse grid" ;
  }
  if ( !cmd.CMD_USE_SPARSE  ) {
    s << " Tensor grid" ;
  }
  if ( cmd.CMD_USE_CUDA  ) {
    s << " CUDA(" << cmd.CMD_USE_CUDA_DEV  << ")" ;
  }
  if ( cmd.CMD_USE_ATOMIC  ) {
    s << " ATOMIC" ;
  }
  if ( cmd.CMD_USE_TRIALS  ) {
    s << " TRIALS(" << cmd.CMD_USE_TRIALS  << ")" ;
  }
  if ( cmd.CMD_USE_BELOS  ) {
    s << " BELOS" ;
  }
  if ( cmd.CMD_USE_MUELU  ) {
    s << " MUELU" ;
  }
  if ( cmd.CMD_USE_MEANBASED  ) {
    s << " MEAN-BASED" ;
  }
  if ( cmd.CMD_VTUNE  ) {
    s << " VTUNE" ;
  }
  if ( cmd.CMD_VERBOSE  ) {
    s << " VERBOSE" ;
  }
  if ( cmd.CMD_PRINT  ) {
    s << " PRINT" ;
  }
  if ( cmd.CMD_SUMMARIZE  ) {
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
   if ( cmd.CMD_USE_THREADS  ) { s << "THREADS , " << cmd.CMD_USE_THREADS  ; }
   else if ( cmd.CMD_USE_OPENMP  ) { s << "OPENMP , " << cmd.CMD_USE_OPENMP  ; }
   else if ( cmd.CMD_USE_CUDA  ) { s << "CUDA" ; }

   if ( cmd.CMD_USE_FIXTURE_QUADRATIC  ) { s << " , QUADRATIC-ELEMENT" ; }
   else { s << " , LINEAR-ELEMENT" ; }

   s << " , FIXTURE , "
     << cmd.CMD_USE_FIXTURE_X << "x"
     << cmd.CMD_USE_FIXTURE_Y << "x"
     << cmd.CMD_USE_FIXTURE_Z ;

   if ( cmd.CMD_USE_ATOMIC ) { s << " , USING ATOMICS" ; }
   if ( cmd.CMD_USE_BELOS  ) { s << " , USING BELOS" ; }
   if ( cmd.CMD_USE_MUELU  ) { s << " , USING MUELU" ; }

   if ( cmd.CMD_USE_UQ ) {
     s << " , KL MEAN , " << cmd.CMD_USE_MEAN ;
     s << " , KL VAR , " << cmd.CMD_USE_VAR ;
     s << " , KL COR , " << cmd.CMD_USE_COR ;
     if ( cmd.CMD_USE_UQ_FAKE ) {
       s << " , UQ FAKE , " << cmd.CMD_USE_UQ_FAKE ;
     }
     else {
       s << " , UQ DIM , " << cmd.CMD_USE_UQ_DIM ;
       s << " , UQ ORDER , " << cmd.CMD_USE_UQ_ORDER;
       if ( cmd.CMD_USE_SPARSE  ) { s << " , USING SPARSE GRID" ; }
       else { s << " , USING TENSOR GRID" ; }
     }
     if ( cmd.CMD_USE_UQ_ENSEMBLE  ) {
       s << " , USING UQ ENSEMBLE , " << cmd.CMD_USE_UQ_ENSEMBLE ;
     }
     if ( cmd.CMD_USE_MEANBASED  ) { s << " , MEAN-BASED PREC" ; }
   }

  }

  std::vector< std::pair<std::string,std::string> > headers;

  headers.push_back(std::make_pair("ELEMS","count"));
  headers.push_back(std::make_pair("NODES","count"));
  if ( cmd.CMD_USE_UQ  ) {
    headers.push_back(std::make_pair("SAMPLES","count"));
    headers.push_back(std::make_pair("NEWTON","iter"));
    headers.push_back(std::make_pair("CG","iter"));
    headers.push_back(std::make_pair("IMPORT/NODE","millisec"));
    headers.push_back(std::make_pair("MATRIX_FILL/NODE","millisec"));
    headers.push_back(std::make_pair("BOUNDARY/NODE","millisec"));
    headers.push_back(std::make_pair("MAT_VEC/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("CG/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("PREC/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("PREC SETUP","millisec"));
    headers.push_back(std::make_pair("CG TOTAL","millisec"));
    headers.push_back(std::make_pair("RESPONSE","mean"));
    headers.push_back(std::make_pair("RESPONSE","std.dev."));
  }
  else {
    headers.push_back(std::make_pair("NEWTON","iter"));
    headers.push_back(std::make_pair("CG","iter"));
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
    headers.push_back(std::make_pair("CG/ITER/ROW","millisec"));
    headers.push_back(std::make_pair("ERROR","ratio"));
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

    s << std::scientific;
    s.precision(3);
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
  s << std::setw(widths[i++]) << perf.global_elem_count << " ,";
  s << std::setw(widths[i++]) << perf.global_node_count << " ,";

  if ( cmd.CMD_USE_UQ  ) {
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
    "/usr/local/intel/vtune_amplifier_xe_2013/bin64/amplxe-cl";
  const std::string output_dir = "./vtune/vtune.";
  cmd << vtune_loc
      << " -collect hotspots -result-dir " << output_dir << p_rank
      << " -target-pid " << my_os_pid << " &";
  if (p_rank == 0)
    std::cout << cmd.str() << std::endl;
  system(cmd.str().c_str());
  system("sleep 10");
}

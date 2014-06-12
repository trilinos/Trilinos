#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include <utility>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>

//----------------------------------------------------------------------------
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Command line processing:

enum clp_return_type {CLP_HELP=0,
      CLP_ERROR,
      CLP_OK};

enum { CMD_USE_THREADS = 0
       , CMD_USE_NUMA
       , CMD_USE_CORE_PER_NUMA
       , CMD_USE_CUDA
       , CMD_USE_OPENMP
       , CMD_USE_CUDA_DEV
       , CMD_USE_FIXTURE_X
       , CMD_USE_FIXTURE_Y
       , CMD_USE_FIXTURE_Z
       , CMD_USE_FIXTURE_BEGIN
       , CMD_USE_FIXTURE_END
       , CMD_USE_FIXTURE_QUADRATIC
       , CMD_USE_ATOMIC
       , CMD_USE_TRIALS
       , CMD_USE_BELOS
       , CMD_USE_MUELU
       , CMD_USE_UQ_ENSEMBLE
       , CMD_USE_UQ_DIM
       , CMD_USE_UQ_ORDER
       , CMD_USE_SPARSE
       , CMD_PRINT
       , CMD_SUMMARIZE
       , CMD_ECHO
       , CMD_ERROR
       , CMD_COUNT };

void print_cmdline( std::ostream & s , const int cmd[] )
{
  if ( cmd[ CMD_USE_THREADS ] ) {
    s << " Threads(" << cmd[ CMD_USE_THREADS ]
      << ") NUMA(" << cmd[ CMD_USE_NUMA ]
      << ") CORE_PER_NUMA(" << cmd[ CMD_USE_CORE_PER_NUMA ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_OPENMP ] ) {
    s << " OpenMP(" << cmd[ CMD_USE_OPENMP ]
      << ") NUMA(" << cmd[ CMD_USE_NUMA ]
      << ") CORE_PER_NUMA(" << cmd[ CMD_USE_CORE_PER_NUMA ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_X ] ) {
    s << " Fixture(" << cmd[ CMD_USE_FIXTURE_X ]
      << "x" << cmd[ CMD_USE_FIXTURE_Y ]
      << "x" << cmd[ CMD_USE_FIXTURE_Z ]
      << ")" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_BEGIN ] ) {
    s << " Fixture( " << cmd[ CMD_USE_FIXTURE_BEGIN ]
      << " .. " << cmd[ CMD_USE_FIXTURE_END ]
      << " )" ;
  }
  if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) {
    s << " Quadratic-Element" ;
  }
  if ( cmd[ CMD_USE_UQ_ENSEMBLE ] ) {
    s << " UQ ensemble" ;
  }
  if ( cmd[ CMD_USE_UQ_DIM ] ) {
    s << " UQ dimension(" << cmd[ CMD_USE_UQ_DIM ] << ")" << ")" ;
  }
  if ( cmd[ CMD_USE_UQ_DIM ] ) {
    s << " UQ order(" << cmd[ CMD_USE_UQ_ORDER ] << ")" << ")" ;
  }
  if ( cmd[ CMD_USE_SPARSE ] ) {
    s << " Sparse grid" ;
  }
  if ( !cmd[ CMD_USE_SPARSE ] ) {
    s << " Tensor grid" ;
  }
  if ( cmd[ CMD_USE_CUDA ] ) {
    s << " CUDA(" << cmd[ CMD_USE_CUDA_DEV ] << ")" ;
  }
  if ( cmd[ CMD_USE_ATOMIC ] ) {
    s << " ATOMIC" ;
  }
  if ( cmd[ CMD_USE_TRIALS ] ) {
    s << " TRIALS(" << cmd[ CMD_USE_TRIALS ] << ")" ;
  }
  if ( cmd[ CMD_USE_BELOS ] ) {
    s << " BELOS" ;
  }
  if ( cmd[ CMD_USE_MUELU ] ) {
    s << " MUELU" ;
  }
  if ( cmd[ CMD_PRINT ] ) {
    s << " PRINT" ;
  }
  if ( cmd[ CMD_SUMMARIZE ] ) {
    s << " SUMMARIZE" ;
  }
  s << std::endl ;
}

std::vector< size_t >
print_headers( std::ostream & s , const int cmd[] , const int comm_rank )
{
  if ( 0 == comm_rank ) {
    if ( cmd[ CMD_USE_THREADS ] ) { s << "THREADS , " << cmd[ CMD_USE_THREADS ] ; }
    else if ( cmd[ CMD_USE_OPENMP ] ) { s << "OPENMP , " << cmd[ CMD_USE_OPENMP ] ; }
    else if ( cmd[ CMD_USE_CUDA ] ) { s << "CUDA" ; }

    if ( cmd[ CMD_USE_FIXTURE_QUADRATIC ] ) { s << " , QUADRATIC-ELEMENT" ; }
    else { s << " , LINEAR-ELEMENT" ; }

    if ( cmd[ CMD_USE_ATOMIC ] ) { s << " , USING ATOMICS" ; }
    if ( cmd[ CMD_USE_BELOS ] ) { s << " , USING BELOS" ; }
    if ( cmd[ CMD_USE_MUELU ] ) { s << " , USING MUELU" ; }
    if ( cmd[ CMD_USE_UQ_ENSEMBLE ] ) { s << " , USING UQ ENSEMBLE" ; }
    if ( cmd[ CMD_USE_UQ_DIM ] ) { s << " , UQ DIM , " << cmd[ CMD_USE_UQ_DIM ]; }
    if ( cmd[ CMD_USE_UQ_ORDER ] ) { s << " , UQ ORDER , " << cmd[ CMD_USE_UQ_ORDER ]; }

    if ( cmd[ CMD_USE_SPARSE ] ) { s << " , USING SPARSE GRID" ; }
    else { s << " , USING TENSOR GRID" ; }
  }

  std::vector< std::pair<std::string,std::string> > headers;

  headers.push_back(std::make_pair("ELEMS","count"));
  headers.push_back(std::make_pair("NODES","count"));
  if ( cmd[ CMD_USE_UQ_DIM ] ) {
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

void print_perf_value( std::ostream & s ,
                       const int cmd[] ,
                       const std::vector<size_t> & widths ,
                       const Kokkos::Example::FENL::Perf & perf )
{
  int i=0;
  s << std::setw(widths[i++]) << perf.global_elem_count << " ,";
  s << std::setw(widths[i++]) << perf.global_node_count << " ,";

  if ( cmd[ CMD_USE_UQ_DIM ] ) {
    s << std::setw(widths[i++]) << perf.uq_count << " ,";
    s << std::setw(widths[i++]) << double(perf.newton_iter_count) / perf.uq_count << " ,";
  s << std::setw(widths[i++]) << double(perf.cg_iter_count) / perf.uq_count << " ,";
    s << std::setw(widths[i++]) << ( perf.import_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( perf.fill_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( perf.bc_time * 1000.0 ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( ( perf.mat_vec_time * 1000.0 ) / perf.cg_iter_count ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( ( perf.cg_iter_time * 1000.0 ) / perf.cg_iter_count ) / (perf.global_node_count*perf.uq_count) << " ,";
    s << std::setw(widths[i++]) << ( ( perf.prec_apply_time * 1000.0 ) / perf.cg_iter_count ) / (perf.global_node_count*perf.uq_count) << " ,";
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
    s << std::setw(widths[i++]) << ( ( perf.cg_iter_time * 1000.0 ) / perf.cg_iter_count ) / perf.global_node_count << " ,";
    s << std::setw(widths[i++]) << perf.error_max << " ,";
  }
    s << std::endl ;
}

clp_return_type parse_cmdline( int argc , char ** argv, int cmdline[],
                               const Teuchos::Comm<int>& comm )
{
  Teuchos::ParameterList params;
  Teuchos::CommandLineProcessor clp(false);
  cmdline[CMD_USE_THREADS] = 0;       clp.setOption("threads",                 cmdline+CMD_USE_THREADS,  "number of pthreads threads");

  cmdline[CMD_USE_OPENMP] = 0;        clp.setOption("openmp",                 cmdline+CMD_USE_OPENMP,  "number of openmp threads");

  cmdline[CMD_USE_NUMA] = 0;          clp.setOption("numa",                    cmdline+CMD_USE_NUMA,  "number of numa nodes");

  cmdline[CMD_USE_CORE_PER_NUMA] = 0; clp.setOption("cores",                   cmdline+CMD_USE_CORE_PER_NUMA,  "cores per numa node");

  bool useCuda = false;               clp.setOption("cuda", "no-cuda",         &useCuda,  "use CUDA");

  bool useCudaDev = false;            clp.setOption("cuda-dev", "no-cuda-dev", &useCudaDev,  "use CUDA dev");

  std::string fixtureSpec="2x2x2";    clp.setOption("fixture",                 &fixtureSpec,  "fixture string: \"XxYxZ\"");
                                      clp.setOption("fixture-x",               cmdline+CMD_USE_FIXTURE_X,  "fixture");
                                      clp.setOption("fixture-y",               cmdline+CMD_USE_FIXTURE_Y,  "fixture");
                                      clp.setOption("fixture-z",               cmdline+CMD_USE_FIXTURE_Z,  "fixture");

  std::string fixtureRange;           clp.setOption("fixture-range",           &fixtureRange,  "fixture range: \"x..y\"");
  cmdline[CMD_USE_FIXTURE_BEGIN]=0;   clp.setOption("fixture-begin",           cmdline+CMD_USE_FIXTURE_BEGIN,  "fixture begin");
  cmdline[CMD_USE_FIXTURE_END]=0;     clp.setOption("fixture-end",             cmdline+CMD_USE_FIXTURE_END,  "fixture end");

  bool useQuadratic = false;          clp.setOption("fixture-quadratic", "no-fixture-quadratic", &useQuadratic,  "quadratic");

  bool useEnsemble = false;           clp.setOption("ensemble", "no-ensemble",  &useEnsemble,  "use ensemble");

  cmdline[CMD_USE_UQ_DIM] = 0;        clp.setOption("uq-dim",                   cmdline+CMD_USE_UQ_DIM,  "UQ dimension");

  cmdline[CMD_USE_UQ_ORDER] = 0;      clp.setOption("uq-order",                 cmdline+CMD_USE_UQ_ORDER,  "UQ order");

  bool useSparse = false;           clp.setOption("sparse", "tensor",  &useSparse,  "use sparse or tensor grid");

  bool useAtomic = false;             clp.setOption("atomic", "no-atomic",      &useAtomic,  "atomic");

  cmdline[CMD_USE_TRIALS] = 1;        clp.setOption("trials",                   cmdline+CMD_USE_TRIALS,  "trials");

  bool useBelos = false;              clp.setOption("belos", "no-belos",        &useBelos,  "use Belos solver");

  bool useMueLu = false;              clp.setOption("muelu", "no-muelu",        &useMueLu,  "use MueLu preconditioner");

  bool doPrint = false;               clp.setOption("print", "no-print",        &doPrint,  "print detailed test output");

  bool doSummarize = false;               clp.setOption("summarize", "no-summarize",        &doSummarize,  "summarize Teuchos timers at end of run");

  bool doDryRun = false;              clp.setOption("echo", "no-echo",          &doDryRun,  "dry-run only");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return CLP_HELP;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return CLP_ERROR;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  // cmdline is of type int*, which means we can't use it directly in CommandLineProcessor::setOptions for bools.
  cmdline[CMD_USE_CUDA]              = useCuda;
  cmdline[CMD_USE_CUDA_DEV]          = useCudaDev;
  if (useCudaDev)
    cmdline[CMD_USE_CUDA] = true;
  cmdline[CMD_USE_FIXTURE_QUADRATIC] = useQuadratic;
  cmdline[CMD_USE_UQ_ENSEMBLE]       = useEnsemble;
  cmdline[CMD_USE_SPARSE]            = useSparse;
  cmdline[CMD_USE_ATOMIC]            = useAtomic;
  cmdline[CMD_USE_BELOS]             = useBelos;
  cmdline[CMD_USE_MUELU]             = useMueLu;
  cmdline[CMD_PRINT]                 = doPrint;
  cmdline[CMD_SUMMARIZE]             = doSummarize;
  sscanf( fixtureSpec.c_str() , "%dx%dx%d" ,
          cmdline + CMD_USE_FIXTURE_X ,
          cmdline + CMD_USE_FIXTURE_Y ,
          cmdline + CMD_USE_FIXTURE_Z );
  sscanf( fixtureRange.c_str(), "%d..%d" ,
          cmdline + CMD_USE_FIXTURE_BEGIN ,
          cmdline + CMD_USE_FIXTURE_END );

  if (doDryRun) {
    print_cmdline( std::cout , cmdline );
    cmdline[ CMD_ECHO ] = 1;
  } else {
    cmdline[ CMD_ECHO ] = 0;
  }
  cmdline[ CMD_ERROR ] = 0 ;

  //--------------------------------------------------------------------------

  comm.broadcast( int(0) , int(CMD_COUNT * sizeof(int)) , (char *) cmdline );

  return CLP_OK;

}

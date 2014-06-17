#include <iostream>
#include <vector>
#include <Teuchos_Comm.hpp>

#include "fenl.hpp"

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
       , CMD_VTUNE
       , CMD_PRINT
       , CMD_SUMMARIZE
       , CMD_ECHO
       , CMD_ERROR
       , CMD_COUNT };

// Parse command line
clp_return_type parse_cmdline( int argc , char ** argv, int cmdline[],
                               const Teuchos::Comm<int>& comm );

// Print command line
void print_cmdline( std::ostream & s , const int cmd[] );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Display performance:

// Print timing headers
std::vector< size_t >
print_headers( std::ostream & s , const int cmd[] , const int comm_rank );

// Print times
void print_perf_value( std::ostream & s ,
                       const int cmd[] ,
                       const std::vector<size_t> & widths ,
                       const Kokkos::Example::FENL::Perf & perf );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Profiling:

// Connect executable to vtune for profiling
void connect_vtune(const int p_rank);

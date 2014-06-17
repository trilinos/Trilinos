#include <iostream>
#include <vector>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Kokkos_DefaultNode.hpp>

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
       , CMD_VERBOSE
       , CMD_SUMMARIZE
       , CMD_ECHO
       , CMD_ERROR
       , CMD_COUNT };

// Parse command line
clp_return_type parse_cmdline( int argc , char ** argv, int cmdline[],
                               const Teuchos::Comm<int>& comm );

// Print command line
void print_cmdline( std::ostream & s , const int cmd[] );

// Create Tpetra node
template <typename NodeType>
Teuchos::RCP<NodeType>
createKokkosNode( const int cmd[] , const int comm_rank ) {
  Teuchos::ParameterList params;
  params.set("Verbose", 0);
  if ( cmd[ CMD_USE_THREADS ] )
    params.set("Num Threads", cmd[CMD_USE_THREADS]);
  else if ( cmd[ CMD_USE_OPENMP ] )
    params.set("Num Threads", cmd[CMD_USE_OPENMP]);
  if ( cmd[ CMD_USE_NUMA ] && cmd[ CMD_USE_CORE_PER_NUMA ] ) {
    params.set("Num NUMA", cmd[ CMD_USE_NUMA ]);
    params.set("Num CoresPerNUMA", cmd[ CMD_USE_CORE_PER_NUMA ]);
  }
  if ( cmd[ CMD_USE_CUDA ] )
    params.set("Device", cmd[ CMD_USE_CUDA_DEV ] );
  Teuchos::RCP<NodeType> node = Teuchos::rcp (new NodeType(params));

  if ( cmd[CMD_VERBOSE] ) {
    typedef typename NodeType::device_type Device;
    if (comm_rank == 0)
      Device::print_configuration(std::cout);
    if ( cmd[ CMD_USE_CUDA ] )
      std::cout << "MPI rank " << comm_rank << " attached to CUDA device "
                << cmd[ CMD_USE_CUDA_DEV ] << std::endl;
  }

  return node;
}

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

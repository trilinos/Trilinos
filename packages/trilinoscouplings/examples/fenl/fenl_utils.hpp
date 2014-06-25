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

struct CMD {
  int CMD_USE_THREADS;
  int CMD_USE_OPENMP;
  int CMD_USE_NUMA;
  int CMD_USE_CORE_PER_NUMA;
  bool CMD_USE_CUDA;
  int CMD_USE_CUDA_DEV;
  int CMD_USE_NGPUS;
  int CMD_USE_FIXTURE_X;
  int CMD_USE_FIXTURE_Y;
  int CMD_USE_FIXTURE_Z;
  int CMD_USE_FIXTURE_BEGIN;
  int CMD_USE_FIXTURE_END;
  bool CMD_USE_FIXTURE_QUADRATIC;
  bool CMD_USE_ATOMIC;
  int CMD_USE_TRIALS;
  bool CMD_USE_BELOS;
  bool CMD_USE_MUELU;
  bool CMD_USE_MEANBASED;
  bool CMD_USE_UQ;
  int CMD_USE_UQ_FAKE;
  int CMD_USE_UQ_DIM;
  int CMD_USE_UQ_ORDER;
  double CMD_USE_MEAN;
  double CMD_USE_VAR;
  double CMD_USE_COR;
  bool CMD_USE_SPARSE;
  int CMD_USE_UQ_ENSEMBLE;
  bool CMD_VTUNE;
  bool CMD_VERBOSE;
  bool CMD_PRINT;
  bool CMD_SUMMARIZE;
  int CMD_ECHO;
  int CMD_ERROR;
  int CMD_COUNT;

  CMD() : CMD_USE_THREADS(0),
          CMD_USE_OPENMP(0),
          CMD_USE_NUMA(0),
          CMD_USE_CORE_PER_NUMA(0),
          CMD_USE_CUDA(false),
          CMD_USE_CUDA_DEV(-1),
          CMD_USE_NGPUS(1),
          CMD_USE_FIXTURE_X(2),
          CMD_USE_FIXTURE_Y(2),
          CMD_USE_FIXTURE_Z(2),
          CMD_USE_FIXTURE_BEGIN(0),
          CMD_USE_FIXTURE_END(0),
          CMD_USE_FIXTURE_QUADRATIC(false),
          CMD_USE_ATOMIC(false),
          CMD_USE_TRIALS(1),
          CMD_USE_BELOS(false),
          CMD_USE_MUELU(false),
          CMD_USE_MEANBASED(false),
          CMD_USE_UQ(false),
          CMD_USE_UQ_FAKE(0),
          CMD_USE_UQ_DIM(3),
          CMD_USE_UQ_ORDER(2),
          CMD_USE_MEAN(1),
          CMD_USE_VAR(0.1),
          CMD_USE_COR(0.25),
          CMD_USE_SPARSE(false),
          CMD_USE_UQ_ENSEMBLE(0),
          CMD_VTUNE(false),
          CMD_VERBOSE(false),
          CMD_PRINT(false),
          CMD_SUMMARIZE(false)
    {}
};

// Print command line
void print_cmdline( std::ostream & s , const CMD & cmd );

// Create Tpetra node
template <typename NodeType>
Teuchos::RCP<NodeType>
createKokkosNode( const CMD & cmd , const Teuchos::Comm<int>& comm ) {
  Teuchos::ParameterList params;
  params.set("Verbose", 0);
  if ( cmd.CMD_USE_THREADS  )
    params.set("Num Threads", cmd.CMD_USE_THREADS);
  else if ( cmd.CMD_USE_OPENMP  )
    params.set("Num Threads", cmd.CMD_USE_OPENMP);
  if ( cmd.CMD_USE_NUMA  && cmd.CMD_USE_CORE_PER_NUMA  ) {
    params.set("Num NUMA", cmd.CMD_USE_NUMA );
    params.set("Num CoresPerNUMA", cmd.CMD_USE_CORE_PER_NUMA );
  }
  if ( cmd.CMD_USE_CUDA  )
    params.set("Device", cmd.CMD_USE_CUDA_DEV  );
  Teuchos::RCP<NodeType> node = Teuchos::rcp (new NodeType(params));

  if ( cmd.CMD_VERBOSE ) {
    typedef typename NodeType::device_type Device;
    if (comm.getRank() == 0)
      Device::print_configuration(std::cout);
    std::cout.flush();
    if ( cmd.CMD_USE_CUDA  ) {
      for (int i=0; i<comm.getSize(); ++i) {
        comm.barrier();
        comm.barrier();
        comm.barrier();
        if ( i == comm.getRank() ) {
          std::cout << "MPI rank " << comm.getRank()
                    << " attached to CUDA device "
                    << cmd.CMD_USE_CUDA_DEV  << std::endl;
          std::cout.flush();
        }
        comm.barrier();
        comm.barrier();
        comm.barrier();
      }
    }
  }

  return node;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Display performance:

// Parse command line
clp_return_type parse_cmdline( int argc , char ** argv, CMD & cmdline,
                               const Teuchos::Comm<int>& comm,
                               const bool uq);


// Print timing headers
std::vector< size_t >
print_headers( std::ostream & s , const CMD & cmd , const int comm_rank );


// Print times
void print_perf_value( std::ostream & s ,
                       const CMD & cmd ,
                       const std::vector<size_t> & widths ,
                       const Kokkos::Example::FENL::Perf & perf );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Profiling:

// Connect executable to vtune for profiling
void connect_vtune(const int p_rank);

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

enum clp_return_type {
  CLP_HELP=0,
  CLP_ERROR,
  CLP_OK
};

enum GroupingType {
  GROUPING_NATURAL=0,
  GROUPING_MAX_ANISOTROPY,
  GROUPING_MORTAN_Z
};

enum SamplingType {
  SAMPLING_STOKHOS=0,
  SAMPLING_TASMANIAN,
  SAMPLING_FILE
};

struct CMD {
  bool USE_SERIAL;
  int USE_THREADS;
  int USE_OPENMP;
  int USE_NUMA;
  int USE_CORE_PER_NUMA;
  bool USE_CUDA;
  int USE_CUDA_DEV;
  int USE_NGPUS;
  int USE_FIXTURE_X;
  int USE_FIXTURE_Y;
  int USE_FIXTURE_Z;
  bool USE_FIXTURE_QUADRATIC;
  bool USE_ATOMIC;
  int USE_TRIALS;
  std::string USE_FENL_XML_FILE;
  bool USE_BELOS;
  bool USE_MUELU;
  bool USE_MEANBASED;
  bool USE_UQ;
  SamplingType USE_UQ_SAMPLING;
  int USE_UQ_FAKE;
  int USE_UQ_DIM;
  int USE_UQ_ORDER;
  int USE_UQ_INIT_LEVEL;
  int USE_UQ_MAX_LEVEL;
  double USE_UQ_TOL;
  double USE_DIFF_COEFF_LINEAR;
  double USE_DIFF_COEFF_CONSTANT;
  double USE_MEAN;
  double USE_VAR;
  double USE_COR;
  bool USE_EXPONENTIAL;
  double USE_EXP_SHIFT;
  double USE_EXP_SCALE;
  bool USE_ISOTROPIC;
  double USE_COEFF_SRC;
  double USE_COEFF_ADV;
  bool USE_SPARSE;
  int USE_UQ_ENSEMBLE;
  GroupingType USE_GROUPING;
  bool VTUNE;
  bool VERBOSE;
  bool PRINT;
  bool PRINT_ITS;
  bool SUMMARIZE;
  int ECHO;
  int ERROR;
  int COUNT;

  CMD() : USE_SERIAL(0),
          USE_THREADS(0),
          USE_OPENMP(0),
          USE_NUMA(0),
          USE_CORE_PER_NUMA(0),
          USE_CUDA(false),
          USE_CUDA_DEV(-1),
          USE_NGPUS(1),
          USE_FIXTURE_X(2),
          USE_FIXTURE_Y(2),
          USE_FIXTURE_Z(2),
          USE_FIXTURE_QUADRATIC(false),
          USE_ATOMIC(false),
          USE_TRIALS(1),
          USE_FENL_XML_FILE("fenl.xml"),
          USE_BELOS(false),
          USE_MUELU(false),
          USE_MEANBASED(false),
          USE_UQ(false),
          USE_UQ_SAMPLING(SAMPLING_STOKHOS),
          USE_UQ_FAKE(0),
          USE_UQ_DIM(3),
          USE_UQ_ORDER(2),
          USE_UQ_INIT_LEVEL(1),
          USE_UQ_MAX_LEVEL(7),
          USE_UQ_TOL(1e-2),
          USE_DIFF_COEFF_LINEAR(0.0),
          USE_DIFF_COEFF_CONSTANT(1.0),
          USE_MEAN(1),
          USE_VAR(0.1),
          USE_COR(0.25),
          USE_EXPONENTIAL(false),
          USE_EXP_SHIFT(1.0),
          USE_EXP_SCALE(1.0),
          USE_ISOTROPIC(true),
          USE_COEFF_SRC(1.0),
          USE_COEFF_ADV(0.0),
          USE_SPARSE(false),
          USE_UQ_ENSEMBLE(0),
          USE_GROUPING(GROUPING_NATURAL),
          VTUNE(false),
          VERBOSE(false),
          PRINT(false),
          PRINT_ITS(false),
          SUMMARIZE(false)
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
  if ( cmd.USE_THREADS  )
    params.set("Num Threads", cmd.USE_THREADS);
  else if ( cmd.USE_OPENMP  )
    params.set("Num Threads", cmd.USE_OPENMP);
  if ( cmd.USE_NUMA  && cmd.USE_CORE_PER_NUMA  ) {
    params.set("Num NUMA", cmd.USE_NUMA );
    params.set("Num CoresPerNUMA", cmd.USE_CORE_PER_NUMA );
  }
  if ( cmd.USE_CUDA  )
    params.set("Device", cmd.USE_CUDA_DEV  );
  Teuchos::RCP<NodeType> node = Teuchos::rcp (new NodeType(params));

  if ( cmd.VERBOSE ) {
    typedef typename NodeType::execution_space Device;
    if (comm.getRank() == 0)
      Device::print_configuration(std::cout);
    std::cout.flush();
    if ( cmd.USE_CUDA  ) {
      for (int i=0; i<comm.getSize(); ++i) {
        comm.barrier();
        comm.barrier();
        comm.barrier();
        if ( i == comm.getRank() ) {
          std::cout << "MPI rank " << comm.getRank()
                    << " attached to CUDA device "
                    << cmd.USE_CUDA_DEV  << std::endl;
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

// Get maximum memory usage
size_t get_max_memory_usage();

// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file block.cpp
    \brief An example of partitioning global ids with Block.
*/

#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

using namespace std;

/*! \example block.cpp
    An example of the use of the Block algorithm to partition data.
    \todo error handling
    \todo write some examples that don't use teuchos
*/

// Zoltan2 is templated.  What data types will we use for
// scalars (coordinate values and weights), for local ids, and
// for global ids?
//
// If Zoltan2 was compiled with explicit instantiation, we will
// use the the library's data types.  These macros are defined
// in Zoltan2_config.h.

#ifdef HAVE_ZOLTAN2_INST_FLOAT_INT_LONG
typedef float scalar_t;
typedef int localId_t;
typedef long globalId_t;
#else
  #ifdef HAVE_ZOLTAN2_INST_DOUBLE_INT_LONG
  typedef double scalar_t;
  typedef int localId_t;
  typedef long globalId_t;
  #else
    #ifdef HAVE_ZOLTAN2_INST_FLOAT_INT_INT
    typedef float scalar_t;
    typedef int localId_t;
    typedef int globalId_t;
    #else
      #ifdef HAVE_ZOLTAN2_INST_DOUBLE_INT_INT
      typedef double scalar_t;
      typedef int localId_t;
      typedef int globalId_t;
      #else
      typedef float scalar_t;
      typedef int localId_t;
      typedef int globalId_t;
      #endif
    #endif
  #endif
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_ZOLTAN2_MPI
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank=0, nprocs=1;
#endif

  ///////////////////////////////////////////////////////////////////////
  // Generate some input data.

  size_t localCount = 40*(rank+1);
  globalId_t *globalIds = new globalId_t [localCount];

  if (rank==0)
    for (int i=0, num=40; i <= nprocs ; i++, num+=40)
      cout << "Rank " << i << " has " << num << " ids." << endl;

  globalId_t offset = 0;
  for (int i=1; i <= rank; i++)
    offset += 40*i;

  for (localId_t i=0; i < localCount; i++)
    globalIds[i] = offset++;
   
  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 input adapter with no weights

  // TODO explain
  typedef Zoltan2::BasicUserTypes<scalar_t, globalId_t, localId_t, globalId_t> myTypes;

  // TODO explain
  typedef Zoltan2::BasicIdentifierInput<myTypes> inputAdapter_t;

  std::vector<const scalar_t *> noWeights;
  std::vector<int> noStrides;

  inputAdapter_t ia(localCount, globalIds, noWeights, noStrides);
   
  ///////////////////////////////////////////////////////////////////////
  // Create parameters for an Block problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("debug_procs", "0");
  params.set("error_check_level", "debug_mode_assertions");

  Teuchos::ParameterList &parParams = params.sublist("partitioning");
  parParams.set("algorithm", "block");
  parParams.set("imbalance_tolerance", 1.1);
  parParams.set("num_global_parts", nprocs);
   
  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 partitioning problem

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params, 
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif
   
  ///////////////////////////////////////////////////////////////////////
  // Solve the problem

  problem.solve();

  ///////////////////////////////////////////////////////////////////////
  // Check the solution.

  if (rank == 0)
    problem.printMetrics(cout);

  if (rank == 0)
    cout << "PASS" << endl;

  delete [] globalIds;
#ifdef HAVE_ZOLTAN2_MPI
  MPI_Finalize();
#endif
}


// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TrilinosCouplings_config.h"
#include "TrilinosCouplings_IntrepidPoissonExample_SolveWithBelos.hpp"
#include "TrilinosCouplings_TpetraIntrepidPoissonExample.hpp"
#include <BelosTpetraAdapter.hpp>

#include "Teuchos_Comm.hpp"

#include <cstdlib>

namespace TrilinosCouplings {
namespace TpetraIntrepidPoissonExample {

void
solveWithBelos (bool& converged,
                int& numItersPerformed,
                const std::string& solverName,
                const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
                const int maxNumIters,
                const int num_steps,
                const Teuchos::RCP<multivector_type>& X,
                const Teuchos::RCP<const sparse_matrix_type>& A,
                const Teuchos::RCP<const multivector_type>& B,
                const Teuchos::RCP<const operator_type>& M_left,
                const Teuchos::RCP<const operator_type>& M_right)
{
  typedef multivector_type MV;
  typedef operator_type OP;

  // Invoke the generic solve routine.
  IntrepidPoissonExample::solveWithBelos<ST, MV, OP> (converged, numItersPerformed,
                                                      solverName, tol, maxNumIters,
                                                      num_steps,
                                                      X, A, B, M_left, M_right);
}

/// \brief Solve the linear system(s) AX=B with Belos by cloning to a new
/// node type.
///
/// This works just like the above solveWithBelos() function except it clones
/// the matrix, LHS, RHS, and preconditioner(s) to the new node type.  This is
/// useful for solving the linear system on, e.g., a GPU.
template <typename CloneNode>
void
cloneAndSolveWithBelos (
  bool& converged,
  int& numItersPerformed,
  const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
  const int maxNumIters,
  const int num_steps,
  const Teuchos::RCP<multivector_type>& X,
  const Teuchos::RCP<const sparse_matrix_type>& A,
  const Teuchos::RCP<const multivector_type>& B,
  const std::string& prec_type,
  const Teuchos::RCP<const operator_type>& M_left=Teuchos::null,
  const Teuchos::RCP<const operator_type>& M_right=Teuchos::null) {

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef Tpetra::CrsMatrix<ST, LO, GO, CloneNode>    clone_sparse_matrix_type;
  typedef Tpetra::Operator<ST, LO, GO, CloneNode>     clone_operator_type;
  typedef Tpetra::MultiVector<ST, LO, GO, CloneNode>  clone_multi_vector_type;
  typedef clone_multi_vector_type MV;
  typedef clone_operator_type OP;

  // Clone Matrix, RHS, LHS
  RCP<ParameterList> plClone = parameterList();
  RCP<clone_sparse_matrix_type> A_clone;
  RCP<clone_multi_vector_type> B_clone, X_clone;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Clone System", clone_system);
    A_clone = rcp(new clone_sparse_matrix_type(A, Teuchos::Copy));
    B_clone = rcp(new clone_multi_vector_type(B, Teuchos::Copy));
    X_clone = rcp(new clone_multi_vector_type(X, Teuchos::Copy));
  }

  // Clone preconditioner(s)
  RCP<const clone_operator_type> M_left_clone, M_right_clone;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Clone Preconditioner", clone_prec);

#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
    if (M_left != Teuchos::null && prec_type == "MueLu") {
      RCP< const MueLu::TpetraOperator<ST,LO,GO,Node> > M_muelu =
        rcp_dynamic_cast<const MueLu::TpetraOperator<ST,LO,GO,Node> >(M_left);
      M_left_clone = rcp(new clone_operator_type(M_muelu));
    }
    if (M_right != Teuchos::null && prec_type == "MueLu") {
      RCP< const MueLu::TpetraOperator<ST,LO,GO,Node> > M_muelu =
        rcp_dynamic_cast<const MueLu::TpetraOperator<ST,LO,GO,Node> >(M_right);
      M_right_clone = rcp(new clone_operator_type(M_muelu));
    }
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      prec_type == "MueLu", std::runtime_error, "Tpetra scaling example: "
      "In order to precondition with MueLu, you must have built Trilinos "
      "with the MueLu package enabled.");
#endif // HAVE_TRILINOSCOUPLINGS_MUELU
  }

  // Solve
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Clone Solve", clone_solve);
    IntrepidPoissonExample::solveWithBelos<ST,MV,OP> (
      converged, numItersPerformed, tol, maxNumIters, num_steps,
      X_clone, A_clone, B_clone, M_left_clone, M_right_clone);
  }

  // Copy X_clone back into X
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Clone Solution", clone_sol);
    RCP<multivector_type> X_host = rcp(new multivector_type(X_clone, Teuchos::Copy));
    X->update(1.0, *X_host, 0.0);
  }
}

void
solveWithBelosGPU (
  bool& converged,
  int& numItersPerformed,
  const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
  const int maxNumIters,
  const int num_steps,
  const int ranks_per_node,
  const int gpu_ranks_per_node,
  const int device_offset,
  const std::string& prec_type,
  const Teuchos::RCP<multivector_type>& X,
  const Teuchos::RCP<const sparse_matrix_type>& A,
  const Teuchos::RCP<const multivector_type>& B,
  const Teuchos::RCP<const operator_type>& M_left,
  const Teuchos::RCP<const operator_type>& M_right) {


  // FIXME (mfh 11 Feb 2015) I am removing
  // KokkosClassic::ThrustGPUNode, so this example won't build.  We
  // need to test this example, though!  We should figure out how to
  // test it in a sensible way without needing ThrustGPUNode.
#if 0 // ifdef HAVE_KOKKOSCLASSIC_THRUST
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Compute whether we are a CPU or GPU rank, and GPU device ID
  // The first gpu_ranks_per_node ranks are associated with GPUs
  // (may not be the best choice performance-wise)
  RCP<const Comm<int> > comm = X->getMap()->getComm();
  int num_ranks = comm->getSize();
  int my_rank = comm->getRank();

  // Compute rank local to node
  int node_rank;
  char *str;
  if ((str = std::getenv("MV2_COMM_WORLD_LOCAL_RANK")) != NULL) {
    node_rank = std::atoi(str);
  }
  else if ((str = std::getenv("OMPI_COMM_WORLD_LOCAL_RANK")) != NULL) {
    node_rank = std::atoi(str);
  }
  else {
    // Default to user telling us the number of ranks_per_node and
    // assuming ranks are laid out consecutively on each node
    int num_node = num_ranks / ranks_per_node;
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_node*ranks_per_node != num_ranks, std::logic_error,
      "ranks_per_node does not evenly divide num_ranks");
    node_rank = num_node == 1 ? my_rank : my_rank % num_node;
  }
  bool gpu_rank = node_rank < gpu_ranks_per_node;
  int num_device; cudaGetDeviceCount(&num_device);
  int device_id = node_rank + device_offset;
  TEUCHOS_TEST_FOR_EXCEPTION(
    gpu_ranks_per_node > num_device, std::logic_error,
    "gpu_ranks_per_node cannot exceed number of GPU devices");
  TEUCHOS_TEST_FOR_EXCEPTION(
    gpu_rank && device_id > num_device, std::logic_error,
    "Invalid device ID " << device_id << ".  You probably are trying" <<
    " to run with too many MPI ranks");

  ParameterList node_params;
  if (gpu_rank) {
    std::cout << "MPI rank " << my_rank
              << ":  Attached to GPU " << device_id
              << std::endl;

#if TPETRA_USE_KOKKOS_DISTOBJECT && defined(KOKKOS_ENABLE_CUDA)
    if (!Kokkos::HostSpace::execution_space::is_initialized())
      Kokkos::HostSpace::execution_space::initialize();
    if (!Kokkos::Cuda::is_initialized())
      Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(device_id) );
#endif
    node_params.set("Verbose", 1);
    node_params.set("Device Number", device_id);
    RCP<KokkosClassic::ThrustGPUNode> gpu_node =
      rcp(new KokkosClassic::ThrustGPUNode(node_params));
    cloneAndSolveWithBelos(
      converged, numItersPerformed, tol, maxNumIters, num_steps,
      gpu_node, X, A, B, prec_type, M_left, M_right);
#if TPETRA_USE_KOKKOS_DISTOBJECT && defined(KOKKOS_ENABLE_CUDA)
    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
#endif
  }
  else {
    // Note for the non-GPU ranks, we still have to clone since new
    // matrices are created which call fillComplete() (which involves
    // communication)
    RCP<Node> node = X->getMap()->getNode();
    cloneAndSolveWithBelos(
      converged, numItersPerformed, tol, maxNumIters, num_steps,
      node, X, A, B, prec_type, M_left, M_right);
  }

#else // don't have ThrustGPUNode

  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "KokkosClassic::ThrustGPUNode has been deprecated "
    "and removed.  Please rewrite this function to use "
    "Tpetra::KokkosCompat::KokkosCudaWrapperNode instead.");

#endif // whether ThrustGPUNode exists
}

} // namespace TpetraIntrepidPoissonExample
} // namespace TrilinosCouplings
